/*  vcfannotate.c -- Annotate and edit VCF/BCF files.

    Copyright (C) 2013-2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include <dlfcn.h>
#include "bcftools.h"
#include "vcmp.h"
#include "filter.h"

typedef struct _plugin_t plugin_t;

/**
 *   Plugin API:
 *   ----------
 *   const char *about(void)
 *      - short description used by 'bcftools plugin -l'
 *
 *   const char *usage(void)
 *      - longer description used by 'bcftools +name -h'
 *
 *   int init(int argc, char **argv, bcf_hdr_t *in_hdr, bcf_hdr_t *out_hdr)
 *      - called once at startup, allows to initialize local variables.
 *      Return 1 to suppress normal VCF/BCF header output, -1 on critical
 *      errors, 0 otherwise.
 *
 *   bcf1_t *process(bcf1_t *rec)
 *      - called for each VCF record, return NULL for no output
 *
 *   void destroy(void)
 *      - called after all lines have been processed to clean up
 */
typedef void (*dl_version_f) (const char **, const char **);
typedef int (*dl_init_f) (int, char **, bcf_hdr_t *, bcf_hdr_t *);
typedef char* (*dl_about_f) (void);
typedef char* (*dl_usage_f) (void);
typedef bcf1_t* (*dl_process_f) (bcf1_t *);
typedef void (*dl_destroy_f) (void);

struct _plugin_t
{
    int argc;
    char *name, **argv;
    dl_version_f version;
    dl_init_f init;
    dl_about_f about;
    dl_usage_f usage;
    dl_process_f process;
    dl_destroy_f destroy;
    void *handle;
};


struct _args_t;

typedef struct _rm_tag_t
{
    char *key;
    int hdr_id;
    void (*handler)(struct _args_t *, bcf1_t *, struct _rm_tag_t *);
}
rm_tag_t;

typedef struct
{
    char **cols;
    int ncols, mcols;
    char **als;
    int nals, mals;
    kstring_t line;
    int rid, start, end;
}
annot_line_t;

typedef struct _annot_col_t
{
    int icol, replace;
    char *hdr_key;
    int (*setter)(struct _args_t *, bcf1_t *, struct _annot_col_t *, void*);
}
annot_col_t;

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int output_type;

    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE

    plugin_t plugin;
    int nplugin_paths;
    char **plugin_paths;

    char **argv, *output_fname, *regions_list, *targets_list;
    int argc, drop_header, verbose;
}
args_t;

char *msprintf(const char *fmt, ...);

static void init_plugin_paths(args_t *args)
{
    if ( args->nplugin_paths!=-1 ) return;

    char *path = getenv("BCFTOOLS_PLUGINS");
    if ( path )
    {
        args->nplugin_paths = 1;
        args->plugin_paths  = (char**) malloc(sizeof(char*));
        char *ss = args->plugin_paths[0] = strdup(path);
        while ( *ss )
        {
            if ( *ss==':' )
            {
                *ss = 0;
                args->plugin_paths = (char**) realloc(args->plugin_paths,sizeof(char*)*(args->nplugin_paths+1));
                args->plugin_paths[args->nplugin_paths] = ss+1;
                args->nplugin_paths++;
            }
            ss++;
        }
    }
    else
        args->nplugin_paths = 0;
}

static void *dlopen_plugin(args_t *args, const char *fname)
{
    init_plugin_paths(args);

    void *handle;
    char *tmp;
    if ( fname[0]!='/' )    // not an absolute path
    {
        int i;
        for (i=0; i<args->nplugin_paths; i++)
        {
            tmp = msprintf("%s/%s.so", args->plugin_paths[i],fname);
            handle = dlopen(tmp, RTLD_NOW); // valgrind complains about unfreed memory, not our problem though
            if ( args->verbose )
            {
                if ( !handle ) fprintf(stderr,"%s:\n\t%s\n", tmp,dlerror());
                else fprintf(stderr,"%s: ok\n", tmp);
            }
            free(tmp);
            if ( handle ) return handle;
        }
    }

    handle = dlopen(fname, RTLD_NOW);
    if ( args->verbose )
    {
        if ( !handle ) fprintf(stderr,"%s:\n\t%s\n", fname,dlerror());
        else fprintf(stderr,"%s: ok\n", fname);
    }

    return handle;
}

static void print_plugin_usage_hint(void)
{
    fprintf(stderr, "\nNo functional bcftools plugins were found");
    if ( !getenv("BCFTOOLS_PLUGINS") )
        fprintf(stderr,". The environment variable BCFTOOLS_PLUGINS is not set.\n\n");
    else
        fprintf(stderr,
                " in BCFTOOLS_PLUGINS=\"%s\".\n\n"
                "- Is the plugin path correct?\n\n"
                "- Are all shared libraries, namely libhts.so, accessible? Verify with\n"
                "   on Mac OS X: `otool -L your/plugin.so` and set DYLD_LIBRARY_PATH if they are not\n"
                "   on Linux:    `ldd your/plugin.so` and set LD_LIBRARY_PATH if they are not\n"
                "\n"
                "- If not installed systemwide, set the environment variable LD_LIBRARY_PATH (linux) or\n"
                "DYLD_LIBRARY_PATH (mac) to include directory where *libhts.so* is located.\n"
                "\n",
                getenv("BCFTOOLS_PLUGINS")
               );
}

static int load_plugin(args_t *args, const char *fname, int exit_on_error, plugin_t *plugin)
{
    plugin->name = strdup(fname);

    plugin->handle = dlopen_plugin(args, fname);
    if ( !plugin->handle )
    {
        if ( exit_on_error )
        {
            print_plugin_usage_hint();
            error("Could not load \"%s\".\n\n", fname);
        }
        return -1;
    }

    dlerror();
    plugin->init = (dl_init_f) dlsym(plugin->handle, "init");
    char *ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }

    plugin->version = (dl_version_f) dlsym(plugin->handle, "version");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }

    plugin->about = (dl_about_f) dlsym(plugin->handle, "about");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }

    plugin->usage = (dl_about_f) dlsym(plugin->handle, "usage");
    ret = dlerror();
    if ( ret )
        plugin->usage = plugin->about;

    plugin->process = (dl_process_f) dlsym(plugin->handle, "process");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }

    plugin->destroy = (dl_destroy_f) dlsym(plugin->handle, "destroy");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }

    return 0;
}

static void init_plugin(args_t *args)
{
    static int warned_bcftools = 0, warned_htslib = 0;

    int ret = args->plugin.init(args->plugin.argc,args->plugin.argv,args->hdr,args->hdr_out);
    if ( ret<0 ) error("The plugin exited with an error: %s\n", args->plugin.name);
    const char *bver, *hver;
    args->plugin.version(&bver, &hver);
    if ( strcmp(bver,bcftools_version()) && !warned_bcftools )
    {
        fprintf(stderr,"WARNING: bcftools version mismatch .. bcftools at %s, the plugin \"%s\" at %s\n", bcftools_version(),args->plugin.name,bver);
        warned_bcftools = 1;
    }
    if ( strcmp(hver,hts_version()) && !warned_htslib )
    {
        fprintf(stderr,"WARNING: htslib version mismatch .. bcftools at %s, the plugin \"%s\" at %s\n", hts_version(),args->plugin.name,hver);
        warned_htslib = 1;
    }
    args->drop_header += ret;
}

static int cmp_plugin_name(const void *p1, const void *p2)
{
    plugin_t *a = (plugin_t*) p1;
    plugin_t *b = (plugin_t*) p2;
    return strcmp(a->name,b->name);
}

static int list_plugins(args_t *args)
{
    plugin_t *plugins = NULL;
    int nplugins = 0, mplugins = 0;

    init_plugin_paths(args);

    kstring_t str = {0,0,0};
    int i;
    for (i=0; i<args->nplugin_paths; i++)
    {
        DIR *dp = opendir(args->plugin_paths[i]);
        if ( dp==NULL ) continue;

        struct dirent *ep;
        while ( (ep=readdir(dp)) )
        {
            int len = strlen(ep->d_name);
            if ( strcasecmp(".so",ep->d_name+len-3) ) continue;
            str.l = 0;
            ksprintf(&str,"%s/%s", args->plugin_paths[i],ep->d_name);
            hts_expand(plugin_t, nplugins+1, mplugins, plugins);
            if ( load_plugin(args, str.s, 0, &plugins[nplugins]) < 0 ) continue;
            nplugins++;
            str.l = 0;
            kputs(ep->d_name, &str);
            int l = str.l - 1;
            while ( l>=0 && str.s[l]!='.' ) l--;
            if ( l>=0 ) str.s[l] = 0;
            free(plugins[nplugins-1].name);
            plugins[nplugins-1].name = strdup(str.s);  // use a short name
        }
        closedir(dp);
    }
    if ( nplugins )
    {
        qsort(plugins, nplugins, sizeof(plugins[0]), cmp_plugin_name);

        for (i=0; i<nplugins; i++)
            printf("\n-- %s --\n%s", plugins[i].name, plugins[i].about());
        printf("\n");
    }
    else
        print_plugin_usage_hint();
    free(str.s);
    return nplugins ? 0 : 1;
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    init_plugin(args);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_plugin");
    if ( !args->drop_header )
    {
        args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
        if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
        bcf_hdr_write(args->out_fh, args->hdr_out);
    }
}

static void destroy_data(args_t *args)
{
    free(args->plugin.name);
    args->plugin.destroy();
    dlclose(args->plugin.handle);
    if ( args->hdr_out ) bcf_hdr_destroy(args->hdr_out);
    if ( args->nplugin_paths>0 )
    {
        free(args->plugin_paths[0]);
        free(args->plugin_paths);
    }
    if ( args->filter )
        filter_destroy(args->filter);
    if (args->out_fh) hts_close(args->out_fh);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Run user defined plugin\n");
    fprintf(stderr, "Usage:   bcftools plugin <name> [OPTIONS] <file>\n");
    fprintf(stderr, "         bcftools +name [OPTIONS] <file>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "VCF input options:\n");
    fprintf(stderr, "   -e, --exclude <expr>        exclude sites for which the expression is true\n");
    fprintf(stderr, "   -i, --include <expr>        select sites for which the expression is true\n");
    fprintf(stderr, "   -r, --regions <region>      restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file <file>   restrict to regions listed in a file\n");
    fprintf(stderr, "   -t, --targets <region>      similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "   -T, --targets-file <file>   similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "VCF output options:\n");
    fprintf(stderr, "   -o, --output <file>         write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <type>    'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "Plugin options:\n");
    fprintf(stderr, "   -h, --help                  list plugin's options\n");
    fprintf(stderr, "   -l, --list-plugins          list available plugins. See BCFTOOLS_PLUGINS environment variable and man page for details\n");
    fprintf(stderr, "   -v, --verbose               print debugging information on plugin failure\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_plugin(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->nplugin_paths = -1;
    int regions_is_file = 0, targets_is_file = 0, plist_only = 0;

    if ( argc==1 ) usage(args);
    char *plugin_name = NULL;
    if ( argv[1][0]!='-' ) { plugin_name = argv[1]; argc--; argv++; }

    static struct option loptions[] =
    {
        {"verbose",0,0,'v'},
        {"help",0,0,'h'},
        {"list-plugins",0,0,'l'},
        {"output",1,0,'o'},
        {"output-type",1,0,'O'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h?o:O:r:R:li:e:v",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': args->verbose = 1; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'l': plist_only = 1; break;
            case '?':
            case 'h': load_plugin(args, plugin_name, 1, &args->plugin); fprintf(stderr,"%s",args->plugin.usage()); return 0; break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( plist_only )  return list_plugins(args);

    char *fname = NULL;
    if ( optind>=argc || argv[optind][0]=='-' )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
        args->plugin.argc = argc - optind + 1;
        args->plugin.argv = argv + optind - 1;
    }
    else
    {
        fname = argv[optind];
        args->plugin.argc = argc - optind;
        args->plugin.argv = argv + optind;
    }
    optind = 0;
    args->plugin.argv[0] = plugin_name;
    load_plugin(args, plugin_name, 1, &args->plugin);

    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
        args->files->collapse |= COLLAPSE_SOME;
    }
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open or the file not indexed: %s\n", fname);

    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }
        line = args->plugin.process(line);
        if ( line ) bcf_write1(args->out_fh, args->hdr_out, line);
    }
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

