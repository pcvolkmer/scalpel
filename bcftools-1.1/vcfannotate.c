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

#define REPLACE_MISSING  0  // replace only missing values
#define REPLACE_ALL      1  // replace both missing and existing values
#define REPLACE_EXISTING 2  // replace only if tgt is not missing
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
    bcf_sr_regions_t *tgts;

    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE

    rm_tag_t *rm;           // tags scheduled for removal
    int nrm;

    vcmp_t *vcmp;           // for matching annotation and VCF lines by allele
    annot_line_t *alines;   // buffered annotation lines
    int nalines, malines;
    int ref_idx, alt_idx, chr_idx, from_idx, to_idx;   // -1 if not present
    annot_col_t *cols;      // column indexes and setters
    int ncols;

    int *sample_map, nsample_map, sample_is_file;   // map[idst] -> isrc
    int ntmpi, mtmpi, ntmpf, mtmpf, ntmps, mtmps;
    int mtmpi2, mtmpf2, mtmps2;
    int mtmpi3, mtmpf3, mtmps3;
    int32_t *tmpi, *tmpi2, *tmpi3;
    float *tmpf, *tmpf2, *tmpf3;
    char *tmps, *tmps2, **tmpp, **tmpp2;

    char **argv, *output_fname, *targets_fname, *regions_list, *header_fname;
    char *remove_annots, *columns, *rename_chrs, *sample_names;
    int argc, drop_header, tgts_is_vcf;
}
args_t;

char *msprintf(const char *fmt, ...);

void remove_id(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_id(args->hdr,line,NULL);
}
void remove_filter(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    if ( !tag->key ) bcf_update_filter(args->hdr, line, NULL, 0);
    else bcf_remove_filter(args->hdr, line, tag->hdr_id, 1);
}
void remove_qual(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_float_set_missing(line->qual);
}
void remove_info(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    // remove all INFO fields
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    int i;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *inf = &line->d.info[i];
        if ( inf->vptr_free )
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        line->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = NULL;
    }
}
void remove_info_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_info(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_format(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    // remove all FORMAT fields except GT
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    int i;
    for (i=0; i<line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];
        const char *key = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id);
        if ( key[0]=='G' && key[1]=='T' && !key[2] ) continue;

        if ( fmt->p_free )
        {
            free(fmt->p - fmt->p_off);
            fmt->p_free = 0;
        }
        line->d.indiv_dirty = 1;
        fmt->p = NULL;
    }
}

static void remove_hdr_lines(bcf_hdr_t *hdr, int type)
{
    int i = 0, nrm = 0;
    while ( i<hdr->nhrec )
    {
        if ( hdr->hrec[i]->type!=type ) { i++; continue; }
        bcf_hrec_t *hrec = hdr->hrec[i];
        if ( type==BCF_HL_FMT )
        {
            // everything except FORMAT/GT
            int id = bcf_hrec_find_key(hrec, "ID");
            if ( id>=0 && !strcmp(hrec->vals[id],"GT") ) { i++; continue; }
        }
        nrm++;
        hdr->nhrec--;
        if ( i < hdr->nhrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
        bcf_hrec_destroy(hrec);
    }
    if ( nrm ) bcf_hdr_sync(hdr);
}

static void init_remove_annots(args_t *args)
{
    kstring_t str = {0,0,0};
    char *ss = args->remove_annots;
    while ( *ss )
    {
        args->nrm++;
        args->rm = (rm_tag_t*) realloc(args->rm,sizeof(rm_tag_t)*args->nrm);
        rm_tag_t *tag = &args->rm[args->nrm-1];
        tag->key = NULL;

        int type = BCF_HL_GEN;
        if ( !strncasecmp("INFO/",ss,5) ) { type = BCF_HL_INFO; ss += 5; }
        else if ( !strncasecmp("INF/",ss,4) ) { type = BCF_HL_INFO; ss += 4; }
        else if ( !strncasecmp("FORMAT/",ss,7) ) { type = BCF_HL_FMT; ss += 7; }
        else if ( !strncasecmp("FMT/",ss,4) ) { type = BCF_HL_FMT; ss += 4; }

        char *se = ss;
        while ( *se && *se!=',' ) se++;
        str.l = 0;
        kputsn(ss, se-ss, &str);

        if ( type!= BCF_HL_GEN )
        {
            int id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,str.s);
            if ( !bcf_hdr_idinfo_exists(args->hdr,type,id) )
            {
                fprintf(stderr,"Warning: The tag \"%s\" not defined in the header\n", str.s);
                args->nrm--;
            }
            else
            {
                tag->key = strdup(str.s);
                if ( type==BCF_HL_INFO ) tag->handler = remove_info_tag;
                else if ( type==BCF_HL_FMT ) tag->handler = remove_format_tag;
                bcf_hdr_remove(args->hdr_out,type,tag->key);
            }
        }
        else if ( !strcasecmp("ID",str.s) ) tag->handler = remove_id;
        else if ( !strncasecmp("FILTER/",str.s,7) )
        {
            tag->handler = remove_filter;
            tag->key = strdup(str.s+7);
            tag->hdr_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, tag->key);
            if ( !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FLT,tag->hdr_id) ) error("Cannot remove %s, not defined in the header.\n", str.s);
            bcf_hdr_remove(args->hdr_out,BCF_HL_FLT,tag->key);
        }
        else if ( !strcasecmp("FILTER",str.s) )
        {
            tag->handler = remove_filter;
            remove_hdr_lines(args->hdr_out,BCF_HL_FLT);
        }
        else if ( !strcasecmp("QUAL",str.s) ) tag->handler = remove_qual;
        else if ( !strcasecmp("INFO",str.s) ) 
        {
            tag->handler = remove_info;
            remove_hdr_lines(args->hdr_out,BCF_HL_INFO);
        }
        else if ( !strcasecmp("FMT",str.s) || !strcasecmp("FORMAT",str.s) )
        {
            tag->handler = remove_format;
            remove_hdr_lines(args->hdr_out,BCF_HL_FMT);
        }
        else if ( str.l )
        {
            if ( str.s[0]=='#' && str.s[1]=='#' )
                bcf_hdr_remove(args->hdr_out,BCF_HL_GEN,str.s+2);
            else
                bcf_hdr_remove(args->hdr_out,BCF_HL_STR,str.s);
            args->nrm--;
        }

        ss = *se ? se+1 : se;
    }
    free(str.s);
    if ( !args->nrm ) error("No matching tag in -x %s\n", args->remove_annots);
}
static void init_header_lines(args_t *args)
{
    htsFile *file = hts_open(args->header_fname, "rb");
    if ( !file ) error("Error reading %s\n", args->header_fname);
    kstring_t str = {0,0,0};
    while ( hts_getline(file, KS_SEP_LINE, &str) > 0 )
    {
        if ( bcf_hdr_append(args->hdr_out,str.s) ) error("Could not parse %s: %s\n", args->header_fname, str.s);
    }
    hts_close(file);
    free(str.s);
}
static int setter_filter(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    hts_expand(int,1,args->mtmpi,args->tmpi);
    args->tmpi[0] = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, tab->cols[col->icol]);
    if ( args->tmpi[0]<0 ) error("The FILTER is not defined in the header: %s\n", tab->cols[col->icol]);
    if ( col->replace!=REPLACE_MISSING )
        bcf_update_filter(args->hdr_out,line,args->tmpi,1);
    else
        bcf_add_filter(args->hdr_out,line,args->tmpi[0]);
    return 0;
}
static int vcf_setter_filter(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);
    hts_expand(int,rec->d.n_flt,args->mtmpi,args->tmpi);
    int i;
    if ( col->replace!=REPLACE_MISSING )
    {
        for (i=0; i<rec->d.n_flt; i++)
        {
            const char *flt = bcf_hdr_int2id(args->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
            args->tmpi[i] = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, flt);
        }
        bcf_update_filter(args->hdr_out,line,args->tmpi,rec->d.n_flt);
        return 0;
    }
    else
    {
        for (i=0; i<rec->d.n_flt; i++)
        {
            const char *flt = bcf_hdr_int2id(args->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
            bcf_add_filter(args->hdr_out,line,bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, flt));
        }
        return 0;
    }
}
static int setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    // possible cases:
    //      IN  ANNOT   OUT     ACHIEVED_BY
    //      x   y       x        -c +ID
    //      x   y       y        -c ID
    //      x   y       x,y      /not supported/
    //      x   .       x        -c +ID, ID
    //      x   .       .        -x ID
    //      .   y       y        -c +ID, -c ID
    //
    annot_line_t *tab = (annot_line_t*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(args->hdr_out,line,tab->cols[col->icol]);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[0]) )
        return bcf_update_id(args->hdr_out,line,tab->cols[col->icol]);
    return 0;
}
static int vcf_setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( rec->d.id && rec->d.id[0]=='.' && !rec->d.id[1] ) return 0;    // don't replace with "."
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(args->hdr_out,line,rec->d.id);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[0]) )
        return bcf_update_id(args->hdr_out,line,rec->d.id);
    return 0;
}
static int setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;   // empty

    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;

    line->qual = strtod(str, &str);
    if ( str == tab->cols[col->icol] )
        error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key,bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
    return 0;
}
static int vcf_setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( bcf_float_is_missing(rec->qual) ) return 0;
    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;
    line->qual = rec->qual;
    return 0;
}
static int setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( str[0]=='1' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key,NULL,1);
    if ( str[0]=='0' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key,NULL,0);
    error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
    return -1;
}
static int vcf_setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int flag = bcf_get_info_flag(args->files->readers[1].header,rec,col->hdr_key,NULL,NULL);
    bcf_update_info_flag(args->hdr_out,line,col->hdr_key,NULL,flag);
    return 0;
}
static int setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(args->hdr, line, col->hdr_key, &args->tmpi, &args->mtmpi);
        if ( ret>0 && args->tmpi[0]!=bcf_int32_missing ) return 0;
    }

    args->ntmpi = 0;
    while ( *end )
    {
        int val = strtol(str, &end, 10);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
        args->ntmpi++;
        hts_expand(int,args->ntmpi,args->mtmpi,args->tmpi);
        args->tmpi[args->ntmpi-1] = val;
        str = end+1;
    }
    return bcf_update_info_int32(args->hdr_out,line,col->hdr_key,args->tmpi,args->ntmpi);
}
static int vcf_setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(args->hdr, line, col->hdr_key, &args->tmpi, &args->mtmpi);
        if ( ret>0 && args->tmpi[0]!=bcf_int32_missing ) return 0;
    }

    bcf1_t *rec = (bcf1_t*) data;
    args->ntmpi = bcf_get_info_int32(args->files->readers[1].header,rec,col->hdr_key,&args->tmpi,&args->mtmpi);
    if ( args->ntmpi >=0 )
        bcf_update_info_int32(args->hdr_out,line,col->hdr_key,args->tmpi,args->ntmpi);
    return 0;
}
static int setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(args->hdr, line, col->hdr_key, &args->tmpf, &args->mtmpf);
        if ( ret>0 && !bcf_float_is_missing(args->tmpf[0]) ) return 0;
    }

    args->ntmpf = 0;
    while ( *end )
    {
        double val = strtod(str, &end);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
        args->ntmpf++;
        hts_expand(float,args->ntmpf,args->mtmpf,args->tmpf);
        args->tmpf[args->ntmpf-1] = val;
        str = end+1;
    }
    return bcf_update_info_float(args->hdr_out,line,col->hdr_key,args->tmpf,args->ntmpf);
}
static int vcf_setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(args->hdr, line, col->hdr_key, &args->tmpf, &args->mtmpf);
        if ( ret>0 && !bcf_float_is_missing(args->tmpf[0]) ) return 0;
    }

    bcf1_t *rec = (bcf1_t*) data;
    args->ntmpf = bcf_get_info_float(args->files->readers[1].header,rec,col->hdr_key,&args->tmpf,&args->mtmpf);
    if ( args->ntmpf >=0 )
        bcf_update_info_float(args->hdr_out,line,col->hdr_key,args->tmpf,args->ntmpf);
    return 0;
}
static int setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(args->hdr,line,col->hdr_key,&args->tmps,&args->mtmps);
        if ( ret>0 && (args->tmps[0]!='.' || args->tmps[1]) ) return 0;
    }
    annot_line_t *tab = (annot_line_t*) data;
    return bcf_update_info_string(args->hdr_out,line,col->hdr_key,tab->cols[col->icol]);
}
static int vcf_setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(args->hdr,line,col->hdr_key,&args->tmps,&args->mtmps);
        if ( ret>0 && (args->tmps[0]!='.' || args->tmps[1]) ) return 0;
    }
    bcf1_t *rec = (bcf1_t*) data;
    args->ntmps = bcf_get_info_string(args->files->readers[1].header,rec,col->hdr_key,&args->tmps,&args->mtmps);
    if ( args->ntmps >=0 )
        bcf_update_info_string(args->hdr_out,line,col->hdr_key,args->tmps);
    return 0;
}
static int vcf_setter_format_gt(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_genotypes(args->files->readers[1].header,rec,&args->tmpi,&args->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi,nsrc);

    int i, j, ndst = bcf_get_genotypes(args->hdr,line,&args->tmpi2,&args->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    nsrc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 )  // field not present in dst file
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi2, args->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *dst = args->tmpi2 + nsrc*i;
            if ( args->sample_map[i]==-1 )
            {
                dst[0] = bcf_gt_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi2,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            int32_t *src = args->tmpi  + nsrc*args->sample_map[i];
            int32_t *dst = args->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && dst[0]==bcf_gt_missing ) continue;
            if ( col->replace==REPLACE_MISSING  && dst[0]!=bcf_gt_missing ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi3, args->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *ori = args->tmpi2 + ndst*i;
            int32_t *dst = args->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( args->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && ori[0]==bcf_gt_missing ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && ori[0]!=bcf_gt_missing ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi3,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int vcf_setter_format_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_int32(args->files->readers[1].header,rec,col->hdr_key,&args->tmpi,&args->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi,nsrc);

    int i, j, ndst = bcf_get_format_int32(args->hdr,line,col->hdr_key,&args->tmpi2,&args->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    nsrc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 )
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi2, args->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *dst = args->tmpi2 + nsrc*i;
            if ( args->sample_map[i]==-1 )
            {
                dst[0] = bcf_int32_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi2,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            int32_t *src = args->tmpi  + nsrc*args->sample_map[i];
            int32_t *dst = args->tmpi2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && dst[0]==bcf_int32_missing ) continue;
            if ( col->replace==REPLACE_MISSING  && dst[0]!=bcf_int32_missing ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi3, args->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *ori = args->tmpi2 + ndst*i;
            int32_t *dst = args->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( args->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && ori[0]==bcf_int32_missing ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && ori[0]!=bcf_int32_missing ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key,args->tmpi3,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int vcf_setter_format_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_float(args->files->readers[1].header,rec,col->hdr_key,&args->tmpf,&args->mtmpf);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf,nsrc);

    int i, j, ndst = bcf_get_format_float(args->hdr,line,col->hdr_key,&args->tmpf2,&args->mtmpf2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    nsrc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 )
    {
        if ( col->replace==REPLACE_EXISTING ) return 0;    // overwrite only if present
        hts_expand(float, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpf2, args->tmpf2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            float *dst = args->tmpf2 + nsrc*i;
            if ( args->sample_map[i]==-1 )
            {
                bcf_float_set_missing(dst[0]);
                for (j=1; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = args->tmpf + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf2,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nsrc )     
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            float *src = args->tmpf  + nsrc*args->sample_map[i];
            float *dst = args->tmpf2 + ndst*i;
            if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) bcf_float_set_vector_end(dst[j]);
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(float, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpf3, args->tmpf3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            float *ori = args->tmpf2 + ndst*i;
            float *dst = args->tmpf3 + nsrc*i;
            int keep_ori = 0;
            if ( args->sample_map[i]==-1 ) keep_ori = 1;
            else if ( col->replace==REPLACE_EXISTING && bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_float_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = args->tmpf + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key,args->tmpf3,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int vcf_setter_format_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    args->tmpp[0] = args->tmps;
    int ret = bcf_get_format_string(args->files->readers[1].header,rec,col->hdr_key,&args->tmpp,&args->mtmps);
    args->tmps = args->tmpp[0]; // tmps might be realloced
    if ( ret==-3 ) return 0;    // the tag is not present
    if ( ret<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_format_string(args->hdr_out,line,col->hdr_key,(const char**)args->tmpp,bcf_hdr_nsamples(args->hdr_out));

    int i;
    args->tmpp2[0] = args->tmps2;
    ret = bcf_get_format_string(args->hdr,line,col->hdr_key,&args->tmpp2,&args->mtmps2);
    args->tmps2 = args->tmpp2[0];   // tmps2 might be realloced

    if ( ret<=0 )   // not present in dst
    {
        hts_expand(char,bcf_hdr_nsamples(args->hdr_out)*2,args->mtmps2,args->tmps2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            args->tmps2[2*i]   = '.';
            args->tmps2[2*i+1] = 0;
            args->tmpp2[i] = args->tmps2+2*i;
        }
    }

    for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
    {
        int isrc = args->sample_map[i];
        if ( isrc==-1 ) continue;
        args->tmpp2[i] = args->tmpp[isrc];
    }
    return bcf_update_format_string(args->hdr_out,line,col->hdr_key,(const char**)args->tmpp2,bcf_hdr_nsamples(args->hdr_out));
}
static void set_samples(args_t *args, bcf_hdr_t *src, bcf_hdr_t *dst, int need_samples)
{
    int i;
    if ( !args->sample_names )
    {
        int nmatch = 0, order_ok = 1;
        for (i=0; i<bcf_hdr_nsamples(src); i++)
        {
            int id = bcf_hdr_id2int(dst, BCF_DT_SAMPLE, src->samples[i]);
            if ( id!=-1 ) 
            {
                nmatch++;
                if ( i!=id ) order_ok = 0;
            }
        }
        if ( bcf_hdr_nsamples(src)==bcf_hdr_nsamples(dst) && nmatch==bcf_hdr_nsamples(src) && order_ok && !need_samples ) 
            return;    // the same samples in both files

        if ( !nmatch ) error("No matching samples found in the source and the destination file\n");
        if ( nmatch!=bcf_hdr_nsamples(src) || nmatch!=bcf_hdr_nsamples(dst) ) fprintf(stderr,"%d sample(s) in common\n", nmatch);

        args->nsample_map = bcf_hdr_nsamples(dst);
        args->sample_map  = (int*) malloc(sizeof(int)*args->nsample_map);
        for (i=0; i<args->nsample_map; i++)
        {
            int id = bcf_hdr_id2int(src, BCF_DT_SAMPLE, dst->samples[i]);
            args->sample_map[i] = id;   // idst -> isrc, -1 if not present
        }
        return;
    }

    args->nsample_map = bcf_hdr_nsamples(dst);
    args->sample_map  = (int*) malloc(sizeof(int)*args->nsample_map);
    for (i=0; i<args->nsample_map; i++) args->sample_map[i] = -1;

    int nsamples = 0;
    char **samples = hts_readlist(args->sample_names, args->sample_is_file, &nsamples);
    for (i=0; i<nsamples; i++)
    {
        int isrc, idst;
        char *ss = samples[i], *se = samples[i];
        while ( *se && !isspace(*se) ) se++;
        if ( !*se ) 
        {
            // only one sample name
            isrc = bcf_hdr_id2int(src, BCF_DT_SAMPLE,ss);
            if ( isrc==-1 ) error("Sample \"%s\" not found in the source file\n", ss);
            idst = bcf_hdr_id2int(dst, BCF_DT_SAMPLE,ss);
            if ( idst==-1 ) error("Sample \"%s\" not found in the destination file\n", ss);
            args->sample_map[idst] = isrc;
            continue;
        }
        *se = 0;
        isrc = bcf_hdr_id2int(src, BCF_DT_SAMPLE,ss);
        if ( isrc==-1 ) error("Sample \"%s\" not found in the source file\n", ss);

        ss = se+1;
        while ( isspace(*ss) ) ss++;
        se = ss;
        while ( *se && !isspace(*se) ) se++;

        idst = bcf_hdr_id2int(dst, BCF_DT_SAMPLE,ss);
        if ( idst==-1 ) error("Sample \"%s\" not found in the destination file\n", ss);

        args->sample_map[idst] = isrc;
    }
    for (i=0; i<nsamples; i++) free(samples[i]);
    free(samples);
}
static char *columns_complement(char *columns, void **skip_info, void **skip_fmt)
{
    kstring_t str = {0,0,0};
    char *ss = columns, *se = ss;
    while ( *ss )
    {
        if ( *se && *se!=',' ) { se++; continue; }
        if ( *ss!='^' )
        {
            if ( str.l ) kputc(',',&str);
            kputsn(ss, se-ss, &str);
            if ( !*se ) break;
            ss = ++se;
            continue;
        }

        if ( !strncasecmp("^INFO/",ss,6) )
        {
            if ( !*skip_info )
            {
                *skip_info = khash_str2int_init();
                if ( str.l ) kputc(',',&str);
                kputs("INFO",&str);
            }
            char tmp = *se; *se = 0;
            khash_str2int_inc(*skip_info, strdup(ss+6));
            *se = tmp;
        }
        else if ( !strncasecmp("^FORMAT/",ss,8) || !strncasecmp("^FMT/",ss,5) )
        {
            int n = !strncasecmp("^FMT/",ss,5) ? 5 : 8;
            if ( !*skip_fmt )
            {
                *skip_fmt = khash_str2int_init();
                if ( str.l ) kputc(',',&str);
                kputs("FORMAT",&str);
            }
            char tmp = *se; *se = 0;
            khash_str2int_inc(*skip_fmt, strdup(ss+n));
            *se = tmp;
        }
        else
        {
            if ( !*skip_info )
            {
                *skip_info = khash_str2int_init();
                if ( str.l ) kputc(',',&str);
                kputs("INFO",&str);
            }
            char tmp = *se; *se = 0;
            khash_str2int_inc(*skip_info, strdup(ss+1));
            *se = tmp;
        }

        if ( !*se ) break;
        ss = ++se;
    }
    free(columns);
    return str.s;
}
static void init_columns(args_t *args)
{
    void *skip_fmt = NULL, *skip_info = NULL;
    if ( args->tgts_is_vcf )
        args->columns = columns_complement(args->columns, &skip_info, &skip_fmt);

    kstring_t str = {0,0,0}, tmp = {0,0,0};
    char *ss = args->columns, *se = ss;
    args->ncols = 0;
    int i = -1, has_fmt_str = 0, force_samples = -1;
    while ( *ss )
    {
        if ( *se && *se!=',' ) { se++; continue; }
        int replace = REPLACE_ALL;
        if ( *ss=='+' ) { replace = REPLACE_MISSING; ss++; }
        else if ( *ss=='-' ) { replace = REPLACE_EXISTING; ss++; }
        i++;
        str.l = 0;
        kputsn(ss, se-ss, &str);
        if ( !strcasecmp("-",str.s) ) ;
        else if ( !strcasecmp("CHROM",str.s) ) args->chr_idx = i;
        else if ( !strcasecmp("POS",str.s) ) args->from_idx = i;
        else if ( !strcasecmp("FROM",str.s) ) args->from_idx = i;
        else if ( !strcasecmp("TO",str.s) ) args->to_idx = i;
        else if ( !strcasecmp("REF",str.s) ) args->ref_idx = i;
        else if ( !strcasecmp("ALT",str.s) ) args->alt_idx = i;
        else if ( !strcasecmp("ID",str.s) )
        {
            if ( replace==REPLACE_EXISTING ) error("todo: -ID\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_id : setter_id;
            col->hdr_key = strdup(str.s);
        }
        else if ( !strcasecmp("FILTER",str.s) )
        {
            if ( replace==REPLACE_EXISTING ) error("todo: -FILTER\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_filter : setter_filter;
            col->hdr_key = strdup(str.s);
            if ( args->tgts_is_vcf )
            {
                bcf_hdr_t *tgts_hdr = args->files->readers[1].header;
                int j;
                for (j=0; j<tgts_hdr->nhrec; j++)
                {
                    bcf_hrec_t *hrec = tgts_hdr->hrec[j];
                    if ( hrec->type!=BCF_HL_FLT ) continue;
                    int k = bcf_hrec_find_key(hrec,"ID");
                    assert( k>=0 ); // this should always be true for valid VCFs
                    tmp.l = 0;
                    bcf_hrec_format(hrec, &tmp);
                    bcf_hdr_append(args->hdr_out, tmp.s);
                }
            }
        }
        else if ( !strcasecmp("QUAL",str.s) )
        {
            if ( replace==REPLACE_EXISTING ) error("todo: -QUAL\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_qual : setter_qual;
            col->hdr_key = strdup(str.s);
        }
        else if ( args->tgts_is_vcf && !strcasecmp("INFO",str.s) ) // All INFO fields
        {
            if ( replace==REPLACE_EXISTING ) error("todo: -INFO/TAG\n");
            bcf_hdr_t *tgts_hdr = args->files->readers[1].header;
            int j;
            for (j=0; j<tgts_hdr->nhrec; j++)
            {
                bcf_hrec_t *hrec = tgts_hdr->hrec[j];
                if ( hrec->type!=BCF_HL_INFO ) continue;
                int k = bcf_hrec_find_key(hrec,"ID");
                assert( k>=0 ); // this should always be true for valid VCFs
                if ( skip_info && khash_str2int_has_key(skip_info,hrec->vals[k]) ) continue;
                tmp.l = 0;
                bcf_hrec_format(hrec, &tmp);
                bcf_hdr_append(args->hdr_out, tmp.s);
                int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, hrec->vals[k]);
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                col->icol = -1;
                col->replace = replace;
                col->hdr_key = strdup(hrec->vals[k]);
                switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id) )
                {
                    case BCF_HT_FLAG:   col->setter = vcf_setter_info_flag; break;
                    case BCF_HT_INT:    col->setter = vcf_setter_info_int; break;
                    case BCF_HT_REAL:   col->setter = vcf_setter_info_real; break;
                    case BCF_HT_STR:    col->setter = vcf_setter_info_str; break;
                    default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id));
                }
            }
        }
        else if ( args->tgts_is_vcf && (!strcasecmp("FORMAT",str.s) || !strcasecmp("FMT",str.s)) ) // All FORMAT fields
        {
            bcf_hdr_t *tgts_hdr = args->files->readers[1].header;
            if ( force_samples<0 ) force_samples = replace;
            if ( force_samples>=0 && replace!=REPLACE_ALL ) force_samples = replace;
            int j;
            for (j=0; j<tgts_hdr->nhrec; j++)
            {
                bcf_hrec_t *hrec = tgts_hdr->hrec[j];
                if ( hrec->type!=BCF_HL_FMT) continue;
                int k = bcf_hrec_find_key(hrec,"ID");
                assert( k>=0 ); // this should always be true for valid VCFs
                if ( skip_fmt && khash_str2int_has_key(skip_fmt,hrec->vals[k]) ) continue;
                tmp.l = 0;
                bcf_hrec_format(hrec, &tmp);
                bcf_hdr_append(args->hdr_out, tmp.s);
                int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, hrec->vals[k]);
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                col->icol = -1;
                col->replace = replace;
                col->hdr_key = strdup(hrec->vals[k]);
                if ( !strcasecmp("GT",col->hdr_key) ) col->setter = vcf_setter_format_gt;
                else
                    switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id) )
                    {
                        case BCF_HT_INT:    col->setter = vcf_setter_format_int; break;
                        case BCF_HT_REAL:   col->setter = vcf_setter_format_real; break;
                        case BCF_HT_STR:    col->setter = vcf_setter_format_str; has_fmt_str = 1; break;
                        default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id));
                    }
            }
        }
        else if ( args->tgts_is_vcf && (!strncasecmp("FORMAT/",str.s, 7) || !strncasecmp("FMT/",str.s,4)) )
        {
            char *key = str.s + (!strncasecmp("FMT/",str.s,4) ? 4 : 7);
            if ( force_samples<0 ) force_samples = replace;
            if ( force_samples>=0 && replace!=REPLACE_ALL ) force_samples = replace;;
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key);
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = -1;
            col->replace = replace;
            col->hdr_key = strdup(key);
            if ( !strcasecmp("GT",key) ) col->setter = vcf_setter_format_gt;
            else
                switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id) )
                {
                    case BCF_HT_INT:    col->setter = vcf_setter_format_int; break;
                    case BCF_HT_REAL:   col->setter = vcf_setter_format_real; break;
                    case BCF_HT_STR:    col->setter = vcf_setter_format_str; has_fmt_str = 1; break;
                    default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id));
                }
        }
        else
        {
            if ( replace==REPLACE_EXISTING ) error("todo: -INFO/TAG\n");
            if ( !strncasecmp("INFO/",str.s,5) ) { memmove(str.s,str.s+5,str.l-4); }
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, str.s);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                if ( args->tgts_is_vcf ) // reading annotations from a VCF, add a new header line
                {
                    bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_INFO, "ID", str.s, NULL);
                    if ( !hrec ) error("The tag \"%s\" is not defined in %s\n", str.s,args->files->readers[1].fname);
                    tmp.l = 0;
                    bcf_hrec_format(hrec, &tmp);
                    bcf_hdr_append(args->hdr_out, tmp.s);
                    hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, str.s);
                }
                else
                    error("The tag \"%s\" is not defined in %s\n", str.s, args->targets_fname);
                assert( bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) );
            }

            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = i;
            col->replace = replace;
            col->hdr_key = strdup(str.s);
            switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                case BCF_HT_FLAG:   col->setter = args->tgts_is_vcf ? vcf_setter_info_flag : setter_info_flag; break;
                case BCF_HT_INT:    col->setter = args->tgts_is_vcf ? vcf_setter_info_int  : setter_info_int; break;
                case BCF_HT_REAL:   col->setter = args->tgts_is_vcf ? vcf_setter_info_real : setter_info_real; break;
                case BCF_HT_STR:    col->setter = args->tgts_is_vcf ? vcf_setter_info_str  : setter_info_str; break;
                default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id));
            }
        }
        if ( !*se ) break;
        ss = ++se;
    }
    free(str.s);
    free(tmp.s);
    if ( args->to_idx==-1 ) args->to_idx = args->from_idx;
    free(args->columns);
    if ( skip_info ) khash_str2int_destroy_free(skip_info);
    if ( skip_fmt ) khash_str2int_destroy_free(skip_fmt);
    if ( has_fmt_str )
    {
        int n = bcf_hdr_nsamples(args->hdr_out) > bcf_hdr_nsamples(args->files->readers[1].header) ? bcf_hdr_nsamples(args->hdr_out) : bcf_hdr_nsamples(args->files->readers[1].header);
        args->tmpp  = (char**)malloc(sizeof(char*)*n);
        args->tmpp2 = (char**)malloc(sizeof(char*)*n);
    }
    if ( force_samples>=0 )
        set_samples(args, args->files->readers[1].header, args->hdr, force_samples==REPLACE_ALL ? 0 : 1);
}

static void rename_chrs(args_t *args, char *fname)
{
    int n, i;
    char **map = hts_readlist(fname, 1, &n);
    if ( !map ) error("Could not read: %s\n", fname);
    for (i=0; i<n; i++)
    {
        char *ss = map[i];
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) error("Could not parse: %s\n", fname);
        *ss = 0;
        int rid = bcf_hdr_name2id(args->hdr_out, map[i]);
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->hdr_out, BCF_HL_CTG, "ID", map[i], NULL);
        if ( !hrec ) continue;  // the sequence not present
        int j = bcf_hrec_find_key(hrec, "ID");
        assert( j>=0 );
        free(hrec->vals[j]);
        ss++;
        while ( *ss && isspace(*ss) ) ss++;
        char *se = ss;
        while ( *se && !isspace(*se) ) se++;
        *se = 0;
        hrec->vals[j] = strdup(ss);
        args->hdr_out->id[BCF_DT_CTG][rid].key = hrec->vals[j];
    }
    for (i=0; i<n; i++) free(map[i]);
    free(map);
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    if ( args->remove_annots ) init_remove_annots(args);
    if ( args->header_fname ) init_header_lines(args);
    if ( args->targets_fname && args->tgts_is_vcf )
    {
        // reading annots from a VCF
        if ( !bcf_sr_add_reader(args->files, args->targets_fname) )
            error("Failed to open or the file not indexed: %s\n", args->targets_fname);
    }
    if ( args->columns ) init_columns(args);
    if ( args->targets_fname && !args->tgts_is_vcf )
    {
        if ( !args->columns ) error("The -c option not given\n");
        if ( args->chr_idx==-1 ) error("The -c CHROM option not given\n");
        if ( args->from_idx==-1 ) error("The -c POS option not given\n");
        if ( args->to_idx==-1 ) args->to_idx = -args->from_idx - 1;

        args->tgts = bcf_sr_regions_init(args->targets_fname,1,args->chr_idx,args->from_idx,args->to_idx);
        if ( !args->tgts ) error("Could not initialize the annotation file: %s\n", args->targets_fname);
        if ( !args->tgts->tbx ) error("Expected tabix-indexed annotation file: %s\n", args->targets_fname);
        args->vcmp = vcmp_init();
    }

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_annotate");
    if ( !args->drop_header )
    {
        if ( args->rename_chrs ) rename_chrs(args, args->rename_chrs);

        args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
        if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
        bcf_hdr_write(args->out_fh, args->hdr_out);
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nrm; i++) free(args->rm[i].key);
    free(args->rm);
    if ( args->hdr_out ) bcf_hdr_destroy(args->hdr_out);
    if (args->vcmp) vcmp_destroy(args->vcmp);
    for (i=0; i<args->ncols; i++)
        free(args->cols[i].hdr_key);
    free(args->cols);
    for (i=0; i<args->malines; i++)
    {
        free(args->alines[i].cols);
        free(args->alines[i].als);
        free(args->alines[i].line.s);
    }
    free(args->alines);
    if ( args->tgts ) bcf_sr_regions_destroy(args->tgts);
    free(args->tmpi);
    free(args->tmpf);
    free(args->tmps);
    free(args->tmpp);
    free(args->tmpi2);
    free(args->tmpf2);
    free(args->tmps2);
    free(args->tmpp2);
    free(args->tmpi3);
    free(args->tmpf3);
    if ( args->filter )
        filter_destroy(args->filter);
    if (args->out_fh) hts_close(args->out_fh);
    free(args->sample_map);
}

static void buffer_annot_lines(args_t *args, bcf1_t *line, int start_pos, int end_pos)
{
    if ( args->nalines && args->alines[0].rid != line->rid ) args->nalines = 0;

    int i = 0;
    while ( i<args->nalines )
    {
        if ( line->pos > args->alines[i].end )
        {
            args->nalines--;
            if ( args->nalines && i<args->nalines )
            {
                annot_line_t tmp = args->alines[i];
                memmove(&args->alines[i],&args->alines[i+1],(args->nalines-i)*sizeof(annot_line_t));
                args->alines[args->nalines] = tmp;
            }
        }
        else i++;
    }

    if ( args->ref_idx==-1 && args->nalines ) return;

    while ( !bcf_sr_regions_overlap(args->tgts, bcf_seqname(args->hdr,line), start_pos,end_pos) )
    {
        args->nalines++;
        hts_expand0(annot_line_t,args->nalines,args->malines,args->alines);
        annot_line_t *tmp = &args->alines[args->nalines-1];
        tmp->rid   = line->rid;
        tmp->start = args->tgts->start;
        tmp->end   = args->tgts->end;
        tmp->line.l = 0;
        kputs(args->tgts->line.s, &tmp->line);
        char *s = tmp->line.s;
        tmp->ncols = 1;
        hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols);
        tmp->cols[0] = s;
        while ( *s )
        {
            if ( *s=='\t' )
            {
                tmp->ncols++;
                hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols);
                tmp->cols[tmp->ncols-1] = s+1;
                *s = 0;
            }
            s++;
        }
        if ( args->ref_idx != -1 )
        {
            assert( args->ref_idx < tmp->ncols );
            assert( args->alt_idx < tmp->ncols );
            s = tmp->cols[args->alt_idx];
            tmp->nals = 1;
            hts_expand(char*,tmp->nals,tmp->mals,tmp->als);
            tmp->als[0] = s;
            while ( *s )
            {
                if ( *s==',' )
                {
                    tmp->nals++;
                    hts_expand(char*,tmp->nals,tmp->mals,tmp->als);
                    tmp->als[tmp->nals-1] = s+1;
                    *s = 0;
                }
                s++;
            }
            int iseq = args->tgts->iseq;
            if ( bcf_sr_regions_next(args->tgts)<0 || args->tgts->iseq!=iseq ) break;
        }
        else break;
    }
}

static void annotate(args_t *args, bcf1_t *line)
{
    int i, j;
    for (i=0; i<args->nrm; i++)
        args->rm[i].handler(args, line, &args->rm[i]);

    if ( args->tgts )
    {
        // Buffer annotation lines. When multiple ALT alleles are present in the
        // annotation file, at least one must match one of the VCF alleles.
        int len = 0;
        bcf_get_variant_types(line);
        for (i=1; i<line->n_allele; i++)
            if ( len > line->d.var[i].n ) len = line->d.var[i].n;
        int end_pos = len<0 ? line->pos - len : line->pos;
        buffer_annot_lines(args, line, line->pos, end_pos);
        for (i=0; i<args->nalines; i++)
        {
            if ( line->pos > args->alines[i].end || end_pos < args->alines[i].start ) continue;
            if ( args->ref_idx != -1 )
            {
                if ( vcmp_set_ref(args->vcmp, line->d.allele[0], args->alines[i].cols[args->ref_idx]) < 0 ) continue;   // refs not compatible
                for (j=0; j<args->alines[i].nals; j++)
                {
                    if ( line->n_allele==1 && args->alines[i].als[j][0]=='.' && args->alines[i].als[j][1]==0 ) break;   // no ALT allele in VCF and annot file has "."
                    if ( vcmp_find_allele(args->vcmp, line->d.allele+1, line->n_allele - 1, args->alines[i].als[j]) >= 0 ) break;
                }
                if ( j==args->alines[i].nals ) continue;    // none of the annot alleles present in VCF's ALT
            }
            break;
        }

        if ( i<args->nalines )
        {
            for (j=0; j<args->ncols; j++)
                if ( args->cols[j].setter(args,line,&args->cols[j],&args->alines[i]) )
                    error("fixme: Could not set %s at %s:%d\n", args->cols[j].hdr_key,bcf_seqname(args->hdr,line),line->pos+1);
        }
    }
    else if ( args->files->nreaders == 2 && bcf_sr_has_line(args->files,1) )
    {
        bcf1_t *aline = bcf_sr_get_line(args->files,1);
        for (j=0; j<args->ncols; j++)
            if ( args->cols[j].setter(args,line,&args->cols[j],aline) )
                error("fixme: Could not set %s at %s:%d\n", args->cols[j].hdr_key,bcf_seqname(args->hdr,line),line->pos+1);
    }
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Annotate and edit VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools annotate [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a, --annotations <file>       VCF file or tabix-indexed file with annotations: CHR\\tPOS[\\tVALUE]+\n");
    fprintf(stderr, "   -c, --columns <list>           list of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details\n");
    fprintf(stderr, "   -e, --exclude <expr>           exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "   -h, --header-lines <file>      lines which should be appended to the VCF header\n");
    fprintf(stderr, "   -i, --include <expr>           select sites for which the expression is true (see man pagee for details)\n");
    fprintf(stderr, "   -o, --output <file>            write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "   -r, --regions <region>         restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file <file>      restrict to regions listed in a file\n");
    fprintf(stderr, "       --rename-chrs <file>       rename sequences according to map file: from\\tto\n");
    fprintf(stderr, "   -s, --samples [^]<list>        comma separated list of samples to annotate (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "   -S, --samples-file [^]<file>   file of samples to annotate (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "   -x, --remove <list>            list of annotations to remove (e.g. ID,INFO/DP,FORMAT/DP,FILTER). See man page for details\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfannotate(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->ref_idx = args->alt_idx = args->chr_idx = args->from_idx = args->to_idx = -1;
    int regions_is_file = 0;

    static struct option loptions[] =
    {
        {"output",1,0,'o'},
        {"output-type",1,0,'O'},
        {"annotations",1,0,'a'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"remove",1,0,'x'},
        {"columns",1,0,'c'},
        {"rename-chrs",1,0,1},
        {"header-lines",1,0,'h'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "h:?o:O:r:R:a:x:c:i:e:S:s:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case 'c': args->columns = strdup(optarg); break;
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
            case 'x': args->remove_annots = optarg; break;
            case 'a': args->targets_fname = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'h': args->header_fname = optarg; break;
            case  1 : args->rename_chrs = optarg; break;
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
    }
    else fname = argv[optind];

    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_fname && hts_file_type(args->targets_fname) & (FT_VCF|FT_BCF) )
    {
        args->tgts_is_vcf = 1;
        args->files->require_index = 1;
        args->files->collapse |= COLLAPSE_SOME;
    }
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open or the file not indexed: %s\n", fname);

    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( line->errcode ) error("Encountered error, cannot proceed. Please check the error output above.\n");
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }
        annotate(args, line);
        bcf_write1(args->out_fh, args->hdr_out, line);
    }
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
