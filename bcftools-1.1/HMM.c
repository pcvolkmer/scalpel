/* The MIT License

   Copyright (c) 2014 Genome Research Ltd.

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
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <htslib/hts.h>
#include "HMM.h"

static inline void multiply_matrix(int n, double *a, double *b, double *dst, double *tmp)
{
    double *out = dst;
    if ( a==dst || b==dst )
        out = tmp;

    int i,j,k;
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            double val = 0;
            for (k=0; k<n; k++) val += MAT(a,n,i,k)*MAT(b,n,k,j);
            MAT(out,n,i,j) = val;
        }
    }
    if ( out!=dst )
        memcpy(dst,out,sizeof(double)*n*n);
}

hmm_t *hmm_init(int nstates, double *tprob, int ntprob)
{
    hmm_t *hmm = (hmm_t*) calloc(1,sizeof(hmm_t));
    hmm->nstates = nstates;
    hmm->curr_tprob = (double*) malloc(sizeof(double)*nstates*nstates);
    hmm->tmp = (double*) malloc(sizeof(double)*nstates*nstates);

    hmm_set_tprob(hmm, tprob, ntprob);

    return hmm;
}

void hmm_set_tprob(hmm_t *hmm, double *tprob, int ntprob)
{
    hmm->ntprob_arr = ntprob;
    if ( ntprob<=0 ) ntprob = 1;

    if ( !hmm->tprob_arr )
        hmm->tprob_arr  = (double*) malloc(sizeof(double)*hmm->nstates*hmm->nstates*ntprob);

    memcpy(hmm->tprob_arr,tprob,sizeof(double)*hmm->nstates*hmm->nstates);

    int i;
    for (i=1; i<ntprob; i++)
        multiply_matrix(hmm->nstates, hmm->tprob_arr, hmm->tprob_arr+(i-1)*hmm->nstates*hmm->nstates, hmm->tprob_arr+i*hmm->nstates*hmm->nstates, hmm->tmp);
}

void hmm_set_tprob_func(hmm_t *hmm, set_tprob_f set_tprob, void *data)
{
    hmm->set_tprob = set_tprob;
    hmm->set_tprob_data = data;
}

static void _set_tprob(hmm_t *hmm, int pos_diff)
{
    assert( pos_diff>=0 );

    int i, n;

    n = hmm->ntprob_arr ? pos_diff % hmm->ntprob_arr : 0;  // n-th precalculated matrix
    memcpy(hmm->curr_tprob, hmm->tprob_arr+n*hmm->nstates*hmm->nstates, sizeof(*hmm->curr_tprob)*hmm->nstates*hmm->nstates);

    if ( hmm->ntprob_arr > 0  )
    {
        n = pos_diff / hmm->ntprob_arr;  // number of full blocks to jump
        for (i=0; i<n; i++)
            multiply_matrix(hmm->nstates, hmm->tprob_arr+(hmm->ntprob_arr-1)*hmm->nstates*hmm->nstates, hmm->curr_tprob, hmm->curr_tprob, hmm->tmp);
    }
}

void hmm_run_viterbi(hmm_t *hmm, int n, double *eprobs, uint32_t *sites)
{
    // Init arrays when run for the first time
    if ( hmm->nsites < n )
    {
        hmm->nsites = n;
        hmm->vpath = (uint8_t*) realloc(hmm->vpath, sizeof(uint8_t)*hmm->nsites*hmm->nstates);
    }
    if ( !hmm->vprob )
    {
        hmm->vprob     = (double*) malloc(sizeof(double)*hmm->nstates);
        hmm->vprob_tmp = (double*) malloc(sizeof(double)*hmm->nstates);
    }


    // Init all states with equal likelihood
    int i,j, nstates = hmm->nstates;
    for (i=0; i<nstates; i++) hmm->vprob[i] = 1./nstates;

    // Run Viterbi
    uint32_t prev_pos = sites[0];
    for (i=0; i<n; i++)
    {
        uint8_t *vpath = &hmm->vpath[i*nstates];
        double *eprob  = &eprobs[i*nstates];

        int pos_diff = sites[i] == prev_pos ? 0 : sites[i] - prev_pos - 1;

        _set_tprob(hmm, pos_diff);
        if ( hmm->set_tprob ) hmm->set_tprob(hmm, prev_pos, sites[i], hmm->set_tprob_data);
        prev_pos = sites[i];

        double vnorm = 0;
        for (j=0; j<nstates; j++)
        {
            double vmax = 0;
            int k, k_vmax = 0;
            for (k=0; k<nstates; k++)
            {
                double pval = hmm->vprob[k] * MAT(hmm->curr_tprob,hmm->nstates,j,k);
                if ( vmax < pval ) { vmax = pval; k_vmax = k; }
            }
            vpath[j] = k_vmax;
            hmm->vprob_tmp[j] = vmax * eprob[j];
            vnorm += hmm->vprob_tmp[j];
        }
        for (j=0; j<nstates; j++) hmm->vprob_tmp[j] /= vnorm;
        double *tmp = hmm->vprob; hmm->vprob = hmm->vprob_tmp; hmm->vprob_tmp = tmp;
    }

    // Find the most likely state
    int iptr = 0;
    for (i=1; i<nstates; i++) 
        if ( hmm->vprob[iptr] < hmm->vprob[i] ) iptr = i;

    // Trace back the Viterbi path, we are reusing vpath for storing the states (vpath[i*nstates])
    for (i=n-1; i>=0; i--)
    {
        assert( iptr<nstates && hmm->vpath[i*nstates + iptr]<nstates );
        iptr = hmm->vpath[i*nstates + iptr];
        hmm->vpath[i*nstates] = iptr;     // reusing the array for different purpose here
    }
}

void hmm_run_fwd_bwd(hmm_t *hmm, int n, double *eprobs, uint32_t *sites)
{
    // Init arrays when run for the first time
    if ( hmm->nsites < n )
    {
        hmm->nsites = n;
        hmm->fwd = (double*) realloc(hmm->fwd, sizeof(double)*(hmm->nsites+1)*hmm->nstates);
    }
    if ( !hmm->fwd )
    {
        hmm->fwd     = (double*) malloc(sizeof(double)*hmm->nstates*(hmm->nsites+1));
        hmm->bwd     = (double*) malloc(sizeof(double)*hmm->nstates);
        hmm->bwd_tmp = (double*) malloc(sizeof(double)*hmm->nstates);
    }


    // Init all states with equal likelihood
    int i,j,k, nstates = hmm->nstates;
    for (i=0; i<nstates; i++) hmm->fwd[i] = 1./hmm->nstates;
    for (i=0; i<nstates; i++) hmm->bwd[i] = 1./hmm->nstates;

    // Run fwd 
    uint32_t prev_pos = sites[0];
    for (i=0; i<n; i++)
    {
        double *fwd_prev = &hmm->fwd[i*nstates];
        double *fwd      = &hmm->fwd[(i+1)*nstates];
        double *eprob    = &eprobs[i*nstates];

        int pos_diff = sites[i] == prev_pos ? 0 : sites[i] - prev_pos - 1;

        _set_tprob(hmm, pos_diff);
        if ( hmm->set_tprob ) hmm->set_tprob(hmm, prev_pos, sites[i], hmm->set_tprob_data);
        prev_pos = sites[i];

        double norm = 0;
        for (j=0; j<nstates; j++)
        {
            double pval = 0;
            for (k=0; k<nstates; k++)
                pval += fwd_prev[k] * MAT(hmm->curr_tprob,hmm->nstates,j,k);
            fwd[j] = pval * eprob[j];
            norm += fwd[j];
        }
        for (j=0; j<nstates; j++) fwd[j] /= norm;
    }

    // Run bwd
    double *bwd = hmm->bwd, *bwd_tmp = hmm->bwd_tmp;
    prev_pos = sites[n-1];
    for (i=0; i<n; i++)
    {
        double *fwd   = &hmm->fwd[(n-i)*nstates];
        double *eprob = &eprobs[(n-i-1)*nstates];
        
        int pos_diff = sites[n-i-1] == prev_pos ? 0 : prev_pos - sites[n-i-1] - 1;

        _set_tprob(hmm, pos_diff);
        if ( hmm->set_tprob ) hmm->set_tprob(hmm, sites[n-i-1], prev_pos, hmm->set_tprob_data);
        prev_pos = sites[n-i-1];

        double norm = 0;
        for (j=0; j<nstates; j++)
        {
            double pval = 0;
            for (k=0; k<nstates; k++)
                pval += bwd[k] * eprob[k] * MAT(hmm->curr_tprob,hmm->nstates,k,j);
            bwd_tmp[j] = pval;
            norm += pval;
        }
        for (j=0; j<nstates; j++)
        {
            bwd_tmp[j] /= norm;
            fwd[j] *= bwd_tmp[j];   // fwd now stores fwd*bwd
        }
        double *tmp = bwd_tmp; bwd_tmp = bwd; bwd = tmp;
    }
}

void hmm_destroy(hmm_t *hmm)
{
    free(hmm->vprob);
    free(hmm->vprob_tmp);
    free(hmm->vpath);
    free(hmm->curr_tprob);
    free(hmm->tmp);
    free(hmm->tprob_arr);
    free(hmm->fwd);
    free(hmm->bwd);
    free(hmm->bwd_tmp);
    free(hmm);
}

