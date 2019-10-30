#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void center_step(void *, void *, void *, void *);
extern void hl_center_step(void *, void *, void *, void *);
/* extern void hl_loc_step(void *, void *, void *, void *); */
extern void norming(void *, void *, void *);
extern void pairdiffc(void *, void *, void *);
extern void pairsumc(void *, void *, void *);
extern void Q2internals(void *, void *, void *);
extern void signed_ranks(void *, void *, void *);
extern void spat_med_step(void *, void *, void *, void *);
extern void spatial_ranks(void *, void *, void *);
extern void sum_of_diff_sign_outers(void *, void *, void *);
extern void sum_of_rank_outers(void *, void *, void *);
extern void sum_of_sign_outers(void *, void *, void *);
extern void symm_huber(void *, void *, void *, void *, void *, void *);
extern void symm_mvtmle(void *, void *, void *, void *, 
void *);
extern void sum_of_diff_sign_select(void *, void *, void *, 
void *);
extern void symm_huber_inc(void *, void *, void *, void *, 
void *, void *, void *);
extern void symm_mvtmle_inc(void *, void *, void *, void *, 
void *, void *);

static const R_CMethodDef CEntries[] = {
    {"center_step",        (DL_FUNC) &center_step,        4},
    {"hl_center_step",     (DL_FUNC) &hl_center_step,     4},
   /* {"hl_loc_step",        (DL_FUNC) &hl_loc_step,        4}, */
    {"norming",            (DL_FUNC) &norming,            3},
    {"pairdiffc",           (DL_FUNC) &pairdiffc,           3},
    {"pairsumc",            (DL_FUNC) &pairsumc,            3},
    {"Q2internals",        (DL_FUNC) &Q2internals,        3},
    {"signed_ranks",       (DL_FUNC) &signed_ranks,       3},
    {"spat_med_step",      (DL_FUNC) &spat_med_step,      4},
    {"spatial_ranks",      (DL_FUNC) &spatial_ranks,      3},
    {"sum_of_diff_sign_outers", (DL_FUNC)  &sum_of_diff_sign_outers, 3},
    {"sum_of_rank_outers", (DL_FUNC) &sum_of_rank_outers, 3},
    {"sum_of_sign_outers", (DL_FUNC) &sum_of_sign_outers, 3},
    {"symm_huber",         (DL_FUNC) &symm_huber,         6},
    {"symm_huber_inc",     (DL_FUNC) &symm_huber_inc,     7},
    {"sum_of_diff_sign_select", (DL_FUNC) &sum_of_diff_sign_select, 4},
    {"symm_mvtmle",    (DL_FUNC) &symm_mvtmle, 5},
    {"symm_mvtmle_inc",    (DL_FUNC) &symm_mvtmle_inc, 6},
    {NULL, NULL, 0}
};

void R_init_SpatialNP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
