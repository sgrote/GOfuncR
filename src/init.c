#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _GOfuncR_binom_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GOfuncR_binom_randset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GOfuncR_conti_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GOfuncR_conti_randset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GOfuncR_hyper_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GOfuncR_hyper_randset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GOfuncR_unlock_environment(SEXP);
extern SEXP _GOfuncR_wilcox_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GOfuncR_wilcox_randset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GOfuncR_binom_category_test",  (DL_FUNC) &_GOfuncR_binom_category_test,  4},
    {"_GOfuncR_binom_randset",        (DL_FUNC) &_GOfuncR_binom_randset,        8},
    {"_GOfuncR_conti_category_test",  (DL_FUNC) &_GOfuncR_conti_category_test,  4},
    {"_GOfuncR_conti_randset",        (DL_FUNC) &_GOfuncR_conti_randset,        8},
    {"_GOfuncR_hyper_category_test",  (DL_FUNC) &_GOfuncR_hyper_category_test,  4},
    {"_GOfuncR_hyper_randset",        (DL_FUNC) &_GOfuncR_hyper_randset,        9},
    {"_GOfuncR_unlock_environment",   (DL_FUNC) &_GOfuncR_unlock_environment,   1},
    {"_GOfuncR_wilcox_category_test", (DL_FUNC) &_GOfuncR_wilcox_category_test, 4},
    {"_GOfuncR_wilcox_randset",       (DL_FUNC) &_GOfuncR_wilcox_randset,       8},
    {NULL, NULL, 0}
};

void R_init_GOfuncR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
