#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GOfuncR_binom_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP GOfuncR_binom_randset(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GOfuncR_conti_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP GOfuncR_conti_randset(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GOfuncR_hyper_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP GOfuncR_hyper_randset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GOfuncR_unlock_environment(SEXP);
extern SEXP GOfuncR_wilcox_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP GOfuncR_wilcox_randset(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GOfuncR_binom_category_test",  (DL_FUNC) &GOfuncR_binom_category_test,  4},
    {"GOfuncR_binom_randset",        (DL_FUNC) &GOfuncR_binom_randset,        5},
    {"GOfuncR_conti_category_test",  (DL_FUNC) &GOfuncR_conti_category_test,  4},
    {"GOfuncR_conti_randset",        (DL_FUNC) &GOfuncR_conti_randset,        5},
    {"GOfuncR_hyper_category_test",  (DL_FUNC) &GOfuncR_hyper_category_test,  4},
    {"GOfuncR_hyper_randset",        (DL_FUNC) &GOfuncR_hyper_randset,        6},
    {"GOfuncR_unlock_environment",   (DL_FUNC) &GOfuncR_unlock_environment,   1},
    {"GOfuncR_wilcox_category_test", (DL_FUNC) &GOfuncR_wilcox_category_test, 4},
    {"GOfuncR_wilcox_randset",       (DL_FUNC) &GOfuncR_wilcox_randset,       5},
    {NULL, NULL, 0}
};

void R_init_GOfuncR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
