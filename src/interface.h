/* NOTE (fixes_j / CRAN): Changed r_oja() parameter types from long* to int*.
 * On Unix, R int is 4 bytes but C long is 8 bytes, causing incompatible
 * pointer type warnings and potential memory corruption in .C() calls.
 */
#pragma once

extern "C" void r_oja(int* rows, int* cols, double* data, double* vec_out, double* mat_out, int* func, double* param1, double* param2, int* param3, int* param4, int* dbg);