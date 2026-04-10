#ifndef SUPERMAG_ERROR_H
#define SUPERMAG_ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

#define SUPERMAG_OK              0
#define SUPERMAG_ERR_NULL_PTR    1
#define SUPERMAG_ERR_INVALID_DIM 2
#define SUPERMAG_ERR_NO_CONVERGE     3
#define SUPERMAG_ERR_ALLOC            4
#define SUPERMAG_ERR_INVALID_MODEL    5
#define SUPERMAG_ERR_NOT_IMPLEMENTED  6

const char* supermag_error_string(int code);

#ifdef __cplusplus
}
#endif

#endif
