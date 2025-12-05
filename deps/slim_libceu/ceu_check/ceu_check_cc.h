/*!
 * @file ceu_check_cc.h
 * @brief Header to get compiler information at compile time
 *
 * @see https://sourceforge.net/p/predef/wiki/Architectures/
 */
#ifndef CEU_CHECK_CC_H
#define CEU_CHECK_CC_H

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @brief Get compiler information.
 *
 * @return Returned buffer, should be freed manually.
 */
char* ceu_check_get_compiler_info(void);

#ifdef __cplusplus
}
#endif
#endif /* CEU_CHECK_CC_H        */
