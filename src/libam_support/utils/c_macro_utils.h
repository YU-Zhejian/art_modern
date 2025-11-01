//
// Created by yuzj on 11/1/25.
//

#ifndef ART_MODERN_LIBAM_SUPPORT_UTILS_C_MACROS_H
#define ART_MODERN_LIBAM_SUPPORT_UTILS_C_MACROS_H

// See: https://www.gnu.org/software/libtool/manual/html_node/C-header-files.html
#undef ART_MODERN_BEGIN_C_DECLS
#undef ART_MODERN_END_C_DECLS

#ifdef __cplusplus
# define ART_MODERN_BEGIN_C_DECLS extern "C" {
# define ART_MODERN_END_C_DECLS }
#else
# define ART_MODERN_BEGIN_C_DECLS /* empty */
# define ART_MODERN_END_C_DECLS /* empty */
#endif

#endif //ART_MODERN_LIBAM_SUPPORT_UTILS_C_MACROS_H
