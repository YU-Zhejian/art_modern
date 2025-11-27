#ifndef CEU_CC_UNKNOWN_H
#define CEU_CC_UNKNOWN_H
#if !defined(CEU_COMPILER_NAME)
#ifdef CEU_COMPILER_IS_EDG
#define CEU_COMPILER_NAME "unknown EDG-based compiler"
#else

#define CEU_COMPILER_NAME "unknown"
#endif
#endif
#endif /* CEU_CC_UNKNOWN_H */
