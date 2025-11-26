// rustc --print target-list

char* ceu_check_get_target_triplet_part_one(void){
#if defined(i386) || defined(__i386) || defined(__i386__) || defined(__i486__) || defined(__i586__) || defined(__i686__) || defined(_M_IX86) || defined(__X86__) || defined(_X86_)
  return "i386";
#elif defined(__x86_64__) || defined(__x86_64) || defined(_M_AMD64) || defined(__amd64__) || defined(__amd64)
    return "x86_64";
#endif
  return "unknown";
}
