#ifndef _PIRS_H
#define _PIRS_H

#define PROGRAM "pirs (Profile-based Illumina pair-end Reads Simulator)"
#define AUTHOR  "Jianying Yuan (BGI-Shenzhen)"
#define CONTACT "yuanjianying@genomics.org.cn"

#include "config.h"

extern void pirs_simulate(int argc, char *argv[]);
extern void pirs_diploid(int argc, char *argv[]);
extern void program_info();

#endif // _PIRS_H
