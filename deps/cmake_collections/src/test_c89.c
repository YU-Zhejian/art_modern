#include "test_c89_globals.h"

int main(int argc, char** argv)
{
    int ok = 0;
    ok |= c89main(argc, argv);
    return ok;
}
