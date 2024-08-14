#include "test_c99_globals.h"

int main(int argc, char** argv)
{
    int ok = 0;
    ok |= c89main(argc, argv);
    ok |= c99main(argc, argv);
    return ok;
}
