#include "test_c11_globals.h"

int main(int argc, char** argv)
{
    int ok = 0;
    ok |= c89main(argc, argv);
    ok |= c99main(argc, argv);
    ok |= c11main();
    return ok;
}
