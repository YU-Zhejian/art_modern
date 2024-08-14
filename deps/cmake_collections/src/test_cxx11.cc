#include "test_cxx11_ax_globals.hh"
#include "test_cxx11_globals.hh"

int main(int argc, char** argv)
{
    int ok = 0;
    cxx98main(argc, argv);
    cxx11main(argc, argv);
    return ok;
}
