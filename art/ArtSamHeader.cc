#include "ArtSamHeader.hh"
#include <iostream>
#include <sstream>

using namespace std;
using namespace labw::art_modern;

std::string ArtSamHeader::printHeader() const
{
    ostringstream fout;

    fout << "@HD\t" << "VN:" << VN << "\tSO:" << SO << endl;
    fout << "@PG\t" << "ID:" << ID << "\tPN:" << PN << "\tCL:" << CL << endl;
    for (size_t i = 0; i < SN.size(); i++) {
        fout << "@SQ\tSN:" << SN[i] << "\tLN:" << LN[i] << endl;
    }
    return fout.str();
}
