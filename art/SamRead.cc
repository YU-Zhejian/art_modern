#include <algorithm>
#include <iostream>

#include "SamRead.hh"
#include "seq_utils.hh"

using namespace std;
using namespace labw::art_modern;

void SamRead::reverse_comp()
{
    seq = revcomp(seq);
    reverse(qual.begin(), qual.end());
}

void SamRead::printRead(ostream& fout) const
{
    fout << qname << "\t" << flag << "\t" << rname << "\t" << pos << "\t" << mapQ
         << "\t" << cigar << "\t" << rNext << "\t" << pNext << "\t" << tLen
         << "\t" << seq << "\t" << qual << endl;
}
