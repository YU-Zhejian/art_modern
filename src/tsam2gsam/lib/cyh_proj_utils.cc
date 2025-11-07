#include "tsam2gsam/lib/cyh_proj_utils.hh"

#include "libam_support/Constants.hh"
#include "libam_support/Dtypes.h"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace labw::art_modern {
namespace {

    static constexpr char code2base[] = "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
                                        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
                                        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
                                        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
                                        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
                                        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
                                        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
                                        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
                                        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
                                        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
                                        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
                                        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
                                        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
                                        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
                                        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
                                        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    /**
     * From HTSLib 1.22 sam_internals.h.
     *
     * Convert a nibble encoded BAM sequence to a string of bases.
     *
     * We do this 2 bp at a time for speed. Equiv to:
     *
     * @code
     * for (i = 0; i < len; i++)
     *     seq[i] = seq_nt16_str[bam_seqi(nib, i)];
     * @endcode
     *
     * @param nib nibble encoded BAM sequence
     * @param seq destination string of bases.
     * @param len length of the sequence.
     */
    void nibble2base_default(const uint8_t* nib, std::string& seq, const int len)
    {
        int i = len / 2;
        const int len2 = i;
        seq[0] = 0;
        for (i = 0; i < len2; i++) {
            // Note size_t cast helps gcc optimiser.
            std::memcpy(&seq[i << 1], &code2base[static_cast<size_t>(nib[i]) * 2], 2);
        }
        if ((i *= 2) < len) {
            seq[i] = seq_nt16_str[bam_seqi(nib, i)];
        }
    }
} // namespace
std::vector<am_cigar_t> merge_cigars(const std::vector<am_cigar_t>& g_cigar)
{
    std::vector<am_cigar_t> g_cigar_cp = g_cigar;

    // Convert terminal insertions of the CIGAR string to soft clipping.

    for (std::size_t i = 0; i < g_cigar_cp.size(); i += 1) {
        if (bam_cigar_op(g_cigar_cp[i]) == BAM_CINS) {
            g_cigar_cp[i] = bam_cigar_gen(bam_cigar_oplen(g_cigar_cp[i]), BAM_CSOFT_CLIP);
        } else {
            break;
        }
    }
    for (std::ptrdiff_t i = g_cigar_cp.size() - 1; i >= 0; i -= 1) {
        if (bam_cigar_op(g_cigar_cp[i]) == BAM_CINS) {
            g_cigar_cp[i] = bam_cigar_gen(bam_cigar_oplen(g_cigar_cp[i]), BAM_CSOFT_CLIP);
        } else {
            break;
        }
    }

    std::vector<am_cigar_t> g_cigar_merged;
    g_cigar_merged.reserve(g_cigar_cp.size());

    am_cigar_ops_t prev_cigar_op = bam_cigar_op(g_cigar_cp.front());
    am_cigar_len_t prev_cigar_length = bam_cigar_oplen(g_cigar_cp.front());
    std::size_t cigar_id = 0;

    // Skip leading zero-lengthen CIGARs.
    for (; cigar_id < g_cigar_cp.size(); cigar_id += 1) {
        if (bam_cigar_oplen(g_cigar_cp[cigar_id]) != 0) {
            prev_cigar_op = bam_cigar_op(g_cigar_cp[cigar_id]);
            prev_cigar_length = bam_cigar_oplen(g_cigar_cp[cigar_id]);
            cigar_id++;
            break;
        }
    }

    am_cigar_len_t this_cigar_length = 0;
    am_cigar_ops_t this_cigar_op = 0;
    for (; cigar_id < g_cigar_cp.size(); cigar_id += 1) {
        this_cigar_length = bam_cigar_oplen(g_cigar_cp[cigar_id]);
        this_cigar_op = bam_cigar_op(g_cigar_cp[cigar_id]);
        if (this_cigar_length == 0) {
            continue; // Skip CIGAR of zero length.
        }
        if (this_cigar_op == prev_cigar_op) {
            prev_cigar_length += this_cigar_length;
        } else {
            g_cigar_merged.emplace_back(bam_cigar_gen(prev_cigar_length, prev_cigar_op));
            prev_cigar_op = this_cigar_op;
            prev_cigar_length = this_cigar_length;
        }
    }
    if (prev_cigar_length != 0) {
        g_cigar_merged.emplace_back(bam_cigar_gen(prev_cigar_length, prev_cigar_op));
    }
    //    if(g_cigar != g_cigar_merged){
    //        std::cerr << "Merging CIGARs: " << labw::art_modern::cigar_arr_to_str(g_cigar) << " -> " <<
    //        labw::art_modern::cigar_arr_to_str(g_cigar_merged) << std::endl;
    //    }
    return g_cigar_merged;
}

std::string bam_qual_to_str(const bam1_t* aln)
{
    if (bam_get_qual(aln)[0] != 0xff) {
        const auto* quals = bam_get_qual(aln);
        std::string retq;
        retq.resize(aln->core.l_qseq);
        for (decltype(aln->core.l_qseq) k = 0; k < aln->core.l_qseq; k++) {
            retq[k] = static_cast<char>(quals[k] + PHRED_OFFSET);
        }
        std::cerr << retq << std::endl;
        return retq;
    }
    return "*";
}

void clear_pe_flag(uint16_t& flag)
{
    flag &= ~BAM_FREAD1;
    flag &= ~BAM_FREAD2;
    flag &= ~BAM_FMREVERSE;
    flag &= ~BAM_FMUNMAP;
    flag &= ~BAM_FPAIRED;
    flag &= ~BAM_FPROPER_PAIR;
}

std::string bam_seq_to_str(const bam1_t* aln)
{
    const auto* seq = bam_get_seq(aln);
    std::string rets;
    rets.resize(aln->core.l_qseq);
    nibble2base_default(seq, rets, aln->core.l_qseq);
    return rets;
}

void bam_set_seq_eqlen(bam1_t* aln, const char* seq)
{
    auto* cp = bam_get_seq(aln);
    int i = 0;
    for (i = 0; i + 1 < aln->core.l_qseq; i += 2) {
        *cp++ = static_cast<uint8_t>(seq_nt16_table[static_cast<unsigned char>(seq[i])] << 4
            | seq_nt16_table[static_cast<unsigned char>(seq[i + 1])]);
    }
    for (; i < aln->core.l_qseq; i++) {
        *cp++ = static_cast<uint8_t>(seq_nt16_table[static_cast<unsigned char>(seq[i])] << 4);
    }
}

const char* get_sam_mode(const char* file_name, const bool write)
{
    std::string const fstr { file_name };
    const bool ends_in_bam = (fstr.find(".bam") != std::string::npos);
    if (write) {
        return ends_in_bam ? "wb" : "w";
    }
    return ends_in_bam ? "rb" : "r";
}

void print_version(const std::string& prog_name)
{
    std::cout << prog_name << " ver. " << TSAM2GSAM_VERSION << " with HTSLib " << hts_version()
              << ", EWAHBoolArray, cgranges, compiled at " << __DATE__ << ", " << __TIME__ << "." << std::endl;
}
} // namespace labw::art_modern
