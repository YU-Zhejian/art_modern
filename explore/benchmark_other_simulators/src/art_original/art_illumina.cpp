/*
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ART -- Artificial Read Transcription, Illumina Q version 
 * Authors: Weichun Huang 2008-2016
 * Contributors: Jason Myers (2011) 
 * License: GPL v3 
 * ############################################################################
 * #    This program is free software: you can redistribute it and/or modify  #
 * #    it under the terms of the GNU General Public License as published by  #
 * #    the Free Software Foundation, either version 3 of the License, or     #
 * #    (at your option) any later version.                                   #
 * #                                                                          #
 * #    This program is distributed in the hope that it will be useful,       #
 * #    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
 * #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
 * #    GNU General Public License for more details.                          #
 * #                                                                          #
 * #    You should have received a copy of the GNU General Public License     #
 * #    along with this program.  If not, see <http://www.gnu.org/licenses/>. #
 * ############################################################################
 * <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*/

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <time.h>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include "art.h"
#include "empdist.h"
#include "readSeqFile.h"
#include "samRead.h"


using namespace std;
double empdist::prob_err[HIGHEST_QUAL];
gsl_rng* art::gsl_R;
int art::gaussain_mean;
double art::gaussain_sigma;

int main(int argc, char* argv[]){
    cout <<endl;
    cout << "    ====================ART===================="<<endl;
    cout << "             ART_Illumina (2008-2016)          "<<endl;
    cout << "          Q Version 2.5.8 (June 6, 2016)       "<<endl;
    cout << "     Contact: Weichun Huang <whduke@gmail.com> "<<endl; 
    cout << "    -------------------------------------------"<<endl<<endl;

    short min_qual_s=0; 
    short max_qual_s=93; 
    const short max_qual=93; 
    char max_q_c = (char)(max_qual_s+33);
    bool mask_n=true; 
    short max_num_n=1; 
    long mask_read_count=0;
    int len_ref_id=250;

    int maxNumIndel = -1; //a negative value means no limit
    //caluate CPUT time
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    // Boolean flags showing the state of the run
    bool err_free_sam=false;
    bool sam_out=false;
    bool is_matepair=false;
    bool arg_success = true;
    bool is_pairend_read=false;
    bool sep_flag = false;
    bool mean_flag = false;
    bool sDev_flag = false;
    bool in_flag = false;
    bool out_flag = false;
    bool len_flag = false;
    bool fold_flag = false;
    bool cc_flag = false;
    bool help_flag = false;
    bool show_flag = false;
    bool qs_flag = false;
    bool rate_flag = false;
    bool first_qual = false;
    bool second_qual = false;
    bool no_ALN = false;
    bool fixed_seed = false;
    unsigned int rand_seed = 0;
    bool amplicon = false;
    bool use_cigarM = false;
//    bool use_org_qscore = false;

    // Command-line Arguments
    string suffID=""; 
    string qual_file1 = "";
    string qual_file2 = "";
    string seqsys="";
    char* seq_file = (char*) "";
    string out_file_prefix = "";
    int read_len  = 0;
    string p_cigar="";
    double x_fold  = 0;
    long read_count = 0;
    string num= "";
    int mean = 0;
    short q_shift_up=0;//2;
    short q_shift_up_2=0;//4;
    double std_dev = 0.0;
    double insRate=0.00009;
    double delRate=0.00011;
    double insRate2 = 0.00015;
    double delRate2 = 0.00023;

    // Get Command-line arguments
    for(int i=1;i < argc;++i){
        char* arg = argv[i];
            if(!strcmp(arg, "--matepair") || !strcmp(arg, "-mp") ){
        	num="1";
        	is_pairend_read=true;
        	is_matepair=true;
	    }
	    else if(!strcmp(arg, "--samout") || !strcmp(arg, "-sam") ){
        	sam_out=true;
	    }
	    else if(!strcmp(arg, "--cigarM") || !strcmp(arg, "-M") ){
        	use_cigarM=true;
	    }
	    else if(!strcmp(arg, "--errfree") || !strcmp(arg, "-ef") ){
        	err_free_sam=true;
        	sam_out=true;
	    }
	    else if(!strcmp(arg, "--maskN") || !strcmp(arg, "-nf") ){ 
		i++;
		max_num_n = atoi(argv[i]);
		if(max_num_n>0){ mask_n=true; }
		else{ mask_n=false; }
	    }
	    else if(!strcmp(arg, "--amplicon") || !strcmp(arg, "-amp") ){
        	amplicon=true;
	    }
	    else if(!strcmp(arg, "--paired") || !strcmp(arg, "-p") ){
        	num="1";
        	is_pairend_read=true;
            } else if(!strcmp(arg, "--help") || !strcmp(arg, "-h")){
		arg_success = false;
		help_flag = true;
		break;
            } else if(!strcmp(arg, "--quiet") || !strcmp(arg, "-q")){
		show_flag = true;
            } else if(!strcmp(arg, "--in") || !strcmp(arg, "-i")){
		i++;
		seq_file = argv[i];
		in_flag = true;
            } else if(!strcmp(arg, "--out") || !strcmp(arg, "-o")){
		i++;
		out_file_prefix = argv[i];
		out_flag = true;
            } else if(!strcmp(arg, "--len") || !strcmp(arg, "-l")){
		i++;
		read_len = atoi(argv[i]);
		p_cigar.append(argv[i]).append("=");
		len_flag = true;
		if(read_len <= 0){
		   cerr << "Fatal Error: The read length must be a positive integer." << endl;
		   arg_success = false;
		}
            } else if(!strcmp(arg, "--maxIndel") || !strcmp(arg, "-k")){
		i++;
		maxNumIndel = atoi(argv[i]);
            } else if(!strcmp(arg, "--insRate") || !strcmp(arg, "-ir")){
		i++;
		insRate = atof(argv[i]);
		rate_flag = true;
            } else if(!strcmp(arg, "--delRate") || !strcmp(arg, "-dr")){
		i++;
		delRate = atof(argv[i]);
		rate_flag = true;
            } else if(!strcmp(arg, "--insRate2") || !strcmp(arg, "-ir2")){
		i++;
		insRate2 = atof(argv[i]);
		rate_flag = true;
            } else if(!strcmp(arg, "--delRate2") || !strcmp(arg, "-dr2")){
		i++;
		delRate2 = atof(argv[i]);
		rate_flag = true;
            } else if(!strcmp(arg, "--fcov") || !strcmp(arg, "-f")){
		i++;
		x_fold = atof(argv[i]);
		fold_flag = true;
		if(x_fold < 0){
		   cerr << "Fatal Error: The fold coverage must be a positive." << endl;
		   arg_success = false;
		}
            } else if(!strcmp(arg, "--rcount") || !strcmp(arg, "-c")){
		i++;
		read_count= atoi(argv[i]);
		cc_flag = true;
		if(read_count< 0){
		   cerr << "Fatal Error: the read count must be a positive integer." << endl;
		   arg_success = false;
		}
            } else if(!strcmp(arg, "--mflen") || !strcmp(arg, "-m")){
		i++;
		mean = atoi(argv[i]);
		mean_flag = true;
		if(mean < 0){
		   cerr << "Input Error: The mean fragment length must be a positive." << endl;
		   arg_success = false;
		}
            } else if(!strcmp(arg, "--minQ") || !strcmp(arg, "-qL")){
		i++;
		min_qual_s = atoi(argv[i]);
		if(min_qual_s<0 || min_qual_s>max_qual){
		   cerr << "Input Error: The minimum quality score must be an integer in [0,"<<max_qual<<"]" << endl;
		   arg_success = false;
		}
            } else if(!strcmp(arg, "--maxQ") || !strcmp(arg, "-qU")){
		i++;
		max_qual_s = atoi(argv[i]);
		if(max_qual_s<=0 || max_qual_s>max_qual){
		   cerr << "Input Error: The quality score must be an integer in [1,"<<max_qual<<"]"<< endl;
		   arg_success = false;
		}
            } else if(!strcmp(arg, "--qShift") || !strcmp(arg, "-qs")){
		i++;
		q_shift_up = atoi(argv[i]);
		qs_flag = true;
            } else if(!strcmp(arg, "--qShift2") || !strcmp(arg, "-qs2")){
		i++;
		q_shift_up_2 = atoi(argv[i]);
		qs_flag = true;
//            } else if(!strcmp(arg, "--qOrig") || !strcmp(arg, "-qo")){ 
//		use_org_qscore = true;
	    } else if(!strcmp(arg, "--sdev") || !strcmp(arg, "-s")){
		i++;
		std_dev = atof(argv[i]);
		sDev_flag = true;
            } else if(!strcmp(arg, "--seqSys") || !strcmp(arg, "-ss")){
		i++;
		seqsys = argv[i];
            } else if(!strcmp(arg, "--qprof1") || !strcmp(arg, "-1")){
		i++;
		qual_file1 = argv[i];
		first_qual = true;
            } else if(!strcmp(arg, "--qprof2") || !strcmp(arg, "-2")){
		i++;
		qual_file2 = argv[i];
		second_qual = true;
            } else if(!strcmp(arg, "--sepProf") || !strcmp(arg, "-sp")){
//		i++;
		sep_flag = true;
            } else if(!strcmp(arg, "--id") || !strcmp(arg, "-d")){
		i++;
		suffID = argv[i];
            } else if(!strcmp(arg, "--noALN") || !strcmp(arg, "-na")){ 
		no_ALN = true;
            } else if(!strcmp(arg, "--rndSeed") || !strcmp(arg, "-rs")){
		i++;
		rand_seed = abs(atoi(argv[i]));
		fixed_seed=true;
            } else {
		arg_success = false;
		cerr << "Fatal Error: " << arg << ", is not a valid parameter." << endl;
		break;
            }
	
    }

    if(no_ALN && !sam_out){
		cerr << "Warning: your simulation will not output any ALN or SAM file with your parameter settings!" << endl;
    }
	
    // Make sure the minimum requirements to run were met if the help tag was not given
    if(help_flag){
	arg_success = false;
    }
    else if (fold_flag && cc_flag){
	cerr << "Error: please use only one of the two parameters: read count (-c) or fold coverage (-f)." << endl << endl;
	exit(1);
    }
    else if(!in_flag || !out_flag || !len_flag || !(fold_flag || cc_flag)){
	arg_success = false;
//	cerr << "Fatal Error: An input-file, output-file prefix, read length, and fold coverage must be specified." << endl << endl;
    }

    if(mean_flag && sDev_flag){
	is_pairend_read = true;
	if(mean>=2000){
	       	is_matepair = true;
	}
	else if(is_matepair){
		cerr<<"Warning: a mate-pair simulation may be not appropriate for DNA fragment size < 2000bp"<<endl;
	}
	num = "1";
    }

    // Make sure the minimum requirements to run were given for a paired end simulation
    if (is_pairend_read){ 
	    //use the both default profiles or provide both profiles
	    if(qual_file1.empty() && !qual_file2.empty() ){
		    cerr<<"Please provide the quality profile of the first read" <<endl;
		    exit(1);
	    }
	    else if(!qual_file1.empty() && qual_file2.empty() ){
		    cerr<<"Please provide the quality profile of the second read" <<endl;
		    exit(1);
	    }
	   
	    if(mean_flag && sDev_flag){
//		    art::ini_read_pair_rand(abs(mean),fabs(std_dev));
		    art::ini_read_pair_rand(abs(mean),fabs(std_dev),rand_seed);
		    if(art::gaussain_mean<=read_len){
			    cerr<<"Fatal Error: The read length must be shorter than the mean fragment length specified." <<endl;
			    exit(1);
		    }
	    } else if (!amplicon){
		    cerr << "Fatal Error: A mean fragment length and a standard deviation must be specified." << endl << endl;
		    exit(1);
	    }
    }

    if(!arg_success){
	cout << "===== USAGE ====="<<endl << endl;
	cout << "art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -o <outfile_prefix>"<<endl;
	cout << "art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -c <num_reads_per_sequence> -o <outfile_prefix>"<<endl;
	cout << "art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -m <mean_fragsize> -s <std_fragsize> -o <outfile_prefix>"<<endl;
	cout << "art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -c <num_reads_per_sequence> -m <mean_fragsize> -s <std_fragsize> -o <outfile_prefix>"<<endl<<endl;
	cout << "===== PARAMETERS =====" << endl << endl;
	cout << "  -1   --qprof1   the first-read quality profile" << endl;
	cout << "  -2   --qprof2   the second-read quality profile" << endl;
	cout << "  -amp --amplicon amplicon sequencing simulation" << endl;
	cout << "  -c   --rcount   number of reads/read pairs to be generated per sequence/amplicon (not be used together with -f/--fcov)" << endl;
	cout << "  -d   --id       the prefix identification tag for read ID" << endl;
	cout << "  -ef  --errfree  indicate to generate the zero sequencing errors SAM file as well the regular one" << endl;
       	cout << "                  NOTE: the reads in the zero-error SAM file have the same alignment positions"<<endl;
       	cout << "                  as those in the regular SAM file, but have no sequencing errors"<<endl;
	cout << "  -f   --fcov     the fold of read coverage to be simulated or number of reads/read pairs generated for each amplicon" << endl;
	cout << "  -h   --help     print out usage information"<<endl;
	cout << "  -i   --in       the filename of input DNA/RNA reference" << endl;
	cout << "  -ir  --insRate  the first-read insertion rate (default: 0.00009)"<< endl;
	cout << "  -ir2 --insRate2 the second-read insertion rate (default: 0.00015)" << endl;
	cout << "  -dr  --delRate  the first-read deletion rate (default:  0.00011)" << endl;
	cout << "  -dr2 --delRate2 the second-read deletion rate (default: 0.00023)" << endl;
	cout << "  -k   --maxIndel the maximum total number of insertion and deletion per read (default: up to read length)" << endl;
	cout << "  -l   --len      the length of reads to be simulated" << endl;
	cout << "  -m   --mflen    the mean size of DNA/RNA fragments for paired-end simulations" << endl;
	cout << "  -mp  --matepair indicate a mate-pair read simulation" << endl;
	cout << "  -M  --cigarM    indicate to use CIGAR 'M' instead of '=/X' for alignment match/mismatch" << endl;
	cout << "  -nf  --maskN    the cutoff frequency of 'N' in a window size of the read length for masking genomic regions"<<endl;
       	cout << "                  NOTE: default: '-nf 1' to mask all regions with 'N'. Use '-nf 0' to turn off masking" <<endl;
	cout << "  -na  --noALN    do not output ALN alignment file" << endl;
	cout << "  -o   --out      the prefix of output filename" << endl;
	cout << "  -p   --paired   indicate a paired-end read simulation or to generate reads from both ends of amplicons" << endl;
       	cout << "                  NOTE: art will automatically switch to a mate-pair simulation if the given mean fragment size >= 2000"<<endl;
	cout << "  -q   --quiet    turn off end of run summary" << endl;
	cout << "  -qL  --minQ     the minimum base quality score" << endl;
	cout << "  -qU  --maxQ     the maxiumum base quality score" << endl;
	cout << "  -qs  --qShift   the amount to shift every first-read quality score by " << endl;
	cout << "  -qs2 --qShift2  the amount to shift every second-read quality score by" << endl;
       	cout << "                  NOTE: For -qs/-qs2 option, a positive number will shift up quality scores (the max is 93) " <<endl;
       	cout << "                  that reduce substitution sequencing errors and a negative number will shift down " <<endl;
       	cout << "                  quality scores that increase sequencing errors. If shifting scores by x, the error"<<endl;
       	cout << "                  rate will be 1/(10^(x/10)) of the default profile." <<endl;
//	cout << "  -qo  --qOrig    indicate to output the original quality scores regardless the original scores were shifted or not" << endl;
	cout << "  -rs  --rndSeed  the seed for random number generator (default: system time in second)"  << endl;
       	cout << "                  NOTE: using a fixed seed to generate two identical datasets from different runs"<<endl;
	cout << "  -s   --sdev     the standard deviation of DNA/RNA fragment size for paired-end simulations." << endl;
	cout << "  -sam --samout   indicate to generate SAM alignment file" << endl;
	cout << "  -sp  --sepProf  indicate to use separate quality profiles for different bases (ATGC)" << endl;
	cout << "  -ss  --seqSys   The name of Illumina sequencing system of the built-in profile used for simulation" << endl;
       	cout << "       NOTE: sequencing system ID names are:"<<endl;
        cout << "            GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)"<<endl;
	cout << "           HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)"<<endl;
	cout << "           HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)"<<endl;
	cout << "           MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)"<<endl;

	cout << "===== NOTES ====="<< endl<< endl;
	cout << "* ART by default selects a built-in quality score profile according to the read length specified for the run." << endl << endl;
	cout << "* For single-end simulation, ART requires input sequence file, output file prefix, read length, and read count/fold coverage." << endl << endl;;
//	cout << "  Example: art --in reference_DNA.fa --out sim1 --len 35 --fcov 2 -sam" << endl << endl;;
	cout << "* For paired-end simulation (except for amplicon sequencing), ART also requires the parameter values of" << endl;
	cout << "  the mean and standard deviation of DNA/RNA fragment lengths" << endl << endl;
//	cout << "  Example: art --paired --in reference_DNA.fa --out sim2 --len 35 --fcov 2 --mflen 200 --sdev 3.5 -sam" << endl << endl;

	cout << "===== EXAMPLES ====="<< endl <<endl;
	cout << " 1) single-end read simulation" <<endl;
	cout << " 	art_illumina -ss HS25 -sam -i reference.fa -l 150 -f 10 -o single_dat" <<endl<<endl;
	cout << " 2) paired-end read simulation" <<endl;
	cout << "       art_illumina -ss HS25 -sam -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_dat" <<endl<<endl;
	cout << " 3) mate-pair read simulation" <<endl;
	cout << "       art_illumina -ss HS10 -sam -i reference.fa -mp -l 100 -f 20 -m 2500 -s 50 -o matepair_dat" <<endl<<endl;
	cout << " 4) amplicon sequencing simulation with 5' end single-end reads " <<endl;
	cout << " 	art_illumina -ss GA2 -amp -sam -na -i amp_reference.fa -l 50 -f 10 -o amplicon_5end_dat" <<endl<<endl;
	cout << " 5) amplicon sequencing simulation with paired-end reads" <<endl;
	cout << "       art_illumina -ss GA2 -amp -p -sam -na -i amp_reference.fa -l 50 -f 10 -o amplicon_pair_dat" <<endl<<endl;
	cout << " 6) amplicon sequencing simulation with matepair reads" <<endl;
	cout << "       art_illumina -ss MSv1 -amp -mp -sam -na -i amp_reference.fa -l 150 -f 10 -o amplicon_mate_dat" <<endl<<endl;
	cout << " 7) generate an extra SAM file with zero-sequencing errors for a paired-end read simulation" <<endl;
	cout << "       art_illumina -ss HSXn -ef -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_twosam_dat" <<endl<<endl;
	cout << " 8) reduce the substitution error rate to one 10th of the default profile"<<endl;
       	cout << "       art_illumina -i reference.fa -qs 10 -qs2 10 -l 50 -f 10 -p -m 500 -s 10 -sam -o reduce_error"<<endl<<endl;
	cout << " 9) turn off the masking of genomic regions with unknown nucleotides 'N'"<<endl;
	cout << "       art_illumina -ss HS20 -nf 0  -sam -i reference.fa -p -l 100 -f 20 -m 200 -s 10 -o paired_nomask" <<endl<<endl;
	cout << " 10) masking genomic regions with >=5 'N's within the read length 50"<<endl;
	cout << "       art_illumina -ss HSXt -nf 5 -sam -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_maskN5" <<endl<<endl;

	exit(1);
    }

    if(max_num_n>read_len){
		cerr << "Error: the cutoff frequency of 'N' for maksing genome region: "<<max_num_n<<" > the read length (" <<read_len<<")"<<endl;
		exit(1);
    }
    if(!fixed_seed){
	    rand_seed=(unsigned int) time(NULL); 
    }
    srand (rand_seed);

//    string seqfasta=out_file_prefix+num+".fa";
//    string qualfasta=out_file_prefix+num+".qual";
    string alnfasta=out_file_prefix+num+".aln";
    string fqfile=out_file_prefix+num+".fq";

    string samfile=out_file_prefix+".sam";
    ofstream SAMFILE;
    string samfile_ef=out_file_prefix+"_errFree.sam";
    ofstream SAMFILE_EF;
    if(sam_out) {
	    SAMFILE.open(samfile.c_str(),ios::binary); 
	    if(!SAMFILE.is_open()) { cerr<<"Can not open output file: "<<samfile<<endl; exit(0); }
	    if(err_free_sam){
		    SAMFILE_EF.open(samfile_ef.c_str(),ios::binary); 
		    if(!SAMFILE_EF.is_open()) { cerr<<"Can not open output file: "<<samfile_ef<<endl; exit(0); }
	    }
    }

    ofstream FQFILE(fqfile.c_str(),ios::binary);
    if(!FQFILE.is_open()) { cerr<<"Can not open output file: "<<fqfile<<endl; exit(0); }
	ofstream ALNFILE;	
	if (!no_ALN) {
	       	ALNFILE.open(alnfasta.c_str(),ios::binary);
	       	if(!ALNFILE.is_open()) { cerr<<"Can not open output file: "<<alnfasta<<endl; exit(0); }
	       	ALNFILE<<"##ART_Illumina\tread_length\t"<<read_len<<endl;
	}	
    empdist::set_err_prob();
    empdist qdist;

    if(qual_file1.empty() && !seqsys.empty()){
	    qdist.setdist(seqsys, sep_flag, read_len);
    }
    else {
	    qdist.setdist(qual_file1, qual_file2, sep_flag, read_len);
	    if (!seqsys.empty()){
		    cout<<"Warnning: ART ignores the specified sequencing system as a customized quality profile provided.";
	    }
    }

    int profile_size=qdist.qual_dist_first.size();
    int profile_size_2=qdist.qual_dist_second.size();
    if(sep_flag){
	    profile_size=qdist.a_qual_dist_first.size();
	    profile_size_2=qdist.a_qual_dist_second.size();
    }

    if(read_len > profile_size){
	if(profile_size == 0){
	    cerr << "Fatal Error: " <<  qual_file1 << ", is not a valid profile." << endl << endl;
	} else {
            cerr<<"Fatal Error: The read length, "<<read_len<<", exceeds the maximum first read profile length, "<<profile_size<< "." <<endl<<endl;
	}
        exit(1);
    }

    if((read_len > profile_size_2)  && is_pairend_read){
	if(profile_size_2 == 0){
	    if (!seqsys.empty())
		    cerr << "Error: the built-in " << seqsys<<"  sequencing system supports only single-end sequencing" << endl << endl;
	    else 
		    cerr << "Fatal Error: " <<  qual_file2 << ", is not a valid profile." << endl << endl;
	} else {
            cerr<<"Fatal Error: The read length, "<<read_len<<", exceeds the maximum second read profile length, " <<profile_size_2 <<"." <<endl<<endl;
	}
        exit(1);
    }

    if(is_pairend_read && profile_size != profile_size_2){
        cerr<<"Fatal Error: The first profile read length, " << qdist.qual_dist_first.size()
	<< ", does not match the second profile read length, " << qdist.qual_dist_second.size() << "."
	<< endl << endl;
	exit(1);
    }

    //increase overal base quality value
    if(!sep_flag){
	    for(size_t i=0; i<qdist.qual_dist_first.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.qual_dist_first[i].begin(); it!=qdist.qual_dist_first[i].end(); it++){
			    if(q_shift_up!=0){
				    if(q_shift_up<0 && (-q_shift_up>it->second)){ it->second=min_qual_s; }
				    else{ it->second+=q_shift_up; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }
	    for(size_t i=0; i<qdist.qual_dist_second.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.qual_dist_second[i].begin(); it!=qdist.qual_dist_second[i].end(); it++){
			    if(q_shift_up_2!=0){
				    if(q_shift_up_2<0 && (-q_shift_up_2>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up_2; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }

    } else{
	    for(size_t i=0; i<qdist.a_qual_dist_first.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.a_qual_dist_first[i].begin(); it!=qdist.a_qual_dist_first[i].end(); it++){
			    if(q_shift_up!=0){
				    if(q_shift_up<0 && (-q_shift_up>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }
	    for(size_t i=0; i<qdist.a_qual_dist_second.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.a_qual_dist_second[i].begin(); it!=qdist.a_qual_dist_second[i].end(); it++){
			    if(q_shift_up_2!=0){
				    if(q_shift_up_2<0 && (-q_shift_up_2>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up_2; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }

	    for(size_t i=0; i<qdist.c_qual_dist_first.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.c_qual_dist_first[i].begin(); it!=qdist.c_qual_dist_first[i].end(); it++){
			    if(q_shift_up!=0){
				    if(q_shift_up<0 && (-q_shift_up>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }
	    for(size_t i=0; i<qdist.c_qual_dist_second.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.c_qual_dist_second[i].begin(); it!=qdist.c_qual_dist_second[i].end(); it++){
			    if(q_shift_up_2!=0){
				    if(q_shift_up_2<0 && (-q_shift_up_2>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up_2; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }

	    for(size_t i=0; i<qdist.g_qual_dist_first.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.g_qual_dist_first[i].begin(); it!=qdist.g_qual_dist_first[i].end(); it++){
			    if(q_shift_up!=0){
				    if(q_shift_up<0 && (-q_shift_up>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }
	    for(size_t i=0; i<qdist.g_qual_dist_second.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.g_qual_dist_second[i].begin(); it!=qdist.g_qual_dist_second[i].end(); it++){
			    if(q_shift_up_2!=0){
				    if(q_shift_up_2<0 && (-q_shift_up_2>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up_2; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }

	    for(size_t i=0; i<qdist.t_qual_dist_first.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.t_qual_dist_first[i].begin(); it!=qdist.t_qual_dist_first[i].end(); it++){
			    if(q_shift_up!=0){
				    if(q_shift_up<0 && (-q_shift_up>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }
	    for(size_t i=0; i<qdist.t_qual_dist_second.size(); i++){
		    for(map<unsigned int, unsigned short>::iterator it=qdist.t_qual_dist_second[i].begin(); it!=qdist.t_qual_dist_second[i].end(); it++){
			    if(q_shift_up_2!=0){
				    if(q_shift_up_2<0 && (-q_shift_up_2>it->second)){ it->second=0; }
				    else{ it->second+=q_shift_up_2; if(it->second>max_qual_s) it->second=max_qual_s; }
			    }
			    if(it->second < min_qual_s) it->second = min_qual_s;
			    if(it->second > max_qual_s) it->second = max_qual_s;
		    }
	    }
   }
    samHeader sH;
    sH.getRefseqID(seq_file);
    sH.ID="01";
    sH.PN="ART_Illumina"; 
    for(int i=0;i < argc; ++i) { sH.CL.append(argv[i]);sH.CL.append(" "); }
    if(!fixed_seed){
	    stringstream tss;
	    tss<<"-rs "<<rand_seed;
	    sH.CL.append(tss.str());
    }
    if (!no_ALN){
	    sH.printAlnHeader(ALNFILE);
    }

    if(sam_out){
	    sH.printHeader(SAMFILE);
	    if(err_free_sam) sH.printHeader(SAMFILE_EF);
    }

    samRead sR;
    string srID;

    vector<short> qual;
    readSeqFile seq_reader(seq_file);
    string id;
    art a_art; 
    seqRead a_read;

//    a_read.set_rate(read_len,insRate,2,a_read.ins_rate);
//    a_read.set_rate(read_len,delRate,2,a_read.del_rate);
    a_read.set_rate(read_len, insRate, a_read.ins_rate, maxNumIndel);
    a_read.set_rate(read_len, delRate, a_read.del_rate, maxNumIndel);
//    void set_rate(int read_len, double p, vector <double>& rate, int max_num=0, double cdf_cutoff=0.999999){
    string aln_read,aln_ref;
    ostringstream osID;
    int num_seq=0;
    string read_id;
//    string seqfasta2="";
//    string qualfasta2="";
    string alnfasta2="";
    string fqfile2="";
    if(is_pairend_read){
	samRead sR2; 
	sR.rNext="=";
	sR2.rNext="=";

//        seqfasta2=out_file_prefix+"2.fa";
//        qualfasta2=out_file_prefix+"2.qual";
        alnfasta2=out_file_prefix+"2.aln";
        fqfile2=out_file_prefix+"2.fq";
//        ofstream SEQFILE2(seqfasta2.c_str(),ios::binary);
//        if(!SEQFILE2.is_open()) { cout<<"can not open output file: "<<seqfasta2<<endl; exit(0); }

//        ofstream QUALFILE2(qualfasta2.c_str(),ios::binary);
//        if(!QUALFILE2.is_open()) { cout<<"can not open output file: "<<qualfasta2<<endl; exit(0); }

        ofstream FQFILE2(fqfile2.c_str(),ios::binary);
        if(!FQFILE2.is_open()) { cerr<<"Can not open output file: "<<fqfile2<<endl; exit(0); }
       	ofstream ALNFILE2;
       	if (!no_ALN) {
	       	ALNFILE2.open(alnfasta2.c_str(),ios::binary);
	       	if(!ALNFILE2.is_open()) { cerr<<"Can not open output file: "<<alnfasta2<<endl; exit(0); }
	       	ALNFILE2<<"##ART_Illumina\tread_length\t"<<read_len<<endl;
	       	sH.printAlnHeader(ALNFILE2);
       	}

        seqRead a_read_2;

//        a_read_2.set_rate(read_len,insRate2,2,a_read_2.ins_rate);
//        a_read_2.set_rate(read_len,delRate2,2,a_read_2.del_rate);

        a_read_2.set_rate(read_len,insRate2,a_read_2.ins_rate, maxNumIndel);
        a_read_2.set_rate(read_len,delRate2,a_read_2.del_rate, maxNumIndel);

        vector<short> qual_2;
        string read_id_2;
        string aln_read_2,aln_ref_2;
        while(seq_reader.next_seq(id,a_art.ref_seq)){ 
	    std::replace(a_art.ref_seq.begin(), a_art.ref_seq.end(), 'U', 'T'); //replace U with T
//            size_t p1=id.find_first_of(' '); if(p1==string::npos) p1=10; size_t p2=id.find_first_of('\t'); if(p2==string::npos) p2=10;            p1=p1<p2?p1:p2; id=id.substr(0,p1); 
            istringstream isID; isID.str(id); isID>>id; id=id.substr(0,len_ref_id); 
            num_seq++;
            a_art.ini_set(read_len);
            if(mask_n){ 
              a_art.mask_n_region(max_num_n);
            }
            //long t_num_read=(unsigned long) a_art.ref_seq.size()/read_len*x_fold;
            long t_num_read= 0;
	    if(cc_flag){
		     t_num_read= 2*read_count;
	    }
	    else{
		    t_num_read= static_cast<unsigned long>(a_art.ref_seq.size()/read_len*x_fold);
		    if(amplicon){ t_num_read = 2*(long) x_fold; }
	    }
            while(t_num_read>0){
//              osID<<num_seq<<fixed<<setfill('0')<< setw(10)<< t_num_read;
                osID<<id<<'-'<<suffID<<t_num_read;
                read_id = osID.str();
                osID.str("");
                a_read.clear();
                a_read_2.clear();
                //a_art.next_read(a_read);
		bool has_reads=true;
		if(is_matepair){
		       	if(amplicon) has_reads=a_art.next_matepair_ampread_indel_cmp(a_read, a_read_2);
		       	else has_reads=a_art.next_pair_read_indel_mate(a_read, a_read_2); 
		}
	       	else{
		       	if(amplicon) has_reads=a_art.next_pair_ampread_indel_cmp(a_read, a_read_2); 
			else has_reads=a_art.next_pair_read_indel_cmp(a_read, a_read_2); 
		}
		if (!has_reads){
			cerr<<"Warning: the reference sequence "<<id<< " (length "<<a_art.ref_seq.size()<<"bps ) is skipped as it < the defined read length ("<< read_len<<" bps)"<<endl;
			break;
		}
                if(mask_n){ 
                  if(a_read.is_plus_strand){
                    size_t bpos2=a_art.ref_seq.size()-a_read_2.bpos-read_len;
                    if(a_art.masked_pos.count(a_read.bpos)>0 || a_art.masked_pos.count(bpos2)>0){
                      t_num_read-=2;
                      continue;
                    }
                  }
                  else{
                    size_t bpos1=a_art.ref_seq.size()-a_read.bpos-read_len;
                    if(a_art.masked_pos.count(bpos1)>0 || a_art.masked_pos.count(a_read_2.bpos)>0){
                      t_num_read-=2;
                      continue;
                    }
                  }
                }
                qual.clear();
                qual_2.clear();

		if(!sep_flag){
		       	qdist.get_read_qual(qual, read_len, true);
		       	qdist.get_read_qual(qual_2,read_len, false);
	       	}
	       	else{
		       	qdist.get_read_qual_1st(a_read.seq_read, qual);
		       	qdist.get_read_qual_2nd(a_read_2.seq_read, qual_2);
	       	}

                a_read.add_error(qual);
                a_read_2.add_error(qual_2);
		srID=read_id;
                read_id_2=read_id+"/2";
                read_id+="/1";

                FQFILE<<"@"<<read_id<<endl<<a_read.seq_read<<endl<<"+"<<endl;
                for(size_t k=0; k<qual.size(); k++){
                    FQFILE<<(char)(qual[k]+33);
                }
                FQFILE<<endl;

		if(! a_read.get_aln(aln_read,aln_ref)){
		       	aln_ref=a_read.seq_ref;
		       	aln_read=a_read.seq_read;
	       	}
		if (!no_ALN) {
		       	ALNFILE<<">"<<id<<"\t"<<read_id<<"-1\t"<<a_read.bpos;
		       	if(a_read.is_plus_strand) ALNFILE<<"\t+\n";
		       	else ALNFILE<<"\t-\n";
		       	ALNFILE<<aln_ref<<endl<<aln_read<<endl;
	       	}

                FQFILE2<<"@"<<read_id_2<<endl<<a_read_2.seq_read<<endl<<"+"<<endl;
                for(size_t k=0; k<qual_2.size(); k++){
                    FQFILE2<<(char)(qual_2[k]+33);
                }
                FQFILE2<<endl;
	       	if(! a_read_2.get_aln(aln_read_2,aln_ref_2)){
		       	aln_ref_2=a_read_2.seq_ref;
		       	aln_read_2=a_read_2.seq_read;
	       	}
	       	if (!no_ALN) { 
			ALNFILE2<<">"<<id<<"\t"<<read_id_2<<"\t"<<a_read_2.bpos;
		       	if(a_read_2.is_plus_strand) ALNFILE2<<"\t+\n";
		       	else ALNFILE2<<"\t-\n";
		       	ALNFILE2<<aln_ref_2<<endl<<aln_read_2<<endl;
	       	}

		if(sam_out){
		       	sR.qname=srID;
		       	sR.rname=id;

		       	sR2.qname=srID;
		       	sR2.rname=id;

			sR.seq=a_read.seq_read;
		       	sR.qual.resize(qual.size());
		       	for(size_t k=0; k<qual.size(); k++){
			       	sR.qual[k]=(char)(qual[k]+33);
		       	}
		
			sR2.seq=a_read_2.seq_read;
		       	sR2.qual.resize(qual_2.size());
		       	for(size_t k=0; k<qual_2.size(); k++){
			       	sR2.qual[k]=(char)(qual_2[k]+33);
		       	}

			sR.flag=0x01 | 0x02 | 0x40;
		       	sR2.flag=0x01 | 0x02 | 0x80;
		       	if(a_read.is_plus_strand){
				       	sR.pos=a_read.bpos+1;
				       	sR.flag =sR.flag | 0x20;
				       	sR2.flag =sR2.flag |0x10;
				       	sR2.pos=a_art.ref_seq.size()-(a_read_2.bpos+read_len-1);
					sR2.reverse_comp();
		       	}
		       	else{
				       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+read_len-1);
				       	sR.reverse_comp();
				       	sR2.pos=a_read_2.bpos+1;
				       	sR.flag = sR.flag | 0x10;
				       	sR2.flag = sR2.flag | 0x20;
		       	}
		       	sR.getCigar(aln_ref,aln_read,use_cigarM);
			sR2.getCigar(aln_ref_2,aln_read_2,use_cigarM);
		       	sR.pNext=sR2.pos;
		       	sR2.pNext=sR.pos;
		       	if(sR2.pos>sR.pos){
			       	sR.tLen=sR2.pos+a_read_2.seq_read.size()-sR.pos;
			       	sR2.tLen=-sR.tLen;
		       	}
		       	else{
			       	sR2.tLen=sR.pos+a_read.seq_read.size()-sR2.pos;
			       	sR.tLen=-sR2.tLen;
		       	} 
			sR.printRead(SAMFILE);
			sR2.printRead(SAMFILE);
			if(err_free_sam){
			       	sR.seq=a_read.seq_ref;
				sR.cigar=p_cigar;
			       	for(size_t k=0; k<sR.qual.size(); k++) sR.qual[k]=max_q_c;
			       	sR2.seq=a_read_2.seq_ref;
				sR2.cigar=p_cigar;
			       	for(size_t k=0; k<sR2.qual.size(); k++) sR2.qual[k]=max_q_c;
			       	if(a_read.is_plus_strand) sR2.reverse_comp();
				else sR.reverse_comp();
			
				//get new error free reads when the original reads having indels 
				if (sR.seq.length()!=read_len){
					sR.seq=a_art.ref_seq.substr(sR.pos-1, read_len);
				}
				if (sR2.seq.length()!=read_len){
					sR2.seq=a_art.ref_seq.substr(sR2.pos-1, read_len);
				}

			       	sR.printRead(SAMFILE_EF);
			       	sR2.printRead(SAMFILE_EF);
			}
		}

                t_num_read-=2;
            }
        }
        FQFILE2.close();
       	if (!no_ALN) ALNFILE2.close();
    }
    else{
        while(seq_reader.next_seq(id,a_art.ref_seq)){
	    std::replace(a_art.ref_seq.begin(), a_art.ref_seq.end(), 'U', 'T'); //replace U with T
            istringstream isID; isID.str(id); isID>>id; id=id.substr(0,len_ref_id); 
            num_seq++;
            a_art.ini_set(read_len);
            long t_num_read= 0;
	    if(cc_flag){
		     t_num_read= read_count;
	    }
	    else{
		    t_num_read= static_cast<unsigned long>(a_art.ref_seq.size()/read_len*x_fold);
		    if(amplicon){ t_num_read = (long) x_fold; }
	    }
            if(mask_n){ 
              a_art.mask_n_region(max_num_n);
            }
            while(t_num_read>0){
                osID<<id<<'-'<<suffID<<t_num_read;
                read_id = osID.str();
                osID.str("");
                a_read.clear();
		if(amplicon){
			a_art.next_ampread_indel(a_read);
		}
		else{ 
			a_art.next_read_indel(a_read);
		}
                if(mask_n){ 
                  if(a_read.is_plus_strand){
                    if(a_art.masked_pos.count(a_read.bpos)>0){
                      t_num_read-=1;
                      continue;
                    }
                  }
                  else{
                    size_t bpos=a_art.ref_seq.size()-a_read.bpos-read_len;
                    if(a_art.masked_pos.count(bpos)>0){
                      t_num_read-=1;
                      continue;
                    }
                  }
                }
 
                qual.clear();
		if(!sep_flag){
		       	qdist.get_read_qual(qual, read_len, true);
	       	}
	       	else{
		       	qdist.get_read_qual_1st(a_read.seq_read, qual);
	       	}
                a_read.add_error(qual);

                FQFILE<<"@"<<read_id<<endl<<a_read.seq_read<<endl<<"+"<<endl;
                for(size_t k=0; k<qual.size(); k++){
                    FQFILE<<(char)(qual[k]+33);
                }
                FQFILE<<endl;
	       	if(! a_read.get_aln(aln_read,aln_ref)){
		       	aln_ref=a_read.seq_ref;
		       	aln_read=a_read.seq_read;
	       	}
	       	if (!no_ALN) {
		       	ALNFILE<<">"<<id<<"\t"<<read_id<<"\t"<<a_read.bpos;
		       	if(a_read.is_plus_strand) ALNFILE<<"\t+\n";
		       	else ALNFILE<<"\t-\n";

		       	ALNFILE<<aln_ref<<endl<<aln_read<<endl;
	       	}

		if(sam_out){
		       	sR.qname=read_id;
		       	sR.rname=id;
		       	sR.flag =0;

			sR.seq=a_read.seq_read;
		       	sR.qual.resize(qual.size());
		       	for(size_t k=0; k<qual.size(); k++){
			       	sR.qual[k]=(char)(qual[k]+33);
		       	}

		       	if(a_read.is_plus_strand){
			       	sR.pos=a_read.bpos+1;
		       	}
		       	else{
			       	sR.flag = 0x10;
			       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+a_read.seq_read.size()-1);
			       	sR.reverse_comp();
		       	}
		       	sR.getCigar(aln_ref,aln_read,use_cigarM);
			sR.printRead(SAMFILE);

			if(err_free_sam){
			       	sR.seq=a_read.seq_ref;
			       	if(!a_read.is_plus_strand) sR.reverse_comp();
				sR.cigar=p_cigar;
			       	for(size_t k=0; k<sR.qual.size(); k++) sR.qual[k]=max_q_c;
				//get new error free reads when the original reads having indels 
				if (sR.seq.length()!=read_len){
					sR.seq=a_art.ref_seq.substr(sR.pos-1, read_len);
				}
			       	sR.printRead(SAMFILE_EF);
			}
		}


                t_num_read--;
            }
        }
    }

    FQFILE.close();
    if (!no_ALN) { ALNFILE.close(); }
    if(sam_out){ SAMFILE.close(); }
    if(err_free_sam){ SAMFILE_EF.close(); }

    if(!is_pairend_read){
	if(amplicon){
	       	cout << "              Amplicon 5'-end sequencing simulation" << endl << endl;
	}
	else{
	       	cout << "                  Single-end Simulation" << endl << endl;
	}
    } else if(is_matepair) {
	if(amplicon){
	       	cout << "             Amplicon matepair sequencing simulation" << endl << endl;
	}
	else{
	       	cout << "                  Matepair-end sequencing simulation" << endl << endl;
	}
    } else {
	if(amplicon){
	       	cout << "             Amplicon paired-end sequencing simulation" << endl << endl;
	}
	else{
	       	cout << "                  Paired-end sequencing simulation" << endl << endl;
	}
    }

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout<< "Total CPU time used: "<< cpu_time_used<<endl<<endl;
    cout<< "The random seed for the run: "<< rand_seed<<endl<<endl;

    if(!show_flag){

	cout << "Parameters used during run" << endl; 
	cout << "\tRead Length:\t" << read_len << endl;
	if(mask_n){
	       	cout << "\tGenome masking 'N' cutoff frequency: \t" <<max_num_n<<" in "<< read_len<<endl;
	}
	else{
	       	cout << "\t'N' genomic regions masking turned off"<< endl;
	}
       	if (amplicon){
	       	if(is_pairend_read)
		       	cout << "\t# Read Pairs per Amplion:  " << x_fold << endl;
	       	else 
			cout << "\t# Reads per Amplion:       " << x_fold << endl;
       	}
       	else{
	       	cout << "\tFold Coverage:            " << x_fold << "X" << endl;
       	}
	if(is_pairend_read && !amplicon){
	    cout << "\tMean Fragment Length:     " << mean << endl;
	    cout << "\tStandard Deviation:       " << std_dev << endl;
	}
	if(rate_flag){
            cout << "\tFirst Insertion Rate:     " << insRate << endl;
	    cout << "\tSecond Insertion Rate:    " << insRate2 << endl;
	    cout << "\tFirst Deletion Rate:      " << delRate << endl;
	    cout << "\tSecond Deletion Rate:     " << delRate2 << endl;
	}
	if(qs_flag){
	    cout << "\tFirst quality shift:      " << q_shift_up << endl;
	    cout << "\tSecond quality shift:     " << q_shift_up_2 << endl;
	}
	if(!sep_flag){
	    cout << "\tProfile Type:             Combined" << endl;
	} else {
	    cout << "\tProfile Type:             Separated" << endl;
	}
	cout << "\tID Tag:                   " << suffID.c_str() << endl << endl;;
	cout << "Quality Profile(s)" << endl; 

	if(is_pairend_read){
	    if(first_qual){
		cout << "\tFirst Read:   " << qual_file1.c_str() <<" (user's profile)" << endl;
	    } else if(!qdist.ssystem.empty()){
		cout << "\tFirst Read:   " <<qdist.ssystem<<" Length "<< profile_size <<" R1"<<" (built-in profile) "<<endl;
	    } else {
		cout << "\tFirst Read:   " <<" EMP" << profile_size <<"R1"<<" (built-in profile) "<<endl;
	    }
	    if(second_qual){
		cout << "\tSecond Read:  " << qual_file2.c_str()<<" (user's profile)"<< endl<< endl;
	    } else if(!qdist.ssystem.empty()){
		cout << "\tFirst Read:   " <<qdist.ssystem<<" Length "<< profile_size <<" R2"<<" (built-in profile) "<<endl<<endl;
	    } else {
		cout << "\tSecond Read:  " <<" EMP" << profile_size_2 <<"R2"<<" (built-in profile) "<<endl<<endl;
	    }
	} else {
	    if(first_qual){
		cout << "\t" << qual_file1.c_str()<<" (user's profile)"<< endl << endl;
	    } else if(!qdist.ssystem.empty()){
		cout << "\tFirst Read:   " <<qdist.ssystem<<" Length "<< profile_size <<" R1"<<" (built-in profile) "<<endl<<endl;
	    } else {
		cout << "\t " <<" EMP" << profile_size <<"R1"<<" (built-in profile) "<<endl<<endl;
	    }
	}
	
	cout << "Output files" << endl << endl;

	if(is_pairend_read){
	    cout << "  FASTQ Sequence Files:" << endl; 
	    cout << "\t the 1st reads: " << fqfile << endl;
	    cout << "\t the 2nd reads: " << fqfile2 << endl << endl;
	    if(!no_ALN){
		    cout << "  ALN Alignment Files:" << endl; 
		    cout << "\t the 1st reads: " << alnfasta << endl;
		    cout << "\t the 2nd reads: " << alnfasta2 << endl << endl;;
	    }
	} else {
	    cout << "  FASTQ Sequence File:" << endl; 
	    cout << "\t" << fqfile << endl << endl;
	    if (!no_ALN){
		    cout << "  ALN Alignment File:" << endl; 
		    cout << "\t" << alnfasta << endl << endl;
	    }
	}
	if(sam_out){
	    cout << "  SAM Alignment File:" << endl; 
	    cout << "\t" << samfile << endl << endl;
	}

	if(err_free_sam){
	    cout << "  Error Free SAM Alignment File:" << endl; 
	    cout << "\t" << samfile_ef << endl << endl;
	}
    }    

   return 0;
}

