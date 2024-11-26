/*
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ART -- Artificial Read Transcription, Illumina Q version 
 * Authors: Weichun Huang 2008-2016
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
#include <time.h>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include "art.h"
#include "empdist.h"
#include "readSeqFile.h"


using namespace std;

bool art::next_read_indel(seqRead& a_read){
    long pos=(long) floor(r_prob()*valid_region); //pos in [0 ..len-1]   
    int slen =a_read.get_indel(read_len);            
//ensure get a fixed read length 
    if((pos+read_len-slen)>ref_seq.length()){ 
	    slen =a_read.get_indel_2(read_len);
    } 
    a_read.is_plus_strand=true;
    if(r_prob()>0.5){
        a_read.is_plus_strand=false;
    }
    if(a_read.is_plus_strand){
        a_read.seq_ref=ref_seq.substr(pos, read_len-slen);
    }
    else{
        a_read.seq_ref=ref_seq_cmp.substr(pos, read_len-slen);
    }
    a_read.bpos=pos;
    a_read.ref2read();
    return true;
}

bool art::next_read(seqRead& a_read){
    long pos=(long) floor(r_prob()*valid_region); //pos in [0 ..len-1]           
    //string seq_ref;
    a_read.is_plus_strand=true;
    if(r_prob()>0.5){
        a_read.is_plus_strand=false;
    }
    if(a_read.is_plus_strand){
        a_read.seq_ref=ref_seq.substr(pos, read_len);
    }
    else{
        a_read.seq_ref=ref_seq_cmp.substr(pos, read_len);
    }
    a_read.bpos=pos;
    a_read.seq_read=a_read.seq_ref;
    return true;
}

//on the same strand of two reads
bool art::next_pair_read_indel(seqRead& read_1, seqRead& read_2){
    int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    while (fragment_len<read_len || fragment_len>ref_seq.length()){
        fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    }
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob());
    long pos_2=pos_1+fragment_len-read_len;
    int slen_1 =read_1.get_indel(read_len);
    int slen_2 =read_2.get_indel(read_len);   

//ensure get a fixed read length 
    if((pos_1+read_len-slen_1)>ref_seq.length()){ 
	    slen_1 =read_1.get_indel_2(read_len);
    } 
    if((pos_2+read_len-slen_2)>ref_seq.length()){ 
	    slen_2 =read_2.get_indel_2(read_len);
    }

    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.seq_ref=ref_seq.substr(pos_1, read_len-slen_1);
        read_2.seq_ref=ref_seq.substr(pos_2, read_len-slen_2);
    }
    else{
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len-slen_1);
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len-slen_2);
    }
    read_1.is_plus_strand=is_plus_strand;
    read_1.bpos=pos_1;
    read_1.ref2read();
    read_2.is_plus_strand=is_plus_strand;
    read_2.bpos=pos_2;
    read_2.ref2read();
    return true;
}

//matepair-end read: the second read is reverse complemenaty strand 
bool art::next_pair_read_indel_mate(seqRead& read_1, seqRead& read_2){
    if(read_len>ref_seq.length()){
	    return false; //ref_seq is too short.
    }
    int fragment_len=1;
    if(gaussain_mean-2*gaussain_sigma>ref_seq.length()){
	    //when reference length < mean-2*std, fragment_len sets to be reference length
	   fragment_len=ref_seq.length();
    }
    else{
	    fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
	    while (fragment_len<read_len || fragment_len>ref_seq.length()){
		    fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
	    }
    }

    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob())+fragment_len-read_len;
    long pos_2=ref_seq.length()-(pos_1+2*read_len-fragment_len);
    int slen_1 =read_1.get_indel(read_len);
    int slen_2 =read_2.get_indel(read_len);   

//ensure get a fixed read length 
    if((pos_1+read_len-slen_1)>ref_seq.length()){ 
	    slen_1 =read_1.get_indel_2(read_len);
    } 
    if((pos_2+read_len-slen_2)>ref_seq.length()){ 
	    slen_2 =read_2.get_indel_2(read_len);
    }

    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.is_plus_strand=true;
        read_1.seq_ref=ref_seq.substr(pos_1, read_len-slen_1); 
        read_2.is_plus_strand=false;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len-slen_2);
    }
    else{
        read_1.is_plus_strand=false;
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len-slen_1);
        read_2.is_plus_strand=true;
        read_2.seq_ref=ref_seq.substr(pos_2, read_len-slen_2);
    }
    read_1.bpos=pos_1;
    read_1.ref2read();
    read_2.bpos=pos_2;
    read_2.ref2read();
    return true;
}


//paired-end read: the second read is reverse complemenaty strand 
bool art::next_pair_read_indel_cmp(seqRead& read_1, seqRead& read_2){
    if(read_len>ref_seq.length()){
	    return false; //ref_seq is too short.
    }
    int fragment_len=1;
    if(gaussain_mean-2*gaussain_sigma>ref_seq.length()){
	    //when reference length < mean-2*std, fragment_len sets to be reference length
	   fragment_len=ref_seq.length();
    }
    else{
	    fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
	    while (fragment_len<read_len || fragment_len>ref_seq.length()){
		    fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
	    }
    }
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob());
    //long pos_2=pos_1+fragment_len-read_len;
    long pos_2=ref_seq.length()-pos_1-fragment_len;
    int slen_1 =read_1.get_indel(read_len);
    int slen_2 =read_2.get_indel(read_len);   

    //ensure get a fixed read length 
    if((pos_1+read_len-slen_1)>ref_seq.length()){ 
	    slen_1 =read_1.get_indel_2(read_len);
    } 
    if((pos_2+read_len-slen_2)>ref_seq.length()){ 
	    slen_2 =read_2.get_indel_2(read_len);
    }

    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.is_plus_strand=true;
        read_1.seq_ref=ref_seq.substr(pos_1, read_len-slen_1); 
        read_2.is_plus_strand=false;
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len-slen_2);
    }
    else{
        read_1.is_plus_strand=false;
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len-slen_1);
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.is_plus_strand=true;
        read_2.seq_ref=ref_seq.substr(pos_2, read_len-slen_2);
    }
    read_1.bpos=pos_1;
    read_1.ref2read();
    read_2.bpos=pos_2;
    read_2.ref2read();
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
    return true;
}


//for amplicon sequencing simulation
// amplicon single end sequencing
bool art::next_ampread_indel(seqRead& a_read){
    if(read_len>ref_seq.length()){
	    return false; //ref_seq is too short.
    }
    long pos=(long) 0; 
    int slen =0; 

    if(read_len==ref_seq.length())
	    slen =a_read.get_indel_2(read_len);            
    else
	    slen =a_read.get_indel(read_len);            
    
    a_read.is_plus_strand=true;
//    if(r_prob()>0.5){
//        a_read.is_plus_strand=false;
//    }
//    if(a_read.is_plus_strand){
    if(slen>=0)
	    a_read.seq_ref=ref_seq.substr(pos, read_len-slen);
    else{
	    a_read.seq_ref=ref_seq.substr(pos, read_len-slen);
	    if(a_read.seq_ref.length()<(read_len-slen)) 
	    slen =a_read.get_indel_2(read_len);            
	    a_read.seq_ref=ref_seq.substr(pos, read_len-slen);
    }
//    }
//    else{
//       a_read.seq_ref=ref_seq_cmp.substr(pos, read_len-slen);
//    }
    a_read.bpos=pos;
    a_read.ref2read();
    return true;
}

//amplicon paired-end reads: the second read is reverse complemenaty strand 
bool art::next_pair_ampread_indel_cmp(seqRead& read_1, seqRead& read_2){
    if(read_len>ref_seq.length()){
	    return false; //ref_seq is too short.
    }
    long pos_1=(long) 0;
    long pos_2=pos_1;
    int slen_1 =read_1.get_indel(read_len);
    int slen_2 =read_2.get_indel(read_len);   
    bool is_plus_strand=true;

    //ensure get a fixed read length 
    if((read_len-slen_1)>ref_seq.length()){ 
	    slen_1 =read_1.get_indel_2(read_len);
    } 
    if((read_len-slen_2)>ref_seq.length()){ 
	    slen_2 =read_2.get_indel_2(read_len);
    }

//    if(r_prob()>0.5){
//        is_plus_strand=false;
//    }
    if(is_plus_strand){
        read_1.is_plus_strand=true;
	
        read_1.seq_ref=ref_seq.substr(pos_1, read_len-slen_1); 
        read_2.is_plus_strand=false;
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len-slen_2);
    }
    else{
        read_1.is_plus_strand=false;
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len-slen_1);
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.is_plus_strand=true;
        read_2.seq_ref=ref_seq.substr(pos_2, read_len-slen_2);
    }
    read_1.bpos=pos_1;
    read_1.ref2read();
    read_2.bpos=pos_2;
    read_2.ref2read();
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
    return true;
}

//amplicon matepaired reads: the second read is reverse complemenaty strand 
bool art::next_matepair_ampread_indel_cmp(seqRead& read_1, seqRead& read_2){
    if(read_len>ref_seq.length()){
	    return false; //ref_seq is too short.
    }
    long pos=(long) ref_seq.length()-read_len;

    int slen_1 =read_1.get_indel(read_len);
    int slen_2 =read_2.get_indel(read_len);   

    long pos_1=pos+slen_1;
    long pos_2=pos+slen_2;
    //ensure no negative position 
    if(pos_1<0 || pos_2 <0){ 
	    slen_1 =read_1.get_indel_2(read_len);
	    slen_2 =read_2.get_indel_2(read_len);   
	    pos_1=pos+slen_1;
	    pos_2=pos+slen_2;
    }

    bool is_plus_strand=true;
//    if(r_prob()>0.5){
//        is_plus_strand=false;
//    }
    if(is_plus_strand){
        read_1.is_plus_strand=true;
        read_1.seq_ref=ref_seq.substr(pos_1, read_len-slen_1); 
        read_2.is_plus_strand=false;
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len-slen_2);
    }
    else{
        read_1.is_plus_strand=false;
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len-slen_1);
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.is_plus_strand=true;
        read_2.seq_ref=ref_seq.substr(pos_2, read_len-slen_2);
    }
    read_1.bpos=pos_1;
    read_1.ref2read();
    read_2.bpos=pos_2;
    read_2.ref2read();
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
    return true;
}





//bool parse_arg(int num, char* arg){
//    bool success=true;
//    int i=1;
//    for(;i<ARGC;++i){
//       	char* pch = ARGV[i];
//       	if( *pch != '-' || *(pch+1) == '\0') break;
//       	while(*++pch && success){
//            switch(*pch){
//                // Version Information
//            case 'v':
//            case 'V':
//                prtVersion = true;
//                break;
//                // Help Information
//            case 'h':
//            case 'H':
//                prtVersion = true;
//                prtUsage = true;
//                break;
//                // Unbuffered Output
//            case 't':
//            case 'T':
//                html_out=false;
//                break;
//            case 'b':
//            case 'B':
//                if(i<ARGC) browser=ARGV[++i]; 
//                else success=PrintErr("Invalid value for option: \"%c\"", *pch);
//                break;
//            case 'o':
//            case 'O':
//                if(i<ARGC) outFile=ARGV[++i]; 
//                else success=PrintErr("Invalid value for option: \"%c\"", *pch);
//                break;
//            case 'c':
//            case 'C':
//                if(i<ARGC) cfgFile=ARGV[++i]; 
//                else success=PrintErr("Invalid value for option: \"%c\"", *pch);
//                break;
//                // Error Reporting
//            default:
//                success=PrintErr("Unreconized switch, \"%c\"", *pch);
//                prtUsage = true;
//                prtVersion = true;
//                break;
//            }
//       	}
//    }
//
//    // Print Version and/or Usage Information and exit
//    if(prtVersion || prtUsage){
//       	if(prtUsage) printUsage();
//       	if(prtVersion) printVer();
//       	return false;
//    }
//
//    if(i>ARGC){ printUsage(); return false;}
//
//    aceFile=ARGV[i];
//    return success;;
//}

//
//ngrs::ngrs(int level, unsigned int tm, double rate, unsigned int branch){
//	const gsl_rng_type *T1, *T2;
//	gsl_rng *r1, *r2;
//	gsl_rng_default_seed=(unsigned int) time(NULL);
////	gsl_rng_env_setup();
//	T1 = gsl_rng_default;
//	T2 = gsl_rng_default;
//	r1 = gsl_rng_alloc(T1);
//	r2 = gsl_rng_alloc(T2);
//
//	depth = level;
//	time_lambda= tm;
//	branch_lambda = branch;
//	max_branch = branch;;
//	real_aveRate = rate;
//
//	node root(0,1,0,0);
//	tree.push_back(root);
//	size_t i = 0;
//	for (; i < tree.size(); i++) {
//		if(tree[i].level>=depth) break;
//		 node newNode=tree[i];
//		 newNode.parent = i;
//		 newNode.level += 1;
//		 double timepass = gsl_ran_exponential(r1, time_lambda);
//		 newNode.time = timepass; 
//		 //gsl_ran_gaussian(const gsl_rng * r, double sigma);
//		 for(size_t num=0; num<branch; num++){
//			 double substitution = gsl_ran_poisson(r2, real_aveRate);
//			 newNode.rate = substitution*timepass/time_lambda/LargeNumber;
////			 cout<<substitution<<"\t"<<timepass<<"\t"<<time_lambda<<endl;
//			 tree.push_back(newNode);
//		 }
//	}
//	gsl_rng_free (r1);
//	gsl_rng_free (r2);
//}
//
//ngrs::ngrs(int level, unsigned int tm, double rate, unsigned int mean_bch, unsigned int max_bch) {
//	const gsl_rng_type *T1, *T2, *T3;
//	gsl_rng *r1, *r2, *r3;
//	gsl_rng_default_seed=(unsigned int)time(NULL);
////	gsl_rng_env_setup();
//
//	T1 = gsl_rng_default;
//	T2 = gsl_rng_default;
//	T3 = gsl_rng_default;
//	r1 = gsl_rng_alloc(T1);
//	r2 = gsl_rng_alloc(T2);
//	r3 = gsl_rng_alloc(T3);
//
//	depth = level;
//	time_lambda= tm;
//	branch_lambda = mean_bch;
//	max_branch = max_bch;
//	real_aveRate = rate;
//
//	size_t k = gsl_ran_poisson(r3, branch_lambda);
//	double timepass = gsl_ran_exponential(r1, time_lambda);
//	double substitution = gsl_ran_poisson(r2, real_aveRate);
//
//	node root(0,1,0,0);
//	tree.push_back(root);
//	size_t i = 0;
//	for (; i < tree.size(); i++) {
//		if(tree[i].level>=depth) break;
//		 size_t k = gsl_ran_poisson(r3, branch_lambda);
//		 if(k==0) continue;
//		 if(k>max_branch) k = max_branch;
//		 node newNode=tree[i];
//		 newNode.parent = i;
//		 newNode.level += 1;
//		 double timepass = gsl_ran_exponential(r1, time_lambda);
//		 newNode.time = timepass; 
//		 //gsl_ran_gaussian(const gsl_rng * r, double sigma);
//		 for(size_t num=0; num<k; num++){
//			 double substitution = gsl_ran_poisson(r2, real_aveRate);
//			 newNode.rate = substitution*timepass/time_lambda/LargeNumber;
//			 tree.push_back(newNode);
//		 }
//	}
//	gsl_rng_free (r1);
//	gsl_rng_free (r2);
//	gsl_rng_free (r3);
//}

//bool art::next_read_indel(seqRead& a_read){
//    long pos=(long) floor(r_prob()*valid_region); //pos in [0 ..len-1]
//    map<int,char,less<int> > indel;    
//    int slen =a_read.get_indel(read_len, indel);            
//    string seq_ref;
//    a_read.is_plus_strand=true;
//    if(r_prob()>0.5){
//        a_read.is_plus_strand=false;
//    }
//    if(is_plus_strand){
//        seq_ref=ref_seq.substr(pos, read_len-slen);
//    }
//    else{
//        seq_ref=ref_seq_cmp.substr(pos, read_len-slen);
//    }
//
//    a_read.bpos=pos;
//
//    //int del_len=(indel.size()-slen)/2;
//    int ins_len=(indel.size()+slen)/2;
//    int aln_len=seq_ref.size()+ins_len;
//
//    string seq_read; seq_read.resize(read_len); 
//    a_read.aln_ref.resize(aln_len);
//    a_read.aln_read.resize(aln_len);
//
//    for(int i=0,j=0, k=0; i<aln_len; i++){
//        if(indel.count(i)==0){
//            a_read.aln_read[i]=seq_ref[k];
//            a_read.aln_ref[i]=seq_ref[k];
//            seq_read[j]=seq_ref[k]; j++; k++; 
//        }
//        else if(indel[i]!='-'){
//            a_read.aln_ref[i]='-';
//            a_read.aln_read[i]=indel[i];
//            seq_read[j]=indel[i]; j++;
//        }
//        else {
//            a_read.aln_ref[i]=seq_ref[k]; k++;
//            a_read.aln_read[i]='-';
//        }
//    }
//
//}
