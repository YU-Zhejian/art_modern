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
#pragma once
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_pow_int.h>

#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <cstdio>
#include <set>


#include "seqRead.h"


using namespace std;

class art{
public:
    //art();
    string ref_seq;
    string ref_seq_cmp;
    int read_len;
    long valid_region;
    void ini_set(int readLen){
        read_len=readLen;
        valid_region=ref_seq.size()-readLen;
        ref_seq_cmp.resize(ref_seq.size());
        size_t size=ref_seq.size();
        for(size_t i=0; i<size; i++){
//            ref_seq[i]=tolower(ref_seq[i]);
            ref_seq[i]=toupper(ref_seq[i]);
            size_t k=size-i-1;
            switch(ref_seq[i]){
                case 'A':
                    ref_seq_cmp[k]='T'; break;
                case 'C':
                    ref_seq_cmp[k]='G'; break;
                case 'G':
                    ref_seq_cmp[k]='C'; break;
                case 'T':
                    ref_seq_cmp[k]='A'; break;
                default:
                    ref_seq_cmp[k]='N';
            }
        }
    }
    bool next_read_indel(seqRead& a_read);
    bool next_pair_read_indel(seqRead& read_1, seqRead& read_2);
    bool next_pair_read_indel_cmp(seqRead& read_1, seqRead& read_2);
    bool next_pair_read_indel_mate(seqRead& read_1, seqRead& read_2);
    bool next_read(seqRead& a_read);

    bool next_ampread_indel(seqRead& a_read);
    bool next_pair_ampread_indel_cmp(seqRead& read_1, seqRead& read_2);
    bool next_matepair_ampread_indel_cmp(seqRead& read_1, seqRead& read_2);

    set <size_t> masked_pos;
    //region with n
    void mask_n_region(int max_num_n){
      masked_pos.clear(); //reset mask position
      if(ref_seq.size()<read_len) return;
      int num_n=0;
      for(int i=0; i<read_len-1; i++){
        if(ref_seq[i]=='N') num_n++;
      }
      for(size_t i=0; i<=ref_seq.size()-read_len; i++){
        if(ref_seq[i+read_len-1]=='N') num_n++;
        if(num_n>=max_num_n){ masked_pos.insert(i); 
        }
        if(ref_seq[i]=='N') num_n--;
      }
    }
    
    static gsl_rng* gsl_R;
    static int gaussain_mean;
    static double gaussain_sigma;

    static void ini_read_pair_rand(int mean, double sigma){
        gsl_rng_default_seed=(unsigned int)time(NULL);
        const gsl_rng_type *rndT=gsl_rng_default;
        gsl_R=gsl_rng_alloc(rndT);
        gaussain_mean=mean;
        gaussain_sigma=sigma;
    }; 

    static void ini_read_pair_rand(int mean, double sigma, unsigned int gseed){
        gsl_rng_default_seed=gseed;
        const gsl_rng_type *rndT=gsl_rng_default;
        gsl_R=gsl_rng_alloc(rndT);
        gaussain_mean=mean;
        gaussain_sigma=sigma;
    }; 


  
    static gsl_rng* gsl_g;
    static gsl_rng* gsl_p;
    static void ini_bias_model(){
        gsl_rng_default_seed=(unsigned int)time(NULL);
        const gsl_rng_type *rndT=gsl_rng_default;
        gsl_g=gsl_rng_alloc(rndT);
        const gsl_rng_type *rndT2=gsl_rng_default;
        gsl_p=gsl_rng_alloc(rndT2);
    }; 
 
    double ave_depth;
    double para_x1; 

    int get_depth_cov(){ 
      double depth_rate=ave_depth+para_x1*ave_depth*gsl_ran_gaussian(gsl_g,1);
      if(depth_rate<0){ depth_rate=0; return 0; }
      return gsl_ran_poisson (gsl_p, depth_rate);
    };

    static void free_gsl(){
      gsl_rng_free(gsl_R);
      gsl_rng_free(gsl_g);
      gsl_rng_free(gsl_p);
    };
    

    bool next_read_indel_bias(seqRead& a_read){
      long pos=a_read.bpos;//(long) floor(r_prob()*valid_region); //pos in [0 ..len-1]   
      int slen=a_read.get_indel(read_len);            
      a_read.is_plus_strand=true;
      if(pos+read_len>ref_seq.length()){ 
        a_read.is_plus_strand=false;
        pos=ref_seq.length()-pos-1;
      }
      else if(r_prob()>0.5){
        a_read.is_plus_strand=false;
        pos=ref_seq.length()-pos-read_len;
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
    };


    bool next_pair_read_indel_bias(seqRead& read_1, seqRead& read_2){
      int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
      while (fragment_len<read_len || fragment_len>ref_seq.length()){
        fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
      } 
      long pos_1=read_1.bpos; //(long) floor((ref_seq.length()-fragment_len)*r_prob());
      long pos_2=pos_1+fragment_len-read_len;
      bool is_plus_strand=true;
      int slen_1 =read_1.get_indel(read_len);
      int slen_2 =read_2.get_indel(read_len);   

      if(pos_1+fragment_len>ref_seq.length()){ 
        pos_2=pos_1-fragment_len+read_len;
        pos_1=ref_seq.length()-pos_1-1;
        pos_2=ref_seq.length()-pos_2-1;
        is_plus_strand=false;
      }
      else if(r_prob()>0.5){
        is_plus_strand=false;
        pos_1=ref_seq.length()-pos_1-read_len;
        pos_2=ref_seq.length()-pos_2-read_len;
        long tm_pos=pos_1;
        pos_1=pos_2;
        pos_2=tm_pos;
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
      //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
      return true;
    };

};
//
//class ngrs{
//private:
//    static const unsigned int LargeNumber=10000;
//    unsigned int max_branch;
//    vector < node > tree;
//    unsigned int depth;
//    unsigned int time_lambda; //the mean time of between special event in thousand year
//    unsigned int branch_lambda;
//    double real_aveRate;
//    map <int, string> subtree;
//public:
//    ngrs(int level, unsigned int tm, double rate, unsigned int branch);
//    ngrs(int level, unsigned int tm, double rate, unsigned int mean_bch, unsigned int max_bch); 
//    string getTree();
//    string getSubTree(int root);
//    size_t node_count(){ return tree.size(); };
//    //~ngrs();
//    //void print(ostream& OUT);
//    //void addLevel(int );
//};


