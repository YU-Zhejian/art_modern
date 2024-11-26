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

#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <cstdio>
#include <set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_sys.h>
//#include <gsl/gsl_machine.h>
//#include <gsl/gsl_precision.h>
//#include <gsl/gsl_nan.h>
//#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_cdf.h>

#include "empdist.h"

using namespace std;

class seqRead{
public:
    seqRead(){
	    int len=35;
//	    double error_1st []={0.0021933, 0.0022350, 0.0022800, 0.0023540, 0.0023960, 0.0024700, 0.0024540, 0.0024000, 0.0022900, 0.0023000, 0.0029220, 0.0049160, 0.0070740, 0.0091220, 0.0119580, 0.0166400, 0.0223840, 0.0244940, 0.0266900, 0.0307500, 0.0325920, 0.0349900, 0.0380880, 0.0430380, 0.0470620, 0.0481820, 0.0505380, 0.0562680, 0.0572500, 0.0564340, 0.0599940, 0.0614700, 0.0733940, 0.0785200, 0.0870333};
	    //qual_1st=27, 27, 26, 26, 26, 26, 26, 26, 26, 26, 25, 23, 22, 20, 19, 18, 17, 16, 16, 15, 15, 15, 14, 14, 13, 13, 13, 12, 12, 12, 12, 12, 11, 11, 11
//	    double error_2nd[]={0.00763, 0.0077725, 0.007898, 0.00779, 0.007644, 0.008728, 0.008476, 0.01096, 0.010886, 0.011398, 0.01085, 0.011326, 0.01012, 0.011252, 0.011654, 0.012712, 0.013806, 0.013734, 0.016438, 0.035104, 0.035396, 0.037918, 0.039818, 0.041352, 0.027216, 0.04486, 0.047014, 0.054538, 0.067542, 0.092088, 0.092308, 0.10725, 0.113132, 0.116895, 0.103563333};
            //qual_2nd=21, 21, 21, 21, 21, 21, 21, 20, 20, 19, 20, 19, 20, 19, 19, 19, 19, 19, 18, 15, 15, 14, 14, 14, 16, 13, 13, 13, 12, 10, 10, 10, 9, 9, 10

            double error_1st []={0.001195, 0.001195, 0.001235, 0.0013425, 0.00146, 0.0014975, 0.001485, 0.0014525, 0.0013275, 0.0011675, 0.0011675, 0.00215, 0.0055125, 0.0097875, 0.011615, 0.01253, 0.0176525, 0.0279425, 0.0326375, 0.026815, 0.0257325, 0.03285, 0.04062, 0.0448075, 0.043745, 0.04701, 0.0504975, 0.051395, 0.05926, 0.0622425, 0.0556075, 0.054135, 0.06112, 0.082305, 0.1124975};
           //qual_1st=29, 29, 29, 29, 28, 28, 28, 28, 29, 29, 29, 27, 23, 20, 19, 19, 18, 16, 15, 16, 16, 15, 14, 13, 14, 13, 13, 13, 12, 12, 13, 13, 12, 11, 9
            double error_2nd[]={0.006395, 0.00664, 0.0070175, 0.007225, 0.006895, 0.006125, 0.007495, 0.00923, 0.011055, 0.0126975, 0.009685, 0.0079975, 0.009095, 0.0100525, 0.0116925, 0.01151, 0.01174, 0.01361, 0.01402, 0.01626, 0.0415125, 0.06215, 0.040825, 0.023065, 0.02608, 0.0300525, 0.0507825, 0.0699575, 0.05991, 0.0638375, 0.10265, 0.1286775, 0.1198675, 0.1104025, 0.11194};
           //qual_2nd=22, 22, 22, 21, 22, 22, 21, 20, 20, 19, 20, 21, 20, 20, 19, 19, 19, 19, 19, 18, 14, 12, 14, 16, 16, 15, 13, 12, 12, 12, 10, 9, 9, 10, 10

	    cal_err_rate_1st.assign(error_1st,error_1st+len);
	    cal_err_rate_2nd.assign(error_2nd,error_2nd+len);
    }

    vector <double> cal_qual_1st;
    vector <double> cal_qual_2nd;
    vector<gsl_rng*> gsl_p_1st;
    vector<gsl_rng*> gsl_p_2nd;
    void ini_ran_qual(){
      gsl_rng_default_seed=(unsigned int)time(NULL);
      const gsl_rng_type *rndT=gsl_rng_default;
      for(size_t k=0; k<cal_err_rate_1st.size(); k++){
        cal_qual_1st.push_back(-10*log10(cal_err_rate_1st[k]));
        cal_qual_2nd.push_back(-10*log10(cal_err_rate_2nd[k]));
        gsl_p_1st.push_back(gsl_rng_alloc(rndT));
        gsl_p_2nd.push_back(gsl_rng_alloc(rndT));
      }
    };

    static const int max_N_qual=10;
    vector <double> cal_err_rate_1st;
    vector <double> cal_err_rate_2nd;
    vector <double> ins_rate; //Bionomial cum_prob gsl_cdf_binomial_Q(unsigned int k, double p, unsigned int n) 
    vector <double> del_rate; //Binomial
    vector<double> sub_rate; //Binomial
    void set_rate(int read_len, double p, int max_num, vector <double>& rate){
        rate.resize(max_num);
	if(max_num>read_len) max_num=read_len;
        for(size_t i=1; i<=max_num; i++){
            rate[i-1]= gsl_cdf_binomial_Q(i, p, read_len);
        }
    };

    //when max_num =-1, no limit on the number of indels 
    //the maxium number of indels is set by cdf_cutoff to save computation time 
    void set_rate(int read_len, double p, vector <double>& rate, int max_num=-1, double cdf_cutoff=0.999999){
        rate.clear();
       	if(max_num==0) return;
	//p ==0 no error
	if(p<0.000000000000000000000000000001) return; //when rate < 10^-30, set it 0
	double tp=gsl_cdf_binomial_Q(0, p, read_len);
	double p_cdf=tp;
        for(size_t i=1;i<read_len;i++){
	    tp=gsl_cdf_binomial_Q(i, p, read_len);
            rate.push_back(tp);
	    if(max_num>0 && (i>=max_num)) break;
	    p_cdf+=tp; 
	    if(p_cdf>=cdf_cutoff) break;
        }
    };

    static char rand_base(){
        short base=(short)ceil(r_prob()*4);
        switch(base){
            case 1:
                return 'A';
            case 2:
                return 'C';
            case 3:
                return 'G';
            case 4:
                return 'T';  
	    default:
                return 'N';  
        }
    };

    //static bool with_indel;
    int get_indel(int read_len);
    //number of deletions <= number of insertions
    int get_indel_2(int read_len);
    map<int,char,less<int> > indel;
    map<int,char> substitution;
    bool is_plus_strand;
    unsigned long bpos; //parent
    string seq_read;
    string seq_ref;
    void clear(){
        indel.clear();
        substitution.clear();
        seq_read.clear();
        seq_ref.clear();
    }
    //string aln_read;
    //string aln_ref;
    bool get_aln(string& aln_read, string& aln_ref){
        if(indel.size()==0) return false;
        map<int,char,less<int> >::iterator it;
        aln_read=seq_read; aln_ref=seq_ref;
        for(it=indel.begin(); it!=indel.end(); it++){
            if(it->second!='-'){
                aln_ref.insert(it->first,1,'-');
            }
            else{
                aln_read.insert(it->first,1,'-');
            }
        }
        return true;
    };

    void ref2read();
    
    //based on based on calibrated position-depended error rates
    int add_calib_error_1st(){
        int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_1st[i]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

   int add_calib_error_2nd(){
        int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_2nd[i]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

    //based on calibrated position-depended error rates
    int add_calib_error_1st(vector<short>&qual){
        int num=0;
        if(qual.size()<seq_read.size()) qual.resize(seq_read.size());
        for(size_t i=0; i<seq_read.size(); i++){
          int q=gsl_ran_poisson(gsl_p_1st[i],cal_qual_1st[i]);
          qual[i]=q;
          if(r_prob()<empdist::prob_err[q]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

   int add_calib_error_2nd(vector<short>&qual){
        int num=0;
        if(qual.size()<seq_read.size()) qual.resize(seq_read.size());
        for(size_t i=0; i<seq_read.size(); i++){
          int q=gsl_ran_poisson(gsl_p_2nd[i],cal_qual_2nd[i]);
          qual[i]=q;
          if(r_prob()<empdist::prob_err[q]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };


    int add_calib_error_1st_stat(vector<short>&err_pos){
        int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_1st[i]){
              err_pos[i]=1;
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

   int add_calib_error_2nd_stat(vector<short>&err_pos){
        int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_2nd[i]){
              err_pos[i]=1;
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };


    //based on empirical dist of quali scores
    int add_error(vector<short>&qual){
        if(qual.size()!=seq_read.size()){
            cerr<<"Error: the number of bases is not equal to the number of quality scores!\n";
	    cerr<< "qual size: " << qual.size() << ",  read len: " << seq_read.size() << endl;
	    exit(1);
            return 0;
        }
        int num=0;
        for(size_t i=0; i<qual.size(); i++){
	    if(seq_read[i]=='N'){
		    //qual[i]=(short) floor(max_N_qual*r_prob()); //reassign N base quality score < max_N_qual   
		    qual[i]=(short)1;
		    continue;
	    }
            if(r_prob()<empdist::prob_err[qual[i]]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

    //based on empirical dist of quali scores
    int add_error_stat(vector<short>&qual, vector<long>&err_pos){
        if(qual.size()!=seq_read.size()){
            cerr<<qual.size()<<"\t"<<seq_read.size()<<endl;
            cerr<<"Error: the number of bases is not equal to the number of quality score!\n";
	    exit(1);
        }
        int num=0;
        for(size_t i=0; i<qual.size(); i++){
            if(r_prob()<empdist::prob_err[qual[i]]){
              err_pos[i]+=1;
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
    };

    //base on bionomial error substitution rate
    int add_error(int read_len){
        int sub_num=0;
        for(int i=(int)sub_rate.size()-1; i>=0; i--){
            if(sub_rate[i]>=r_prob()){
                sub_num=i+1;
                for(int j=i; j>=0;){
                    int pos=(int) floor(read_len*r_prob());
                    if(substitution.count(pos)==0){
                        char achar=seq_read[pos];
                        while(seq_read[pos]==achar){ achar=rand_base(); }
                        substitution[pos]=achar;
                        seq_read[pos]=achar;
                        j--;
                    }
                }
                break;
            }
        }
        return sub_num;
    };

    //vector<char> aln_read;
    //vector<char> aln_ref;
};
