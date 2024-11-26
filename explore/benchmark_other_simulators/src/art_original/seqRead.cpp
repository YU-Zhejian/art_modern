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

#include "seqRead.h"

int seqRead::get_indel(int read_len){
    indel.clear();
    //if(ins_rate.size()>=read_len) {cerr<<"fatal error\n";  exit(1)};
    int ins_len=0, del_len=0;
    //deletion
    for(int i=(int)del_rate.size()-1; i>=0; i--){
        if(del_rate[i]>=r_prob()){
            del_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor((read_len-1)*r_prob()); //invalid deletion positions: 0 or read_len-1
                if(pos==0) continue;
                if(indel.count(pos)==0){
                    indel[pos]='-';
                    j--;
                }
            }
            break;
        }
    }

    for(int i=ins_rate.size()-1; i>=0; i--){

	if((read_len-del_len-ins_len)<(i+1)) continue; //ensure that enough unchanged position for mutation

        if(ins_rate[i]>=r_prob()){
            ins_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor(r_prob()*read_len);
                if(indel.count(pos)==0){
                    short base=(short)ceil(r_prob()*4);
                    switch(base){
                                case 1:
                                    indel[pos]='A';   break;
                                case 2:
                                    indel[pos]='C';   break;
                                case 3:
                                    indel[pos]='G';   break;
                                case 4:
                                    indel[pos]='T';  
                    }
                    j--;
                }
            }
            break;
        }
    }
    return (ins_len-del_len);
};

//number of deletions <= number of insertions
int seqRead::get_indel_2(int read_len){
    indel.clear();
    //if(ins_rate.size()>=read_len) {cerr<<"fatal error\n";  exit(1)};
    int ins_len=0, del_len=0;

    for(int i=ins_rate.size()-1; i>=0; i--){
        if(ins_rate[i]>=r_prob()){
            ins_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor(r_prob()*read_len);
                if(indel.count(pos)==0){
                    short base=(short)ceil(r_prob()*4);
                    switch(base){
                                case 1:
                                    indel[pos]='A';   break;
                                case 2:
                                    indel[pos]='C';   break;
                                case 3:
                                    indel[pos]='G';   break;
                                case 4:
                                    indel[pos]='T';  
                    }
                    j--;
                }
            }
            break;
        }
    }

    //deletion
    for(int i=(int)del_rate.size()-1; i>=0; i--){
	if(del_len==ins_len) break;

	if((read_len-del_len-ins_len)<(i+1)) continue; //ensure that enough unchanged position for mutation

        if(del_rate[i]>=r_prob()){
            del_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor((read_len-1)*r_prob()); //invalid deletion positions: 0 or read_len-1
                if(pos==0) continue;
                if(indel.count(pos)==0){
                    indel[pos]='-';
                    j--;
                }
            }
            break;
        }
    }
    return (ins_len-del_len);
};


void seqRead::ref2read(){
    if(indel.size()==0){
        seq_read=seq_ref;
        return;
    }
    seq_read.clear();
    int k=0;
    for(size_t i=0; i<seq_ref.size();){
        //cout<<i<<"\t"<<k<<endl;
        if(indel.count(k)==0){
            seq_read.push_back(seq_ref[i]); i++; k++; 
        }
        else if(indel[k]=='-'){
            i++;k++;
        }
        else{
            seq_read.push_back(indel[k]); k++;
        }
    }
    while(indel.count(k)>0){
        seq_read.push_back(indel[k]);
        k++;
    }
}

