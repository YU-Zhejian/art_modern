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

#include <sstream>
#include "samRead.h"


void samHeader::getRefseqID(char* seqfile){
       	SN.clear();
       	LN.clear();
       	readSeqFile seq_reader(seqfile); 
	string id, seq;
       	while(seq_reader.next_seq(id,seq)){ 
		SN.push_back(id);
	       	LN.push_back(seq.size());
	       	seq.clear();
       	}
//       	seq_reader.~readSeqFile();
}

void samHeader::printHeader(ostream& fout){
       	if(SN.size()!=LN.size()){
	       	cerr<<"fatal error in determining the number of sequences"<<endl;
	       	exit(1);
       	}
       	fout<<"@HD\t"<<"VN:"<<VN<<"\tSO:"<<SO<<endl;
       	for(size_t i=0; i<SN.size(); i++){
	       	fout<<"@SQ\t"<<"SN:"<<SN[i]<<"\tLN:"<<LN[i]<<endl;
       	}
       	fout<<"@PG\t"<<"ID:"<<ID<<"\tPN:"<<PN<<"\tCL:"<<CL<<endl;
}

void samHeader::printAlnHeader(ostream& fout){
       	if(SN.size()!=LN.size()){
	       	cerr<<"fatal error in determining the number of sequences"<<endl;
	       	exit(1);
       	}
       	fout<<"@CM\t"<<CL<<endl;
       	for(size_t i=0; i<SN.size(); i++){
	       	fout<<"@SQ\t"<<SN[i]<<"\t"<<LN[i]<<endl;
       	}
       	fout<<"##Header End"<<endl;
}

//make sure flag set before call getCigar
void samRead::getCigar(string & aln_ref, string& aln_read, bool use_M){
	if(aln_ref.length()!=aln_read.length()){
		cerr<<"fatal error: wrong alignment"<<endl;
		exit(1);
	}
//	vector<char> cType; 
//	vector<int> len; 
	cigar="";
	
	char t, t2;
	int k=0;
	if(flag & 0x10){ //reverse complement SEQ 
		int ins_len=0, del_len=0;
		for(int i=aln_ref.length()-1; i>=0; i--){
		       	if(aln_ref[i]==aln_read[i]){
			       	t='=';
				if (use_M) t='M';
		       	}
		       	else if(aln_ref[i]=='-'){
			       	t='I';
				ins_len++;
		       	}
		       	else if(aln_read[i]=='-'){
			       	t='D';
				del_len++;
		       	}
		       	else{
			       	t='X';
				if (use_M) t='M';
		       	}
		       	if(t!=t2 && k>0){
			       	//			cType.push_back(t2);
				//			len.push_back(k);
				ostringstream oss;
			       	oss<<k<<t2;
			       	cigar.append(oss.str());
			       	k=0;
		       	}
		       	k++; t2=t;
	       	}
		pos-=del_len-ins_len; //adjust start position
	}
	else{ 
		for(int i=0; i<aln_ref.length(); i++){
		       	if(aln_ref[i]==aln_read[i]){
			       	t='=';
				if (use_M) t='M';
		       	}
		       	else if(aln_ref[i]=='-'){
			       	t='I';
		       	}
		       	else if(aln_read[i]=='-'){
			       	t='D';
		       	}
		       	else{
			       	t='X';
				if (use_M) t='M';
		       	}
		       	if(t!=t2 && k>0){
			       	//			cType.push_back(t2);
				//			len.push_back(k);
				ostringstream oss;
			       	oss<<k<<t2;
			       	cigar.append(oss.str());
			       	k=0;
		       	}
		       	k++; t2=t;
	       	}
	}
       	ostringstream oss;
       	oss<<k<<t2;
	cigar.append(oss.str());
}

void samRead::reverse_comp(){
       	reverse(seq.begin(), seq.end());  
	reverse(qual.begin(), qual.end());  
	for(int i=0; i<seq.length(); i++){
	       	if (seq[i] == 'A') seq[i] = 'T';
	       	else if(seq[i] == 'T') seq[i] = 'A';
	       	else if(seq[i] == 'C') seq[i] = 'G';
	       	else if(seq[i] == 'G') seq[i] = 'C';
	       	else seq[i] = 'N';
       	}
}

void samRead::printRead(ostream& fout){
       	fout<<qname<<"\t"<<flag<<"\t"<<rname<<"\t"<<pos<<"\t"<<mapQ<<"\t"<<cigar<<"\t"<<rNext<<"\t"<<pNext<<"\t"<<tLen<<"\t"<<seq<<"\t"<<qual<<endl;
}
