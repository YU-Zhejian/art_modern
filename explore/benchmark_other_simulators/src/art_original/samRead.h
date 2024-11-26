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
#include <string>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <cstdio>
#include <set>
#include <algorithm>
#include "readSeqFile.h"

using namespace std;



class samHeader{
	public:
		samHeader(){ 
			VN="1.4";
		       	SO="unsorted";	
		}

	       	//@HD
		string VN; //=1.4
	       	string SO; //=unsorted	

		//@SQ
		vector<string> SN;
	       	vector<long> LN;
		
	       	//@PG
		string ID;
	       	string PN;
	       	string CL;

		void getRefseqID(char* seqfile);
		void printHeader(ostream& fout);
		void printAlnHeader(ostream& fout);
};

class samRead{

	public:
		samRead(){
			flag=0;
			mapQ=99;
		       	rNext="*";
		       	pNext=0;
			tLen=0;
		};
		string qname;
		int flag;
		string rname;
		long pos;
		short mapQ;
		string cigar;
		string rNext;
		long pNext;
		int tLen;
		string seq;
		string qual;
	       	void reverse_comp();
		void getCigar(string & aln_ref, string& aln_read, bool use_M=false);
		void printRead(ostream& fout);
};
