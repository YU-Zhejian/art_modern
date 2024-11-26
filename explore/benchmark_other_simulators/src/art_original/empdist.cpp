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

#include <cmath>
#include <cstdlib>
#include <ctime>
#include "empdist.h"

empdist::empdist(){
       	//map <string,string> mid2name;
	mid2name["GA1"]="Genome Analyzer I";
	mid2name["GA2"]="Genome Analyzer II";
	mid2name["HS10"]="HiSeq 1000";
	mid2name["HS20"]="HiSeq 2000";
	mid2name["HS25"]="HiSeq 2500"; 
	mid2name["MSv1"]="MiSeq v1";
	mid2name["MSv3"]="MiSeq v3";
	mid2name["MinS"]="MiniSeq TruSeq";
	mid2name["NS50"]="NextSeq 500 v2";
	mid2name["HSXn"]="HiSeqX v2.5 PCR free"; 
	mid2name["HSXt"]="HiSeqX v2.5 TruSeq"; 
};

bool empdist::setdist(string sequencer,  bool sep_qual, int length){
       	map <string, map<int, string> > m2q1, m2q2;
       	map <string, map<int, string> >::iterator it1, it2;
       	map<int, string>::iterator qt1,qt2;

	m2q1["GA1"][36]=QUAL_DIST_ONE36;
	m2q2["GA1"][36]=QUAL_DIST_TWO36;
	m2q1["GA1"][44]=QUAL_DIST_ONE44;
	m2q2["GA1"][44]=QUAL_DIST_TWO44;
	m2q1["GA2"][50]=QUAL_DIST_ONE50;
	m2q2["GA2"][50]=QUAL_DIST_TWO50;
	m2q1["GA2"][75]=QUAL_DIST_ONE75;
	m2q2["GA2"][75]=QUAL_DIST_TWO75;
	m2q1["HS10"][100]=QUAL_DIST_ONE100;
	m2q2["HS10"][100]=QUAL_DIST_TWO100;
	m2q1["MSv1"][250]=QUAL_DIST_ONE250;
	m2q2["MSv1"][250]=QUAL_DIST_TWO250;
	m2q1["HS20"][100]=HiSeq2000L100R1;
	m2q2["HS20"][100]=HiSeq2000L100R2;
	m2q1["HS25"][125]=HiSeq2500L125R1;
	m2q2["HS25"][125]=HiSeq2500L125R2;
	m2q1["HS25"][150]=HiSeq2500L150R1;
	m2q2["HS25"][150]=HiSeq2500L150R2; 
//MiniSeq sequences only single end
	m2q1["MinS"][50] = MiniSeqTruSeq_L50;
	m2q1["MSv3"][250]= MiSeqv3_L250_R1;
	m2q2["MSv3"][250]= MiSeqv3_L250_R2;  
	m2q1["NS50"][75] = NextSeq500v2_L75_R1;
	m2q2["NS50"][75] = NextSeq500v2_L75_R2;
	m2q1["HSXn"][150]= HiSeqXPCRfree_L150_R1;
	m2q2["HSXn"][150]= HiSeqXPCRfree_L150_R2;
	m2q1["HSXt"][150]= HiSeqXtruSeq_L150_R1;
	m2q2["HSXt"][150]= HiSeqXtruSeq_L150_R2;

	comb_sym = ".";
       	a_sym = "A";
       	t_sym = "T";
       	g_sym = "G";
       	c_sym = "C";
       	n_sym = "N"; 
	sep_qual_flag = sep_qual;
	it1=m2q1.find(sequencer);
	it2=m2q2.find(sequencer);
	if(it1==m2q1.end()){
	       	cerr<<"ART does not has a built-in profile for the given sequencing system: "<<sequencer<<endl; 
	       	cerr<<"List of all built-in profiles are:"<<endl<<endl; 
		short i=0;
		for (it1=m2q1.begin(); it1!=m2q1.end(); ++it1){
		       	for (qt1=it1->second.begin(); qt1!=it1->second.end(); ++qt1){
			       	i++;
			       	cerr <<i<<") "<<it1->first<<":"<< mid2name[it1->first] <<" length: "<<qt1->first<<endl;
			}
		}
		cerr<<endl;
		exit(1);
	}
	else if (it1->second.rbegin()->first < length){
	       	cerr<<"Error: the requested read length > the max length ("<<it1->second.rbegin()->first<<"bp) of any built-in profile for "<<sequencer; 
		exit(1);
	}
	ssystem=mid2name[it1->first];

       	qt1=it1->second.lower_bound(length); 
        istringstream  distss;
       	distss.str(qt1->second);
       	read_emp_dist(distss, true);

//if has only single end profile
	if(it2==m2q2.end()) return true;

       	qt2=it2->second.lower_bound(length); 
       	distss.clear();
       	distss.str(qt2->second);
       	read_emp_dist(distss, false);
	return true;
}

bool empdist::setdist(string file_name, string file_name1, bool sep_qual, int length){

    comb_sym = ".";
    a_sym = "A";
    t_sym = "T";
    g_sym = "G";
    c_sym = "C";
    n_sym = "N";

    sep_qual_flag = sep_qual;

    if(file_name.empty() && file_name1.empty()){
        istringstream  distss;
	if(length <= 36){
            distss.str(QUAL_DIST_ONE36);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(QUAL_DIST_TWO36);
            read_emp_dist(distss, false); 
	    ssystem=mid2name["GA1"];
	} else if(length <= 44){
            distss.str(QUAL_DIST_ONE44);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(QUAL_DIST_TWO44);
            read_emp_dist(distss, false);
	    ssystem=mid2name["GA1"];
	} else if(length <= 50){
            distss.str(QUAL_DIST_ONE50);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(QUAL_DIST_TWO50);
            read_emp_dist(distss, false);
	    ssystem=mid2name["GA2"];
	} else if(length <= 75){
            distss.str(QUAL_DIST_ONE75);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(QUAL_DIST_TWO75);
            read_emp_dist(distss, false);
	    ssystem=mid2name["GA2"];
	} else if(length <= 100){
            distss.str(HiSeq2000L100R1);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(HiSeq2000L100R2);
            read_emp_dist(distss, false);
	    ssystem=mid2name["HS20"];
	} else if(length <= 125){
            distss.str(HiSeq2500L125R1);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(HiSeq2500L125R2);
            read_emp_dist(distss, false);
	    ssystem=mid2name["HS25"];
	} else if(length <= 150){
            distss.str(HiSeq2500L150R1);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(HiSeq2500L150R2);
            read_emp_dist(distss, false);
	    ssystem=mid2name["HS25"];
	} else if(length <= 250){
            distss.str(QUAL_DIST_ONE250);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(QUAL_DIST_TWO250);
            read_emp_dist(distss, false);
	    ssystem=mid2name["MSv1"];
	}else {
		cerr<<"No read quality profile can generate "<< length << "bp reads!"<<endl; 
		exit(1);
	}
    } else if(file_name.empty()){
        istringstream  distss;
	if(length <= 36){
            distss.str(QUAL_DIST_ONE36);
            read_emp_dist(distss, true);
	} else if(length <= 44){
            distss.str(QUAL_DIST_ONE44);
            read_emp_dist(distss, true);
	} else if(length <= 50){
            distss.str(QUAL_DIST_ONE50);
            read_emp_dist(distss, true);
	} else if(length <= 75){
            distss.str(QUAL_DIST_ONE75);
            read_emp_dist(distss, true);
	} else if(length <= 100){
            distss.str(HiSeq2000L100R1);
            read_emp_dist(distss, true);
	} else if(length <= 125){
            distss.str(HiSeq2500L125R1);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(HiSeq2500L125R2);
            read_emp_dist(distss, false);
	} else if(length <= 150){
            distss.str(HiSeq2500L150R1);
            read_emp_dist(distss, true);
            distss.clear();
            distss.str(HiSeq2500L150R2);
            read_emp_dist(distss, false);
	} else {
            distss.str(QUAL_DIST_ONE250);
            read_emp_dist(distss, true);
	}

        read_emp_dist(file_name1, false);
    } else if(file_name1.empty()){
        read_emp_dist(file_name, true);

        istringstream  distss;
	if(length <= 36){
            distss.str(QUAL_DIST_TWO36);
            read_emp_dist(distss, false);
	} else if(length <= 44){
            distss.str(QUAL_DIST_TWO44);
            read_emp_dist(distss, false);
	} else if(length <= 50){
            distss.str(QUAL_DIST_TWO50);
            read_emp_dist(distss, false);
	} else if(length <= 75){
            distss.str(QUAL_DIST_TWO75);
            read_emp_dist(distss, false);
	} else if(length <= 100){
            distss.str(HiSeq2000L100R2);
            read_emp_dist(distss, false);
	} else if(length <= 125){
            distss.str(HiSeq2500L125R2);
            read_emp_dist(distss, false);
	} else if(length <= 150){
            distss.str(HiSeq2500L150R2);
            read_emp_dist(distss, false);
	} else {
            distss.str(QUAL_DIST_TWO250);
            read_emp_dist(distss, false);
	}
    } else {
        read_emp_dist(file_name, true);
        read_emp_dist(file_name1, false);
    }
//    srand ( (unsigned int)time(NULL) );
   return true;
}

//generate quali vector from dist of one read from pair-end [default first read]
bool empdist::get_read_qual(vector<short>& read_qual, int len, bool first){
    if(first){
        if(len>(int)qual_dist_first.size()){
            cerr<<"Error: Maximum read length allowed is:"<<qual_dist_first.size() << endl; 
            return false;
        }
        return get_read_qual(qual_dist_first, read_qual, len);
    }
    else{
        if(len>(int)qual_dist_second.size()){
            cerr<<"Error: Maximum length allowed for 2nd read is:"<<qual_dist_first.size() << endl; 
            return false;
        }
        return get_read_qual(qual_dist_second, read_qual, len);
    }
}

//generate quali vector from dist of first read from pair-end
bool empdist::get_read_qual(vector<short>& read_qual, bool first){
    if(first){
        return get_read_qual(qual_dist_first, read_qual, (int)qual_dist_first.size());
    }
    else{
        return get_read_qual(qual_dist_second, read_qual, (int)qual_dist_second.size());
    }
}

bool empdist::get_read_qual(vector< map <unsigned int, unsigned short> >&qual_dist, vector<short>& read_qual, int len){
    if((int)qual_dist.size()<len) return false;
    unsigned int cumCC;
    map <unsigned int, unsigned short>::iterator it;
    for(int i=0; i<len; i++){
        cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
        it=qual_dist[i].lower_bound(cumCC);
        //if(it!=qual_dist[i].end()){
        read_qual.push_back(it->second);
        //}
    }
    return true;
}

bool empdist::get_read_qual_1st(string seq, vector<short>& read_qual){
    int len=seq.size();
    if(len==0) return false;
    if(a_qual_dist_first.size()<len || t_qual_dist_first.size()<len || g_qual_dist_first.size()<len || c_qual_dist_first.size()<len ){
	    cerr<<"Error: The required read length exceeds the length of the read quality profile"<<endl;
	    exit(1);
    }

    unsigned int cumCC;
    map <unsigned int, unsigned short>::iterator it;
    for(int i=0; i<len; i++){
        cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
	if(seq[i]=='A'){
	       	it=a_qual_dist_first[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else if(seq[i]=='C'){
	       	it=c_qual_dist_first[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else if(seq[i]=='G'){
	       	it=g_qual_dist_first[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else if(seq[i]=='T'){
	       	it=g_qual_dist_first[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else{
		//return random quality less than 10 
	       	read_qual.push_back((short) r_prob()*10); 
	}
    }
    return true;
}

bool empdist::get_read_qual_2nd(string seq, vector<short>& read_qual){
    int len=seq.size();
    if(len==0) return false;
    if(a_qual_dist_second.size()<len || t_qual_dist_second.size()<len || g_qual_dist_second.size()<len || c_qual_dist_second.size()<len ){
	    cerr<<"Error: The required read length exceeds the length of the read quality profile"<<endl;
	    exit(1);
    }

    unsigned int cumCC;
    map <unsigned int, unsigned short>::iterator it;
    for(int i=0; i<len; i++){
        cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
	if(seq[i]=='A'){
	       	it=a_qual_dist_second[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else if(seq[i]=='C'){
	       	it=c_qual_dist_second[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else if(seq[i]=='G'){
	       	it=g_qual_dist_second[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else if(seq[i]=='T'){
	       	it=g_qual_dist_second[i].lower_bound(cumCC);
	       	read_qual.push_back(it->second);
	}
	else{
		//return random quality less than 10 
	       	read_qual.push_back((short) r_prob()*10); 
	}
    }
    return true;
}

empdist::~empdist(){
    qual_dist_first.clear();
    qual_dist_second.clear();
}

bool empdist::read_emp_dist(string infile, bool is_first){
    ifstream distss(infile.c_str());
    if(!distss){
        cerr<<"Fatal Error: Cannot open the distribution file: "<<infile<< endl;
        return false;
    }
    return read_emp_dist(distss, is_first);
}

bool empdist::read_emp_dist(istream& input, bool is_first){
    int linenum=0;
    int read_pos;
    char alt_read_pos;

    while(!input.eof()){
        bool a_flag = false;
        bool t_flag = false;
        bool g_flag = false;
        bool c_flag = false;

        string aLine;
        getline(input, aLine);
        if(aLine.length()==0) continue;
        if(aLine[0]=='#') continue;
		
	if( !sep_qual_flag && (aLine[0] != comb_sym[0])){ 
	    continue;
	} else if( sep_qual_flag && (aLine[0] == a_sym[0])){ 
	    a_flag = true;
	} else if( sep_qual_flag && (aLine[0] == t_sym[0])){ 
	    t_flag = true;
	} else if( sep_qual_flag && (aLine[0] == g_sym[0])){ 
	    g_flag = true;
	} else if( sep_qual_flag && (aLine[0] == c_sym[0])){ 
	    c_flag = true;
	} else if( sep_qual_flag && (aLine[0] == n_sym[0])){
	   continue;
	} else if( sep_qual_flag && (aLine[0] == comb_sym[0])){
	   continue;
	} 

        istringstream ss(aLine);
	ss>>alt_read_pos;
	ss>>read_pos;

        if(read_pos!=linenum){
	    linenum = 0;
	    if(read_pos!=linenum){
                cerr<<"Fatal error (1): Wrong format of input distribution."<<endl;
                exit(1);
	    }
        }

        unsigned short t_int;
        vector<unsigned short> qual;

        while (ss >> t_int){ qual.push_back(t_int); }

        getline(input, aLine);
        ss.clear();ss.str(aLine);
	ss>>alt_read_pos;
	ss>>read_pos;

        if(read_pos!=linenum){
            cerr<<"Fatal error (2): Wrong format of input distribution."<<endl;
            exit(2);
        }

        unsigned long t_uint;
        vector<unsigned long> count;

        while (ss >> t_uint){  count.push_back(t_uint); }

        if(count.size()!=qual.size()){
            cerr<<"Fatal error (3): Wrong format of input distribution."<<endl;
            exit(3);
        }

        double denom=count[count.size()-1]/(double)max_dist_number;
        map<unsigned int, unsigned short> dist;

        for(size_t i=0; i<count.size(); i++){
            unsigned int cc=(unsigned int)ceil(count[i]/denom);           
            dist[cc]=qual[i];
        }
        if(dist.size()>0){
            linenum++;

	    if(!sep_qual_flag && is_first){
		qual_dist_first.push_back(dist);
	    } else if(!sep_qual_flag && !is_first){
		qual_dist_second.push_back(dist);
	    } else if(a_flag && is_first){
		a_qual_dist_first.push_back(dist);
	    } else if(a_flag && !is_first){
		a_qual_dist_second.push_back(dist);
	    } else if(t_flag && is_first){
		t_qual_dist_first.push_back(dist);
	    } else if(t_flag && !is_first){
		t_qual_dist_second.push_back(dist);
	    } else if(g_flag && is_first){
		g_qual_dist_first.push_back(dist);
	    } else if(g_flag && !is_first){
		g_qual_dist_second.push_back(dist);
	    } else if(c_flag && is_first){
		c_qual_dist_first.push_back(dist);
	    } else if(c_flag && !is_first){
		c_qual_dist_second.push_back(dist);
	    } else {
		cerr << "Unexpected Error: Profile was not read in correctly." << endl;
		exit(5);	
	    }
        }
    }
    //max_len_read=linenum;
    if(linenum==0) return false;
    return true;
}

bool empdist::read_emp_dist(istream& input, map <unsigned int, unsigned int>& dist){
    while(!input.eof()){
        string aLine;
        getline(input, aLine);
        if(aLine.length()==0) continue;
        if(aLine[0]=='#') continue;
        istringstream ss(aLine);
        unsigned int t_int;
        vector<unsigned int> value;
        while (ss >> t_int){ value.push_back(t_int); }
        getline(input, aLine);
        ss.clear();ss.str(aLine);
        vector<unsigned int> count;
        while (ss >> t_int){  count.push_back(t_int); }
        if(count.size()!=value.size()){
            cerr<<"Fatal error (4): Wrong format of input distribution."<<endl;
            exit(4);
        }
        double denom=count[count.size()-1]/(double)max_dist_number;
        for(size_t i=0; i<count.size(); i++){
            unsigned int cc=(unsigned int)ceil(count[i]/denom);           
            dist[cc]=value[i];
        }
        break;
    }
    if(dist.size()==0) return false;
    return true;
}

bool empdist::read_emp_dist(string infile, map <unsigned int, unsigned int>& dist){
    ifstream distss(infile.c_str());
    if(!distss){
        cerr<<"Cannot open the distribution file "<<infile<<endl;
        return false;
    }
    return read_emp_dist(distss, dist);
}
