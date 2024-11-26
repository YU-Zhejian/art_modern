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

#include <vector>
#include <fstream>
#include <string>
using namespace std;

class readSeqFile{
  private:
      char* fileName;
      ifstream infile;
  public:
      readSeqFile(char* file_name);
      readSeqFile();
      ~readSeqFile();
      //close opened file and open another file
      bool reSetFile(char* file_name);
  	  //read next num_next seq in fasta format
      int next_seq(vector<string>& nextid,vector<string>& nextseq, int num_next);
      int next_seq(vector<string>& nextid,vector<string>& nextseq);
      int next_seq(string& geneID,string& geneSEQ);
      //set pointer to begining
      bool restart();
 
};
