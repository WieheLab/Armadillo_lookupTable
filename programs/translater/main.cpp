#include <math.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include "utilities.hpp"
#include <boost/algorithm/string.hpp>

using namespace std;

struct multicolumn_sort
{
  multicolumn_sort(int column){this->column=column;}
  bool operator () (vector<string> &a, vector<string> &b)
  {
    return atof(a[column].c_str())>atof(b[column].c_str());
  }
  int column;
};

//functions
string sys_call(string);
string aa_convert_code(string);

int main(int argc, char *argv[])
{  
   if (argc <2){cout << "USAGE: ./translater -i [fasta] -f [1 (default),2,3] -s [flag  to suppress stop codons]\n"; exit(1);}
   //ARG HANDLING
   int i=0, frame=1;
   string fasta_filename="";
   bool printStopsFlag=true;

   while(i<argc)
     {
       string arg=argv[i];
       string next_arg;
       if (i<argc-1){next_arg=argv[i+1];}else{next_arg="";}
       
       if ((arg.substr(0,1)=="-")&&(next_arg.substr(0,1)=="-")){cerr << "incorrectly formatted cmdline\n"; exit(1);}
  
       if (arg == "-i")
	 {
	   fasta_filename=next_arg;
	 }
       if (arg == "-f")
	 {
	   frame=atoi(next_arg.c_str());
	 }
	if(arg == "-s" || arg == "-nostops")
	{
	   printStopsFlag=false;
	}
     
       i++;
     }
   
   //read fasta
   map<string,string> sequences;
   vector<string> sequence_names;
   //read_fasta_file(fasta_filename, sequences, sequence_names);
   
   //translate
   map<string,string> dna_to_aa_tranx_map;
   get_aa_tranx_map(dna_to_aa_tranx_map);


   ///OPEN SEQUENCE FILE AND PUT INTO MAP
  ifstream file(fasta_filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "ERROR reading fasta file: could not open \"" << fasta_filename << "\"...exiting...\n"; exit(1);}

  string file_str;
  string sequence_string, name;
  int name_count=0, seq_count=0;

  while (!getline(file, file_str).eof())
    {
      chomp(file_str);
      if (file_str.substr(0,1) == ">")
        {
          name=file_str.substr(1,file_str.size()-1);
          //sequence_names.push_back(name);
          //name_count++;
        }
      else
        {
          sequence_string=file_str.substr(0,file_str.size());
	  if(sequence_string.size()<3)
		continue;
          boost::to_upper(sequence_string);
          //sequence_hash[name]+=sequence_string;
          //seq_count++;
  	 string aa_sequence;
	 translate_dna_to_aa(sequence_string,aa_sequence,frame,dna_to_aa_tranx_map);
	 int wcardPOS=aa_sequence.find('*');

	if(!printStopsFlag && wcardPOS>=0)
	 {
	   //cerr << "found wildcard "<<wcardPOS<<"\n";
	   //cerr << ">" << name<<"\n"<<aa_sequence<<"\n";
	   continue;
	 }	
	else	
	  cout << ">" << name<<"\n"<<aa_sequence<<"\n";

	//cout << sequence_string<<"\n";
        }
    }



   /*
   int size=sequence_names.size();
   for(int i=0; i<sequence_names.size(); i++)
     {
       print_pct_progress(i, size, 2);
       string aa_sequence;
       translate_dna_to_aa(sequences[sequence_names[i]],aa_sequence,frame,dna_to_aa_tranx_map);
       cout << ">" << sequence_names[i] << "\n" << aa_sequence << "\n"; 
     }	*/

    return 0;
}

string sys_call(string cmd)
{
  string data;
  FILE * stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  cmd.append(" 2>&1");
  
  stream = popen(cmd.c_str(), "r");
  if (stream)
    {
      while (!feof(stream))
	if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
      pclose(stream);
    }
  return data;
}

string aa_convert_code(string R)
{
  if (R.length()==3)
    {
      if (R=="ALA"){return "A";}
      if (R=="CYS"){return "C";}
      if (R=="ASP"){return "D";}
      if (R=="GLU"){return "E";}
      if (R=="PHE"){return "F";}
      if (R=="GLY"){return "G";}
      if (R=="HIS"){return "H";}
      if (R=="ILE"){return "I";}
      if (R=="LYS"){return "K";}
      if (R=="LEU"){return "L";}
      if (R=="MET"){return "M";}
      if (R=="ASN"){return "N";}
      if (R=="PRO"){return "P";}
      if (R=="GLN"){return "Q";}
      if (R=="ARG"){return "R";}
      if (R=="SER"){return "S";}
      if (R=="THR"){return "T";}
      if (R=="VAL"){return "V";}
      if (R=="TRP"){return "W";}
      if (R=="TYR"){return "Y";}
      cerr << "ERROR, unrecognized amino acid: '"<< R << "'. Exiting!\n"; exit(1);
    }
  if (R.length()==1)
    {
      if (R=="A"){return "ALA";}
      if (R=="C"){return "CYS";}
      if (R=="D"){return "ASP";}
      if (R=="E"){return "GLU";}
      if (R=="F"){return "PHE";}
      if (R=="G"){return "GLY";}
      if (R=="H"){return "HIS";}
      if (R=="I"){return "ILE";}
      if (R=="K"){return "LYS";}
      if (R=="L"){return "LEU";}
      if (R=="M"){return "MET";}
      if (R=="N"){return "ASN";}
      if (R=="P"){return "PRO";}
      if (R=="Q"){return "GLN";}
      if (R=="R"){return "ARG";}
      if (R=="S"){return "SER";}
      if (R=="T"){return "THR";}
      if (R=="V"){return "VAL";}
      if (R=="W"){return "TRP";}
      if (R=="Y"){return "TYR";}
      cerr << "ERROR, unrecognized amino acid: '" << R << "'. Exiting!\n"; exit(1);
    }
}
