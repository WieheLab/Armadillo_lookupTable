#include <math.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <climits>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>

using namespace std;

///functions
void chomp(string &);
void read_fasta_file(string, map<string, string> &, vector<string> &);
void tokenize(const string&, std::vector<string>&, const string& );
int match_score(char, char, double, double);
template <class T> 
void print_vector_of_vector(vector<vector<T> > &); 
double Max(double,double,double,int &);
double Max(double,double,int &);
void pairwise_align_sequences_semiglobal_w_affine_gap(map<char,map<char,int> > &, string,string,double,double,string &, string &,double &);
void pairwise_align_sequences_global_w_affine_gap(map<char,map<char,int> > &, string,string,double,double,string &, string &,double &);
vector < vector <string> > read_delimited_file(string, string);
void dna_alignment_summary_stats(string,string,int &,int &,int &,int &);

///globals

int main(int argc, char *argv[])
{  
  if (argc <2){cout << "USAGE: ./align_all_to_first -i [fasta file]  -go [gap open] -ge [gap extend] -matrix [scoring matrix file] -retain_ins_to_file [insertion file] -align_type [default=global,semi-global]\n"; exit(1);}
  int i=0;
  string score_matrix_file, sequence_filename, align_type="global", insertion_out_file="";
  double gap_open=-2, gap_extend=-0.5, match_weight=1,mismatch_weight=-1;
   while(i<argc)
     {
       string arg=argv[i];
       string next_arg;
       if (i<argc-1){next_arg=argv[i+1];}else{next_arg="";}
       
       if ((arg.substr(0,1)=="-")&&(next_arg.substr(0,1)=="-")){cerr << "incorrectly formatted cmdline\n"; exit(1);}
       if (arg == "-i")
	 {
	   sequence_filename=next_arg;
	 }
       if (arg == "-matrix")
	 {
	   score_matrix_file=next_arg;
	   if (score_matrix_file == "EDNAFULL.txt")
	     {
	       gap_open=-2;
	       gap_extend=-1;
	     }
	   else if(score_matrix_file == "BLOSUM62.txt")
	     {
	       gap_open=-11;
	       gap_extend=-1;
	     }
	 }
       if (arg == "-align_type")
	 {
	   align_type=next_arg;
	 }
       if (arg == "-retain_ins_to_file")
	 {
	   insertion_out_file=next_arg;
	 }
       if (arg == "-DNA" || arg == "-nt")
	 {
	   gap_open=-2;
	   gap_extend=-1;
	 }
       if (arg == "-AA" || arg == "-protien")
	 {
	   gap_open=-11;
	   gap_extend=-1;
	 }
       i++;
     }

   i=0;
   while(i<argc)
     {
       string arg=argv[i];
       string next_arg;
	 if (i<argc-1){next_arg=argv[i+1];}else{next_arg="";}
        if (arg == "-go")
	 {
	   gap_open=-1*atof(next_arg.c_str());
	 }
       if (arg == "-ge")
	 {
	   gap_extend=-1*atof(next_arg.c_str());
	 }
       i++;
     }
   
   ///Read scoring matrix file
   vector<vector<string> > file_contents=read_delimited_file(score_matrix_file,"\t");
   map<char,map<char,int> > scoring_matrix;
   vector<char> col_header;
   for(int j=0; j<file_contents[0].size(); j++){col_header.push_back(file_contents[0][j][0]);}
   for(int i=1; i<file_contents.size(); i++)
     {
       for(int j=1; j<file_contents[i].size(); j++)
	 {
	   scoring_matrix[file_contents[i][0][0]][col_header[j-1]]=atoi(file_contents[i][j].c_str());
	 }
     }
   
   ///Read fasta
   map<string,string> sequence_data;
   vector<string> sequence_names;
   read_fasta_file(sequence_filename, sequence_data, sequence_names);
   string ref_sequence=sequence_data[sequence_names[0]];
   //cout << ">" << sequence_names[0] << "\n" << sequence_data[sequence_names[0]] << "\n"; 
   ///Align all sequences to the first, keeping length of first sequence, thus insertion sequences are ignored, but reported
   int insertion_count=0, deletion_count=0;
   cout << ">" << sequence_names[0] << "\n" << ref_sequence << "\n";
   cerr << "---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|\n"; 
   ofstream output_insertions;
   output_insertions.open(insertion_out_file.c_str());
   output_insertions.close();
   if (insertion_out_file!="")
     {
       output_insertions.open(insertion_out_file.c_str(), std::ofstream::out | std::ofstream::app);
       if (!output_insertions.is_open()){cerr << "can't open output file " << insertion_out_file << "\n"; exit(1);}
     }
   for(int i=1; i<sequence_names.size(); i++)
     {
       if (sequence_names.size()>100){if (i%(sequence_names.size()/100)==0){cerr << ".";}}
       string alignment_sequence1="",alignment_sequence2="";
       double score=0;
       if (align_type=="global"){
	 pairwise_align_sequences_global_w_affine_gap(scoring_matrix, ref_sequence,sequence_data[sequence_names[i]], gap_open, gap_extend,alignment_sequence1, alignment_sequence2, score);}
       else if (align_type=="semi-global"){
	 pairwise_align_sequences_semiglobal_w_affine_gap(scoring_matrix, ref_sequence,sequence_data[sequence_names[i]], gap_open, gap_extend,alignment_sequence1, alignment_sequence2, score);}
       else{cerr << "unrecongized option for alignment type...exiting\n"; exit(1);}

       ///trim off any gap overhang
       string alignment_sequence1_trimmed, alignment_sequence2_trimmed;
       alignment_sequence1_trimmed=alignment_sequence1;
       alignment_sequence2_trimmed=alignment_sequence2;

       int start_seq=0, end_seq=alignment_sequence2.size()-1;
    
       for(int j=0; j<alignment_sequence2.size(); j++)
	 {
	   //   cout << "here: <" << alignment_sequence1[j] << ">\n"; 
	   if (alignment_sequence1[j]!='-'){start_seq=j;break;}
	   
	 }
       bool end_with_gaps=true;
       for(int j=alignment_sequence2.size()-1; j>=0; j--)
	 {
	   if (alignment_sequence1[j]!='-'){end_seq=j; break;}
	 }

       if (end_seq < alignment_sequence2.size()-1){alignment_sequence2_trimmed.erase(end_seq+1); alignment_sequence1_trimmed.erase(end_seq+1);}
       if (start_seq >0){alignment_sequence2_trimmed.erase(0,start_seq); alignment_sequence1_trimmed.erase(0,start_seq);}
      
       ///check for insertions and discard
       if (alignment_sequence1_trimmed.find("-") != string::npos)
	 {
	   insertion_count++; /*cerr << "insertion in " << sequence_names[i] << ". Sequence discarded from alignment\n";*/ 
	   if (insertion_out_file!="")
	     {
	       output_insertions << /*">" << sequence_names[0] << "\n" << alignment_sequence1_trimmed << "\n" <<*/ ">" << sequence_names[i] << "\n" << alignment_sequence2_trimmed << "\n"; 
	     }
	   continue; 
	 }
       if (alignment_sequence2_trimmed.find("-") != string::npos)
	 {
	   deletion_count++; //cerr << "deletion in " << sequence_names[i] << "\n"; 
	 }
       cout << ">" << sequence_names[i] << "\n" << alignment_sequence2_trimmed << "\n";
     }
   cerr << "\n"; 
   cerr << "num insertions: " << insertion_count << "\n"; 
   cerr << "num deletions: " << deletion_count << "\n"; 
   return 0;
}


double Max(double a, double b, int &dir)
{
  if (a>=b){dir=1;return a;}else{dir=2;return b;}
}
double Max(double a, double b, double c, int &dir)
{
  //cout << "MAX: " << a  << " " << b << " " << c << "\n"; 
  if ((a>=b) && (a>=c)){dir=1;return a;}
  if ((b>=a) && (b>=c)){dir=2;return b;}
  else{dir=3; return c;}
}

int match_score(char a, char b, double match_weight, double mismatch_weight)
{
  if (a==b){return match_weight;}
  else{return mismatch_weight;}
  
}
template <class T>
void print_vector_of_vector(vector<vector<T> > &M)
{
  for(int j=0; j<M[0].size(); j++){cout << "\t" << j;}cout << "\n"; 
  for(int i=0; i<M.size(); i++)
    {
      cout << i;
      for(int j=0; j<M[i].size(); j++)
	{
	  cout << "\t" << M[i][j];
	  //if (j!=M[i].size()-1){cout << "\t";}
	}
	cout << "\n"; 
    }
}

void read_fasta_file(string filename, map<string, string> &sequence_hash, vector<string> &sequence_names)
{
   ///OPEN SEQUENCE FILE AND PUT INTO MAP
  ifstream file(filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "could not open " << filename << " ...exiting...\n"; exit(1);}

  string file_str;
  string sequence_string, name;
  int name_count=0, seq_count=0;
  bool duplicate_name=false;
  while (!getline(file, file_str).eof())
    {
      chomp(file_str);
     
      if (file_str.substr(0,1) == ">")
	{
	  name=file_str.substr(1,file_str.size()-1);
	  
	  if (sequence_hash.find(name)!=sequence_hash.end()){cerr << "duplicate name in fasta file: " << name << " ...skipping\n"; duplicate_name=true; continue;}else{duplicate_name=false;} 
	  sequence_names.push_back(name);
	  name_count++;
	}
      else
	{
	  sequence_string=file_str.substr(0,file_str.size());
	  boost::to_upper(sequence_string);
	  if (duplicate_name==false)
	    {
	      sequence_hash[name]+=sequence_string;
	      seq_count++;
	    }
	}
    }
  if (name_count != seq_count){cerr << "WARNING: The number of sequence names does not match the number of sequences\n";}
}

void chomp(string &str)
{
  std::string whitespaces ("\n\r");
  
  std::size_t found = str.find_last_not_of(whitespaces);
  if (found!=std::string::npos)
    str.erase(found+1);
  else
    str.clear();            // str is all whitespace
}



void tokenize(const string& str, std::vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
vector < vector <string> > read_delimited_file(string filename, string delimiter)
{
  vector < vector <string> > file_contents;

    ifstream file(filename.c_str(), ios::in );
    if (!file.is_open()) {cerr << "could not open " << filename << " ...exiting...\n"; exit(1);}
    string line_str;
    while (getline(file, line_str))
    {
      chomp(line_str);
      vector <string> strs;
        tokenize(line_str, strs, delimiter);
        file_contents.push_back(strs);
    }
    file.close();
    if (file_contents.size()==0){cerr << "Problem reading file: " << filename << " ...exiting\n"; exit(1);}
    return file_contents;
}



void pairwise_align_sequences_semiglobal_w_affine_gap(map<char, map<char, int> > &scoring_matrix, string sequence1,string sequence2,double gap_open,double gap_extend,string &alignment_sequence1, string &alignment_sequence2,double &score)
{

  ///ERROR: not properly handling the end of alignments when sequences are equal
  // cerr << "ERROR NOT WORKING PROPERLY\n"; exit(1);

  ///initialize
  alignment_sequence1="";
  alignment_sequence2="";

  ///convert sequences to vectors of chars
  vector<char> sequence1_vctr(sequence1.begin(),sequence1.end());
  vector<char> sequence2_vctr(sequence2.begin(),sequence2.end());
  sequence1_vctr.insert(sequence1_vctr.begin(),'-'); ///start with gap
  sequence2_vctr.insert(sequence2_vctr.begin(),'-'); ///start with gap
  
  ///initialize matrices
  vector<vector<double> > M(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<double> > A(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<double> > B(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<char> > M_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));
  vector<vector<char> > A_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));
  vector<vector<char> > B_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));

  ///set boundaries
  double boundary=-INT_MAX;
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      M[i][0]=0;
      B[i][0]=boundary;
    }
  for(int j=1; j<sequence2_vctr.size(); j++)
    {
      M[0][j]=0;
      A[0][j]=boundary; 
    }
  
  ///fill in pointer matrices
  M[0][0]=0; 
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      for(int j=1; j<sequence2_vctr.size(); j++)
	{
	  
	   int num;
	   A[i][j]=Max(M[i-1][j]+gap_open+gap_extend, A[i-1][j]+gap_extend, B[i-1][j]+gap_open+gap_extend,num);
	   if (num==1){A_ptr[i][j]='M';}else if (num==2){A_ptr[i][j]='A';}else if (num==3){A_ptr[i][j]='B';}

	   B[i][j]=Max(M[i][j-1]+gap_open+gap_extend,B[i][j-1]+gap_extend,A[i][j-1]+gap_open+gap_extend,num);
	   if (num==1){B_ptr[i][j]='M';}else if (num==2){B_ptr[i][j]='B';}else if (num==3){B_ptr[i][j]='A';}

	   M[i][j]=scoring_matrix[sequence1_vctr[i]][sequence2_vctr[j]]+Max(M[i-1][j-1],A[i-1][j-1],B[i-1][j-1], num);
	   if (num==1){M_ptr[i][j]='M';}else if (num==2){M_ptr[i][j]='A';}else if (num==3){M_ptr[i][j]='B';}

	 }
     }
  
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      M_ptr[i][0]='A';
      A_ptr[i][0]='A';//new
      B_ptr[i][0]='A';//new
      
    }
  for(int i=1; i<sequence2_vctr.size(); i++)
    {
      M_ptr[0][i]='B';
      B_ptr[0][i]='B';//new
      A_ptr[0][i]='B';//new
    }
  
   string status="continue";
   int s1=sequence1_vctr.size()-1, s2=sequence2_vctr.size()-1;
   vector<vector< char> > align(sequence1_vctr.size()+sequence2_vctr.size(),vector<char>(2));
   int counter=0;
   int num;
  
   ///semi global alignment starts with max of last row and last column, not corner cell
   double max_last_val=0;
   int maxi,maxj;
   for(int i=0; i<sequence1_vctr.size(); i++)
     {
       if (M[i][sequence2_vctr.size()-1]>max_last_val){max_last_val=M[i][sequence2_vctr.size()-1]; maxi=i;maxj=sequence2_vctr.size()-1;}
     }
   for(int j=0; j<sequence2_vctr.size(); j++)
     {
       if (M[sequence1_vctr.size()-1][j]>max_last_val){max_last_val=M[sequence1_vctr.size()-1][j]; maxi=sequence1_vctr.size()-1;maxj=j;}
     }
   ///record score as maximum of last row and last col
   score=max_last_val;
    
   // exit(1);
   ///go from last cell to maxi/maxj
   if (s1==maxi){for(int j=s2; j>maxj; j--){align[counter][0]='-';align[counter][1]=sequence2_vctr[j];counter++;}}
   if (s2==maxj){for(int i=s1; i>maxi; i--){align[counter][0]=sequence1_vctr[i];align[counter][1]='-';counter++;}}
   s1=maxi;
   s2=maxj;
   //cout << s1 << "\t" << s2 << "\n"; 
   char current_matrix='M';
   int iter=0;

   while(status=="continue")
     {
       if (++iter>100000){cerr << "infinite loop\n"; exit(1);}
       if (s1==0)
	 {
	   for(int i=s2; i>0; i--){align[counter][0]='-'; align[counter][1]=sequence2_vctr[i];counter++;}
	   status="end";
	   break;
	 }
       if (s2==0)
	 {
	   for(int i=s1; i>0; i--){align[counter][1]='-'; align[counter][0]=sequence1_vctr[i];counter++;}
	   status="end";
	   break;
	 }

       /*
       if ((s1==0) && (s2!=0))//base condition, top row, go from here to 0,0 filling with gaps
	 {
	   for(int j=s2; j>=0; j--){cout << s1 << "," << j << "\n"; align[counter][0]='-'; counter++;}
	   break;
	 }
     
       if ((s1!=0) && (s2==0))//base condition, first col, go from here to 0,0 filling with gaps
	 {
	   for(int j=s1; j>=0; j--){align[counter][1]='-'; counter++;}
	   break;
	   }
       */
       if (current_matrix=='M')
	 {
	   current_matrix=M_ptr[s1][s2]; //go to this matrix next
	   ///if in M, go diagonal, if going diagonal it's a match
	   align[counter][0]=sequence1_vctr[s1]; 
	   align[counter][1]=sequence2_vctr[s2];
	   s1--; 
	   s2--; 
	 }
       else if (current_matrix=='A')
	 {
	   current_matrix=A_ptr[s1][s2];
	   ///if in A, go up, if going up, it's a del relative to first seq
	   align[counter][0]=sequence1_vctr[s1]; 
	   align[counter][1]='-';
	   s1--; 
	 }
       else if (current_matrix=='B')
	 {
	   current_matrix=B_ptr[s1][s2];
	   //if in B, go left, if going left, it's an ins relative to first seq 
	   align[counter][0]='-'; 
	   align[counter][1]=sequence2_vctr[s2];
	   s2--; 
	 }
       else {cerr << "should not reach here\n"; exit(1);}
       // cout << "--> " << s1 << "," << s2 << " in " << current_matrix << "\n"; 
       ///<< align[counter][0] << "|" << align[counter][1] << "\n"; 
       if (current_matrix=='x'){cerr << "should not reach x\n"; exit(1);}
       counter++;
       
     }

   for(int i=counter-1; i>=0; i--)
     {
       alignment_sequence1+=align[i][0];
     }
   for(int i=counter-1; i>=0; i--)
     {
       alignment_sequence2+=align[i][1];
     }

   return;
}

void pairwise_align_sequences_global_w_affine_gap(map<char, map<char, int> > &scoring_matrix, string sequence1,string sequence2,double gap_open,double gap_extend,string &alignment_sequence1, string &alignment_sequence2,double &score)
{
  ///initialize
  alignment_sequence1="";
  alignment_sequence2="";

  ///convert sequences to vectors of chars
  vector<char> sequence1_vctr(sequence1.begin(),sequence1.end());
  vector<char> sequence2_vctr(sequence2.begin(),sequence2.end());
  sequence1_vctr.insert(sequence1_vctr.begin(),'-'); ///start with gap
  sequence2_vctr.insert(sequence2_vctr.begin(),'-'); ///start with gap
  
  sequence1_vctr.insert(sequence1_vctr.end(),'-'); ///pad last with gap
  sequence2_vctr.insert(sequence2_vctr.end(),'-'); ///pad last with gap

  
  ///initialize matrices
  vector<vector<double> > M(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<double> > A(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<double> > B(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<char> > M_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));
  vector<vector<char> > A_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));
  vector<vector<char> > B_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));
  //getchar();
  
  ///set boundaries
  double boundary=-INT_MAX;
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      M[i][0]=boundary;
      B[i][0]=boundary;
      A[i][0]=gap_open+(i*gap_extend);
    }
  for(int j=1; j<sequence2_vctr.size(); j++)
    {
      M[0][j]=boundary;
      A[0][j]=boundary; 
      B[0][j]=gap_open+(j*gap_extend);
    }
  
  ///fill in pointer matrices
  M[0][0]=0; 
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      for(int j=1; j<sequence2_vctr.size(); j++)
	{
	  
	   int num;
	   A[i][j]=Max(M[i-1][j]+gap_open+gap_extend, A[i-1][j]+gap_extend, B[i-1][j]+gap_open+gap_extend,num);
	   if (num==1){A_ptr[i][j]='M';}else if (num==2){A_ptr[i][j]='A';}else if (num==3){A_ptr[i][j]='B';}

	   B[i][j]=Max(M[i][j-1]+gap_open+gap_extend,B[i][j-1]+gap_extend,A[i][j-1]+gap_open+gap_extend,num);
	   if (num==1){B_ptr[i][j]='M';}else if (num==2){B_ptr[i][j]='B';}else if (num==3){B_ptr[i][j]='A';}

	   M[i][j]=scoring_matrix[sequence1_vctr[i]][sequence2_vctr[j]]+Max(M[i-1][j-1],A[i-1][j-1],B[i-1][j-1], num);
	   if (num==1){M_ptr[i][j]='M';}else if (num==2){M_ptr[i][j]='A';}else if (num==3){M_ptr[i][j]='B';}

	 }
     }
  //  cout << "M:\n"; 
  //print_vector_of_vector(M);
  
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      M_ptr[i][0]='A';
      A_ptr[i][0]='A';//new
      B_ptr[i][0]='A';//new
      
    }
  for(int i=1; i<sequence2_vctr.size(); i++)
    {
      M_ptr[0][i]='B';
      B_ptr[0][i]='B';//new
      A_ptr[0][i]='B';//new
    }
  
   ///global alignment starts at corner cell
   score=M[sequence1_vctr.size()-1][sequence2_vctr.size()-1];

   ///do traceback
   char current_matrix='M';
   int iter=0;
   string status="continue";
   int s1=sequence1_vctr.size()-1, s2=sequence2_vctr.size()-1;
   vector<vector< char> > align(sequence1_vctr.size()+sequence2_vctr.size(),vector<char>(2));
   int counter=0;

   while(status=="continue")
     {

       if (++iter>2000){cerr << "infinite loop\n"; exit(1);}
       if ((s1==0)&&(s2==0)){status="end";break;}
       if (current_matrix=='M')
	 {
	   current_matrix=M_ptr[s1][s2]; //go to this matrix next
	   ///if in M, go diagonal, if going diagonal it's a match
	   align[counter][0]=sequence1_vctr[s1]; 
	   align[counter][1]=sequence2_vctr[s2];
	   s1--; 
	   s2--; 
	 }
       else if (current_matrix=='A')
	 {
	   current_matrix=A_ptr[s1][s2];
	   ///if in A, go up, if going up, it's a del relative to first seq
	   align[counter][0]=sequence1_vctr[s1]; 
	   align[counter][1]='-';
	   s1--; 
	 }
       else if (current_matrix=='B')
	 {
	   current_matrix=B_ptr[s1][s2];
	   //if in B, go left, if going left, it's an ins relative to first seq 
	   align[counter][0]='-'; 
	   align[counter][1]=sequence2_vctr[s2];
	   s2--; 
	 }
       else {cerr << "should not reach here\n"; exit(1);}
       if (current_matrix=='x'){cerr << "should not reach x\n"; exit(1);}
       counter++;
       
     }

   for(int i=counter-1; i>=0; i--)
     {
       alignment_sequence1+=align[i][0];
     }
   for(int i=counter-1; i>=0; i--)
     {
       alignment_sequence2+=align[i][1];
     }
  
   return;
}

void dna_alignment_summary_stats(string alignment_sequence1,string alignment_sequence2,int &mismatch_count, int &mismatch_wo_gaps_count, int &num_ins, int &num_del)
{
  mismatch_count=0; mismatch_wo_gaps_count=0; num_ins=0; num_del=0; 
  if (alignment_sequence1.length() != alignment_sequence2.length()){cerr << "incorrect alignment, unequal lengths, exiting...\n"; exit(1);}
  for(int i=0; i<alignment_sequence1.length(); i++)
    {
      if (alignment_sequence1[i] != alignment_sequence2[i])
	{
	  mismatch_count++;
	  if ((alignment_sequence1[i] != '-') && (alignment_sequence2[i] != '-')){mismatch_wo_gaps_count++;}
	}
      if (alignment_sequence1[i] == '-'){num_ins++;}
      if (alignment_sequence2[i] == '-'){num_del++;}
    }
  return;
}
