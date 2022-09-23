#include <math.h>
#include <cstdlib>
#include <iostream>
#include <random>
#include <iostream>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>

using namespace std;

//classes
class S5F_mut
{
public:
  string fivemer, group, subst_group;
  double score,score25,score75;
  map<char, double> substitutions;

  ~S5F_mut(){};
  S5F_mut(){};
  S5F_mut(string _fivemer,double _score, string _group, double _score25, double _score75)
  {
    fivemer=_fivemer;
    score=_score;
    group=_group;
    score25=_score25;
    score75=_score75;
  }
  
};

class Seq //nucleotide level sequence object (later could have each base point to a aa object)
{
public:
  int aa_num;
  string aa;
  string base;
  double S5F_mut_score;
  double simulated_aa_positional_frequency;

  ~Seq(){};
  Seq(){};
  Seq(string _base, int _aa_num,string _aa,double _S5F_mut_score)
  {
    base=_base;
    aa_num=_aa_num;
    aa=_aa;
    S5F_mut_score=_S5F_mut_score;
  }
};

///functions
void chomp(string &);
void tokenize(const string&, std::vector<string>&, const string&);
vector < vector <string> > read_delimited_file(string, string);
void read_fasta_file(string, map<string, string> &, vector<string> &);
void load_S5F_files(string,string, map<string,S5F_mut> &);
void translate_dna_to_aa(string &, string &, int, map<string,string> &);
void get_aa_tranx_map(map<string,string> &);
void process_fasta_sequence_to_seq_vector(string &,vector<Seq> &, map<string,string> &, map<string,S5F_mut> &);
void number_of_mutations_two_seqs(string &, string &, int &);
void simulate_S5F_mutation(string , int &, map<string,S5F_mut> &, mt19937 &, uniform_real_distribution<double> &, bool, vector<string> &, bool);
bool dna_sequence_has_stop_codon_in_reading_frame(string);
void simulate_clonal_lineage(string, vector<double> &,  map <string, S5F_mut> &, std::mt19937 &, std::uniform_real_distribution<double> &, vector<string> &);
void make_frequency_table_from_MSA(vector<string> &, map<int,map<char,double> > &);
void print_frequency_table(map<int,map<char,double> > &, string);
void print_transfac(map<int,map<char,double> > &, string);
void mean_std(vector<double> &, double &, double &);
void simulate_sequence_to_mutation_frequency(string, double, map<string,S5F_mut> &, mt19937 &, uniform_real_distribution<double> &, bool, vector<string> &, bool);

///templated functions
template <typename Type>
void vector2D_to_3D(vector<vector<Type> > &, int, vector<vector<vector<Type> > > &);
template <typename Type>
void vector1D_to_2D(vector<Type> &, int, vector<vector<Type> > &);
template <typename Type>
string convert_to_string(Type);


///TODO:
///1. fix X problem
///TODO
//1. Correct for clone if we do bootstrap in here
//2. Compute boostrapped mean frquency table correctly. 


int main(int argc, char *argv[])
{  
  if (argc <2){cout << "USAGE:positional_aa_matrix_makerc -i [aligned fasta] -o [frequency table] -transfac [transfac formatted]\n"; exit(1);}

  ///get cmdline args
  int i=0;
  string fasta_filename, output_filename, transfac_filename;
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
    
       if (arg == "-o")
	 {
	   output_filename=next_arg;
	 }
       if (arg == "-transfac")
	 {
	   transfac_filename=next_arg;
	 }
       i++;
     }

  
   ///read input sequence alignment
   map <string, string> sequences;
   vector <string> sequence_names;   
   cerr << "reading fasta file\n"; 
   read_fasta_file(fasta_filename, sequences, sequence_names);
   cerr << "done\n"; 

   //amino acids
   vector<char> aas={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'};

   string ref_aa_sequence=sequences[sequence_names[0]];
   vector<string> MSA_sequences;
   for(int i=1; i<sequence_names.size(); i++){MSA_sequences.push_back(sequences[sequence_names[i]]);}
   map<int,map<char,double> > freq_table;
   make_frequency_table_from_MSA(MSA_sequences, freq_table);

   //print out table
   print_frequency_table(freq_table, output_filename);
   print_transfac(freq_table,transfac_filename);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                FUNCTION DEFINITIONS
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simulate_sequence_to_mutation_frequency_N_times(string sequence, double mutation_frequency, int max_iter, map<string,S5F_mut> &S5F_model, mt19937 &gen, uniform_real_distribution<double> &dis, bool kill_stop_seqs, vector<string> &mature_mutant_sequences, map<int, map<char,double> > &mature_mutant_positional_aa_freqs, bool dont_count_remutation)
{
  /*EMPTY AS OF NOW*/
  
  return;
}

void print_transfac(map<int,map<char,double> > &frequency_table, string filename)
{
  vector<char> aas={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
  ostringstream output;
  output << "ID BLANK\nBF BLANK\nP0";
  for(int i=0; i<aas.size(); i++){output << "\t" << aas[i];} 
  output << "\n"; 
  
  for(int i=1; i<=frequency_table.size(); i++)
    {
      output << setfill('0') << setw(2) << i;
      double total=0;
      int X_counter=0;
      for(int j=0; j<aas.size(); j++)
	{
	  int count=(int) round(frequency_table[i][aas[j]]*100000);
	  output << "\t" << count;
	}
      output << "\n";
    }

  if (filename.empty()){cout << "TRANSFAC:\n" <<  output.str();}
  else
    {
      ofstream output_stream;
      output_stream.open(filename.c_str());
      if (!output_stream.is_open()){cerr<< "can't open output file " << filename << "\n"; exit(1);}
      output_stream << output.str();
      output_stream.close();
    }
  return;
}

void print_frequency_table(map<int,map<char,double> > &frequency_table, string filename)
{
  vector<char> aas={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'};
  
  ostringstream output;
  output.precision(10);

  for(int i=0; i<aas.size(); i++){output << "\t" << aas[i];} 
  output << "\n"; 
  for(int i=0; i<frequency_table.size(); i++)
    {
       output << i+1;
       double total=0;
       for(int j=0; j<aas.size(); j++)
	 {
	   output << "\t" << setprecision(10) << frequency_table[i+1][aas[j]];
	   total+=frequency_table[i+1][aas[j]];
	 }
       output << "\t" << total << "\n"; 
    }
  if (filename.empty()){cout << "FREQEUNCY TABLE:\n" << output.str();}
  else
    {
      ofstream output_stream;
      output_stream.open(filename.c_str());
      if (!output_stream.is_open()){cerr<< "can't open output file " << filename << "\n"; exit(1);}
      output_stream << output.str();
      output_stream.close();
    }
  return;
}
void make_frequency_table_from_MSA(vector<string> &sequences, map<int,map<char,double> > &frequency_table)
{
  int N=sequences.size();

  for(int i=0; i<sequences.size(); i++)
    {
      for(int j=0; j<sequences[i].length(); j++)
	{
	  frequency_table[j+1][sequences[i].at(j)]+=1/(double) N;
	}
    }
  return;
}



void simulate_clonal_lineage(string starting_sequence, vector<double> &clone_mut_freqs,  map <string, S5F_mut> &S5F_5mers, std::mt19937 &gen, std::uniform_real_distribution<double> &dis, vector<string> &mature_clonal_members)
{
  int sequence_length=starting_sequence.length();
  vector<int> observed_mutation_counts;
  for(int i=0; i<clone_mut_freqs.size(); i++)
    {
      int mutation_count=(int) round(clone_mut_freqs[i]*sequence_length);
      observed_mutation_counts.push_back(mutation_count);
    }
  
  //sort observed mutation counts
  sort(observed_mutation_counts.begin(), observed_mutation_counts.end(), std::greater<int>());

  //Simulate maturation 

  if (dna_sequence_has_stop_codon_in_reading_frame(starting_sequence)){cerr << "germline has stop codon...exiting\n"; exit(1);}

 
   //cerr << "Simulating maturation of clonal lineage...\n"; 

   vector<string> ancestor_sequences;
  
   string ref_sequence=starting_sequence;
   ancestor_sequences.push_back(ref_sequence);

   for(int i=0; i<observed_mutation_counts.size(); i++)
     {
       if (i==0)
	 {//start from UCA
	   vector<string> mutant_sequences;
	   if (observed_mutation_counts[i] == 0)
	     {
	       mature_clonal_members.push_back(starting_sequence);
	     }
	   else
	     {
	       simulate_S5F_mutation(starting_sequence, observed_mutation_counts[i], S5F_5mers, gen, dis,true, mutant_sequences, true);
	    
	       for(int j=0; j<mutant_sequences.size(); j++)
		 {
		   ancestor_sequences.push_back(mutant_sequences[j]); 
		 }
	       mature_clonal_members.push_back(mutant_sequences[mutant_sequences.size()-1]);
	     }
	   // int m;
	   // number_of_mutations_two_seqs(ref_sequence,mutant_sequences[mutant_sequences.size()-1], m);
	   //	   cerr << "exp: " << observed_mutation_counts[i] << "\tobs: " << m << "\n"; 
	 }
       else //start from random ancestor
	 {
	   //randomly choose ancestor and mutate from there
	   std::uniform_int_distribution<int> distribution(0,ancestor_sequences.size()-1);
	   int rand_index = distribution(gen);  
	   string random_ancestor_sequence=ancestor_sequences[rand_index];
	   int ancestor_mutation_count;
	   number_of_mutations_two_seqs(ref_sequence,random_ancestor_sequence, ancestor_mutation_count);
	   //repeat draw if the ancestor is more mutated than the observed
	   while (ancestor_mutation_count > observed_mutation_counts[i])
	     {
	       rand_index=distribution(gen);
	       random_ancestor_sequence=ancestor_sequences[rand_index];
	       number_of_mutations_two_seqs(ref_sequence,random_ancestor_sequence, ancestor_mutation_count);
	     }
	   
	   vector<string> mutant_sequences;
	 
	   int mut=(observed_mutation_counts[i]-ancestor_mutation_count);
	   if (mut==0) //use ancestor sequence if no further mutations are needed
	     {
	       mature_clonal_members.push_back(random_ancestor_sequence);
	     }
	   else
	     {
	       //cerr << "MUT: " << mut << "\n"; 
	       simulate_S5F_mutation(random_ancestor_sequence,mut, S5F_5mers, gen, dis, true, mutant_sequences, true);
	       for(int j=0; j<mutant_sequences.size(); j++)
		 {
		   ancestor_sequences.push_back(mutant_sequences[j]);
		 }
	       assert(mutant_sequences.size() > 0);
	       mature_clonal_members.push_back(mutant_sequences[mutant_sequences.size()-1]);
	       //	       cerr << "SIZE: " << mutant_sequences.size() << "\n"; 
	     }
	   
	 }
       
     }

}

void simulate_S5F_mutation(string sequence, int &num_mutations, map<string,S5F_mut> &S5F_model, mt19937 &gen, uniform_real_distribution<double> &dis, bool kill_stop_seqs, vector<string> &mutant_sequences, bool dont_count_remutation)
{

  string ref_sequence=sequence;
  ///iterate num_mutations times
  int accrued_mutations=0;
  while (accrued_mutations<num_mutations)
    {
      //cerr << "A: " << accrued_mutations << "\t" << num_mutations << "\n"; 
      ///get mutability scores  
      vector<double> mut_scores;
      mut_scores.push_back(1);///first two positions set to neutral
      mut_scores.push_back(1);///
      double sum_mut_scores=2.0;
      for(int i=2; i<sequence.length()-2; i++)
	{
	  string fivemer=sequence.substr(i-2,5);
	  double mut_score=S5F_model[fivemer].score;
	  mut_scores.push_back(mut_score);
	  sum_mut_scores+=mut_score;
	}
      
      mut_scores.push_back(1);///last two positions set to neutral
      mut_scores.push_back(1);///
      sum_mut_scores+=2.0;
      // cerr << accrued_mutations << "\tsum of mut scores: " << sum_mut_scores << "\n"; 
      
      ///convert mutability scores to probability of mutating position  
      vector<double> mut_probability_ladder(sequence.length(),0); ///cumulative distribution of probabilities
      mut_probability_ladder[0]=mut_scores[0]/(double) sum_mut_scores;
      //cerr << "0" << "\t" << sequence[0] << "\t" << (mut_scores[0]/(double) sum_mut_scores) << "\t" << mut_probability_ladder[0] << "\n"; 
      for(int i=1; i<sequence.length(); i++)
	{
	  mut_probability_ladder[i]=mut_probability_ladder[i-1]+(mut_scores[i]/(double) sum_mut_scores); 
	  //cerr << i << "\t" << sequence[i] << "\t" << (mut_scores[i]/(double) sum_mut_scores) << "\t" << mut_probability_ladder[i] << "\n"; 
	}
      
      ///draw position randomly according to probability ladder
      double R=dis(gen);
      int mutate_position_i;
      for(int i=0; i<mut_probability_ladder.size(); i++)  ///OPTIMIZE: combine the two loops into one
	{
	  if (R<mut_probability_ladder[i])
	    {
	      mutate_position_i=i;
	      break;
	    }
	}
      //cerr << accrued_mutations << "\t" << R << "\t" << mutate_position_i << "\t" << sequence[mutate_position_i] << "\t" << mut_scores[mutate_position_i] << "\t"; 

      ///mutate position according to substitution model
      map<char, double> substitution_probs;
      
      if ((mutate_position_i<=1) || (mutate_position_i>=sequence.length()-2))///edge cases
	{
	  substitution_probs['A']=1.0/3.0; substitution_probs['C']=1.0/3.0; substitution_probs['G']=1.0/3.0; substitution_probs['T']=1.0/3.0; 
	  substitution_probs[sequence[mutate_position_i]]=0;
	}
      else
	{
	  string fivemer_to_mutate=sequence.substr(mutate_position_i-2,5);
	  substitution_probs=S5F_model[fivemer_to_mutate].substitutions;
	}
      
      //cerr << "A: " << substitution_probs['A'] << " C: " << substitution_probs['C'] << " G: " << substitution_probs['G'] << " T: " <<  substitution_probs['T'] << "\t"; 
      double R2=dis(gen);
      double cuml=0;
      char base_to_mutate_to='X';
      for(map<char,double>::iterator it = substitution_probs.begin(); it != substitution_probs.end(); ++it)
	{
	  if (R2<(cuml+it->second))
	    {
	      base_to_mutate_to=it->first;
	      break;
	    }
	  cuml+=it->second;
	}
      //cerr << accrued_mutations << "\t" << mutate_position_i << "\t" << sequence[mutate_position_i];
      string sequence_copy=sequence;
      sequence[mutate_position_i]=base_to_mutate_to;

      //cerr << R2 << "\t" << base_to_mutate_to << "\n"; 
      //cerr << "->" << base_to_mutate_to << "\n"; 
      ///store 
      if (dna_sequence_has_stop_codon_in_reading_frame(sequence))
	{sequence=sequence_copy; continue;}///discard if hits a stop codon and start from prev sequence
      else
	{
	  mutant_sequences.push_back(sequence);
	}

      if (dont_count_remutation)//mutations at the same position don't count towards total number of mutations (obs mut freq)
	{
	  number_of_mutations_two_seqs(ref_sequence, sequence, accrued_mutations);
	}
      else
	{
	  accrued_mutations++;
	}
    }

  return;
}


void number_of_mutations_two_seqs(string &s1, string &s2, int &mutation_count)
{
  ///Assumes sequences are already properly aligned
  assert(s1.length()==s2.length());

  mutation_count=0;
  for(int i=0; i<s1.size(); i++)
    {
      if ((s1[i] == '-') || (s2[i] == '-')){continue;} ///Not counting gaps as mutations currently
      if (s1[i] != s2[i]) {mutation_count++;}
    }
  return;
}


void process_fasta_sequence_to_seq_vector(string &sequence,vector<Seq> &seq_vector, map<string,string> &dna_to_aa_map, map<string,S5F_mut> &S5F_5mers)
{
  seq_vector.clear();
  string aa="X";
  int aa_counter=0;
  for(int i=0; i<sequence.length(); i++)
    {
      Seq temp;
      temp.base=sequence.substr(i,1);
      //get amino acid
      if (i%3==0)
	{
	  //how to deal with gaps -> gaps should have to be 3mers in order for alignment to be valid between two fxnl sequences, check for this
	  string codon=sequence.substr(i,3);
	  if (dna_to_aa_map.find(codon) != dna_to_aa_map.end())
	    aa=dna_to_aa_map[codon];
	  else
	    aa="X";
	  aa_counter++;
	}
      //get S5F fivemer mutability score
      double mut_score;
      if ((i<2)|| (i>sequence.length()-2))
	{
	  mut_score=-1;
	}
      else
	{
	   string fivemer=sequence.substr(i-2,5);
	   mut_score=S5F_5mers[fivemer].score;
	}
      temp.aa=aa;
      temp.aa_num=aa_counter;
      temp.S5F_mut_score=mut_score;
      seq_vector.push_back(temp);
    }
  return;
}


template <typename Type>
void vector1D_to_2D(vector<Type> &vector1D, int interval, vector<vector<Type> > &vector2D)
{
  vector2D.clear();///start fresh
  assert(interval<=vector1D.size());
  vector<Type> row;
  for(int i=0; i<vector1D.size(); i++)
    {
      row.push_back(vector1D[i]);
      if (((i+1)%interval==0)||(i==vector1D.size()-1))
	{
	  vector2D.push_back(row);
	  row.clear();
	}
    }
  return;
}

template <typename Type>
void vector2D_to_3D(vector<vector<Type> > &vector2D, int interval, vector<vector<vector<Type> > > &vector3D)
{
  vector3D.clear();
  // cerr << vector2D.size() << " by " << vector2D[0].size() << "\n"; 
  //slice each row by interval, and stack into one big 2D vector
  vector<vector<Type> > all_sliced_rows;
  int max_num_splits=0;
  for(int i=0; i<vector2D.size(); i++)
    {
      vector<vector<Type> > dim2;
      vector1D_to_2D(vector2D[i],interval,dim2);
      for(int j=0; j<dim2.size(); j++)
	{
	  all_sliced_rows.push_back(dim2[j]);
	}
      if (dim2.size() > max_num_splits){max_num_splits=dim2.size();}
    }
  //cerr << "asr: " <<  all_sliced_rows.size() << "\n"; 
  //re-apportion the slices in the right order to a 3D vector
  for(int i=0; i<max_num_splits; i++)
    {
      //cerr << i << "\n";
      vector<vector<Type> > temp;
      for(int j=0; j<vector2D.size(); j++)
	{
	  temp.push_back(all_sliced_rows[i+(j*max_num_splits)]);
	  //cerr << " " << j << "\n"; 
	}
      vector3D.push_back(temp);
    }
  //cerr << vector3D.size() << " by " << vector3D[0].size() << " by " << vector3D[0][0].size() << "\n"; 
  return;
}

void translate_dna_to_aa(string &dna, string &aa, int reading_frame, map<string,string> &dna_to_aa_tranx_map)
{
  for(int i=reading_frame-1; i<dna.size(); i+=3)
    {
      string codon=dna.substr(i,3);
      boost::to_upper(codon);
      if (dna_to_aa_tranx_map.find(codon)!=dna_to_aa_tranx_map.end())
	{aa+=dna_to_aa_tranx_map[codon];}
      else
	{aa+="X";}

    }
  return;
}


void load_S5F_files(string mutability_filename, string substitution_filename, map<string, S5F_mut> &S5F_5mers)
{
  //Open mutability CSV and read in
  ifstream file1(mutability_filename.c_str(), std::ios::in );
  if (!file1.is_open()) {cerr << "could not open " << mutability_filename << " ...exiting...\n"; exit(1);}

  string file1_str;
  int counter=0;
  while (!getline(file1, file1_str).eof())
    {
      if (counter==0){counter++; continue;}//skip header line

      chomp(file1_str);
      vector<string> tokens; 
      tokenize(file1_str, tokens," ");
      for(int i=0; i<tokens.size(); i++){boost::replace_all(tokens[i],"\"","");}
      S5F_mut temp(tokens[0],atof(tokens[1].c_str()),tokens[2],atof(tokens[3].c_str()),atof(tokens[4].c_str()));
      S5F_5mers[tokens[0]]=temp;
      counter++;
    }
  
  //Open substitution CSV and read in
  ifstream file2(substitution_filename.c_str(), std::ios::in );
  if (!file2.is_open()) {cerr << "could not open " << mutability_filename << " ...exiting...\n"; exit(1);}

  string file2_str;
  counter=0;
  while (!getline(file2, file2_str).eof())
    {
      if (counter==0){counter++; continue;}//skip header line

      chomp(file2_str);
      vector<string> tokens; 
      tokenize(file2_str, tokens," ");
      for(int i=0; i<tokens.size(); i++){boost::replace_all(tokens[i],"\"","");}
      if (S5F_5mers.find(tokens[0]) != S5F_5mers.end())
	{
	  S5F_5mers[tokens[0]].substitutions['A']=atof(tokens[1].c_str());
	  S5F_5mers[tokens[0]].substitutions['C']=atof(tokens[2].c_str());
	  S5F_5mers[tokens[0]].substitutions['G']=atof(tokens[3].c_str());
	  S5F_5mers[tokens[0]].substitutions['T']=atof(tokens[4].c_str());
	  S5F_5mers[tokens[0]].subst_group=tokens[5];
	}
      else{cerr << "ERROR parsing substitution.csv.  No fivemer " << tokens[0] << " found in mutability scores\n"; exit(1);}
      counter++;
    }
  
  return;
}
void read_fasta_file(string filename, map<string, string> &sequence_hash, vector<string> &sequence_names)
{
   ///OPEN SEQUENCE FILE AND PUT INTO MAP
  ifstream file(filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "could not open " << filename << " ...exiting...\n"; exit(1);}

  string file_str;
  string sequence_string, name;
  int name_count=0, seq_count=0;
  map<string, int> duplicate_names;
  while (!getline(file, file_str).eof())
    {
      chomp(file_str);
      if (file_str.substr(0,1) == ">")
	{
	  name=file_str.substr(1,file_str.size());
	  sequence_names.push_back(name);
	  if (duplicate_names.count(name)==0){duplicate_names[name]=1;}else{cerr << "ERROR: duplicate names of sequences are not allowed (" << name << ")...exiting.\n"; }
	  name_count++;
	}
      else
	{
	  sequence_string=file_str.substr(0,file_str.size());
	  boost::to_upper(sequence_string);
	  sequence_hash[name]+=sequence_string;
	  seq_count++;
	}
    }
  return;
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

void get_aa_tranx_map(map<string,string> &dna_to_aa_tranx_map)
{
  ///ALA
  dna_to_aa_tranx_map["GCT"]="A";  
  dna_to_aa_tranx_map["GCC"]="A";  
  dna_to_aa_tranx_map["GCA"]="A";  
  dna_to_aa_tranx_map["GCG"]="A";  
  ///ARG
  dna_to_aa_tranx_map["CGT"]="R";  
  dna_to_aa_tranx_map["CGC"]="R";  
  dna_to_aa_tranx_map["CGA"]="R";  
  dna_to_aa_tranx_map["CGG"]="R"; 
  dna_to_aa_tranx_map["AGA"]="R"; 
  dna_to_aa_tranx_map["AGG"]="R"; 
  ///ASN
  dna_to_aa_tranx_map["AAT"]="N"; 
  dna_to_aa_tranx_map["AAC"]="N"; 
  ///ASP
  dna_to_aa_tranx_map["GAT"]="D"; 
  dna_to_aa_tranx_map["GAC"]="D"; 
  ///CYS
  dna_to_aa_tranx_map["TGT"]="C"; 
  dna_to_aa_tranx_map["TGC"]="C"; 
  ///GLN
  dna_to_aa_tranx_map["CAA"]="Q"; 
  dna_to_aa_tranx_map["CAG"]="Q"; 
  ///GLU
  dna_to_aa_tranx_map["GAA"]="E"; 
  dna_to_aa_tranx_map["GAG"]="E"; 
  ///GLY
  dna_to_aa_tranx_map["GGT"]="G"; 
  dna_to_aa_tranx_map["GGC"]="G"; 
  dna_to_aa_tranx_map["GGA"]="G"; 
  dna_to_aa_tranx_map["GGG"]="G"; 
  ///HIS
  dna_to_aa_tranx_map["CAT"]="H"; 
  dna_to_aa_tranx_map["CAC"]="H"; 
  ///ILE
  dna_to_aa_tranx_map["ATT"]="I"; 
  dna_to_aa_tranx_map["ATC"]="I"; 
  dna_to_aa_tranx_map["ATA"]="I"; 
  ///LEU
  dna_to_aa_tranx_map["TTA"]="L"; 
  dna_to_aa_tranx_map["TTG"]="L"; 
  dna_to_aa_tranx_map["CTT"]="L"; 
  dna_to_aa_tranx_map["CTC"]="L"; 
  dna_to_aa_tranx_map["CTA"]="L"; 
  dna_to_aa_tranx_map["CTG"]="L"; 
  ///LYS
  dna_to_aa_tranx_map["AAA"]="K"; 
  dna_to_aa_tranx_map["AAG"]="K"; 
  ///MET
  dna_to_aa_tranx_map["ATG"]="M"; 
  ///PHE
  dna_to_aa_tranx_map["TTT"]="F"; 
  dna_to_aa_tranx_map["TTC"]="F"; 
  ///PRO
  dna_to_aa_tranx_map["CCT"]="P"; 
  dna_to_aa_tranx_map["CCC"]="P"; 
  dna_to_aa_tranx_map["CCA"]="P"; 
  dna_to_aa_tranx_map["CCG"]="P"; 
  //SER
  dna_to_aa_tranx_map["TCT"]="S"; 
  dna_to_aa_tranx_map["TCC"]="S"; 
  dna_to_aa_tranx_map["TCA"]="S"; 
  dna_to_aa_tranx_map["TCG"]="S"; 
  dna_to_aa_tranx_map["AGT"]="S"; 
  dna_to_aa_tranx_map["AGC"]="S"; 
  //THR
  dna_to_aa_tranx_map["ACT"]="T"; 
  dna_to_aa_tranx_map["ACC"]="T"; 
  dna_to_aa_tranx_map["ACA"]="T"; 
  dna_to_aa_tranx_map["ACG"]="T"; 
  ///TRP
  dna_to_aa_tranx_map["TGG"]="W"; 
  //TYR
  dna_to_aa_tranx_map["TAT"]="Y"; 
  dna_to_aa_tranx_map["TAC"]="Y"; 
  ///VAL
  dna_to_aa_tranx_map["GTT"]="V"; 
  dna_to_aa_tranx_map["GTC"]="V"; 
  dna_to_aa_tranx_map["GTA"]="V"; 
  dna_to_aa_tranx_map["GTG"]="V"; 
  ///STOP
  dna_to_aa_tranx_map["TAA"]="*"; 
  dna_to_aa_tranx_map["TGA"]="*"; 
  dna_to_aa_tranx_map["TAG"]="*"; 
  ///GAP
  dna_to_aa_tranx_map["---"]="-";
 
}


template <typename Type>
string convert_to_string(Type t)
{
  ostringstream convert;
  convert << t;
  return convert.str();
}

bool dna_sequence_has_stop_codon_in_reading_frame(string sequence)
{
  for(int i=0; i<sequence.length(); i+=3)
    {
      if (sequence[i]=='T')
	{
	  if (sequence[i+1] == 'G')
	    {
	      if (sequence[i+2] == 'A'){ return true;}
	    }
	  if (sequence[i+1] == 'A')
	    {
	      if ((sequence[i+2] == 'A') || (sequence[i+2] == 'G')){return true;}
	    }
	}
    }
  return false;
}

vector < vector <string> > read_delimited_file(string filename, string delimiter)
{
  vector < vector <string> > file_contents;

    ifstream file(filename.c_str(), ios::in );

    string line_str;
    while (!getline(file, line_str).eof())
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

void mean_std(vector<double> &list, double &mean, double &std)
{
  int n=list.size();
  double total=0;
  for(int i=0; i<n; i++)
    {
      total+=list[i];
    }
  mean=total/(double) n;
  double SSE=0;
  for(int i=0; i<n; i++)
    {
      SSE+=pow((list[i]-mean),2);
    }
  
  std=sqrt(SSE/((double) (n-1)));
}
