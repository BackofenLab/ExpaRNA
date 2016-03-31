#ifndef FIND_PAIR_MCS_H_
#define FIND_PAIR_MCS_H_

#include<string>
#include<vector>
#include<algorithm>
#include<utility>
#include"ExpaRNAdata.h"
#include"Util_String.h"
#include"ExpaRNAOptions.h"

//--------------------------------------------------------------------------
// class PatternMapValue
// value-type of the pattern-hash-Maps
//--------------------------------------------------------------------------
class EPMFinding
{

   public:
						EPMFinding(Molecule&  myMol1,Molecule&  myMol2,PatternPairMap& ensMap, ExpaRNAOptions& myOptions);
	  virtual			~EPMFinding() {};
	  void				start();

   private:

      const Molecule&   mol1,mol2;
      PatternPairMap&	mcsPatterns;
      ExpaRNAOptions&	mcsOptions;

      vector<int>    	check_str(string my_str);
      void           	resizePatternMapCount(unsigned int count);

      //----------------------------------------------------------------------------
      // definitions for the mcs-algorithm
      //----------------------------------------------------------------------------
      void           	calculateMCS();
      void           	calculatePatterns(int minPatSize);

      vector<vector<int> > Mnb; // not bond-Matrix
      vector<vector<int> > Mb ; // bond-Matrix
      vector<vector<int> > Ml ; // loop-Matrix
      vector<vector<int> > Mu ; // used-Matrix/considered positions:
                             //                -1 : not considered
                             //                -2 : considered
                             //                -3 : Mu[i][j]!=0, but taken by overlapping case
      string seq1;
      string seq2;
      string str1;
      string str2;
      vector<int> bplist1;
      vector<int> bplist2;

      int max(int a,int b);
      int get_left_loop_pos(int pos,int offset,int seqnum);
      int get_right_loop_pos(int pos,int offset,int seqnum);
      int get_loop_size(int pos,int seqnum);
      int common_loop_right_max_extended(int starti,int startj);
      int common_loop_left_max_extended(int starti,int startj);
      int get_num_matchings(int i,int j,int size,bool modifying);
      void calc_size_inner_loop(int i,int j);
      void calc_size_bond(int i,int i_,int j,int j_);
      void calc_size_not_bond(int i,int i_,int j,int j_);
      void find_pat_seq(int i,int j,int length,intVec& pattern1,intVec& pattern2);
      void find_pat_bond(int i,int j,intVec& pattern1,intVec& pattern2);
      void print_M(); //helper to print the arrays
};

#endif /*FIND_PAIR_MCS_H_*/
