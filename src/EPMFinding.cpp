#include "EPMFinding.h"

EPMFinding::EPMFinding(Molecule&  myMol1,Molecule&  myMol2,PatternPairMap& ensMap, ExpaRNAOptions& myOptions)
   :mol1(myMol1),mol2(myMol2),mcsPatterns(ensMap),mcsOptions(myOptions)
{
   seq1="U"+myMol1.getSequence()+"V";
   seq2="X"+myMol2.getSequence()+"Y";
   str1="("+myMol1.getStructure()+")";
   str2="("+myMol2.getStructure()+")";
   /*
   seq1="UAAAAAACCCCCCCCCCCCUUUUUUUV";
   seq2="XGGGGAACCUUUCUCUCUCUUAAAAAY";
   str1="(....((............)).....)";
   str2="(....((............)).....)";
   */
   bplist1=check_str(str1);
   bplist2=check_str(str2);

   Mnb.resize(seq1.size());
   Mb.resize(seq1.size());
   Ml.resize(seq1.size());
   Mu.resize(seq1.size());

   for(unsigned int i=0;i<Mnb.size();i++)
   {
      Mnb[i].resize(seq2.size(),-1);
      Mb[i].resize(seq2.size(),-1);
      Ml[i].resize(seq2.size(),-1);
      Mu[i].resize(seq2.size(),-1);
   }
}

vector<int> EPMFinding::check_str(string my_str)
{
	vector<int> temp_bplist;
	stack<int> structure;
	temp_bplist.resize(my_str.size());
    for(unsigned int i=0;i<my_str.size();i++)
    		if (my_str[i]=='(')
        		structure.push(i);
        else
        		if (my_str[i]==')')
          	{
	    			if (structure.empty())
	      		{
					cerr << "...structure not correct. Please use a correct structure!" << endl;
					exit(1);
	      		}
            		temp_bplist[i]=structure.top(); // index basepair partner
            		temp_bplist[structure.top()]=i; // index basepair partner
            		structure.pop();
          	}
        		else
	  			temp_bplist[i]=-1;      // not paired
	if (!structure.empty())
	{
		cerr << "...structure not correct. Please use a correct structure!" << endl;
	  	exit(1);
	}
	return (temp_bplist);
}

void  EPMFinding::resizePatternMapCount(unsigned int count)
{
   mcsPatterns.makeOrderedMap();
   if (mcsOptions.verboseOut) cout << "Hold only the biggest " << count << " EPMs...";
   PatternPairMap::orderedMapITER I = mcsPatterns.getOrderedMap2().begin();
   unsigned int i=1;
   PatternPairMap::patListTYPE tempList;
   while ((i<=count) && (I!=mcsPatterns.getOrderedMap().end()))
      {
         I++;
         i++;
      }
   mcsPatterns.getOrderedMap2().erase(I,mcsPatterns.getOrderedMap2().end());
   mcsPatterns.updateFromMap();
   if (mcsOptions.verboseOut) cout << "done!" <<endl;
}

void EPMFinding::start()
{
   cout << endl << " -->> starting EPM finding algorithm" << endl;
   calculateMCS();
   cout << "Count EPMs with minimal size gamma = " << mcsOptions.minPatternSize << "...";
   calculatePatterns(mcsOptions.minPatternSize);

   if (mcsOptions.maxPatternNumber > 0) //&& (mcsPatterns.size() > maxCount))
      resizePatternMapCount(mcsOptions.maxPatternNumber);
   if (mcsOptions.verboseOut) cout << "...resting " << mcsPatterns.size() << " EPMs." << endl;
   cout << " -->> end EPM finding algorithm" << endl;
}

void EPMFinding::calculatePatterns(int minPatSize)
{
   SinglePattern pattern1,pattern2;
   intVec pat1vec,pat2vec;
   int count = 0;
   string patId;

   for(unsigned int i=0;i<Ml.size();i++)
      for(unsigned int j=0;j<Ml[i].size();j++)
         {
    	  if ((Ml[i][j]>0 && ( Mu[i][j]==-2 || ( Mb[i][j]>0 && Mu[i][j]>=0 ))) || ((Mnb[i][j]>=minPatSize && Ml[i+1][j+1]<=0) && ((Mu[i][j]==-2) || (Mu[i][j]>=0))) )
          {
        	if ((Ml[i][j]>=minPatSize)||(Mnb[i][j]>=minPatSize))
            {
               pat1vec.clear();  	pat2vec.clear();
               find_pat_seq(i,j,common_loop_right_max_extended(i,j),pat1vec,pat2vec);
               sort(pat1vec.begin(),pat1vec.end());
               sort(pat2vec.begin(),pat2vec.end());
               // new bug found!!! hack:
               if (pat1vec.size() >= minPatSize) {
            	   count ++;
            	   patId = "pat_" + biu::Util_String::int2str(count);
            	   pattern1 = SinglePattern(mol1.getId(),mol1.getActualStructureId(),patId,pat1vec);
            	   pattern2 = SinglePattern(mol2.getId(),mol2.getActualStructureId(),patId,pat2vec);
            	   mcsPatterns.add(patId,pat1vec.size(),pattern1,pattern2);
               }//else{
            	   //cout << "size:" << pat1vec.size() << " Ml:" << Ml[i][j] << " Mnb:" << Mnb[i][j] << "Mu:" << Mu[i][j] << "Mb:" << Mb[i][j]<< endl;
               //}
            }
         }
      }
   cout << " Finished! Found: " << count << endl;
}

void EPMFinding::calculateMCS()
{
	if (mcsOptions.verboseOut) if (mcsOptions.verboseOut) cout << "Calculating..." << flush;
	// from inner to outer loops
	for(unsigned int i=0;i<bplist1.size();i++)
    		if ( bplist1[i]!=-1 && bplist1[i]<(int)i )           // basepair : (bplist1[i],i)
      		for(unsigned int j=0;j<bplist2.size();j++)
				if ( bplist2[j]!=-1 && bplist2[j]<(int)j )   // basepair : (bplist2[j],j)
	  			{
	    				calc_size_inner_loop(bplist1[i]+1,bplist2[j]+1);
	    				if ( seq1[bplist1[i]]==seq2[bplist2[j]] && seq1[i]==seq2[j] )
	      				calc_size_bond(bplist1[i],i,bplist2[j],j);
	    				else
	      				calc_size_not_bond(bplist1[i],i,bplist2[j],j);
	  			}
	if (mcsOptions.verboseOut) cout << " Finished!" << endl;
	//print_M();
}


//---------------------------------------------------------------------------------------
// mcs algorithm
//---------------------------------------------------------------------------------------

void EPMFinding::print_M()
{
  int m=0;
  int i_=-1;
  int j_=-1;
  cout << "Mnb :" << endl;
  for(unsigned int i=0;i<Mnb.size();i++)
	  cout << seq1[i] << "\t" ;
  cout << endl;
  for(unsigned int i=0;i<Mnb.size();i++)
    {
      for(unsigned int j=0;j<Mnb[i].size();j++)
	{
	  if (m<Mnb[i][j])
	    {
	      m=Mnb[i][j];
	      i_=i;
	      j_=j;
	    }
	  cout << Mnb[i][j] << "\t" ;
	}
      cout << endl;
    }

  cout << "Mb :" << endl;
  for(unsigned int i=0;i<Mb.size();i++)
    {
      for(unsigned int j=0;j<Mb[i].size();j++)
	{
	  if (m<Mb[i][j])
	    {
	      m=Mb[i][j];
	      i_=i;
	      j_=j;
	    }
      	  cout << Mb[i][j] << "\t" ;
	}
      cout << endl;
    }

  cout << "Ml :" << endl;
  for(unsigned int i=0;i<Ml.size();i++)
    {
      for(unsigned int j=0;j<Ml[i].size();j++)
	{
	  if (m<Ml[i][j])
	    {
	      m=Ml[i][j];
	      i_=i;
	      j_=j;
	    }
	  cout << Ml[i][j] << "\t" ;
	}
      cout << endl;
    }

  cout << "Mu :" << endl;
  for(unsigned int i=0;i<Mu.size();i++)
    {
      for(unsigned int j=0;j<Mu[i].size();j++)
	cout << Mu[i][j] << "\t" ;
      cout << endl;
    }

  cout << "max (" << i_ << "," << j_ << ") : " << m << endl;
}

int EPMFinding::max(int a,int b)
{
  return (a>b?a:b);
}

int EPMFinding::get_left_loop_pos(int pos,int offset,int seqnum)
{
  if(seqnum==1)
    {
      while (offset>0)
        {
          if (str1[pos]==')')
            {
              pos=bplist1[pos];
            }
          else
            pos--;
          offset--;
        }
    }
  else if(seqnum==2)
    {
      while (offset>0)
        {
          if (str2[pos]==')')
            {
              pos=bplist2[pos];
            }
          else
            pos--;
          offset--;
        }
    }

  return pos;
}
//wandert offset-positionen im loop entlang, neue bonds werden unten übergangen
//und zählen nicht als schritt, nur ungebundene basen zählen
//rückgabe der position der offset-base,kann bond-base sein
int EPMFinding::get_right_loop_pos(int pos,int offset,int seqnum)
{
  if(seqnum==1)
    {
      while (offset>0)
        {
          if (str1[pos]=='(')
            {
              pos=bplist1[pos];
            }
          else
            pos++;
          offset--;
        }
    }
  else if(seqnum==2)
    {
      while (offset>0)
        {
          if (str2[pos]=='(')
            {
              pos=bplist2[pos];
            }
          else
            pos++;
          offset--;
        }
    }
  return pos;
}
//?aufruf immer innerhalb schließender bindung?
//wandert seq nach rechts ab, neue bonds werden übersprungen, keine wanderung im stacking
//solange bis schließende base gefunden wird, size ist größe ohne diese!
int EPMFinding::get_loop_size(int pos,int seqnum)
{
  int size=0;
  if (seqnum==1)
    while (str1[pos]!=')' && pos<(int)seq1.size())
      {
        if (str1[pos]=='(')
          {
            pos=bplist1[pos]+1;
            size+=2;
          }
        else
          {
            pos++;
            size++;
          }
      }
  else if (seqnum==2)
    while (str2[pos]!=')' && pos<(int)seq2.size())
      {
        if (str2[pos]=='(')
          {
            pos=bplist2[pos]+1;
            size+=2;
          }
        else
          {
            pos++;
            size++;
          }
      }
  return size;
}
//berechnet maximalen loopmatch ab einer position
//wandert loop entlang bis mismatch(seq&str)
//wenn bond: wandert außen entlang,wenn bondmatch dann +2 für größe des matches
int EPMFinding::common_loop_right_max_extended(int starti,int startj)
{
 //   cout << "proc common_loop_right_max_extended " << starti << " " << startj << endl;
  int i=starti;
  int j=startj;
  int size=0;

  while (seq1[i]==seq2[j] && str1[i]==str2[j] && i<(int)seq1.size() && j<(int)seq2.size())
    {
      if (str1[i]==')' && bplist1[i]<starti && i>starti && j>startj)
        return size;
      if (str1[i]=='(')
        {
          i=bplist1[i];
          j=bplist2[j];
        }
      else
        {
          i++;
          j++;
        }
      size++;
      // cout << "proc common_loop_right_max_extended size: " << size << " i: " << i <<  " j: " << j << endl;
    }
  return size;
}

int EPMFinding::common_loop_left_max_extended(int starti,int startj)
{
	//cout << "proc common_loop_left_max_extended  i: " << starti <<  " j: " << startj << endl;
  int i=starti;
  int j=startj;
  int size=0;

  while (seq1[i]==seq2[j] && str1[i]==str2[j] && i>=0 && j>=0)
    {
      if (str1[i]=='(' && bplist1[i]>starti)
        return size;
      if (str1[i]==')')
        {
          i=bplist1[i];
          j=bplist2[j];
        }
      else
        {
          i--;
          j--;
        }

      size++;
      //cout << "proc common_loop_left_max_extended size: " << size << " i: " << i <<  " j: " << j << endl;
    }
  return size;
}

int EPMFinding::get_num_matchings(int i,int j,int size,bool modifying)
{
  // modifying changes Mu values in inner loops
	// cout << "proc get_num_matchngs " << i << " " << j << " " << size << endl;

  int v=0;
  //  Ml[i][j]=0;
  while (seq1[i]==seq2[j] && str1[i]==str2[j] && size> 0 && i<(int)seq1.size() && j<(int)seq2.size())
    {
      if (str1[i]=='(')
        {
          if (modifying)
	    {
	      Mu[i][j]=-2;
	      Mu[bplist1[i]][bplist2[j]]=-2;
	    }
          if (Mb[i][j]==-1)
            return v+Mnb[i][j];
          v+=Mb[i][j];
          i=bplist1[i]+1;
          j=bplist2[j]+1;
          size-=2;
          //cout << "mmmmm21 " << i << " " << j << " " << v << endl;
        }
      else if ( str1[i]==')' )
        {
          if (modifying)
            Mu[i][j]=-2;
          if ( Mnb[i][j]>0 )
            v+=Mnb[i][j];
          i++;
          j++;
          size--;
          //cout << "mmmmm22 " << i << " " << j << " v: " << v << "size: " << size << endl;
        }
      else
        {
          if (modifying)
            Mu[i][j]=-2;
          v++;
          i++;
          j++;
          size--;
          //cout << "mmmmm23 " << i << " " << j << " " << v << endl;
        }

    }
  return v;
}

void EPMFinding::calc_size_inner_loop(int i,int j)
{
  int loop1_size=get_loop_size(i,1);
  int loop2_size=get_loop_size(j,2);
  //cout << "proc calc_size_inner_loop i " << i << " j " << j << " loop1: " << loop1_size << " loop2: " << loop2_size<< endl;

  for(int k=0;k<loop1_size;k++)
    for(int l=0;l<loop2_size;l++)
      {
        int k_= get_right_loop_pos(i,k,1);
        int l_= get_right_loop_pos(j,l,2);
        if(Mu[k_][l_]==-1)
          {
            int size=common_loop_right_max_extended(k_,l_);
            Ml[k_][l_]=get_num_matchings(k_,l_,size,true);
            Mu[k_][l_]=-2;
            //cout << "int : " << size << "  " << k_ << "  " << l_ << "  " << Ml[k_][l_] << endl;
              }
      }
}

void EPMFinding::calc_size_bond(int i,int i_,int j,int j_)
{
 // cout << "proc calc_size_bond " << i << " " << i_ << " " << j << " " << j_ << endl;
  int loop1_size    = get_loop_size(i+1,1) ;
  int loop2_size    = get_loop_size(j+1,2) ;
  int min_loop_size = loop1_size           ;

  if (min_loop_size>loop2_size){
    min_loop_size=loop2_size;
  }

  int loop_r = common_loop_right_max_extended(i+1,j+1);
  int loop_l = common_loop_left_max_extended(i_-1,j_-1);
  //cout << "loop1_size: " << loop1_size << " loop2_size: " << loop2_size << " loop_r: " << loop_r << " loop_l: " << loop_l << endl;
  int max=0,km=-2,ki=0,kj=0,im=0,jm=0,i_m=0,j_m=0;
  int tmp;

  if (loop_r+loop_l>min_loop_size)    // overlapping case
    {
      for(int k=min_loop_size-loop_l;k<=loop_r;k++)
        {
          im=get_right_loop_pos(i+1,loop1_size-(min_loop_size-k),1);
          jm=get_right_loop_pos(j+1,loop2_size-(min_loop_size-k),2);
          int split=get_num_matchings(i+1,j+1,k,false)+get_num_matchings(im,jm,min_loop_size-k,false);

          if (max<split)
            {
              max=split;
              km=k;
	      ki=im;
	      kj=jm;
            }
	    // cout << "k: " << k << " im : " << im << " jm : " << jm << " split : " << split << endl;
        }
      tmp = get_num_matchings(i+1,j+1,km,true)+get_num_matchings(im,jm,min_loop_size-km,true);
      Mu[i+1][j+1] = km;
      Mu[ki][kj]   = -3;
    }
  else                                 // non-overlapping case
    {
      im=get_right_loop_pos(i+1,loop1_size-loop_l,1);	//orig
      jm=get_right_loop_pos(j+1,loop2_size-loop_l,2);	//orig

      //i_m=get_left_loop_pos(i_-1,loop1_size-loop_l,1);
      //j_m=get_left_loop_pos(j_-1,loop2_size-loop_l,2);

      //int max1=get_num_matchings(i+1,j+1,loop_r,true) + get_num_matchings(i_m,j_m,loop_l,true);

      max=get_num_matchings(i+1,j+1,loop_r,true) + get_num_matchings(im,jm,loop_l,true);

      // cout << "i_m" << i_m << " j_m: " << j_m << " max1: " << max1 << endl;
      Mu[im][jm]   = -3;
      Mu[i+1][j+1] = -3;
      // cout << "non-overlapping " << i << " " << i_ << " " << j << " " << j_ << endl;
    }

  Mb[i][j]=max+2;


/*   cout << "min_loop_size : " << min_loop_size << endl;
   cout << "loop_r (" << i+1 << "," << j+1 << ") : " << loop_r << endl;
   cout << "loop_l (" << i_-1 << "," << j_-1 << ") : " << loop_l << endl;

   cout << "loop1 size (" << i+1 << ") : " << get_loop_size(i+1,1) << endl;
   cout << "loop2 size (" << j+1 << ") : " << get_loop_size(j+1,2) << endl;
   cout << "Mb[" << i << "][" << j << "] = " << Mb[i][j] << endl;
 */
}

void EPMFinding::calc_size_not_bond(int i,int i_,int j,int j_)
{
  //  cout << "proc calc_size_not_bond " << i << " " << i_ << " " << j << " " << j_ << endl;

  if ( seq1[i]==seq2[j] && Ml[i+1][j+1]>=0 )
    {
      Mnb[i  ][j  ] = Ml[i+1][j+1]+1 ;
      Mu [i+1][j+1] = -3;
      Mnb[i_ ][j_ ] = 0;
      //  cout << "calc_size_not_bond matches1 " << i << " and " << j << " : " << " " << Mnb[i][j] << endl;
    }
  else if ( seq1[i_]==seq2[j_] )
    {
      int size=common_loop_left_max_extended(i_-1,j_-1);
      if (size>0)
	{
	  int loop1_pos=get_left_loop_pos(i_-1,size-1,1);
	  int loop2_pos=get_left_loop_pos(j_-1,size-1,2);
	  Mnb[i_][j_] = Ml[loop1_pos][loop2_pos]+1;
	  Mu [loop1_pos][loop2_pos] = -3;
	  Mnb[i ][j ] = 0;
	 //      cout << "calc_size_not_bond matches2 " << i_ << " and " << j_ << " : " << size << " " << loop1_pos << " " << loop2_pos << " " << Mnb[i_][j_] << endl;
	}
      else
	{
	  Mnb[i_][j_] = 1;
	  Mnb[i ][j ] = 0;
	}
    }
  else
    {
      Mnb[i][j] = Mnb[i_][j_] = 0;
    }
}


void EPMFinding::find_pat_seq(int i,int j,int length,intVec& pattern1, intVec& pattern2)
{
	//    cout << "find_pat_seq: i->" << i << " j->" << j << " length : " << length << " seq " << endl;
	//cout << "pat 1: " << biu::Util_String::intvec2str(pattern1,":") << "  pat2: " << biu::Util_String::intvec2str(pattern2,":") << endl;
  if (i<0 || i>=(int)seq1.size() || j<0 || j>=(int)seq2.size())
    return;

  int im,jm;

  for(int l=0;l<length;l++)
    {
      im=get_right_loop_pos(i,l,1);
      jm=get_right_loop_pos(j,l,2);
      //     cout << "find_pat_seq:  l : " << l << " im : " << im << " jm : " << jm << endl;
      if ( str1[im]=='.' && str2[jm]=='.' && seq1[im]==seq2[jm] )
	{
	  pattern1.push_back((unsigned int)(im));
	  pattern2.push_back((unsigned int)(jm));
	  //cout << "find_pat_seq str= . pat1: " << biu::Util_String::intvec2str(pattern1,":") << "  pat2: " << biu::Util_String::intvec2str(pattern2,":") << endl;
	}
      else if ( str1[im]=='(' && str2[jm]=='(')
	{
    	  //cout << "find_pat_seq str = ( before pat1: " << biu::Util_String::intvec2str(pattern1,":") << "  pat2: " << biu::Util_String::intvec2str(pattern2,":") << endl;
      find_pat_bond(im,jm,pattern1,pattern2);
      //cout << "find_pat_seq str = ( after pat1: " << biu::Util_String::intvec2str(pattern1,":") << "  pat2: " << biu::Util_String::intvec2str(pattern2,":") << endl;
      l++;
	}
      else if ( str1[im]==')' && str2[jm]==')')
	{
    	  //cout << "find_pat_seq str = ) before pat1: " << biu::Util_String::intvec2str(pattern1,":") << "  pat2: " << biu::Util_String::intvec2str(pattern2,":") << endl;
    	  find_pat_bond(im,jm,pattern1,pattern2);
    	  //cout << "find_pat_seq str = ) after pat1: " << biu::Util_String::intvec2str(pattern1,":") << "  pat2: " << biu::Util_String::intvec2str(pattern2,":") << endl;
	}
    }
}

void EPMFinding::find_pat_bond(int i,int j,intVec& pattern1,intVec& pattern2)
{
  //cout << "find_pat_bond: i->" << i << " j->" << j << " bond " << endl;
  if (i<0 || i>=(int)seq1.size() || j<0 || j>=(int)seq2.size())
    return;
  if (str1[i]==')' && seq1[i]==seq2[j])        // right bond
    {
      //cout << "find_pat_bond: bne " << bplist1[i] << " " << i << " " << bplist2[j] << " " << j << endl;
      pattern1.push_back(i);
      pattern2.push_back(j);
      if (str1[i-1]==str2[j-1] && seq1[i-1]==seq2[j-1] )
	{
	  int size=common_loop_left_max_extended(i-1,j-1);
	  find_pat_seq(get_left_loop_pos(i-1,size-1,1),get_left_loop_pos(j-1,size-1,2),size,pattern1,pattern2);
	}
    }
  else if (str1[i]=='(' && seq1[i]==seq2[j] && Mnb[i][j]>0)        // left bond
    {
     // cout << "find_pat_bond: bne " << i << " " << bplist1[i] << " " << j << bplist2[j] << endl;
      pattern1.push_back(i);
      pattern2.push_back(j);
      if (str1[i+1]==str2[j+1] && seq1[i+1]==seq2[j+1] )
	{
	  int size=common_loop_right_max_extended(i+1,j+1);
	  find_pat_seq(i+1,j+1,size,pattern1,pattern2);
	}
    }
  else if (Mb[i][j]!=-1)         // bonds are equal
    {
      //cout << "find_pat_bond: beq " << i << " " << bplist1[i] << " " << j << " " << bplist2[j] << endl;
      pattern1.push_back(i);
      pattern1.push_back(bplist1[i]);
      pattern2.push_back(j);
      pattern2.push_back(bplist2[j]);

      if (Mu[i+1][j+1]>=0)          // overlapping
	{
	  //cout << "find_pat_bond: bond round " << Mu[i+1][j+1] << endl;
	  find_pat_seq(i+1,j+1,Mu[i+1][j+1],pattern1,pattern2);

	  int loop1_size=get_loop_size(i+1,1);
	  int loop2_size=get_loop_size(j+1,2);
	  int min_loop_size=loop1_size;
	  if (min_loop_size>loop2_size)
	    min_loop_size=loop2_size;

	  int size=min_loop_size-Mu[i+1][j+1];

	  int pos1=get_left_loop_pos(bplist1[i]-1,size-1,1);
	  int pos2=get_left_loop_pos(bplist2[j]-1,size-1,2);

	/*  cout << "loop1_size  : " << loop1_size << endl;
	  cout << "loop2_size  : " << loop2_size << endl;
	  cout << "minloopsize : " << min_loop_size << endl;
	  cout << "size        : " << size << endl;
	  cout << "pos1        : " << pos1 << endl;
	  cout << "pos2        : " << pos2 << endl;
	  cout << "i,j         : " << i << "," << j << endl;
	  cout << "Mu(pos1,pos2):" << Mu[pos1][pos2] << endl;
	  cout << "Mu(i+1,j+i):" << Mu[i+1][j+1] << endl;*/
	  // before if (Mu[pos1][pos2]!=-3)  !!!
	  if (Mu[pos1][pos2]!=3)/*((Mu[pos1][pos2]==-3)&&(Ml[pos1][pos2]>0))*/
	  {  find_pat_seq(pos1,pos2,size,pattern1,pattern2);}

	}
      else  if ((Mu[i+1][j+1]==-3))              // non-overlapping
	{
	  //cout << "find_pat_bond: inner loop cut - i" << i << " j" << j << endl;
	  int loop_r=common_loop_right_max_extended(i+1,j+1);
	  int loop_l=common_loop_left_max_extended(bplist1[i]-1,bplist2[j]-1);
	  find_pat_seq(i+1,j+1,loop_r,pattern1,pattern2);

	  int loop1_size=get_loop_size(i+1,1);
	  int loop2_size=get_loop_size(j+1,2);
	  int min_loop_size=loop1_size;
	  if (min_loop_size>loop2_size)
	    min_loop_size=loop2_size;

	  int pos1=get_left_loop_pos(bplist1[i]-1,loop_l-1,1);
	  int pos2=get_left_loop_pos(bplist2[j]-1,loop_l-1,2);

	  if (loop_r+loop_l<=min_loop_size)
	    find_pat_seq(pos1,pos2,loop_l,pattern1,pattern2);
	}
      /*
      else
	{
	  find_pat_seq(i+1,j+1,2,pattern1,pattern2);
	  }*/
    }
}

