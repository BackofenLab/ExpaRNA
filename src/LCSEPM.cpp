#include "LCSEPM.h"
#include "Util_String.h"

LCSEPM::~LCSEPM()
{
}

void    LCSEPM::calculateLCSEPM()
{
	struct timeval *Tps, *Tpf;
	struct timezone *Tzp;
	Tps = (struct timeval*) malloc(sizeof(struct timeval));
	Tpf = (struct timeval*) malloc(sizeof(struct timeval));
	Tzp = 0;

	gettimeofday (Tps, Tzp);

	LCSEPMscore = 0;
	LCSEPMsize  = 0;

	preProcessing();

	gettimeofday (Tpf, Tzp);
	long double time = ((Tpf->tv_sec-Tps->tv_sec)*1000000 + Tpf->tv_usec-Tps->tv_usec)/1000000.0;
	if (lcsOptions.verboseOut) cout << "   time for precomputing secs: " << setprecision(4) << time << endl;

	calculateHoles();

	vector < vector<int> > last_vec;
	int i,k;
	i = 1;
	k = 1;
	LCSEPMscore = D_rec(i,mol1.getLength(),k,mol2.getLength(),last_vec,false);
	cout << "    Score LCS-EPM: "<< LCSEPMscore <<endl;
	calculateTraceback(i,mol1.getLength(),k,mol2.getLength(),last_vec);
	LCSEPMsize = matchedEPMs.getMapBases();
	cout << "    #EPMs: "<< matchedEPMs.size() << " / matched Bases: "<< LCSEPMsize <<endl;
}

void    LCSEPM::calculatePatternBoundaries(PatternPair*   myPair)
{
   const vector<unsigned int>& myPatStr1 = myPair->getFirstPat().getPat();
   const vector<unsigned int>& myPatStr2 = myPair->getSecPat().getPat();

   myPair->resetBounds();

   //if (lcsOptions.k_crossing == false){
	   // determine all boundaries
	   // between two pattern positions have to be space for another pattern
	   for (unsigned int k=1;k < (myPatStr1.size());++k)
	   {
	   if ((myPatStr1[k]-minEPMSize > myPatStr1[k-1])
    		  &&(myPatStr2[k]-minEPMSize > myPatStr2[k-1]))
		  {
    	  myPair->addInsideBounds(std::make_pair(make_pair(myPatStr1[k-1],myPatStr1[k]),make_pair(myPatStr2[k-1],myPatStr2[k])));
		  }
	   }
	   // insert global min/max of the pattern
	   myPair->setOutsideBounds(make_pair(make_pair(myPatStr1.front(),myPatStr1.back()),make_pair(myPatStr2.front(),myPatStr2.back())));
 /*  }else{
	   // EXPERIMENTAL - NOT FINISHED!
	   vector<intPair> mol1Holes,mol2Holes;
	   for (unsigned int k=1;k < (myPatStr1.size());++k)
	   {
		   if (myPatStr1[k]-minEPMSize > myPatStr1[k-1]){
			   mol1Holes.push_back(make_pair(myPatStr1[k-1],myPatStr1[k]));
		   }
		   if (myPatStr2[k]-minEPMSize > myPatStr2[k-1]){
			   mol2Holes.push_back(make_pair(myPatStr2[k-1],myPatStr2[k]));
		   }
	   }
	   for (int i=0;i<mol1Holes.size();++i){
		   for (int j=i;j<mol2Holes.size();++j){
			   myPair->addInsideBounds(make_pair(mol1Holes[i],mol2Holes[j]));
		   }
	   }
   }*/
 }

void LCSEPM::preProcessing()
{
    EPM_Table.resize(mol1.getLength()+1);
    for (unsigned int i = 0; i < EPM_Table.size();++i)
    	EPM_Table[i].resize(mol2.getLength()+1);

    for (PatternPairMap::patListCITER myPair = patterns.getList().begin(); myPair != patterns.getList().end(); ++myPair)
    {
        calculatePatternBoundaries(*myPair);

        (*myPair)->initEPMScore(lcsOptions.EPMscoring);

        EPM_Table[(*myPair)->getOutsideBounds().first.second][(*myPair)->getOutsideBounds().second.second] = (*myPair);

        for(IntPPairCITER h = (*myPair)->getInsideBounds().begin(); h != (*myPair)->getInsideBounds().end(); ++h)
        {
            HoleKeyPTR myHoleKey = HoleKeyPTR(new HoleKey());
            int holeSize = (*h).first.second - (*h).first.first - 1;
            myHoleKey->bounds = (*h);
            myHoleKey->pattern = (*myPair);
            holeOrdering.insert(make_pair(holeSize,myHoleKey));
        }
    }
}


int LCSEPM::max3(int& a, int& b, int& c)
{
	int tmp = max(a,b);
	return (max(tmp,c));
}

int LCSEPM::D_rec(const int& i,const  int& j,const int& k,const int& l,vector < vector<int> >& D_h,const bool debug)
{

	int 					score_EPM;
	int 					pos_before_EPM_Str1;
	int 					pos_before_EPM_Str2;

	D_h.clear();
	D_h.resize(j - i + 2);
	for (unsigned int a = 0; a < D_h.size();++a)
		D_h[a].resize(l - k + 2,0);

	for(unsigned int j_1 = 1; j_1 < (j-i+2); ++j_1)
		for (unsigned int l_2 = 1; l_2 < (l-k+2); ++l_2)
		{
			if (EPM_Table[i + j_1-1][k + l_2-1] == NULL)
			{
				D_h[j_1][l_2] = max(D_h[j_1-1][l_2],D_h[j_1][l_2-1]);
			}
			else
			{
				if (debug) {
					if (EPM_Table[i + j_1-1][k + l_2-1]->getScore()==3){
					SinglePattern my1 = EPM_Table[i + j_1-1][k + l_2-1]->getFirstPat();
					my1.print();
					my1 = EPM_Table[i + j_1-1][k + l_2-1]->getSecPat();
					my1.print();}
				}
				pos_before_EPM_Str1 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().first.first ) - i;
				pos_before_EPM_Str2 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().second.first ) - k;
				if ((pos_before_EPM_Str1 < 0)||(pos_before_EPM_Str2 <0))
					score_EPM = 0;
				else
					score_EPM = D_h[pos_before_EPM_Str1][pos_before_EPM_Str2] + EPM_Table[i + j_1-1][k + l_2-1]->getScore();

				D_h[j_1][l_2] = max3(score_EPM,D_h[j_1-1][l_2],D_h[j_1][l_2-1]);
			}
		}
	return (D_h[j - i + 1][l - k + 1]);
}

void LCSEPM::calculateHoles()
{
	vector < vector<int> > vec;
	for (HoleMapCITER t = holeOrdering.begin();t != holeOrdering.end();++t)
    {
		bool deb=false;
		(*t).second->pattern->setEPMScore(	(*t).second->pattern->getScore() +
											D_rec((*t).second->bounds.first.first+1,(*t).second->bounds.first.second-1,\
											(*t).second->bounds.second.first+1,(*t).second->bounds.second.second-1,vec,deb));
    }
}

void LCSEPM::calculateTraceback(const int i,const  int j,const int k,const int l,vector < vector<int> > holeVec)
{
	int j_1 = holeVec.size()-1;
	int l_2 = holeVec[0].size()-1;
	while ((j_1 >= 1)&&(l_2 >= 1))
	{
		if (holeVec[j_1 - 1][l_2] == holeVec[j_1][l_2])
			--j_1;
		else
			if (holeVec[j_1][l_2 - 1] == holeVec[j_1][l_2])
				--l_2;
			else
			{
				vector < vector<int> > tmpHoleVec;
				matchedEPMs.add(EPM_Table[i + j_1 - 1][k + l_2 - 1]);
				for(IntPPairCITER h = EPM_Table[i + j_1 - 1][k + l_2 - 1]->getInsideBounds().begin(); h != EPM_Table[i + j_1-1][k + l_2-1]->getInsideBounds().end(); ++h)
				{
					int sc = D_rec((*h).first.first+1,(*h).first.second-1,(*h).second.first+1,(*h).second.second-1,tmpHoleVec,false);
					if (sc != 0)
						calculateTraceback((*h).first.first+1,(*h).first.second-1,(*h).second.first+1,(*h).second.second-1,tmpHoleVec);
				}
				int s1 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().first.first ) - i;
				int s2 = (EPM_Table[i + j_1-1][k + l_2-1]->getOutsideBounds().second.first) - k;
				j_1 = s1;
				l_2 = s2;
			}
	}
}
