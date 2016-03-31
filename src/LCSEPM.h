#ifndef MCSERPAP_H_
#define MCSERPAP_H_

#include"ExpaRNAEnsemble.h"
#include"ExpaRNAOptions.h"
#include<map>
#include<algorithm>
#include<iostream>

class LCSEPM
{
    public:

					LCSEPM(	const		Molecule&       myMol1,
    						const 		Molecule&       myMol2,
    						const 		PatternPairMap& myPatterns,
										PatternPairMap& myLCSEPM,
										int&		    mySize,
										int&			myScore,
    						const		ExpaRNAOptions& myOptions)

    	                                   :matchedEPMs(myLCSEPM),
    	                                	mol1(myMol1),mol2(myMol2),
    	                                    patterns(myPatterns),
    	                                    LCSEPMsize(mySize),
    	                                    LCSEPMscore(myScore),
    	                                    minEPMSize(myOptions.minPatternSize),
    	                                    lcsOptions(myOptions) {};

        virtual 		~LCSEPM();

        		void    calculateLCSEPM();

	private:

		struct	HoleKeyS;
        struct	HoleKeyS
        {
        	PatternPairMap::SelfValuePTR    pattern;
            intPPair                        bounds;
        };

        typedef     HoleKeyS                            HoleKey;
        typedef     HoleKey*                            HoleKeyPTR;
        typedef     multimap<int,HoleKeyPTR>            HoleOrderingMapTYPE;
        typedef     HoleOrderingMapTYPE::const_iterator HoleMapCITER;

        void    preProcessing				();
        void    calculateHoles				();
	    void    calculatePatternBoundaries	(PatternPair* myPair);
        void 	calculateTraceback			(const int i,const int j,const int k,const int l,vector < vector<int> > holeVec);
        int 	D_rec						(const int& i,const  int& j,const int& k,const int& l,vector < vector<int> >& D_h,const bool debug);
        int 	max3						(int& a, int& b, int& c);

				vector< vector<PatternPairMap::SelfValuePTR> >	EPM_Table;
				PatternPairMap&									matchedEPMs;
				HoleOrderingMapTYPE    					 		holeOrdering;
        const 	Molecule&					               		mol1,mol2;
        const 	PatternPairMap&         						patterns;
				int& 											LCSEPMsize;
				int& 											LCSEPMscore;
		const	int                     						minEPMSize;
				ExpaRNAOptions									lcsOptions;
};
#endif /*MCSERPAP_H_*/
