#ifndef MCS_IO_H_
#define MCS_IO_H_

#include<iomanip>
#include<sys/time.h>
#include<fstream>
#include"ExpaRNAinput.h"
#include"ExpaRNAdata.h"
#include"EPMFinding.h"
#include"LCSEPM.h"

//**************************************************************************
//**************************************************************************


   typedef __gnu_cxx::hash_map<string,Molecule,StringHash,StringEq>           MoleculeMAP;
   typedef __gnu_cxx::hash_map<string,Molecule,StringHash,StringEq>::iterator MolMapITER;

   //typedef unordered_map<string,Molecule,StringHash,StringEq>           MoleculeMAP;
   //typedef unordered_map<string,Molecule,StringHash,StringEq>::iterator MolMapITER;


class ExpaRNAEnsemble
{
public:

                  	ExpaRNAEnsemble()
                  	{
                  		load_input_data(ensOptions.cmdLineValues);
                  	};

                  	ExpaRNAEnsemble(const ExpaRNAEnsemble& myEns)
											: ensemblePatterns(myEns.ensemblePatterns),
											ensembleMolecules(myEns.ensembleMolecules),
											ensembleLCSEPM(myEns.ensembleLCSEPM),
											sizeLCSEPM(myEns.sizeLCSEPM),
											scoreLCSEPM(myEns.scoreLCSEPM) {};

	virtual			~ExpaRNAEnsemble() {};

   			void	startALL();

   	static ExpaRNAOptions	  					ensOptions;

protected:

			PatternPairMap 						ensemblePatterns;
			MoleculeMAP    						ensembleMolecules;
			PatternPairMap						ensembleLCSEPM;

			int 								sizeLCSEPM;
			int 								scoreLCSEPM;

			void        clearEnsemble();
			void		load_input_data(const vector<string>& mydata);
			void        startLCSEPM();
};

#endif /*MCS_IO_H_*/
