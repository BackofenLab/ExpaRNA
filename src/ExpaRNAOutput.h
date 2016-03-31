/*
 * ExpaRNAOutput.h
 *
 *  Created on: 10.03.2009
 *      Author: heyne
 */

#ifndef EXPARNAOUTPUT_H_
#define EXPARNAOUTPUT_H_


#include "ExpaRNAEnsemble.h"
#include <sys/stat.h>


class ExpaRNAOutput: protected ExpaRNAEnsemble {

public:
	ExpaRNAOutput	();
	ExpaRNAOutput	(const ExpaRNAEnsemble& myEns);

	virtual ~ExpaRNAOutput();

	void 	LCSEPMtoPS();
	void 	output_Fasta();
	void	output_LCSEPM();
	void	output_allEPMs();
	void	output_Clustal();
	void	output_LocARNA();
	void	compareAlign();

private:
			Molecule 	myMol1;
			Molecule 	myMol2;

	void			MapToPS(const string&			myId,
					const int&      		mySize,
					const PatternPairMap& 	myMap);
};

#endif /* EXPARNAOUTPUT_H_ */
