/*
 * ExpaRNAinput.h
 *
 *  Created on: 13.03.2009
 *      Author: heyne
 */

#ifndef EXPARNAINPUT_H_
#define EXPARNAINPUT_H_

#include"ExpaRNAdata.h"

class ExpaRNAinput {

public:

	ExpaRNAinput();
	virtual ~ExpaRNAinput();

	string 		getFASTALine	(string myStr,int start,char myBreak);
	void 		readFASTA		(string myFASTAfile, Molecule& myMol1, Molecule& myMol2);

};

#endif /* EXPARNAINPUT_H_ */
