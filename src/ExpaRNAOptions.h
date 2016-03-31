#ifndef MCSOPTIONS_H_
#define MCSOPTIONS_H_

#include <vector>
#include <string>

class ExpaRNAOptions
{
 public:
	ExpaRNAOptions()
		{
			VERSION_ExpaRNA				="0.8.3";
			EPMscoring				= BY_SIZE;
			minPatternSize 				= 2;
			maxPatternNumber 			= 0;
			calcLCSEPM				= true;
			alignColor				= false;
  			verboseOut				= false;
  			alignOut				= false;
  			locarnaOut				= false;
  			epmOut					= false;
  			allOut					= false;
  			psOut 					= true;
  			ignoreGaps				= true;
  			k_crossing				= false;
  			align_file				= "LCSEPM_align.aln";
  			locarna_file				= "LCSEPM_LocARNA_input.fa";
  			epm_file				= "LCSEPM.epm";
  			all_file				= "allEPM.epm";
  			out_dir					= "./ExpaRNA-results";
  			cmdLineValues.clear();
		};

	ExpaRNAOptions(const ExpaRNAOptions& myOptions): cmdLineValues(myOptions.cmdLineValues),
										mol1Id(myOptions.mol1Id),
										mol2Id(myOptions.mol2Id),
										VERSION_ExpaRNA(myOptions.VERSION_ExpaRNA),
										EPMscoring(myOptions.EPMscoring),
										minPatternSize(myOptions.minPatternSize),
										maxPatternNumber(myOptions.maxPatternNumber),
										calcLCSEPM(myOptions.calcLCSEPM),
										alignColor(myOptions.alignColor),
										verboseOut(myOptions.verboseOut),
										alignOut(myOptions.alignOut),
										locarnaOut(myOptions.locarnaOut),
										epmOut(myOptions.epmOut),
										allOut(myOptions.allOut),
										psOut(myOptions.psOut),
										ignoreGaps(myOptions.ignoreGaps),
										k_crossing(myOptions.k_crossing),
										align_file(myOptions.align_file),
										locarna_file(myOptions.locarna_file),
										epm_file(myOptions.epm_file),
										all_file(myOptions.all_file),
										out_dir(myOptions.out_dir) {};




	virtual ~ExpaRNAOptions();

		std::vector<std::string>	cmdLineValues;

		std::string			mol1Id;
		std::string 			mol2Id;
		std::string			VERSION_ExpaRNA;

		enum 				EPMscoringType	{BY_SIZE, BY_QUADSIZE };

		EPMscoringType			EPMscoring;
		int				minPatternSize;
		int				maxPatternNumber;
  		bool				calcLCSEPM;
  		bool				alignColor;
  		bool				verboseOut;
  		bool				alignOut;
  		bool				locarnaOut;
  		bool				epmOut;
  		bool				allOut;
  		bool				psOut;
  		bool				ignoreGaps;
  		bool				k_crossing;
  		std::string			align_file;
  		std::string			locarna_file;
  		std::string			epm_file;
  		std::string			all_file;
  		std::string			out_dir;
};

#endif /*MCSOPTIONS_H_*/
