///////////////////////////////////////////////////////////////////////////////////////////
// program which finds the longest common subsequence of exact pattern matchings (LCS-EPM)
// EPM-finding algorithm by Rolf Backofen/Sven Siebert (2007)
// LCS-EPM algorithm  by Steffen Heyne/Sebastian Will/Rolf Backofen (2007-2009)
///////////////////////////////////////////////////////////////////////////////////////////

#include"ExpaRNAEnsemble.h"
#include "ExpaRNAOptions.h"

///////////////////////////////////////////////////////////////////////////////////////////
// help screen
///////////////////////////////////////////////////////////////////////////////////////////
void printHelp(std::string myVersion)
{
  cerr << endl;
  cerr << "###############################################################################" << endl;
  cerr << "    ExpaRNA - Exact Pattern Alignment of RNA  - version " << myVersion << endl;
  cerr << "           (c) copyright by Steffen Heyne - 2008-2009" << endl << endl;
  cerr << "  The Longest Common Subsequence of Exact Pattern Matchings (LCS-EPM) algorithm" << endl;
  cerr << "                http://www.bioinf.uni-freiburg.de/Software"<<endl;
  cerr << "###############################################################################"<<endl<<endl;
  cerr << "Usage : 1) ExpaRNA <options> *.fa"                                    << endl;
  cerr << "        2) ExpaRNA <options> *.fa *.align"                            << endl << endl;
  cerr << "   1) two sequences given in fasta format"				   	<< endl;
  cerr << "      if no structure data is available, the sequence is folded"     	<< endl;
  cerr << "      to an mfe structure with Vienna RNAfold."                       	<< endl <<endl;
  cerr << "   2) like 1) and *.align contains reference alignment (use option -A, see below)"	<< endl <<endl;
  cerr << "EPM Selection Options:"	                                                  << endl << endl;
  cerr << "  -s#   : # minmal size of an EPM (gamma)"   					<< endl;
  cerr << "          2 is default and returns all EPMs with at least 2 nucleotides" << endl << endl;;
  cerr << "  -n#   : # determines the maximal number of returned EPMs"       		<< endl;
  cerr << "          0 is default an returns all EPMs"                       		<< endl;
  cerr << "          1 is the largest EPM and so on..."          			<< endl << endl;
  cerr << "EPM Scoring Options:"                                                    << endl << endl;
  cerr << "  -t#   : 1: initial EPM score = EPM size (default)"						<< endl;
  cerr << "        : 2: initial EPM score = (EPM size)^2 (prefers larger patterns in LCS-EPM)"<< endl<<endl;
  cerr << "Input treatment:"                                               			<< endl << endl;
  cerr << "  -g    : Do NOT ignore gaps in input sequences"	                		<< endl << endl;
  cerr << "Output Options:"            						<< endl << endl;
  cerr << "  -d#   : write all output to dir <#>"                		<< endl << endl;
  cerr << "  -o    : write LCS-EPM into file 'LCSEPM_align.aln' as alignment"                		<< endl << endl;
  cerr << "  -i    : write LCS-EPM into file 'LCSEPM_LocARNA_input.fa' as constraint input for LocARNA"   << endl << endl;
  cerr << "  -e    : write LCS-EPM into file 'LCSEPM.epm' as single EPMs"   					<< endl << endl;
  cerr << "  -a    : write all EPMs into file 'allEPM.epm' (depends on -s/ -n)"				<< endl << endl;
  cerr << "  -p    : DO NOT save LCS-EPM as colored postscript file"				<< endl << endl;
  cerr << "  -h    : this info"                                                  	<< endl << endl;
  cerr << "Further Options:"                                                         	<< endl << endl;
  //cerr << "  -L    : Calculate the maximal common arrangment of EPMs (LCS-EPM)"  	<< endl << endl;
  cerr << "  -A    : determine similarity/differences of LCS-EPM with reference alignment"  << endl;
  cerr << "          via two colored postscript files" 						<< endl;
  cerr << "  -v    : verbose output"  						<< endl << endl;

}

///////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////
int main(int argc,char* argv[])
{
  ///////////////////////////////////////////////////////////////////////////////
  // read command line options
  ///////////////////////////////////////////////////////////////////////////////

  int    		opt;
  static ExpaRNAOptions	myOptions;

  while ((opt=getopt(argc,argv,"s:n:Apoieat:d:hvg"))!=-1)
   switch(opt)
   {
      case 's':
            myOptions.minPatternSize = atoi(optarg);
            if (myOptions.minPatternSize<2)
            {
               cerr << endl << " Minimal pattern size not valid!" << endl << endl;
               exit(1);
            }
         break;
      case 'n':
		     myOptions.maxPatternNumber = atoi(optarg);
			    if (myOptions.maxPatternNumber<0)
             {
               cerr << endl << " The number of returned patterns is not valid!" << endl << endl;
	    			exit(1);
             }
         break;
//      case 'L':
//	      myOptions.calcLCSEPM = true;
 //    	break;
      case 'A':
	      myOptions.alignColor = true;
              break;
      case 'p':
	      myOptions.psOut = false;
         break;
      case 'o':
	      //myOptions.align_file 	= optarg;
	      myOptions.alignOut	= true;
         break;
      case 'i':
	      //myOptions.locarna_file 	= optarg;
	      myOptions.locarnaOut 	= true;
         break;
      case 'e':
	      //myOptions.epm_file = optarg;
              myOptions.epmOut	 = true;
         break;
      case 'a':
	      //myOptions.all_file = optarg;
              myOptions.allOut	 = true;
              break;
      case 'd':
	      myOptions.out_dir = optarg;
	      break;
      case 't':
			  if (atoi(optarg)==1){
				  myOptions.EPMscoring = ExpaRNAOptions::BY_SIZE;}
			  else if (atoi(optarg)==2){
  				  myOptions.EPMscoring = ExpaRNAOptions::BY_QUADSIZE;}
			  else { printHelp(myOptions.VERSION_ExpaRNA);
					 exit(1); }
			  break;
      case 'g':
    	  myOptions.ignoreGaps = false;
    	  break;
      case 'h':
    	  printHelp(myOptions.VERSION_ExpaRNA);
    	  exit(0);
    	 break;
      case 'v':
    	  myOptions.verboseOut = true;
    	 break;
      case '?' :
    	  cerr << "Unrecognized option." << endl;
    	  printHelp(myOptions.VERSION_ExpaRNA);
    	  exit(1);
    	 break;
  	}

  	///////////////////////////////////////////////////////////////////////////////
  	// Input recognition
  	///////////////////////////////////////////////////////////////////////////////

	for(;optind<argc;optind++)
		myOptions.cmdLineValues.push_back(argv[optind]);

  	int size=myOptions.cmdLineValues.size();

  	if ( (size<1) || (size>2) )
   	{
    	cerr << "Please use correct input format!" << endl;
    	printHelp(myOptions.VERSION_ExpaRNA);
		exit(1);
   }
  	///////////////////////////////////////////////////////////////////////////////
   //  load data into the main interface and generate the two molecules
   ///////////////////////////////////////////////////////////////////////////////

   ExpaRNAEnsemble::ensOptions = myOptions;
   ExpaRNAEnsemble myNewEnsemble;

   ///////////////////////////////////////////////////////////////////////////////
   //  start algorithms with the provided data & options
   ///////////////////////////////////////////////////////////////////////////////
   myNewEnsemble.startALL();

   return(0);
};
