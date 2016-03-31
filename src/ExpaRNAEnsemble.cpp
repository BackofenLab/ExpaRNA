#include "ExpaRNAEnsemble.h"
#include <dirent.h>
#include"ExpaRNAOutput.h"

// static members

ExpaRNAOptions			ExpaRNAEnsemble::ensOptions;

// other members

void ExpaRNAEnsemble::clearEnsemble()
{
   ensembleMolecules.clear();
   ensemblePatterns = PatternPairMap();
}

void ExpaRNAEnsemble::load_input_data(const vector<string>& my_input)
{
	cout << "###############################################################################"<<endl;
	cout << "    ExpaRNA - Exact Pattern Alignment of RNA  - version " << ensOptions.VERSION_ExpaRNA << endl;
	cout << "           (c) copyright by Steffen Heyne - 2008-2009"<<endl;
	cout << "###############################################################################"<<endl<<endl;
	cout << " --> Creating virtual molecules with provided data..." << endl;

	clearEnsemble();
	Molecule myMol1,myMol2;
	ExpaRNAinput inputReader;

	inputReader.readFASTA(my_input[0],myMol1,myMol2);

	if (ensOptions.ignoreGaps) {
		myMol1.GapMapping();
		myMol2.GapMapping();
	}

	cout << "finished!" << endl << endl;
	if (ensOptions.verboseOut) myMol1.printInfo(false);
	if (ensOptions.verboseOut) myMol2.printInfo(false);

	ensOptions.mol1Id = myMol1.getId();
	ensOptions.mol2Id = myMol2.getId();

	ensembleMolecules.insert(make_pair(myMol1.getId(),myMol1));
	ensembleMolecules.insert(make_pair(myMol2.getId(),myMol2));
}

void ExpaRNAEnsemble::startALL()
{
   struct timeval *Tps, *Tpf;
   struct timezone *Tzp;
   Tps = (struct timeval*) malloc(sizeof(struct timeval));
   Tpf = (struct timeval*) malloc(sizeof(struct timeval));
   Tzp = 0;

   // pattern finding
   gettimeofday (Tps, Tzp);
   EPMFinding PairMCS(ensembleMolecules.find(ensOptions.mol1Id)->second,ensembleMolecules.find(ensOptions.mol2Id)->second,ensemblePatterns,ensOptions);
   PairMCS.start();
   gettimeofday (Tpf, Tzp);
   long double time = ((Tpf->tv_sec-Tps->tv_sec)*1000000 + Tpf->tv_usec-Tps->tv_usec)/1000000.0;
   cout << "time to find EPMs in secs: " << setprecision(4) << time << endl;
   if (ensOptions.verboseOut) cout << "The EPM library contains " << ensemblePatterns.size() << " EPMs." << endl;

   sizeLCSEPM = 0;
   scoreLCSEPM = 0;

   //LCSEPM
   if (ensOptions.calcLCSEPM)
	   startLCSEPM();

   // start all output
   cout << endl << " -->> start writing output" << endl;

   // chdir to output
   //string dirname = "output";
   //chdir(const_cast<char*>(dirname.c_str()));

   ExpaRNAOutput outputWriter(*this);

   // postscript output
   if (ensOptions.psOut)
   {
	   cout << "    -->> Writing LCS-EPM as Postscript file...";
	   outputWriter.LCSEPMtoPS();
	   cout << "  LCS-EPM output done!"<<endl;
   }

   // output LCS-EPM as constraint input for LocARNA
   if (ensOptions.locarnaOut)
	   outputWriter.output_LocARNA();

   // output LCS-EPM as alignment
   if (ensOptions.alignOut)
	   outputWriter.output_Clustal();

   // output LCS-EPM in plain file
   if (ensOptions.epmOut)
	   outputWriter.output_LCSEPM();

   // output all EPM in plain file
   if (ensOptions.allOut)
	   outputWriter.output_allEPMs();

   // extract Alignment matchings for Coloring
   if (ensOptions.alignColor)
	   outputWriter.compareAlign();

   cout << " -->> end writing output" << endl;
}

void ExpaRNAEnsemble::startLCSEPM()
{
	struct timeval *Tps, *Tpf;
	struct timezone *Tzp;
	Tps = (struct timeval*) malloc(sizeof(struct timeval));
	Tpf = (struct timeval*) malloc(sizeof(struct timeval));
	Tzp = 0;

	gettimeofday (Tps, Tzp);
	cout << endl << " -->> starting LCS-EPM algorithm" << endl;
	if (ensOptions.verboseOut) cout << "   EPMs for LCS-EPM: " << ensemblePatterns.size() << endl;
	LCSEPM myLCSEPM(	ensembleMolecules.find(ensOptions.mol1Id)->second,
						ensembleMolecules.find(ensOptions.mol2Id)->second,
						ensemblePatterns,
						ensembleLCSEPM,
						sizeLCSEPM,
						scoreLCSEPM,
						ensOptions);

	myLCSEPM.calculateLCSEPM();

	gettimeofday (Tpf, Tzp);
	long double time = ((Tpf->tv_sec-Tps->tv_sec)*1000000 + Tpf->tv_usec-Tps->tv_usec)/1000000.0;
	cout << "   time for LCS-EPM in secs: " << setprecision(4) << time << endl;
	cout << " -->> end LCS-EPM algorithm" << endl;
}
