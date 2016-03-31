/*
 * ExpaRNAOutput.cpp
 *
 *  Created on: 10.03.2009
 *      Author: heyne
 */

#include "ExpaRNAOutput.h"
#include <iostream>
#include <sstream>
#include <ostream>
#include <fstream>
#include <sys/stat.h>

extern "C"
{
	#include <ViennaRNA/fold_vars.h>
	int    PS_rna_plot(char *string, char *structure, char *file);
	int    PS_rna_plot_a(char *string, char *structure, char *file, char *pre, char *post);
}
ExpaRNAOutput::ExpaRNAOutput () : ExpaRNAEnsemble()
{
	int temp = umask(0);
	//cout << "dir: " << ensOptions.out_dir << endl;

	mkdir(ensOptions.out_dir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
}

ExpaRNAOutput::ExpaRNAOutput(const ExpaRNAEnsemble& myEns) : ExpaRNAEnsemble(myEns)
{
	int temp = umask(0);
	//cout << "dir: " << ensOptions.out_dir << endl;

	mkdir(ensOptions.out_dir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

	myMol1 = ensembleMolecules.find(ensOptions.mol1Id)->second;
	myMol2 = ensembleMolecules.find(ensOptions.mol2Id)->second;

	this->output_Fasta();
}

ExpaRNAOutput::~ExpaRNAOutput() {
	ensembleMolecules.clear();
}


void ExpaRNAOutput::output_Fasta() {

	string outname = ensOptions.out_dir + "/" + "ExpaRNA_input.fa";
	ofstream outfile (outname.c_str());

	string tmpId1 = myMol1.getActualStructureId();
	string tmpId2 = myMol2.getActualStructureId();

	if (myMol1.hasGapMapping()){
		myMol1.setWithGapData();
	}
	if (myMol2.hasGapMapping()){
		myMol2.setWithGapData();
	}

	outfile << ">" << myMol1.getId() << endl;
	outfile << myMol1.getSequence() << endl;
	outfile << myMol1.getStructure() << endl;
	outfile << ">" << myMol2.getId() << endl;
	outfile << myMol2.getSequence() << endl;
	outfile << myMol2.getStructure();
	outfile.close();
}

void ExpaRNAOutput::LCSEPMtoPS()
{
	MapToPS("LCSEPM",sizeLCSEPM,ensembleLCSEPM);
}

void ExpaRNAOutput::MapToPS(	const string&	  myId,
				const int&      mySize,
				const PatternPairMap& myMap)
{
   string func_str="\
   /drawpattern {\n\
      /Panz pattern length def\n\
      0 1 pattern length 1 sub {\n\
         /i exch def\n\
         pattern i get\n\
         newpath\n\
         {\n\
            1 Panz div i mul 0 add 1 1 sethsbcolor\n\
            coor exch 1 sub get aload pop fsize 2.1 div 0 360 arc\n\
            fill\n\
         } forall\n\
      } for\n\
   } bind def\n\
   \n\
   /pattern [\n";
   string clus1_str,clus2_str;
   ID myID;

  /* DIR *dirPTR;
   dirPTR=opendir("./postscript");
   if (dirPTR =NULL)
      cout << "No such dir!"<<endl;*/

   string label1Str,label2Str;

   for (unsigned int i=1;i<=myMol1.getSequence().length();++i)
   {
      if (i % 50 == 0)
         label1Str += biu::Util_String::int2str(i) + " 0.5 0.5 (" + biu::Util_String::int2str(i) + ") Label\n";
   }

   for (unsigned int i=1;i<=myMol2.getSequence().length();++i)
   {
      if (i % 50 == 0)
         label2Str += biu::Util_String::int2str(i) + " 0.5 0.5 (" + biu::Util_String::int2str(i) + ") Label\n";
   }

   for (PatternPairMap::patListCITER i=myMap.getList().begin();i != myMap.getList().end();i++)
   {
      intVec tmpvec1=(*i)->getFirstPat().getPat();
      if (myMol1.hasGapMapping()){
    	  tmpvec1 = myMol1.ReMapPattern(tmpvec1);
	  }
      clus1_str+="["+biu::Util_String::intvec2str(tmpvec1," ")+"]\n";

      intVec tmpvec2=(*i)->getSecPat().getPat();
      if (myMol2.hasGapMapping()){
    	  tmpvec2 = myMol2.ReMapPattern(tmpvec2);
      }
      clus2_str+="["+biu::Util_String::intvec2str(tmpvec2," ")+"]\n";
   }
   clus1_str+="] def\n\n";
   clus2_str+="] def\n\n";
   clus1_str=func_str+clus1_str;
   clus2_str=func_str+clus2_str;

   string tmpId1 = myMol1.getActualStructureId();
   myMol1.setWithGapData();
   string psfilename = ensOptions.out_dir + "/" + myId + "_RNA1_"+biu::Util_String::int2str(mySize) + "_" + myMol1.getActualStructureId()+".ps";
   PS_rna_plot_a(const_cast<char*>(myMol1.getSequence().c_str()),
                 const_cast<char*>(myMol1.getStructure().c_str()),
                 const_cast<char*>(psfilename.c_str()),
                 const_cast<char*>(clus1_str.c_str()),
                 const_cast<char*>(string("drawpattern\ndrawbases\n"+label1Str).c_str()));

   string tmpId2 = myMol2.getActualStructureId();
   myMol2.setWithGapData();
   psfilename = ensOptions.out_dir + "/" + myId + "_RNA2_"+biu::Util_String::int2str(mySize) + "_" + myMol2.getActualStructureId()+".ps";
   PS_rna_plot_a(const_cast<char*>(myMol2.getSequence().c_str()),
                 const_cast<char*>(myMol2.getStructure().c_str()),
                 const_cast<char*>(psfilename.c_str()),
                 const_cast<char*>(clus2_str.c_str()),
                 const_cast<char*>(string("drawpattern\ndrawbases\n"+label2Str).c_str()));

   myMol1.setActualStructureId(tmpId1);
   myMol2.setActualStructureId(tmpId2);
}

void	ExpaRNAOutput::output_LocARNA()
{
	// extract matching edges (pairs of positions) from LCS-EPM
	vector<intPair> matchingsLCSEPM;
	intVec positionsSeq1LCSEPM;
	intVec positionsSeq2LCSEPM;

	for (PatternPairMap::patListCITER i=ensembleLCSEPM.getList().begin();i != ensembleLCSEPM.getList().end();i++)
	{
		positionsSeq1LCSEPM.insert(positionsSeq1LCSEPM.end(),(*i)->getFirstPat().getPat().begin(),(*i)->getFirstPat().getPat().end());
		positionsSeq2LCSEPM.insert(positionsSeq2LCSEPM.end(),(*i)->getSecPat().getPat().begin(),(*i)->getSecPat().getPat().end());
		//SinglePattern my1 = (*i)->getFirstPat();
		//my1.print();
		//my1 = (*i)->getSecPat();
		//my1.print();
	}

	sort(positionsSeq1LCSEPM.begin(),positionsSeq1LCSEPM.end());
	sort(positionsSeq2LCSEPM.begin(),positionsSeq2LCSEPM.end());;

	for (unsigned int i=0;i<positionsSeq1LCSEPM.size();++i)
	{
		matchingsLCSEPM.push_back(make_pair(positionsSeq1LCSEPM[i],positionsSeq2LCSEPM[i]));
	}
	string outname = ensOptions.out_dir + "/" + ensOptions.locarna_file;
	ofstream outLocARNAfile (outname.c_str());

	int last_edge_seq1,last_edge_seq2;
	last_edge_seq1=0;
	last_edge_seq2=0;

	string seq1_0,seq2_0,seq1_1,seq1_2,seq1_3,seq2_1,seq2_2,seq2_3;
	//int edge = 100;
	int edge = 1;
	for (vector<intPair>::iterator i_edge = matchingsLCSEPM.begin(); i_edge != matchingsLCSEPM.end();++i_edge)
	{
		//cout << "first: " << (*i_edge).first << " second: " << (*i_edge).second << endl;

		for (int i=last_edge_seq1+1;i<(*i_edge).first;++i)
		{
			seq1_0.push_back('.');
			seq1_1.push_back('.');
			seq1_2.push_back('.');
			seq1_3.push_back('.');
		}

		for (int j=last_edge_seq2+1;j<(*i_edge).second;++j)
		{
			seq2_0.push_back('.');
			seq2_1.push_back('.');
			seq2_2.push_back('.');
			seq2_3.push_back('.');
		}

		ostringstream edge_st_;
		edge_st_ << edge;
		string edge_st;
		edge_st = edge_st_.str();
		const char *c_str_edge = edge_st.c_str();

		seq1_0.push_back('0'+(edge%10000)/1000);
		seq1_1.push_back('0'+(edge%1000)/100);
		seq1_2.push_back('0'+(edge%100)/10);
		seq1_3.push_back('0'+(edge%10));

//		seq1_1.push_back(c_str_edge[0]);
//		seq1_2.push_back(c_str_edge[1]);
//		seq1_3.push_back(c_str_edge[2]);


//		seq2_1.push_back(c_str_edge[0]);
//		seq2_2.push_back(c_str_edge[1]);
//		seq2_3.push_back(c_str_edge[2]);

		seq2_0.push_back('0'+(edge%10000)/1000);
		seq2_1.push_back('0'+(edge%1000)/100);
		seq2_2.push_back('0'+(edge%100)/10);
		seq2_3.push_back('0'+(edge%10));


		++edge;

		last_edge_seq1= (*i_edge).first;
		last_edge_seq2 = (*i_edge).second;
	}

	// end stuff
	for (int i=last_edge_seq1+1;i<=myMol1.getLength();++i)
	{
		seq1_0.push_back('.');
		seq1_1.push_back('.');
		seq1_2.push_back('.');
		seq1_3.push_back('.');
	}

	for (int j=last_edge_seq2+1;j<=myMol2.getLength();++j)
	{
		seq2_0.push_back('.');
		seq2_1.push_back('.');
		seq2_2.push_back('.');
		seq2_3.push_back('.');
	}

	seq1_0 += " #1";
	seq1_1 += " #2";
	seq1_2 += " #3";
	seq1_3 += " #4";

	seq2_0 += " #1";
	seq2_1 += " #2";
	seq2_2 += " #3";
	seq2_3 += " #4";

	outLocARNAfile << ">" << myMol1.getId() << endl << myMol1.getSequence() << endl;
	outLocARNAfile << seq1_0 << "\n" << seq1_1 << endl << seq1_2 << endl << seq1_3 << endl;
	outLocARNAfile << ">" << myMol2.getId() << endl << myMol2.getSequence() << endl;
	outLocARNAfile << seq2_0 << "\n" << seq2_1 << endl << seq2_2 << endl << seq2_3 << endl << endl;

	outLocARNAfile.close();
}

void	ExpaRNAOutput::output_LCSEPM()
{
	// extract matching edges (pairs of positions) from LCS-EPM
	vector<intPair> matchingsLCSEPM;
	intVec positionsSeq1LCSEPM;
	intVec positionsSeq2LCSEPM;
	string outname = ensOptions.out_dir + "/" + ensOptions.epm_file;
	ofstream outfile (outname.c_str());

	string tmpId1 = myMol1.getActualStructureId();
	string tmpId2 = myMol2.getActualStructureId();

	if (myMol1.hasGapMapping()){
		myMol1.setWithGapData();
	}
	if (myMol2.hasGapMapping()){
		myMol2.setWithGapData();
	}

	outfile << ">" << myMol1.getId() << endl;
	outfile << myMol1.getSequence() << endl;
	outfile << myMol1.getStructure() << endl;
	outfile << ">" << myMol2.getId() << endl;
	outfile << myMol2.getSequence() << endl;
	outfile << myMol2.getStructure() << endl;
	outfile << endl;
	outfile << "score: " << this->scoreLCSEPM << endl;
	outfile << "bases: " << this->sizeLCSEPM << endl;
	outfile << "EPMs : " << this->ensembleLCSEPM.size() << endl;
	outfile << endl;
	outfile << "<ID> <SEQ> <STRUCT> <EPM_POS_RNA1> <EPM_POS_RNA2> <SCORE_EPM>" << endl;
	outfile << endl;

	myMol1.setActualStructureId(tmpId1);
	myMol2.setActualStructureId(tmpId2);

	for (PatternPairMap::patListCITER i=ensembleLCSEPM.getList().begin();i != ensembleLCSEPM.getList().end();i++)
	{
		positionsSeq1LCSEPM.assign((*i)->getFirstPat().getPat().begin(),(*i)->getFirstPat().getPat().end());
		positionsSeq2LCSEPM.assign((*i)->getSecPat().getPat().begin(),(*i)->getSecPat().getPat().end());

		string structure;
		string sequence;

		for (unsigned int j=0;j<positionsSeq1LCSEPM.size();++j)
		{
			structure += myMol1.getStructure(positionsSeq1LCSEPM[j]);
			sequence  += myMol1.getBase(positionsSeq1LCSEPM[j]);
		}

		if (myMol1.hasGapMapping()){
			positionsSeq1LCSEPM = myMol1.ReMapPattern(positionsSeq1LCSEPM);
		}

		if (myMol2.hasGapMapping()){
			positionsSeq2LCSEPM = myMol2.ReMapPattern(positionsSeq2LCSEPM);
		}

		outfile << (*i)->getId() << " "<< sequence << " " << structure;
		outfile << " " << biu::Util_String::intvec2str(positionsSeq1LCSEPM,":");
		outfile << " " << biu::Util_String::intvec2str(positionsSeq2LCSEPM,":");
		outfile << " " << (*i)->getEPMScore();
		outfile << endl;
	}

	outfile.close();
}


void	ExpaRNAOutput::output_allEPMs()
{
	intVec positionsSeq1LCSEPM;
	intVec positionsSeq2LCSEPM;

	string outname = ensOptions.out_dir + "/" + ensOptions.all_file;
	ofstream allfile (outname.c_str());

	string tmpId1 = myMol1.getActualStructureId();
	string tmpId2 = myMol2.getActualStructureId();

	if (myMol1.hasGapMapping()){
		myMol1.setWithGapData();
	}
	if (myMol2.hasGapMapping()){
		myMol2.setWithGapData();
	}

	allfile << ">" << myMol1.getId() << endl;
	allfile << myMol1.getSequence() << endl;
	allfile << myMol1.getStructure() << endl;
	allfile << ">" << myMol2.getId() << endl;
	allfile << myMol2.getSequence() << endl;
	allfile << myMol2.getStructure() << endl;
	allfile << endl;
	allfile << "<ID> <SEQ> <STRUCT> <EPM_POS_RNA1> <EPM_POS_RNA2>" << endl;
	allfile << endl;

	myMol1.setActualStructureId(tmpId1);
	myMol2.setActualStructureId(tmpId2);

	ensemblePatterns.makeOrderedMap();

	for (PatternPairMap::orderedMapITER i = ensemblePatterns.getOrderedMap2().begin();i!=ensemblePatterns.getOrderedMap().end();++i)
	{
		positionsSeq1LCSEPM.assign((*i).second->getFirstPat().getPat().begin(),(*i).second->getFirstPat().getPat().end());
		positionsSeq2LCSEPM.assign((*i).second->getSecPat().getPat().begin(),(*i).second->getSecPat().getPat().end());

		string structure;
		string sequence;

		for (unsigned int j=0;j<positionsSeq1LCSEPM.size();++j)
		{
			structure += myMol1.getStructure(positionsSeq1LCSEPM[j]);
			sequence  += myMol1.getBase(positionsSeq1LCSEPM[j]);
		}

		if (myMol1.hasGapMapping()){
			positionsSeq1LCSEPM = myMol1.ReMapPattern(positionsSeq1LCSEPM);
		}
		if (myMol2.hasGapMapping()){
			positionsSeq2LCSEPM = myMol2.ReMapPattern(positionsSeq2LCSEPM);
		}

		allfile << (*i).second->getId() << " " << sequence << " " << structure;
		allfile << " " << biu::Util_String::intvec2str(positionsSeq1LCSEPM,":");
		allfile << " " << biu::Util_String::intvec2str(positionsSeq2LCSEPM,":");
		allfile << endl;
	}
	allfile.close();
}

void	ExpaRNAOutput::output_Clustal()
{
	// extract matching edges (pairs of positions) from LCS-EPM
	vector<intPair> matchingsLCSEPM;
	intVec positionsSeq1LCSEPM;
	intVec positionsSeq2LCSEPM;

	for (PatternPairMap::patListCITER i=ensembleLCSEPM.getList().begin();i != ensembleLCSEPM.getList().end();i++)
	{
		positionsSeq1LCSEPM.insert(positionsSeq1LCSEPM.end(),(*i)->getFirstPat().getPat().begin(),(*i)->getFirstPat().getPat().end());
		positionsSeq2LCSEPM.insert(positionsSeq2LCSEPM.end(),(*i)->getSecPat().getPat().begin(),(*i)->getSecPat().getPat().end());
	}

	sort(positionsSeq1LCSEPM.begin(),positionsSeq1LCSEPM.end());
	sort(positionsSeq2LCSEPM.begin(),positionsSeq2LCSEPM.end());;

	for (unsigned int i=0;i<positionsSeq1LCSEPM.size();++i)
	{
		matchingsLCSEPM.push_back(make_pair(positionsSeq1LCSEPM[i],positionsSeq2LCSEPM[i]));
	}

	string outname = ensOptions.out_dir + "/" + ensOptions.align_file;
	ofstream outfile (outname.c_str());

	string seq1_aln,seq2_aln,seq1_aln_str,seq2_aln_str;

	int last_edge_seq1,last_edge_seq2;
	last_edge_seq1=0;
	last_edge_seq2=0;

	for (vector<intPair>::iterator i_edge = matchingsLCSEPM.begin(); i_edge != matchingsLCSEPM.end();++i_edge)
	{
		for (int i=last_edge_seq1+1;i<(*i_edge).first;++i)
		{
			seq1_aln.push_back(myMol1.getBase(i));
			seq2_aln.push_back('-');
			seq1_aln_str.push_back(myMol1.getStructure(i));
			seq2_aln_str.push_back('-');
		}
		for (int j=last_edge_seq2+1;j<(*i_edge).second;++j)
		{
			seq1_aln.push_back('-');
			seq2_aln.push_back(myMol2.getBase(j));
			seq1_aln_str.push_back('-');
			seq2_aln_str.push_back(myMol2.getStructure(j));
		}

		seq1_aln.push_back(myMol1.getBase((*i_edge).first));
		seq2_aln.push_back(myMol2.getBase((*i_edge).second));

		seq1_aln_str.push_back(myMol1.getStructure((*i_edge).first));
		seq2_aln_str.push_back(myMol2.getStructure((*i_edge).second));

		last_edge_seq1= (*i_edge).first;
		last_edge_seq2 = (*i_edge).second;
	}

	// for the part after the last edge
	for (int i=last_edge_seq1+1;i<=myMol1.getLength();++i)
	{
				seq1_aln.push_back(myMol1.getBase(i));
				seq2_aln.push_back('-');
				seq1_aln_str.push_back(myMol1.getStructure(i));
				seq2_aln_str.push_back('-');
			}
	for (int j=last_edge_seq2+1;j<=myMol2.getLength();++j)
	{
		seq1_aln.push_back('-');
		seq2_aln.push_back(myMol2.getBase(j));
		seq1_aln_str.push_back('-');
		seq2_aln_str.push_back(myMol2.getStructure(j));
	}

	outfile << "CLUSTAL W (1.83) multiple sequence alignment --- expaRNA 0.7.2 - exact pattern Alignment of RNA --- Score: "<< matchingsLCSEPM.size() << endl <<endl;

	string tmp1 = myMol1.getId() +"      ";
	string tmp2 = myMol2.getId() +"      ";
	if (tmp1.length() < tmp2.length())
		tmp1.resize(tmp2.length(),' ');
	else
		if (tmp2.length() < tmp1.length())
			tmp2.resize(tmp1.length(),' ');
	string tmp3;
	tmp3.resize(tmp1.length(),' ');
	outfile << tmp3 << seq1_aln_str <<endl;
	outfile << tmp1 << seq1_aln << endl;
	outfile << tmp2 << seq2_aln << endl;
	outfile << tmp3 << seq2_aln_str << endl << endl;
	outfile.close();
}

void	ExpaRNAOutput::compareAlign()
{
	std::ifstream fileAlign(ensOptions.cmdLineValues[1].c_str());

	string strRead;
	string alignSeq1,alignSeq2;

	if(fileAlign)
	{
		while(fileAlign)
		{
			strRead.clear();
			while(fileAlign)
			{
				getline(fileAlign,strRead);
				string::size_type loc = strRead.find_first_of("AUTCGautcg-.");
				if( loc != string::npos )
				{
					alignSeq1 += strRead;
					break;
				}
			}
			strRead.clear();
			while(fileAlign)
			{
				getline(fileAlign,strRead);
				string::size_type loc = strRead.find_first_of("AUTCGautcg-.");
				if( loc != string::npos )
				{
					alignSeq2 += strRead;
					break;
				}
			}
		}

		if (alignSeq1.size() != alignSeq2.size() )
		{
			cerr << "Alginment format is wrong. Maybe different sizes!"<<endl;
			exit(1);
		}

		intVec matchingsSeq1Align,matchingsSeq2Align;
		vector<intPair> matchingsAlign;
		int posSeq1,posSeq2;
		posSeq1=1;
		posSeq2=1;

		for (unsigned int i=0;i<alignSeq1.size();++i)
		{
			if ( (alignSeq1[i] == alignSeq2[i]) && (alignSeq2[i]!='-') &&
			     (myMol1.getStructure(posSeq1)== myMol2.getStructure(posSeq2)) &&
			     (myMol1.getBase(posSeq1) == (myMol2.getBase(posSeq2))) )
			{
				matchingsSeq1Align.push_back(posSeq1);
				matchingsSeq2Align.push_back(posSeq2);
				matchingsAlign.push_back(make_pair(posSeq1,posSeq2));

			}
			if (alignSeq1[i] != '-')
				++posSeq1;
			if (alignSeq2[i] != '-')
				++posSeq2;
		}

		cout << "align1:"<< alignSeq1.size()<<endl<< alignSeq1 << endl << endl;
		cout << "align2:"<< alignSeq2.size()<<endl<< alignSeq2 << endl << endl;

		for (unsigned int i=0;i<matchingsSeq1Align.size();++i)
		{
			cout << "match: "<<matchingsSeq1Align[i]<<" with " << matchingsSeq2Align[i]<<endl;
		}

		PatternPairMap myAlignMap;
/*		int act=1;
		int stepSize = 300;

		while (act<matchingsSeq1Align.size())
		{
			string actStr = biu::Util_String::int2str(act);
			int i=1;
			intVec partAlign1,partAlign2;

			while(act<matchingsSeq1Align.size() && (i<=stepSize))
			{
				partAlign1.push_back(matchingsSeq1Align[act]);
				partAlign2.push_back(matchingsSeq2Align[act]);
				++i;++act;
			}
			SinglePattern myAlign1(myMol1.getId(),myMol1.getActualStructureId(),"align"+actStr,partAlign1);
			SinglePattern myAlign2(myMol2.getId(),myMol2.getActualStructureId(),"align"+actStr,partAlign2);

			myAlignMap.add("align"+actStr,matchingsSeq1Align.size(),myAlign1,myAlign2);
		}*/

		intVec matchesSeq1Both;
		intVec matchesSeq2Both;
		intVec matchesSeq1ERPAPOnly;
		intVec matchesSeq2ERPAPOnly;
		intVec matchesSeq1AlignOnly;
		intVec matchesSeq2AlignOnly;
		intVec matchesSeq1ERPAP;
		intVec matchesSeq2ERPAP;
		intVec matchesSeq1ERPAPRest;
		intVec matchesSeq2ERPAPRest;
		intVec matchesSeq1AlignRest;
		intVec matchesSeq2AlignRest;
		intVec matchesSeq1ERPAPAlign;
		intVec matchesSeq2ERPAPAlign;


		// extract matching edges (pairs of positions) from LCS-EPM
		vector<intPair> matchesERPAP;

		for (PatternPairMap::patListCITER i=ensembleLCSEPM.getList().begin();i != ensembleLCSEPM.getList().end();i++)
		{
		    matchesSeq1ERPAP.insert(matchesSeq1ERPAP.end(),(*i)->getFirstPat().getPat().begin(),(*i)->getFirstPat().getPat().end());
		    matchesSeq2ERPAP.insert(matchesSeq2ERPAP.end(),(*i)->getSecPat().getPat().begin(),(*i)->getSecPat().getPat().end());
		}
		sort(matchesSeq1ERPAP.begin(),matchesSeq1ERPAP.end());
		sort(matchesSeq2ERPAP.begin(),matchesSeq2ERPAP.end());;

		for (unsigned int i=0;i<matchesSeq1ERPAP.size();++i)
		{
			matchesERPAP.push_back(make_pair(matchesSeq1ERPAP[i],matchesSeq2ERPAP[i]));
		}


		// extract intersecting edges from LCS-EPM with given alignment (from input file)
		vector<int> v(matchesSeq1ERPAP.size());
		vector<int> w(matchesSeq2ERPAP.size());
		vector<intPair> u(matchesSeq1ERPAP.size());
		vector<int>::iterator it1;
		vector<int>::iterator it2;
		vector<intPair>::iterator it3;

		it1 = set_intersection(matchesSeq1ERPAP.begin(),matchesSeq1ERPAP.end(),matchingsSeq1Align.begin(),matchingsSeq1Align.end(),v.begin());
		it2 = set_intersection(matchesSeq2ERPAP.begin(),matchesSeq2ERPAP.end(),matchingsSeq2Align.begin(),matchingsSeq2Align.end(),w.begin());
		it3 = set_intersection(matchesERPAP.begin(),matchesERPAP.end(),matchingsAlign.begin(),matchingsAlign.end(),u.begin());

		//matchesSeq1Both.insert(matchesSeq1Both.end(),v.begin(),it1);
		//matchesSeq2Both.insert(matchesSeq2Both.end(),w.begin(),it2);

		for (vector<intPair>::iterator i=u.begin();i != it3;++i)
		{
			matchesSeq1Both.push_back(i->first);
			matchesSeq2Both.push_back(i->second);
		}
		// rest von ERPAP
		it1 = set_difference(matchesSeq1ERPAP.begin(),matchesSeq1ERPAP.end(),matchesSeq1Both.begin(),matchesSeq1Both.end(),v.begin());
		it2 = set_difference(matchesSeq2ERPAP.begin(),matchesSeq2ERPAP.end(),matchesSeq2Both.begin(),matchesSeq2Both.end(),w.begin());

		matchesSeq1ERPAPRest.insert(matchesSeq1ERPAPRest.end(),v.begin(),it1);
		matchesSeq2ERPAPRest.insert(matchesSeq2ERPAPRest.end(),w.begin(),it2);

		//rest des alignments
		it1 = set_difference(matchingsSeq1Align.begin(),matchingsSeq1Align.end(),matchesSeq1Both.begin(),matchesSeq1Both.end(),v.begin());
		it2 = set_difference(matchingsSeq2Align.begin(),matchingsSeq2Align.end(),matchesSeq2Both.begin(),matchesSeq2Both.end(),w.begin());

		matchesSeq1AlignRest.insert(matchesSeq1AlignRest.end(),v.begin(),it1);
		matchesSeq2AlignRest.insert(matchesSeq2AlignRest.end(),w.begin(),it2);

		// schnittmenge align Rest / ERPAP Rest
		it1 = set_intersection(matchesSeq1ERPAPRest.begin(),matchesSeq1ERPAPRest.end(),matchesSeq1AlignRest.begin(),matchesSeq1AlignRest.end(),v.begin());
		it2 = set_intersection(matchesSeq2ERPAPRest.begin(),matchesSeq2ERPAPRest.end(),matchesSeq2AlignRest.begin(),matchesSeq2AlignRest.end(),w.begin());

		matchesSeq1ERPAPAlign.insert(matchesSeq1ERPAPAlign.end(),v.begin(),it1);
		matchesSeq2ERPAPAlign.insert(matchesSeq2ERPAPAlign.end(),w.begin(),it2);

		// unique rest ERPAP
		it1 = set_difference(matchesSeq1ERPAPRest.begin(),matchesSeq1ERPAPRest.end(),matchesSeq1ERPAPAlign.begin(),matchesSeq1ERPAPAlign.end(),v.begin());
		it2 = set_difference(matchesSeq2ERPAPRest.begin(),matchesSeq2ERPAPRest.end(),matchesSeq2ERPAPAlign.begin(),matchesSeq2ERPAPAlign.end(),w.begin());

		matchesSeq1ERPAPOnly.insert(matchesSeq1ERPAPOnly.end(),v.begin(),it1);
		matchesSeq2ERPAPOnly.insert(matchesSeq2ERPAPOnly.end(),w.begin(),it2);

		// unique rest alignment
		it1 = set_difference(matchesSeq1AlignRest.begin(),matchesSeq1AlignRest.end(),matchesSeq1ERPAPAlign.begin(),matchesSeq1ERPAPAlign.end(),v.begin());
		it2 = set_difference(matchesSeq2AlignRest.begin(),matchesSeq2AlignRest.end(),matchesSeq2ERPAPAlign.begin(),matchesSeq2ERPAPAlign.end(),w.begin());

		matchesSeq1AlignOnly.insert(matchesSeq1AlignOnly.end(),v.begin(),it1);
		matchesSeq2AlignOnly.insert(matchesSeq2AlignOnly.end(),w.begin(),it2);

		SinglePattern myAlign1,myAlign2;

		myAlign1 = SinglePattern(myMol1.getId(),myMol1.getActualStructureId(),"align+ERPAP",matchesSeq1AlignOnly);
		myAlign2 = SinglePattern(myMol2.getId(),myMol2.getActualStructureId(),"align+ERPAP",matchesSeq2AlignOnly);
		myAlignMap.add("align+ERPAP",matchingsSeq1Align.size(),myAlign1,myAlign2);

		myAlign1 = SinglePattern(myMol1.getId(),myMol1.getActualStructureId(),"align+ERPAP",matchesSeq1ERPAPOnly);
		myAlign2 = SinglePattern(myMol2.getId(),myMol2.getActualStructureId(),"align+ERPAP",matchesSeq2ERPAPOnly);
		myAlignMap.add("align+ERPAP",matchingsSeq1Align.size(),myAlign1,myAlign2);

		myAlign1 = SinglePattern(myMol1.getId(),myMol1.getActualStructureId(),"align+ERPAP",matchesSeq1ERPAPAlign);
		myAlign2 = SinglePattern(myMol2.getId(),myMol2.getActualStructureId(),"align+ERPAP",matchesSeq2ERPAPAlign);
		myAlignMap.add("align+ERPAP",matchingsSeq1Align.size(),myAlign1,myAlign2);

		myAlign1 = SinglePattern(myMol1.getId(),myMol1.getActualStructureId(),"align+ERPAP",matchesSeq1Both);
		myAlign2 = SinglePattern(myMol2.getId(),myMol2.getActualStructureId(),"align+ERPAP",matchesSeq2Both);
		myAlignMap.add("align+ERPAP",matchingsSeq1Align.size(),myAlign1,myAlign2);

		for (unsigned int i=0;i<matchesSeq1ERPAPOnly.size();++i)
			cout << "sec: " << matchesSeq1ERPAPOnly[i]<<endl;

		MapToPS("align",matchingsSeq1Align.size(),myAlignMap);

		cout << "# exact matches within reference alignment  : " << matchingsSeq1Align.size() << endl;
		cout << "# exact matches LCS-ERP and alignment (seq1): " << matchesSeq1Both.size() << endl;
		cout << "# exact matches LCS-ERP and alignment (seq2): " << matchesSeq2Both.size() << endl;
		cout << "# matches in seq1 with different position in seq2: " << matchesSeq1ERPAPAlign.size() << endl;
		cout << "# matches in seq2 with different position in seq1: " << matchesSeq2ERPAPAlign.size() << endl;
		cout << "# matches unique to LCS-ERP in seq1  : " << matchesSeq1ERPAPOnly.size() << endl;
		cout << "# matches unique to LCS-ERP in seq2  : " << matchesSeq2ERPAPOnly.size() << endl;
		cout << "# matches unique to alignment in seq1: " << matchesSeq1AlignOnly.size() << endl;
		cout << "# matches unique to alignment in seq2: " << matchesSeq2AlignOnly.size() << endl;
	}
}
