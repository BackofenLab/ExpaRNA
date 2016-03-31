/*
 * ExpaRNAinput.cpp
 *
 *  Created on: 13.03.2009
 *      Author: heyne
 */

#include "ExpaRNAinput.h"
#include <fstream>

// for getrusage()
#include <sys/resource.h>
#include <sys/types.h>
// for gettimeofday()
#include <sys/time.h>
// for setprecision
#include <iomanip>

extern "C"
{
	#include <ViennaRNA/fold_vars.h>
	float  fold(const char *sequence, char *structure);
}

ExpaRNAinput::ExpaRNAinput() {
	// TODO Auto-generated constructor stub

}

ExpaRNAinput::~ExpaRNAinput() {
	// TODO Auto-generated destructor stub
}

//return string, start from pos start, up to 'break' character, remove all escape chars
std::string ExpaRNAinput::getFASTALine(std::string myStr,int start,char myBreak)
{
	int tmppos=start;
	string retStr;
	retStr.clear();
	while ((tmppos<myStr.size()) && isprint(myStr[tmppos]))
	{
		if (myStr[tmppos] == myBreak) break;
		retStr.push_back(myStr[tmppos]);
		++tmppos;
	}
	return retStr;
}

// read input file in FASTA style, ATTENTION: seq/str only in one line each!!!
void ExpaRNAinput::readFASTA(string myFASTAfile, Molecule& myMol1, Molecule& myMol2)
{
	std::ifstream fa_file(myFASTAfile.c_str(),ios::in);

	string line,seq1,seq1_name,seq1_str,seq1_strname,seq2,seq2_name,seq2_str,seq2_strname;

	if (fa_file)
	{
		getline(fa_file,line);

		if (line[0]=='>')
		{
			seq1_name = getFASTALine(line,1,' ');
		}
		else
		{
			cerr <<" Wrong fasta format!" << endl; exit(1);
		}

		getline(fa_file,line);

		seq1 = getFASTALine(line,0,' ');

		getline(fa_file,line);
		if (line[0] != '>')
		{
			seq1_str = getFASTALine(line,0,' ');
			getline(fa_file,line);
			if (line[0]=='>')
			{
				seq2_name = getFASTALine(line,1,' ');
			}
			else
			{
				cerr << seq1_str << line << " wrong fasta format!" << endl; exit(1);
			}
		}
		else
		{
			if (line[0] =='>')
			{
				seq2_name = getFASTALine(line,1,' ');
			}
			else
			{
				cerr <<" wrong fasta format!" << endl; exit(1);
			}
		}
		fa_file >> line;
		seq2 = getFASTALine(line,0,' ');
		line="";
		if (!fa_file.eof())
		{
			fa_file >> line;
			if (line[0] != '>')
			{
				seq2_str = getFASTALine(line,0,' ');
			}
			else
			{
				cerr <<" wrong fasta format!" << endl; exit(1);
			}
		}
		fa_file.close();

		// check which structure is not present and fold with vienna RNA package
		dangles = 2; // vienna package variable for dangling ends, 2 standard for partition function folding

		struct timeval tp;
		struct rusage ruse;

		gettimeofday( &tp, NULL );
		double start_fold = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;

		getrusage( RUSAGE_SELF, &ruse );
		double startR_fold = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;

		if (seq1_str.length()==0)
		{
		    cout << "...fold seq 1..." << endl;

		    int len1=seq1.length();
		    char *c_str1 = new char[len1+1];
		    float my = fold(seq1.c_str(), c_str1);
		    seq1_str = (string)c_str1;
		    seq1_strname = seq1_name + "_mfe";
		}
		else
			seq1_strname=seq1_name+"_str";

		if (seq2_str.length()==0)
		{
			cout << "...fold seq 2..." << endl;
			int len2=seq2.length();
			char *c_str2 = new char[len2+1];
			float my = fold(seq2.c_str(), c_str2);
			seq2_str = (string)c_str2;
			seq2_strname=seq2_name+"_mfe";
		}
		else
			seq2_strname=seq2_name+"_str";

		gettimeofday( &tp, NULL );
		double end_fold = static_cast<double>( tp.tv_sec ) + static_cast<double>( tp.tv_usec )/1E6;

		getrusage( RUSAGE_SELF, &ruse );
		double endR_fold = static_cast<double>( ruse.ru_utime.tv_sec ) + static_cast<double>( ruse.ru_utime.tv_usec )/1E6;
		cout << endl << "time_wall folding = " << setprecision(3) << end_fold - start_fold << " sec" << endl;
		cout << "time_cpu folding = " << setprecision(3) << endR_fold - startR_fold << " sec" << endl << endl;

		myMol1 = Molecule(seq1_name,seq1,seq1_strname,seq1_str);
		myMol2 = Molecule(seq2_name,seq2,seq2_strname,seq2_str);
	}
	else
    {
       cerr << "Could not open File:" << myFASTAfile << endl;
       exit(1);
    }
	fa_file.close();
}
