// Util_String.cpp: Implementierung der Klasse Util_String.
/*

	<author>	Martin Mann					</author>
	<created>	9.10.2005					</created>

	<info>									</info>
*/
//
//////////////////////////////////////////////////////////////////////

#include "Util_String.h"


namespace biu {

///////////////////////////////////////////////////////////////////
// wandelt int in string um
///////////////////////////////////////////////////////////////////
std::string Util_String::int2str(const int& number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}


///////////////////////////////////////////////////////////////////
// wandelt int vector in String um. Trennzeichen ist delim.
///////////////////////////////////////////////////////////////////
std::string Util_String::intvec2str(const std::vector<int>& V, const std::string delim){
    std::ostringstream oss;
    copy(V.begin(), V.end(), std::ostream_iterator<int>(oss, delim.c_str()));
    std::string tmpstr;
    tmpstr = oss.str();
    if (tmpstr.length()>0) tmpstr.erase(tmpstr.end()-1);
    return tmpstr;
}

std::string Util_String::intvec2str(const std::vector<unsigned int>& V, const std::string delim){
    std::ostringstream oss;
    copy(V.begin(), V.end(), std::ostream_iterator<unsigned int>(oss, delim.c_str()));
    std::string tmpstr;
    tmpstr = oss.str();
    if (tmpstr.length()>0) tmpstr.erase(tmpstr.end()-1);
    return tmpstr;
}

///////////////////////////////////////////////////////////////////
// wandelt integer string in int um
///////////////////////////////////////////////////////////////////
int Util_String::str2int(const std::string& numString) {
	if (numString.empty())
		return 0;
	std::istringstream iss(numString);
	int retInt=0;
	iss >> retInt;
	return retInt;
}


///////////////////////////////////////////////////////////////////
// liefert anzahl der vorkommen von <c> in <str>
///////////////////////////////////////////////////////////////////
int Util_String::countChar( const std::string &str, const char c) {
	std::string::size_type  i=0;
	int ret = 0;
	for ( ;i<str.size();i++)
		if (str[i] == c)
			ret++;
	return ret;
}


///////////////////////////////////////////////////////////////////
// liefert laenge der laengsten wiederholung von <c> in <str>
///////////////////////////////////////////////////////////////////
int Util_String::maxSubseq( const std::string &str, const char c) {
	int m=0, act=0;
	std::string::size_type i=str.find_first_not_of(c);
	for(;i<str.find_last_not_of(c);i++) {
//		std::cerr <<i <<" = " <<str[i] ;
		if (str[i]==c)
			act++;
		else
			act=0;
		if (act > m) m = act;
//		std::cerr <<" act = " <<act <<" max = " << m <<"\n";
	}
	return m;
}


///////////////////////////////////////////////////////////////////
// konvertiert den string in grossbuchstaben
///////////////////////////////////////////////////////////////////
std::string Util_String::str2upperCase(std::string &str) {
	std::string ret(str);
	for (std::string::size_type i=0; i<ret.size(); i++)
		if (ret[i]>96 && ret[i]<123)
			ret[i] = ret[i]-'a'+'A';
	return ret;
}


///////////////////////////////////////////////////////////////////
// prueft ob string nur aus elementen des alphabets besteht
///////////////////////////////////////////////////////////////////
bool Util_String::isAlphStr(const std::string str, const std::string alph) {
	bool isOK = true;
	for (std::string::size_type i=0; i< str.length() && isOK; i++) {
		isOK = alph.find(str[i]) != std::string::npos;	// str[i] exists in alphabet
	}
	return isOK;
}

///////////////////////////////////////////////////////////////////
// entfernt fuehrende und endende "\n" und " " und "\t"
///////////////////////////////////////////////////////////////////
std::string Util_String::chompStr(const std::string &str) {
	std::string ret(str);

	ret = ret.substr(ret.find_first_not_of("\n"));
	ret = ret.substr(0,ret.find_last_not_of("\n")+1);
	ret = ret.substr(ret.find_first_not_of("\r"));		// carriage return
	ret = ret.substr(0,ret.find_last_not_of("\r")+1);
	ret = ret.substr(ret.find_first_not_of(" "));
	ret = ret.substr(0,ret.find_last_not_of(" ")+1);
	ret = ret.substr(ret.find_first_not_of("\t"));
	ret = ret.substr(0,ret.find_last_not_of("\t")+1);
	return ret;
}

} // namespace biu
