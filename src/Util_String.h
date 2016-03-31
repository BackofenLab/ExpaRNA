// Util_String.h:
/*

	<author>	Martin Mann					</author>
	<created>	9.10.2005					</created>

	<info>									</info>
*/
//
//////////////////////////////////////////////////////////////////////

#if !defined(IT_UTIL_STRING_H__INCLUDED)
#define IT_UTIL_STRING_H__INCLUDED

#include <sstream>
#include <string>
#include <vector>
#include <iterator>

namespace biu {

	/////////////////////////////////
	// string utility klasse
	/////////////////////////////////
	class Util_String
	{
	public:
		Util_String()
		{}
		virtual ~Util_String()
		{}


			// typecast
		static std::string int2str(const int& number);
		static int str2int(const std::string& numString);
			// wandelt container vector int in string um mit Trennzeichen
		static std::string intvec2str(const std::vector<int>& V, const std::string delim);
		static std::string intvec2str(const std::vector<unsigned int>& V, const std::string delim);

			// liefert anzahl der vorkommen von <c> in <str>
		static int countChar( const std::string &str, const char c);

			// liefert laenge der laengsten wiederholung von <c> in <str>
		static int maxSubseq( const std::string &str, const char c);

			// konvertiert den string in grossbuchstaben
		static std::string str2upperCase(std::string &str);

			// prueft ob string nur aus elementen des alphabets besteht
		static bool isAlphStr(const std::string str, const std::string alph) ;

			// entfernt fuehrende und endende "\n" und " " und "\t"
		static std::string chompStr(const std::string &str);


	};


} // namspace biu

#endif // !define(IT_UTIL_STRING_H__INCLUDED)
