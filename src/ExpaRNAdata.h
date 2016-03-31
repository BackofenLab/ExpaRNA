#ifndef DATA_SINGLEPATTERN_H_
#define DATA_SINGLEPATTERN_H_

#include <vector>
#include <string>
#include <stack>
#include <ext/hash_map>
//#include <unordered_map> // fur new C++ standard compiler option: -std=gnu++0x
#include <list>
#include <map>
#include <deque>
#include <iostream>
#include <utility>
#include "Util_String.h"
#include <stdlib.h>
#include <stdio.h>
#include "ExpaRNAOptions.h"

using namespace std::rel_ops;
using namespace std;

    typedef vector<unsigned int>                    intVec;
    typedef vector<unsigned int>::const_iterator	intVecCITER;
    typedef vector<unsigned int>::iterator			intVecITER;
    typedef pair<unsigned int,unsigned int>         intPair;
	typedef pair<intPair, intPair>                  intPPair;
	typedef vector<intPPair>::const_iterator        IntPPairCITER;

//--------------------------------------------------------------------------
// class ID
//--------------------------------------------------------------------------
class ID
{
   public:

      ID();
      ID(const string& myMolRef,const string& myStrRef,const string& myId);
      ID(const ID &myid);

      virtual ~ID();

      const string& getMolRef()   const;
      const string& getStrRef()   const;
      const string& getId()       const;

   private:

      std::string molRef;
      std::string strRef;
      std::string id;
};

class IDHash
{
   public:
      size_t operator()(const ID &myID) const
      {
         unsigned long hash = 5381;
         std::string myStr = myID.getId() + myID.getMolRef() + myID.getStrRef();

         for (unsigned int i = 0; i < myStr.length(); i++)
         {
            hash = ((hash << 5) + hash) + myStr[i]; // hash * 33 + cStateString[i]
         }

         return hash;
   }
};

class IDEq {
public:
   bool operator()(const ID &a,const ID &b) const
   {
      return (a.getId() == b.getId());
   }
};

class StringHash
{
   public:
      size_t operator()(const string &myStr) const
      {
         unsigned long hash = 5381;

         for (unsigned int i = 0; i < myStr.length(); i++)
         {
            hash = ((hash << 5) + hash) + myStr[i]; // hash * 33 + cStateString[i]
         }

         return hash;
   }
};

class StringEq {
public:
   bool operator()(const string &a,const string &b) const
   {
      return (a == b);
   }
};

//--------------------------------------------------------------------------
// class SinglePattern
//     is set of positions (pattern) in one RNA structure
//--------------------------------------------------------------------------
class SinglePattern
{
public:

      SinglePattern();
      SinglePattern(const SinglePattern& mySingP);
      SinglePattern(const string& myMolRef,const string& myStrRef,const string& myId,const intVec& mySinglePattern);

	virtual ~SinglePattern();

   const ID&            getID()  const;
   const intVec&        getPat() const;
         void           print();

private:

      ID          Id;
      intVec         pattern;
};

//--------------------------------------------------------------------------
// class PatternPair
//    is able to manage an EPM, consists of 2 singlepatterns, one in each RNA
//--------------------------------------------------------------------------
class PatternPair
   {
      public:
      PatternPair(){};
      PatternPair(const PatternPair &myPair)
                     :id(myPair.id),size(myPair.size),first(myPair.first),second(myPair.second),
                      score(myPair.score), EPMscore(myPair.EPMscore)
      { };

      PatternPair(const string& myId,const int& mySize,const SinglePattern& myFirstPat,const SinglePattern& mySecPat)
                  : id(myId),size(mySize),first(myFirstPat),second(mySecPat)
      {  score = 0;  };

      virtual ~PatternPair()
      {
    	insideBounds.clear();
      };

      const string& 		getId() const;
      const int& 			getSize() const;
      const SinglePattern& 	getFirstPat() const;
      const	SinglePattern& 	getSecPat() const;
			void			resetBounds();
			void			setOutsideBounds(intPPair myPPair);
			intPPair 		getOutsideBounds();
			void			addInsideBounds(intPPair myPPair);
      const vector<intPPair>& getInsideBounds();
			void			setEPMScore(int myScore);
			void			initEPMScore(ExpaRNAOptions::EPMscoringType EPMscoring);
			int 			getScore();
			int 			getEPMScore();

      private:
         string         	id;
         int            	size;
         SinglePattern  	first;
         SinglePattern  	second;

         int				score;
         int				EPMscore;
         vector<intPPair>   insideBounds;
         intPPair           outsideBounds;
   };

//--------------------------------------------------------------------------
// class PatternPairMap
//    manage a set of EPMs (PatternPair)
//--------------------------------------------------------------------------
class PatternPairMap
{
   public:
	  typedef  PatternPair                                  selfValueTYPE;
	  typedef  PatternPair*				               		SelfValuePTR;

      typedef  multimap<int,SelfValuePTR,greater<int> >     orderedMapTYPE;
      typedef  orderedMapTYPE::const_iterator               orderedMapCITER;
      typedef  orderedMapTYPE::iterator                     orderedMapITER;
      typedef  list<SelfValuePTR>                           patListTYPE;
      typedef  patListTYPE::iterator                        patListITER;
      typedef  patListTYPE::const_iterator                  patListCITER;
      typedef  __gnu_cxx::hash_map<string,SelfValuePTR,StringHash,StringEq> PatternIdMapTYPE;
      //typedef  unordered_map<string,SelfValuePTR,StringHash,StringEq> PatternIdMapTYPE;


         PatternPairMap();
         PatternPairMap(const PatternPairMap& myPairMap)
                           :patternList(myPairMap.patternList),
                            patternOrderedMap(myPairMap.patternOrderedMap),
                            idMap(myPairMap.idMap)  {};

      virtual ~PatternPairMap();

               void              add( const string& id,
                                      const int& mysize,
                                      const SinglePattern& first,
                                      const SinglePattern& second);
               void              add(const SelfValuePTR value);
               void              makeOrderedMap();
               void              updateFromMap();
      const    PatternPair&      getPatternPair(const string& id)const;
      const    SelfValuePTR      getPatternPairPTR(const string& id)const;
      const    patListTYPE&      getList() const;
      const    orderedMapTYPE&   getOrderedMap() const;
               orderedMapTYPE&   getOrderedMap2();
      const    int               size()   const;
			   int  			 getMapBases();

   private:

     patListTYPE        patternList;
     orderedMapTYPE     patternOrderedMap;
     PatternIdMapTYPE   idMap;
};

//--------------------------------------------------------------------------
// class Molecule
//       stores all information for one molecule/sequence
//--------------------------------------------------------------------------
class Molecule
{

   public:
							Molecule();
							Molecule(const Molecule& myMol);
							Molecule(const string& myMolId,
                                    const string& mySeq,
                                    const string& myStrId,
                                    const string& myStr);
      virtual             	~Molecule();

      const string&        	getId()                                      const;
      const string&        	getSequence()                                const;
      const char           	getBase(const unsigned int& position)        const;
      const string&        	getStructure()                               const;
      const string&        	getActualStructureId()                       const;
      const char           	getStructure(const unsigned int& position)   const;
      const vector<int>&   	getBasePairs()                               const;
      const int            	getBasePartner(const unsigned int& position) const;
      const int            	getLength()                                  const;

            void           	setActualStructureId(const string& strId);
            void           	printInfo(const bool& shortview);
			void			GapMapping();
			intVec			ReMapPattern(intVec myPattern);
	  const bool			hasGapMapping();
	  const string			getWithGapId();
			void 			setWithGapData();

   private:

      typedef __gnu_cxx::hash_map <string,string,StringHash,StringEq>      MAPTYPE_Structures;
      typedef __gnu_cxx::hash_map <string,vector<int>,StringHash,StringEq> MAPTYPE_BasePairs;
	  //typedef unordered_map <string,string,StringHash,StringEq>      MAPTYPE_Structures;
	  //typedef unordered_map <string,vector<int>,StringHash,StringEq> MAPTYPE_BasePairs;

		    bool 				 gapMapping;

            string               id;
            string               id_withGaps;
            string               actual_structure;
            MAPTYPE_Structures   structures;
            MAPTYPE_Structures   sequences;
            MAPTYPE_BasePairs    basePairs;
            intVec 				 mapPos;

            void                 makeBasePairs(string myStr);
};

#endif /*DATA_SINGLEPATTERN_H_*/
