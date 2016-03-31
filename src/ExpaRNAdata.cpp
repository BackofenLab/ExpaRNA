#include "ExpaRNAdata.h"
#include <math.h>
//--------------------------------------------------------------------------
// class ID for Identification of a pattern
//--------------------------------------------------------------------------
ID::ID()
{
}

ID::ID(const ID& myid)
         :molRef(myid.molRef),strRef(myid.strRef),id(myid.id)
{
}

ID::ID(const string& myMolRef,const string& myStrRef,const string& myId)
      :molRef(myMolRef),strRef(myStrRef),id(myId)
{
}

ID::~ID()
{
}

const string& ID::getMolRef() const
{
   return molRef;
}

const string& ID::getStrRef() const
{
   return strRef;
}

const string& ID::getId() const
{
   return id;
}


//--------------------------------------------------------------------------
// a pattern consist of elements
//--------------------------------------------------------------------------
SinglePattern::SinglePattern()
                              : Id(), pattern()
{
}

SinglePattern::SinglePattern(const SinglePattern &mySingP)
                              : Id(mySingP.Id),pattern(mySingP.pattern)
{
}


SinglePattern::~SinglePattern()
{
	pattern.clear();
}

SinglePattern::SinglePattern(const string& myMolRef,const string& myStrRef,const string& myId,const intVec& mySinglePattern)
   : Id(myMolRef,myStrRef,myId),pattern(mySinglePattern)
   {
   }

const ID& SinglePattern::getID() const
   {
      return Id;
   }

const intVec& SinglePattern::getPat() const
   {
      return pattern;
   }

void SinglePattern::print()
{
   cout << "SinglePatternObject Info:" << endl;
   cout << "  Pattern " << biu::Util_String::intvec2str(getPat(),":")<< endl;
   cout << " ID MolRef: " << Id.getMolRef() << endl;
   cout << " ID StrRef: " << Id.getStrRef() << endl;
   cout << " ID id: " << Id.getId() << endl;
   cout << "---------------------------------" << endl;
}

//--------------------------------------------------------------------------
// class PatternPair
//    is able to manage an EPM, consists of 2 singlepatterns, one in each RNA
//--------------------------------------------------------------------------
const string& PatternPair::getId() const
{
   return id;
};

const int& PatternPair::getSize() const
{
   return size;
};

const SinglePattern& PatternPair::getFirstPat() const
{
   return first;
};

const SinglePattern& PatternPair::getSecPat() const
{
   return second;
};

void PatternPair::resetBounds()
{
	  insideBounds.clear();
}

void PatternPair::setOutsideBounds(intPPair myPPair)
{
	  outsideBounds = myPPair;
};

intPPair PatternPair::getOutsideBounds()
{
	  return outsideBounds;
};

void PatternPair::addInsideBounds(intPPair myPPair)
{
	  insideBounds.push_back(myPPair);
};

const vector<intPPair>& PatternPair::getInsideBounds()
{
	  return insideBounds;
};

void PatternPair::setEPMScore(int myScore)
{
	  score = myScore;
};

void PatternPair::initEPMScore(ExpaRNAOptions::EPMscoringType EPMscoring)
{
	switch (EPMscoring)
	{
		case ExpaRNAOptions::BY_SIZE:

			this->EPMscore = this->getSize();

		   break;
		case ExpaRNAOptions::BY_QUADSIZE:

			//this->EPMscore =  this->getSize()+floor(pow(this->getSize(),1.6));
			this->EPMscore =  this->getSize() * this->getSize();
		break;
	}
	this->score = EPMscore;
};

int PatternPair::getScore()
      {
    	  return score;
      };

int PatternPair::getEPMScore()
      {
    	  return EPMscore;
      };



//--------------------------------------------------------------------------
// class PatternPairMap
//    is able to manage a set of PatternPairs(EPMs), each with 2 SinglePatterns
//--------------------------------------------------------------------------
PatternPairMap::PatternPairMap()
{
   idMap.clear();
   patternList.clear();
   patternOrderedMap.clear();
}

PatternPairMap::~PatternPairMap()
{
	idMap.clear();
	patternList.clear();
	patternOrderedMap.clear();
}

void PatternPairMap::add(const string& id,
                         const int& mysize,
                         const SinglePattern& first,
                         const SinglePattern& second)
{
   SelfValuePTR myP = SelfValuePTR(new PatternPair(id,mysize,first,second));
   patternList.push_back(myP);
   idMap.insert(make_pair(id,myP));
}

void PatternPairMap::add(const SelfValuePTR value)
{
   SelfValuePTR myP = SelfValuePTR(new PatternPair(*value));
   patternList.push_back(myP);
   idMap.insert(make_pair(value->getId(),myP));
}

void  PatternPairMap::makeOrderedMap()
{
   patternOrderedMap.clear();
   for(patListITER i = patternList.begin();i!=patternList.end();i++)
   {
      patternOrderedMap.insert(make_pair((*i)->getSize(),*i));
   }
}

void PatternPairMap::updateFromMap()
{
   if (!patternOrderedMap.empty())
   {
      idMap.clear();
      patternList.clear();
      for (orderedMapITER i=patternOrderedMap.begin();i!=patternOrderedMap.end();i++)
      {
         add(i->second);
      }
   }
}
const PatternPair& PatternPairMap::getPatternPair(const string& id)const
{
   return *(idMap.find(id)->second);
}

const    PatternPairMap::SelfValuePTR  PatternPairMap::getPatternPairPTR(const string& id)const
{
   return (idMap.find(id)->second);
}

const PatternPairMap::patListTYPE& PatternPairMap::getList()const
{
   return patternList;
}
const PatternPairMap::orderedMapTYPE& PatternPairMap::getOrderedMap() const
{
   return patternOrderedMap;
}

PatternPairMap::orderedMapTYPE& PatternPairMap::getOrderedMap2()
{
   return patternOrderedMap;
}

const int PatternPairMap::size() const
{
   return idMap.size();
}

int  PatternPairMap::getMapBases()
{
   int bases = 0;
   for(patListITER i = patternList.begin();i!=patternList.end();i++)
   {
      bases += (*i)->getSize();
   }
   return bases;
}
//--------------------------------------------------------------------------
// class Molecule
//       stores all information for one molecule/sequence
//--------------------------------------------------------------------------
Molecule::Molecule()
{
	gapMapping = false;
}

Molecule::Molecule(const string& myMolId,
                   const string& mySeq,
                   const string& myStrId,
                   const string& myStr)
                   :id(myMolId)
{
	this->gapMapping = false;
	this->sequences.insert(make_pair(myStrId,mySeq));
	this->structures.insert(make_pair(myStrId,myStr));
	this->actual_structure = myStrId;
	this->id_withGaps = actual_structure;
	this->makeBasePairs(getStructure());
}

Molecule::Molecule(const Molecule& myMol)
                     :id(myMol.id),
                      sequences(myMol.sequences),
                      structures(myMol.structures),
                      basePairs(myMol.basePairs),
                      actual_structure(myMol.actual_structure),
                      gapMapping(myMol.gapMapping),
                      id_withGaps(myMol.id_withGaps),
                      mapPos(myMol.mapPos)
{
}

Molecule::~Molecule()
{
	sequences.clear();
	structures.clear();
	basePairs.clear();
}

const string& Molecule::getId() const
   {
      return id;
   }

const string& Molecule::getSequence() const
   {
      return (sequences.find(actual_structure))->second;
   }


const string& Molecule::getStructure() const
   {
         return (structures.find(actual_structure))->second;
   }

const string& Molecule::getActualStructureId() const
   {
         return actual_structure;
   }

void  Molecule::setActualStructureId(const string& strId)
{
   if (structures.find(strId) != structures.end())
      actual_structure = strId;
}
const char Molecule::getStructure(const unsigned int& position) const
   {
      if ((position>0)&&(position<=sequences.find(actual_structure)->second.length()))
         return (structures.find(actual_structure)->second[position-1]);
      else return '#';
   }

const char Molecule::getBase(const unsigned int& position) const
   {
      if ((position>0)&&(position<=sequences.find(actual_structure)->second.length()))
         return (sequences.find(actual_structure)->second[position-1]);
      else return '#';
   }

const vector<int>& Molecule::getBasePairs() const
   {
         return ((basePairs.find(actual_structure))->second);
   }

const int Molecule::getBasePartner(const unsigned int& position) const
   {
      if ((position>0)&&(position<=sequences.find(actual_structure)->second.length())&&(basePairs.find(actual_structure)->second[position-1]!=-1))

         return ((basePairs.find(actual_structure))->second[position-1])+1;
      else return -1;
   }

const int   Molecule::getLength() const

{
   return sequences.find(actual_structure)->second.length();
}

void Molecule::printInfo(const bool& shortview)
{
   cout << "Molecule Information:" << endl;
   cout << "-----------------------" << endl;
   cout << "Molecule/Sequence Id: " << getId() << endl;
   cout << "Structure Id: " << getActualStructureId() << endl;
   cout << "actual Structure: " << actual_structure << endl;
   cout << "Sequence length: " << getSequence().length() << endl;
   if (!shortview)
   {
      cout << "Sequence: " << getSequence() << endl;
      cout << "Structure: " << getStructure() << endl;
      cout << "BasePairs: " << biu::Util_String::intvec2str(getBasePairs(),",") << endl;
   }
   cout << endl;
}

void Molecule::makeBasePairs(string myStr)
{
   vector<int> temp_bplist;
   stack<int> structure;

   if (myStr.size() != this->getLength())
   {
	   cerr << "...structure not correct. Seq/Struct of different lengths!" << endl;
	   actual_structure = "nostructure";
       exit(1);
   }

   temp_bplist.resize(myStr.size());
    for(unsigned int i=0;i<myStr.size();i++)
         if (myStr[i]=='(')
            structure.push(i);
        else
            if (myStr[i]==')')
            {
               if (structure.empty())
               {
                  cerr << "...structure not correct. Please use a correct structure!" << endl;
                  exit(1);
               }
                  temp_bplist[i]=structure.top(); // index basepair partner
                  temp_bplist[structure.top()]=i; // index basepair partner
                  structure.pop();
            }
            else
            temp_bplist[i]=-1;      // not paired
   if (!structure.empty())
   {
      cerr << "...structure not correct. Please use a correct structure!" << endl;
      actual_structure = "nostructure";
      exit(1);
   }
   basePairs.insert(make_pair(actual_structure,temp_bplist));
}

const bool Molecule::hasGapMapping()
{
	return gapMapping;
}

const string Molecule::getWithGapId()
{
	return id_withGaps;
}

void Molecule::setWithGapData()
{
	actual_structure = id_withGaps;
}

void Molecule::GapMapping()
{
	this->gapMapping = true;

	string tmpSeq,tmpStr;

	for (int i=0;i<getSequence().length();++i){
		if ( getSequence()[i] != '-' ) {
			tmpSeq.push_back(getSequence()[i]);
			tmpStr.push_back(getStructure(i+1));
			this->mapPos.push_back(i);
		}
	}

	this->id_withGaps = getActualStructureId();

	string newId = id_withGaps+"_$without_gaps$";

	this->sequences.insert(make_pair(newId,tmpSeq));
	this->structures.insert(make_pair(newId,tmpStr));
	this->setActualStructureId(newId);

	this->makeBasePairs(getStructure());
}

intVec Molecule::ReMapPattern(intVec myPattern)
{
	intVec newPos;
	for (int i=0;i<myPattern.size();++i)
	{
		newPos.push_back(mapPos[myPattern[i]-1]+1);
	}
	return newPos;

}
