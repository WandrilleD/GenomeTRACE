/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/*

This file contains various functions used by DeCo

Created the: 02-03-2016
by: Wandrille Duchemin

Last modified the: 11-01-2018
by: Wandrille Duchemin

*/


#include "DeCoUtils.h"

/*
read the reconcilied tree in a recphyloXML file and create a GeneFamily for each of them

Takes:
	- GeneFamilyList (vector <GeneFamily *> ) : gene families 
	- string fileName : name of the file

*/
void ReadRecPhyLoXMLFile(  vector <GeneFamily *> * GeneFamilyList, string fileName , MySpeciesTree * Stree, bool verbose, bool superverbose)
{
	ifstream fileStream(fileName.c_str());
	if( !fileStream.is_open() ) 
	{
		throw Exception("ReconciledTree::ReconciledTree : could not open reconciled tree file : " + fileName);
		exit(1);
	}

	string RecGeneTreeTag = "recGeneTree";
	string startTag = "clade";

	while( goToNextTag(fileStream,RecGeneTreeTag) )  //goes to the next RecGeneTreeTag -> next reconciled gene tree
	{
		if( goToNextTag(fileStream,startTag) ) // go to clade -> the root of the reconciled gene tree
		{
			ReconciledTree Rtree(fileStream, Stree, superverbose);
			GeneFamilyList->push_back(new GeneFamily(Rtree, verbose, superverbose));
		}
		else
		{ // no clade ? pb but ignore...
			if(verbose)
				cout << "found a "<< RecGeneTreeTag << " without any " << startTag << " in the file " << fileName << " -> ignoring that tree."<< endl;
		}
	}

}


/* return false if error, true if not */
bool StrToDouble(string s, double & d) 
{
	istringstream ss(s);

	return static_cast<bool>( ss >> d );
}








void AddToFile(string filename, string line)
{
	ofstream ofs;
	ofs.open(filename.c_str(),ofstream::out | ofstream::app);
	ofs << line << endl;
	ofs.close();	
}


/*
Makes a map linking each leaf name to a genefamily ID

Takes:
 - vector <GeneFamily *> * GeneFamilyList 

Returns:
	map <string, int>	
*/
map <string, int> makeLeafToGFMap(vector <GeneFamily *> * GeneFamilyList)
{
	map <string, int> LeafToGFMap;


	for(unsigned i = 0 ; i < GeneFamilyList->size(); i++)
	{
		vector < string > leafNames = GeneFamilyList->at(i)->getleafNames();
		for(unsigned j = 0 ; j < leafNames.size(); j++)
			LeafToGFMap[leafNames[j]] = i;
	}
	return LeafToGFMap;
}


/*
Adds the leaf of the reconciled tree to the map

Takes:
	- ReconciledTree * rtree
	- map <string,int> &LeafToSpMap

*/
void fillLeafToSpMap(ReconciledTree * rtree, map <string,int> &LeafToSpMap)
{
	vector <int> leaves = rtree->getLeavesId();

	for(unsigned i = 0 ; i < leaves.size(); i++)
	{
		if( rtree->isRealLeaf(leaves[i]) )
		{
			LeafToSpMap[ rtree->getNodeName(leaves[i]) ] = rtree->getNodeSpecies(leaves[i]) ;
		}
	}
}


/**
 * Print the tree root at node to the file.
 */
void printNewick( MyGeneNode *node, string fileName ) 
{
	string newickStr = TreeTemplateTools::nodeToParenthesis( *node );	
	ios_base::openmode mode = ios::out | ios::binary;
	mode |= ios::app; // append
	ofstream out( fileName.c_str(), mode );
	out << newickStr << ";" << endl;
	out.close();
}
//species version
void printNewick( MySpeciesNode *node, string fileName ) 
{
	string newickStr = TreeTemplateTools::nodeToParenthesis( *node );	
	ios_base::openmode mode = ios::out | ios::binary;
	mode |= ios::app; // append
	ofstream out( fileName.c_str(), mode );
	out << newickStr << ";" << endl;
	out.close();
}




/*
Takes:
	- GeneFamilyList (vector <GeneFamily *> ) : gene families 
	- TopoWeight (double) [default = 1] : weight of the topology part of the score

Returns:
	(double): score associated with topologies
*/
double computeTopoScore( vector <GeneFamily *> GeneFamilyList, double TopoWeight )
{
	double s = 0;
	for(unsigned i = 0; i< GeneFamilyList.size(); i++)
	{
		s -= log( GeneFamilyList[i]->getNormalizedTreeLikelihood());
	}
	s *= TopoWeight;
	return s;
}

/*
Takes:
	- GeneFamilyList (vector <GeneFamily *> ) : gene families 
	- ReconWeight (double) [default = 1] : weight of the reconciliation part of the score

Returns:
	(double): score associated with reconciliations
*/
double computeReconciliationScore(vector <GeneFamily *> GeneFamilyList, double ReconWeight )
{
	double s = 0;
	for(unsigned i = 0; i< GeneFamilyList.size(); i++)
	{
		s +=  GeneFamilyList[i]->getRecScore();
	}
	s *= ReconWeight;
	return s;
}

/*
Takes:
	- RefinedEquivalenceClasses (vector <EquivalenceClass> *) : equivalence classes
	- AGainCost (double) : cost of a single gain
	- ABreakCost (double) : cost of a single adjBreak
	- AdjWeight (double) [default = 1 ] : weight of the adjacency part of the score

Returns:
	(double): score associated with adjacency histories
*/
double computeAdjacenciesScore(vector <EquivalenceClass> * RefinedEquivalenceClasses, double AGainCost, double ABreakCost, double AdjWeight )
{
	double s = 0;
	for(unsigned i = 0; i< RefinedEquivalenceClasses->size(); i++)
	{
		s +=  RefinedEquivalenceClasses->at(i).getNbAdjGain() * AGainCost;
		s +=  RefinedEquivalenceClasses->at(i).getNbAdjBreak() * ABreakCost;
	}
	s *= AdjWeight;
	return s;	
}

/*
Takes:
	- ECFams (vector <EquivalenceClassFamily > *) : equivalence class families
	- AGainCost (double) : cost of a single gain
	- ABreakCost (double) : cost of a single adjBreak
	- AdjWeight (double) [default = 1 ] : weight of the adjacency part of the score

Returns:
	(double): score associated with adjacency histories
*/
double computeAdjacenciesScore(vector <EquivalenceClassFamily> * ECFams, double AGainCost, double ABreakCost, double AdjWeight )
{
	double s = 0;
	for(unsigned i = 0; i< ECFams->size(); i++)
	{
		s +=  ECFams->at(i).getNbAdjGain() * AGainCost;
		s +=  ECFams->at(i).getNbAdjBreak() * ABreakCost;
	}
	s *= AdjWeight;
	return s;	
}

/*
Takes:
	- Lcoevent (vector <CoEvent> ) : list of coevents accross all gene families
	- DupCost (double) : cost of a single duplication
	- LossCost (double) : codt of a single loss
	- HGTCost (double) : cost of a single HGT
	- ReconWeight (double) [default = 1] : weight of the reconciliation part of the score

Returns:
	(double): score balancing when account for coevents
*/
double computeCoEventScore(vector <CoEvent> Lcoevent, double DupCost, double LossCost, double HGTCost, double ReconWeight  )
{
	double s = 0;
	for(unsigned i = 0; i< Lcoevent.size(); i++)
	{
		s += Lcoevent[i].getLinearScore(DupCost,HGTCost,LossCost);
	}
	s *= ReconWeight;
	return s;
}

/*
Takes:
	- GeneFamilyList (vector <GeneFamily *> ) : gene families 
	- RefinedEquivalenceClasses (vector <EquivalenceClass> *) : equivalence classes
	- Lcoevent (vector <CoEvent> ) : list of coevents accross all gene families
	- AGainCost (double) : cost of a single gain
	- ABreakCost (double) : cost of a single adjBreak
	- DupCost (double) : cost of a single duplication
	- LossCost (double) : codt of a single loss
	- HGTCost (double) : cost of a single HGT
	- TopoWeight (double) [default = 1] : weight of the topology part of the score
	- ReconWeight (double) [default = 1] : weight of the reconciliation part of the score
	- AdjWeight (double) [default = 1 ] : weight of the adjacency part of the score

Returns:
	(double): score of the whole system
*/
double computeSystemScore( vector <GeneFamily *> GeneFamilyList, vector <EquivalenceClass> *RefinedEquivalenceClasses, vector <CoEvent> Lcoevent, double AGainCost, double ABreakCost, double DupCost, double LossCost, double HGTCost, double TopoWeight , double ReconWeight , double AdjWeight  )
{
	double s = 0;

	s += computeTopoScore( GeneFamilyList,  TopoWeight);
	s += computeReconciliationScore( GeneFamilyList, ReconWeight );
	s += computeAdjacenciesScore( RefinedEquivalenceClasses, AGainCost, ABreakCost, AdjWeight);
	s += computeCoEventScore(Lcoevent, DupCost, LossCost, HGTCost,  ReconWeight );
	return s;	
}

/*
Takes:
	- GeneFamilyList (vector <GeneFamily *> ) : gene families 
	- ECFams (vector <EquivalenceClassFamily > *) : equivalence class families
	- Lcoevent (vector <CoEvent> ) : list of coevents accross all gene families
	- AGainCost (double) : cost of a single gain
	- ABreakCost (double) : cost of a single adjBreak
	- DupCost (double) : cost of a single duplication
	- LossCost (double) : codt of a single loss
	- HGTCost (double) : cost of a single HGT
	- TopoWeight (double) [default = 1] : weight of the topology part of the score
	- ReconWeight (double) [default = 1] : weight of the reconciliation part of the score
	- AdjWeight (double) [default = 1 ] : weight of the adjacency part of the score

Returns:
	(double): score of the whole system
*/
double computeSystemScore( vector <GeneFamily *> GeneFamilyList, vector <EquivalenceClassFamily> * ECFams, vector <CoEvent> Lcoevent, double AGainCost, double ABreakCost, double DupCost, double LossCost, double HGTCost, double TopoWeight , double ReconWeight , double AdjWeight )
{
	double s = 0;

	s += computeTopoScore( GeneFamilyList,  TopoWeight);
	s += computeReconciliationScore( GeneFamilyList, ReconWeight );
	s += computeAdjacenciesScore( ECFams, AGainCost, ABreakCost, AdjWeight);
	s += computeCoEventScore(Lcoevent, DupCost, LossCost, HGTCost,  ReconWeight );
	return s;	
}

/*
Takes:
	- adjacencies (vector< pair <string,string > >)
	- GeneFamilyList (vector <GeneFamily *> *) : RECONCILIATION MUST BE COMPLETE
	- DeCoLTrefine (bool) : if true, DeCoLT refining method will used even if there are no HGT in the trees
	- Verbose (bool)
	- SuperVerbose (bool)

Returns:
 - (vector < EquivalenceClass >) : a vector of refined Equivalence classes
*/
vector < EquivalenceClass > CreateEquivalenceClasses(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList, bool DeCoLTrefine,bool Verbose, bool SuperVerbose)
{
	vector <EquivalenceClass> EquivalenceClasses;
	map < int, map <int , int > > Gfam1ToGfam2ToIndex; // to easily find the already created EquivalenceClasses
	map < int, map <int , int > >::iterator it1;
	map <int , int >::iterator it2;
	for(unsigned i; i < adjacencies.size(); i++)
	{
		pair <string, string> currentadj = adjacencies[i];
		//Finding which Gfams are concerned.
		int gfam1 = -1;
		int gfam2 = -1;
		for(unsigned j = 0; j < GeneFamilyList->size(); j++)
		{
			GeneFamily * GF = GeneFamilyList->at(j);
			if(gfam1 ==-1)
				if(GF->hasLeaf(currentadj.first))
					gfam1 = j;
			if(gfam2 == -1)
				if(GF->hasLeaf(currentadj.second))
					gfam2 = j;
			if((gfam1 !=-1) && (gfam2 != -1))//both gfam were found
				break;
		}
		if((gfam1 ==-1) || (gfam2 == -1))// at least one of the leaves was not found. Do nothing
		{
			if(Verbose)
				cout << "The adjacency " << currentadj.first << " " << currentadj.second << " had a member which did not correspond to any gene family. Ignoring that adjacency." << endl;
			continue;
		}
		//case were gfam1 and gfam2 were found
		if(SuperVerbose)
			cout << "Adj " << currentadj.first << " - " << currentadj.second << " families " << gfam1 << " " << gfam2 << endl;
		if(gfam1 > gfam2) // we will invert the values so that gfam1 is always the smallest id
		{
			int junk;
			string junkstr;
			junk = gfam1;
			gfam1 = gfam2;
			gfam2 = junk;
			junkstr = currentadj.first;
			currentadj.first = currentadj.second;
			currentadj.second = junkstr;
		}
		int ECindex;
		it1 = Gfam1ToGfam2ToIndex.find(gfam1);
		if(it1 == Gfam1ToGfam2ToIndex.end()) //gfam1 not found -> creating the map
		{
			map <int , int > newmap;
			Gfam1ToGfam2ToIndex.insert(pair < int , map <int , int > > (gfam1,newmap) );
		}
		it2 = Gfam1ToGfam2ToIndex[gfam1].find(gfam2);
		if(it2 == Gfam1ToGfam2ToIndex[gfam1].end()) //gfam2 could not be found in association with gfam1 -> we have to create a new EquivalenceClass instance
		{
			if(SuperVerbose)
				cout << "Creating adjacency class for " << gfam1 << " " << gfam2 << endl;
			EquivalenceClasses.push_back( EquivalenceClass(gfam1,gfam2));
			Gfam1ToGfam2ToIndex[gfam1][gfam2] = EquivalenceClasses.size() -1 ;
			ECindex =  EquivalenceClasses.size() -1 ;
		}
		else
			ECindex	= Gfam1ToGfam2ToIndex[gfam1][gfam2];
		bool success = EquivalenceClasses[ECindex].CheckAddAdj(currentadj.first,currentadj.second,gfam1,gfam2);//adding the adj to the correct Equivalence class
		if(SuperVerbose)
			if(!success)
				cout << "Could not add adj to the Equivalence Class " << ECindex << " " << EquivalenceClasses[ECindex].getGfamily1() << " - " << EquivalenceClasses[ECindex].getGfamily2() << endl;
	}
	if(Verbose)
		cout << "First round of Equivalence Class construction: " << EquivalenceClasses.size() << " Classes built." << endl;
	//Second round: refining equivalence classes
	vector <EquivalenceClass> RefinedEquivalenceClasses;
	for(unsigned i = 0; i <  EquivalenceClasses.size(); i++)
	{
		int gfam1 = EquivalenceClasses[i].getGfamily1();
		int gfam2 = EquivalenceClasses[i].getGfamily2();
		ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
		ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();
		if(SuperVerbose)
		{
			cout << "refining equivalence class " << i; 
			EquivalenceClasses[i].printMe(false);
		}
			
		vector<EquivalenceClass *> newEqClass = EquivalenceClasses[i].refineEqClass(Rtree1, Rtree2, DeCoLTrefine, SuperVerbose); 
		for(unsigned j =0; j < newEqClass.size(); j++)
			RefinedEquivalenceClasses.push_back(*newEqClass[j]);
	}
	//cleaning
	EquivalenceClasses.clear();
	Gfam1ToGfam2ToIndex.clear();
	if(Verbose)
		cout << "Second round of Equivalence Class construction: " << RefinedEquivalenceClasses.size() << " classes after refinement." << endl;
	return RefinedEquivalenceClasses;
}






/*
Will create an Equivalence class for each pair of gene families (including with themselves). Refinement considering all possible adjacencies (including non observed ones)

Takes:
	- adjacencies (vector< pair <string,string > >)
	- GeneFamilyList (vector <GeneFamily *> *): RECONCILIATION MUST BE COMPLETE
	- Verbose (bool)
	- SuperVerbose (bool)

Returns:
 - (vector < EquivalenceClass >) : a vector of refined Equivalence classes
*/
vector < EquivalenceClass > CreateAllPairEquivalenceClasses(vector< pair <string,string > > adjacencies, vector <GeneFamily *> * GeneFamilyList,bool Verbose, bool SuperVerbose)
{

	vector <EquivalenceClass> EquivalenceClasses;
	map < int, map <int , int > > Gfam1ToGfam2ToIndex; // to easily find the already created EquivalenceClasses

	for(unsigned i = 0; i < GeneFamilyList->size(); i++)
	{
		for(unsigned j = i ; j < GeneFamilyList->size() ; j++)
		{
			EquivalenceClasses.push_back(EquivalenceClass(i,j));
			Gfam1ToGfam2ToIndex[i][j] = EquivalenceClasses.size() -1 ;

			EquivalenceClasses.back().setAncestor( 0, GeneFamilyList->at(i)->getRecTree()->getRootId());
			EquivalenceClasses.back().setAncestor( 1, GeneFamilyList->at(j)->getRecTree()->getRootId());
		}
	}


	map < int, map <int , int > >::iterator it1;
	map <int , int >::iterator it2;
	for(unsigned i; i < adjacencies.size(); i++)
	{
		pair <string, string> currentadj = adjacencies[i];
		//Finding which Gfams are concerned.
		int gfam1 = -1;
		int gfam2 = -1;
		for(unsigned j = 0; j < GeneFamilyList->size(); j++)
		{
			GeneFamily * GF = GeneFamilyList->at(j);
			if(gfam1 ==-1)
				if(GF->hasLeaf(currentadj.first))
					gfam1 = j;
			if(gfam2 == -1)
				if(GF->hasLeaf(currentadj.second))
					gfam2 = j;
			if((gfam1 !=-1) && (gfam2 != -1))//both gfam were found
				break;
		}
		if((gfam1 ==-1) || (gfam2 == -1))// at least one of the leaves was not found. Do nothing
		{
			if(Verbose)
				cout << "The adjacency " << currentadj.first << " " << currentadj.second << " had a member which did not correspond to any gene family. Ignoring that adjacency." << endl;
			continue;
		}
		//case were gfam1 and gfam2 were found
		if(SuperVerbose)
			cout << "Adj " << currentadj.first << " - " << currentadj.second << " families " << gfam1 << " " << gfam2 << endl;
		if(gfam1 > gfam2) // we will invert the values so that gfam1 is always the smallest id
		{
			int junk;
			string junkstr;
			junk = gfam1;
			gfam1 = gfam2;
			gfam2 = junk;
			junkstr = currentadj.first;
			currentadj.first = currentadj.second;
			currentadj.second = junkstr;
		}
		int ECindex	= Gfam1ToGfam2ToIndex[gfam1][gfam2];

		bool success = EquivalenceClasses[ECindex].CheckAddAdj(currentadj.first,currentadj.second,gfam1,gfam2);//adding the adj to the correct Equivalence class
		if(SuperVerbose)
			if(!success)
				cout << "Could not add adj to the Equivalence Class " << ECindex << " " << EquivalenceClasses[ECindex].getGfamily1() << " - " << EquivalenceClasses[ECindex].getGfamily2() << endl;
	}

	//Second round: refining equivalence classes
	vector <EquivalenceClass> RefinedEquivalenceClasses;
	for(unsigned i = 0; i <  EquivalenceClasses.size(); i++)
	{
		int gfam1 = EquivalenceClasses[i].getGfamily1();
		int gfam2 = EquivalenceClasses[i].getGfamily2();
		ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
		ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();
		if(SuperVerbose)
		{
			cout << "refining equivalence class " << i; 
			EquivalenceClasses[i].printMe(false);
		}
			
		vector<EquivalenceClass *> newEqClass = EquivalenceClasses[i].refineEqClassWhole(Rtree1, Rtree2, SuperVerbose); 
		for(unsigned j =0; j < newEqClass.size(); j++)
			RefinedEquivalenceClasses.push_back(*newEqClass[j]);
	}
	//cleaning
	EquivalenceClasses.clear();
	Gfam1ToGfam2ToIndex.clear();
	if(Verbose)
		cout << "Second round of Equivalence Class construction: " << RefinedEquivalenceClasses.size() << " classes after refinement." << endl;
	return RefinedEquivalenceClasses;

}

/*
Takes:
	- adjacencies (vector< pair <string,string > >)
	- GeneFamilyList (vector <GeneFamily *> *)
	- DeCoLTrefine (bool) : if true, DeCoLT refining method will used even if there are no HGT in the trees
	- useWholeClass (bool): if true, all possible adjacencies (including non observed ones), will be used for refinement
	- Verbose (bool)
	- SuperVerbose (bool)

Returns:
 - (vector < EquivalenceClassFamily > * ) : a pointer to a vector of refined Equivalence class families
*/
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies, 
																	vector <GeneFamily *> * GeneFamilyList, 
																	bool DeCoLTrefine, bool useWholeClass, 
																	bool Verbose, bool SuperVerbose)
{

	// 1. creating all EquivalenceClassFamilies

	vector <EquivalenceClassFamily> * ECFams = new vector <EquivalenceClassFamily>;
	map < int, map <int , int > > Gfam1ToGfam2ToIndex; // to easily find the already created EquivalenceClasses

	for(unsigned i = 0; i < GeneFamilyList->size(); i++)
	{
		for(unsigned j = i ; j < GeneFamilyList->size() ; j++)
		{
			ECFams->push_back(EquivalenceClassFamily(i,j));
			Gfam1ToGfam2ToIndex[i][j] = ECFams->size() -1 ;

			//EquivalenceClasses.back().setAncestor( 0, GeneFamilyList->at(i)->getRecTree()->getRootId());
			//EquivalenceClasses.back().setAncestor( 1, GeneFamilyList->at(j)->getRecTree()->getRootId());
		}
	}

	cout << "treating "<< adjacencies.size() << "adjacencies."<<endl;

	//2. going through all leaves an adding them to the corresponding EquivalenceClassFamily

	map < int, map <int , int > >::iterator it1;
	map <int , int >::iterator it2;

	for(unsigned i= 0; i < adjacencies.size(); i++)
	{
		pair <string, string> currentadj = adjacencies[i];
		//Finding which Gfams are concerned.
		int gfam1 = -1;
		int gfam2 = -1;
		for(unsigned j = 0; j < GeneFamilyList->size(); j++)
		{
			GeneFamily * GF = GeneFamilyList->at(j);
			if(gfam1 ==-1)
				if(GF->hasLeaf(currentadj.first))
					gfam1 = j;
			if(gfam2 == -1)
				if(GF->hasLeaf(currentadj.second))
					gfam2 = j;
			if((gfam1 !=-1) && (gfam2 != -1))//both gfam were found
				break;
		}
		if((gfam1 ==-1) || (gfam2 == -1))// at least one of the leaves was not found. Do nothing
		{
			if(Verbose)
				cout << "The adjacency " << currentadj.first << " " << currentadj.second << " had a member which did not correspond to any gene family. Ignoring that adjacency." << endl;
			continue;
		}
		//case were gfam1 and gfam2 were found
		if(SuperVerbose)
			cout << "Adj " << currentadj.first << " - " << currentadj.second << " families " << gfam1 << " " << gfam2 << endl;
		if(gfam1 > gfam2) // we will invert the values so that gfam1 is always the smallest id
		{
			int junk;
			string junkstr;
			junk = gfam1;
			gfam1 = gfam2;
			gfam2 = junk;
			junkstr = currentadj.first;
			currentadj.first = currentadj.second;
			currentadj.second = junkstr;
		}
		int ECindex	= Gfam1ToGfam2ToIndex[gfam1][gfam2];

		bool success = ECFams->at(ECindex).CheckAddAdj(currentadj.first,currentadj.second,gfam1,gfam2);//adding the adj to the correct Equivalence class
		if(SuperVerbose)
			if(!success)
				cout << "Could not add adj to the Equivalence Class " << ECindex << " " << ECFams->at(ECindex).getGfamily1() << " - " << ECFams->at(ECindex).getGfamily2() << endl;

		//cout << i << endl;
	}


	//2 bis: delete all equivalence class family with 0 adjacencies if !useWholeClass"
	if(!useWholeClass)
	{
		for(int i = ECFams->size() -1 ; i >= 0; i--)
		{

			if(ECFams->at(i).getNbAdj() == 0 ) // no adjs
			{

 				ECFams->erase(ECFams->begin()+i); // deleting that ECFamily
			}

		}
	}


	//3. refining equivalence classes

	cout << "refining."<<endl;
	
	for(unsigned i = 0; i <  ECFams->size(); i++)
	{
		int gfam1 = ECFams->at(i).getGfamily1();
		int gfam2 = ECFams->at(i).getGfamily2();
		ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
		ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();

		if(Verbose)
		{
			cout << "refining equivalence class " << i << endl; 

		}

		ECFams->at(i).refine(Rtree1, Rtree2, useWholeClass, DeCoLTrefine, SuperVerbose); 


	}
	//cleaning

	Gfam1ToGfam2ToIndex.clear();
	

	return ECFams;

}


/*
=> SAME AS FUNCTION ABOVE BUT TAKES 2 MORE PARAMS (USEFUL for ARt-DeCo algorithm): 
	- map < string, int > speciesChrN, string: species name | int: chromosome number
	- map < int,vector <float> > speciesC0C1, int: species Id | vector <float>: value of the 6 different costs => c0/1(0,0), c0/1(1/0,0/1), c0/1(1,1)

Returns:
 - (vector < EquivalenceClassFamily > * ) : a pointer to a vector of refined Equivalence class families
*/
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies,
																 vector <GeneFamily *> * GeneFamilyList, 
																 MySpeciesTree *speciesTree, 
																 map < string, int > speciesChrNb, 
																 map < int,vector <float> > &speciesC0C1, 
																 map<int,string> &species_id_name, 
																 map<int, map<string,int> > &speGeneAdjNb, float Break, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose)
{

	// 1. creating all EquivalenceClassFamilies

	if(Verbose)
	{
		cout<<endl<<"### Creating Equivalence Class Families ###"<<endl;
		cout<<"1/ Creating all EquivalenceClassFamilies... "<<flush;
	}

	vector <EquivalenceClassFamily> * ECFams = new vector <EquivalenceClassFamily>;
	map < int, map <int , int > > Gfam1ToGfam2ToIndex; // to easily find the already created EquivalenceClasses

	for(unsigned i = 0; i < GeneFamilyList->size(); i++)
	{
		for(unsigned j = i ; j < GeneFamilyList->size() ; j++)  // "j = i" or "j = (i+1)"" ???
		{
			ECFams->push_back(EquivalenceClassFamily(i,j));
			Gfam1ToGfam2ToIndex[i][j] = ECFams->size() -1 ;

			//EquivalenceClasses.back().setAncestor( 0, GeneFamilyList->at(i)->getRecTree()->getRootId());
			//EquivalenceClasses.back().setAncestor( 1, GeneFamilyList->at(j)->getRecTree()->getRootId());
		}
	}


	//2. going through all leaves and adding them to the corresponding EquivalenceClassFamily

	if(Verbose)
	{
		cout<<"DONE"<<endl;
		cout<<"2/ Going through all leaves and adding them to the corresponding EquivalenceClassFamily... "<<flush;
	}
    // vector<int> spe_ID;
	for(unsigned i; i < adjacencies.size(); i++)
	{
		int spe1=-1;
		int spe2=-1;
		pair <string, string> currentadj = adjacencies[i];
		//Finding which Gfams are concerned.
		int gfam1 = -1;
		int gfam2 = -1;
		for(unsigned j = 0; j < GeneFamilyList->size(); j++)
		{
			GeneFamily * GF = GeneFamilyList->at(j);
			if(gfam1 ==-1)
				if(GF->hasLeaf(currentadj.first)){
					gfam1 = j;
					spe1=GF->getRecTree()->getNodeSpecies(GF->getRecTree()->getNode(currentadj.first)->getId());
					// vector<int>::iterator it;
					// it = find (spe_ID.begin(), spe_ID.end(), spe1);
					// if (it == spe_ID.end())
					// {
					// 	cout<<currentadj.first<<" -> "<<spe1<<endl;
					// 	spe_ID.push_back(spe1);
					// }
				}
			if(gfam2 == -1)
				if(GF->hasLeaf(currentadj.second)){
					gfam2 = j;
					spe2=GF->getRecTree()->getNodeSpecies(GF->getRecTree()->getNode(currentadj.second)->getId());
				}
			if((gfam1 !=-1) && (gfam2 != -1))//both gfam were found
				break;
		}
		if((gfam1 ==-1) || (gfam2 == -1))// at least one of the leaves was not found. Do nothing
		{
			if(Verbose)
				cout << "WARNING: The adjacency " << currentadj.first << " " << currentadj.second << " had a member which did not correspond to any gene family. Ignoring that adjacency." << endl;
			continue;
		}
		else
		{
			if(spe1!=spe2)
			{
				cout<<"ERROR: For adjacency: "<< currentadj.first << " " << currentadj.second << ": the genes are on two different species!!! (Respectively: "<< spe1 << " and " << spe2 << endl;
				exit(EXIT_FAILURE);
			}


			if(speGeneAdjNb[spe1].find(currentadj.first) != speGeneAdjNb[spe1].end())
				speGeneAdjNb[spe1][currentadj.first]+=1;
			else
				speGeneAdjNb[spe1][currentadj.first]=1;

			if(speGeneAdjNb[spe2].find(currentadj.second) != speGeneAdjNb[spe2].end())
				speGeneAdjNb[spe2][currentadj.second]+=1;
			else
				speGeneAdjNb[spe2][currentadj.second]=1;
		}





		//case were gfam1 and gfam2 were found
		if(SuperVerbose)
			cout << "Adj " << currentadj.first << " - " << currentadj.second << " families " << gfam1 << " " << gfam2 << endl;
		if(gfam1 > gfam2) // we will invert the values so that gfam1 is always the smallest id
		{
			int junk;
			string junkstr;
			junk = gfam1;
			gfam1 = gfam2;
			gfam2 = junk;
			junkstr = currentadj.first;
			currentadj.first = currentadj.second;
			currentadj.second = junkstr;
		}
		int ECindex	= Gfam1ToGfam2ToIndex[gfam1][gfam2];

		bool success = ECFams->at(ECindex).CheckAddAdj(currentadj.first,currentadj.second,gfam1,gfam2);//adding the adj to the correct Equivalence class
		if(SuperVerbose)
			if(!success)
				cout << "Could not add adj to the Equivalence Class " << ECindex << " " << ECFams->at(ECindex).getGfamily1() << " - " << ECFams->at(ECindex).getGfamily2() << endl;
	}

	//2 bis: delete all equivalence class family with 0 adjacencies if !useWholeClass"

	if(Verbose){
		cout<<"DONE"<<endl;
		cout<<"2bis/ Delete all equivalence class family with 0 adjacencies if !useWholeClass... "<<flush;
	}

	if(!useWholeClass)
	{
		for(int i = ECFams->size() -1 ; i >= 0; i--)
		{
			if(SuperVerbose)
				cout<<"\tProcessing equivalence class family "<< i <<endl;
			if(ECFams->at(i).getNbAdj() == 0 ) // no adjs
			{
 				ECFams->erase(ECFams->begin()+i); // deleting that ECFamily
			}
		}
	}


	//3. refining equivalence classes

	if(Verbose)
	{
		cout<<"DONE"<<endl;
		cout<<endl<<"Species name and ID association:"<<endl;
	}

	//WARNING: Use species tree that can be changed !!!
	vector<MySpeciesNode*> extantSpecies = speciesTree->getLeaves();
	for(int j =0; j < extantSpecies.size() ; j++)
	{
		string name=extantSpecies[j]->getName();
		int id=speciesTree->getRPO(extantSpecies[j]->getId());

		if(Verbose)
			cout<<"\t"<<name<< " -> " << id <<endl;

		species_id_name[id]=name;
	}
	extantSpecies.clear();

	if(Verbose)
	{
		cout<<"3/ Refining equivalence classes"<<endl;
	}
	
	for(unsigned i = 0; i <  ECFams->size(); i++)
	{
		int gfam1 = ECFams->at(i).getGfamily1();
		int gfam2 = ECFams->at(i).getGfamily2();
		ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
		ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();

        // MAYBE IT'S BETTER TO BROWSE ALL GENE TREES IN AN OTHER LOOP CAUSE CURRENTLY, SOME TREES WILL BE BROWSED SEVERAL TIMES...
		vector<int> leaves1 = Rtree1->getLeavesId();
		for(int j =0; j< leaves1.size() ; j++)
		{
			if(Rtree1->isRealLeaf(leaves1[j]))
			{
				// cout<<"Leaf "<<j<<endl;
				int spe1=Rtree1->getNodeSpecies(Rtree1->getNode(leaves1[j])->getId());
				map <string , int >::iterator itStrInt1;
				itStrInt1=speGeneAdjNb[spe1].find(Rtree1->getNodeName(leaves1[j]));
				if( itStrInt1 == speGeneAdjNb[spe1].end())
				{
					speGeneAdjNb[spe1][Rtree1->getNodeName(leaves1[j])]=0;
				}
			}
		}

		vector<int> leaves2 = Rtree2->getLeavesId();
		for(int j =0; j< leaves2.size() ; j++)
		{
			if(Rtree2->isRealLeaf(leaves2[j]))
			{
				// cout<<"Leaf "<<j<<endl;
				int spe2=Rtree2->getNodeSpecies(Rtree2->getNode(leaves2[j])->getId());
				map <string , int >::iterator itStrInt2;
				itStrInt2=speGeneAdjNb[spe2].find(Rtree2->getNodeName(leaves2[j]));
				if( itStrInt2 == speGeneAdjNb[spe2].end())
				{
					speGeneAdjNb[spe2][Rtree2->getNodeName(leaves2[j])]=0;
				}
			}
		}

		if(SuperVerbose)
			cout << "\tRefining equivalence class " << i << endl;

		ECFams->at(i).refine(Rtree1, Rtree2, useWholeClass, DeCoLTrefine, SuperVerbose);
	}


	if(Verbose){
		cout<<"4/ Compute C0/C1 for case Ext/Ext with  stats for all species: "<<endl;
	}



    //4. Compute C0/C1 for case Ext/Ext with  stats for all species

	map < int, map <string , int > >::iterator it1;
	for(it1=speGeneAdjNb.begin();it1!=speGeneAdjNb.end();it1++)
	{
		int spe=(*it1).first;
		int nbGene0Adj=0;
		int nbGene1Adj=0;
		int nbGene2Adj=0;
		int nbGene3Adj=0;
		map <string , int >::iterator it2 = (*it1).second.begin();
		while(it2 != (*it1).second.end())
		{
			if((*it2).second == 0)
			{
				nbGene0Adj+=1;
				++it2;
			}
			else if((*it2).second == 1)
			{
				nbGene1Adj+=1;
				++it2;
			}
			else if((*it2).second == 2)
			{
				nbGene2Adj+=1;
				speGeneAdjNb[spe].erase(it2++);
			}
			else
			{
				nbGene3Adj+=1;
				speGeneAdjNb[spe].erase(it2++);
			}
		}

		string species=species_id_name[spe];
		if(Verbose)
			cout<<"\tProcessing species "<<species<<"("<<spe<<"):"<<endl;

        if(nbGene3Adj>0)
        	cout<<"Genome assembly of species "<<species<<" is not linear or circular then computation of P(v1~v2) won't be exact (For more information see article on ARt-DeCo algorithm)"<<endl;

		int chr=speciesChrNb[species];
		int ctg=nbGene0Adj+(nbGene1Adj/2);

		if (Verbose)
			cout<<"\t\ta/ Compute P(v1~v2)"<<endl;

		//Compute P(v1~v2) for rho(v1~v2)=1
		float Proba=float(ctg-chr)/float(2*ctg*(ctg-1));

		if (Verbose)
			cout<<"\t\tb/ Compute base_log"<<endl;

		//Compute base log to compute c0P(v1~v2) & c1P(v1~v2)      
		float base_log;
		// Case where genome assembly is complete and chr==contig
		if (Proba==0)
			base_log=2.0;
		else
			base_log=ceil(pow((1.0-Proba)/Proba,1.0/Break))*1.05;	// Computation of base_log

		if(Verbose)
			cout<<"\t\t\tbase_log = "<<base_log<<endl;

		if (Verbose)
			cout<<"\t\tc/ Compute C0 and C1 for the 3 different cases (0 vs 0, 0/1 vs 1/0 and 1 vs 1)"<<endl;


		float c0P_11=-log(1.0-Proba)/log(base_log);
		if(Verbose)
			cout<<"\t\t\tc0P_11 = "<<c0P_11<<endl;
		float c0P_01=-log(1.0-(2*Proba))/log(base_log);
		if(Verbose)
			cout<<"\t\t\tc0P_01 = "<<c0P_01<<endl;
		float c0P_00=-log(1.0-(4*Proba))/log(base_log);
		if(Verbose)
			cout<<"\t\t\tc0P_00 = "<<c0P_00<<endl;
		float c1P_11=-log(Proba)/log(base_log);
		if(Verbose)
			cout<<"\t\t\tc1P_11 = "<<c1P_11<<endl;
		float c1P_01=-log(2*Proba)/log(base_log);
		if(Verbose)
			cout<<"\t\t\tc1P_01 = "<<c1P_01<<endl;
		float c1P_00=-log(4*Proba)/log(base_log);
		if(Verbose)
			cout<<"\t\t\tc1P_00 = "<<c1P_00<<endl;

		speciesC0C1[spe].push_back(c0P_11);
		speciesC0C1[spe].push_back(c0P_01);
		speciesC0C1[spe].push_back(c0P_00);
		speciesC0C1[spe].push_back(c1P_11);
		speciesC0C1[spe].push_back(c1P_01);
		speciesC0C1[spe].push_back(c1P_00);

		if(Verbose)
		{
			cout<<"\t\t=> Species "<<species<<" has "<<nbGene0Adj+nbGene1Adj+nbGene2Adj+nbGene3Adj<<" genes:"<<endl;
			cout<<"\t\t\t"<<nbGene0Adj<<" genes with 0 adjacency"<<endl;
			cout<<"\t\t\t"<<nbGene1Adj<<" genes with 1 adjacency"<<endl;
			cout<<"\t\t\t"<<nbGene2Adj<<" genes with 2 adjacencies"<<endl;
			cout<<"\t\t\t"<<nbGene3Adj<<" genes with more than 2 adjacencies"<<endl;
		}
	}

	if(Verbose){
		cout<<"DONE"<<endl;
	}

	//cleaning
	Gfam1ToGfam2ToIndex.clear();

	return ECFams;

}


/*
Takes:
	- adjacencies (vector< pair <string,string > >)
	- map <string,int> LeafToGFMap : leaf names to gene family ids
 	- int nbGfam : number of gene families
	- DeCoLTrefine (bool) : if true, DeCoLT refining method will used even if there are no HGT in the trees
	- useWholeClass (bool): if true, all possible adjacencies (including non observed ones), will be used for refinement
	- Verbose (bool)
	- SuperVerbose (bool)

Returns:
 - (vector < EquivalenceClassFamily > * ) : a pointer to a vector of refined Equivalence class families
*/
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > &adjacencies, 
																map <string,int> &LeafToGFMap, 
																int nbGfam, 
																bool DeCoLTrefine, 
																bool useWholeClass ,bool Verbose, bool SuperVerbose)
{



	// 1. creating  EquivalenceClassFamilies containiners

	vector <EquivalenceClassFamily> * ECFams = new vector <EquivalenceClassFamily>;
	map < int, map <int , int > > Gfam1ToGfam2ToIndex; // to easily find the already created EquivalenceClasses

	//cout << "treating "<< adjacencies.size() << "adjacencies."<<endl;

	//2. going through all leaves an adding them to the corresponding EquivalenceClassFamily

	map < int, map <int , int > >::iterator it1;
	map <int , int >::iterator it2;

	for(unsigned i= 0; i < adjacencies.size(); i++)
	{
		pair <string, string> currentadj = adjacencies[i];
		//Finding which Gfams are concerned.
		int gfam1 = -1;
		int gfam2 = -1;

		map<string,int>::iterator it;

		it = LeafToGFMap.find(currentadj.first);
		if(it != LeafToGFMap.end())
			gfam1 = it->second;

		it = LeafToGFMap.find(currentadj.second);
		if(it != LeafToGFMap.end())
			gfam2 = it->second;


		if((gfam1 ==-1) || (gfam2 == -1))// at least one of the leaves was not found. Do nothing
		{
			if(Verbose)
				cout << "The adjacency " << currentadj.first << " " << currentadj.second << " had a member which did not correspond to any gene family. Ignoring that adjacency." << endl;
			continue;
		}
		//case were gfam1 and gfam2 were found
		if(SuperVerbose)
			cout << "Adj "<< i << " : " << currentadj.first << " - " << currentadj.second << " families " << gfam1 << " " << gfam2 << endl;

		if(gfam1 > gfam2) // we will invert the values so that gfam1 is always the smallest id
		{
			int junk;
			string junkstr;
			junk = gfam1;
			gfam1 = gfam2;
			gfam2 = junk;
			junkstr = currentadj.first;
			currentadj.first = currentadj.second;
			currentadj.second = junkstr;
		}

		bool toCreate = false;
		if(Gfam1ToGfam2ToIndex.find(gfam1) == Gfam1ToGfam2ToIndex.end())
		{
			toCreate = true;
		}
		else if(Gfam1ToGfam2ToIndex[gfam1].find(gfam2) == Gfam1ToGfam2ToIndex[gfam1].end())
		{
			toCreate = true;
		}

		int ECindex;
		if(toCreate)
		{
			ECFams->push_back(EquivalenceClassFamily(gfam1,gfam2));
			Gfam1ToGfam2ToIndex[gfam1][gfam2] = ECFams->size() -1 ;
			ECindex = ECFams->size() -1;
		}
		else
		{
			ECindex	= Gfam1ToGfam2ToIndex[gfam1][gfam2];
		}

		bool success = ECFams->at(ECindex).CheckAddAdj(currentadj.first,currentadj.second,gfam1,gfam2);//adding the adj to the correct Equivalence class
		if(SuperVerbose)
			if(!success)
				cout << "Could not add adj to the Equivalence Class " << ECindex << " " << ECFams->at(ECindex).getGfamily1() << " - " << ECFams->at(ECindex).getGfamily2() << endl;
	}
	if(SuperVerbose)
		cout << "2bis"<< endl;

	//2 bis: delete all equivalence class family with 0 adjacencies if !useWholeClass"
//	if(!useWholeClass)
//	{
//		for(int i = ECFams->size() -1 ; i >= 0; i--)
//		{
//
//			if(ECFams->at(i).getNbAdj() == 0 ) // no adjs
//			{
//				if(SuperVerbose)
//					cout << "clearing ECF " << i << endl;
//
// 				ECFams->erase(ECFams->begin()+i); // deleting that ECFamily
//			}
//
//		}
//	}
	if(useWholeClass)
	{
		for(unsigned gfam1 = 0; gfam1 < nbGfam; gfam1++)
		{
			for(unsigned gfam2 = gfam1 ; gfam2 < nbGfam ; gfam2++)
			{
				bool toCreate = false;
				if(Gfam1ToGfam2ToIndex.find(gfam1) == Gfam1ToGfam2ToIndex.end())
				{
					toCreate = true;
				}
				else if(Gfam1ToGfam2ToIndex[gfam1].find(gfam2) == Gfam1ToGfam2ToIndex[gfam1].end())
				{
					toCreate = true;
				}
		
				
				if(toCreate)
				{
					ECFams->push_back(EquivalenceClassFamily(gfam1,gfam2));
					Gfam1ToGfam2ToIndex[gfam1][gfam2] = ECFams->size() -1 ;
				}

			}
		}

	}

	//cleaning

	Gfam1ToGfam2ToIndex.clear();
	

	return ECFams;

}

/*
Takes:
	- adjacencies (vector< pair <string,string > >)
	- vector < pair <int,int> > &adjacenciesOrientations : the int indicated which extremity of the gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
	- map <string,int> LeafToGFMap : leaf names to gene family ids
 	- int nbGfam : number of gene families
	- DeCoLTrefine (bool) : if true, DeCoLT refining method will used even if there are no HGT in the trees
	- useWholeClass (bool): if true, all possible adjacencies (including non observed ones), will be used for refinement
	- Verbose (bool)
	- SuperVerbose (bool)

Returns:
 - (vector < EquivalenceClassFamily > * ) : a pointer to a vector of refined Equivalence class families
*/
vector < EquivalenceClassFamily > * CreateEquivalenceClassFamilies(vector< pair <string,string > > &adjacencies, 
																vector < pair <int,int> >  &adjacenciesOrientations, map <string,int> &LeafToGFMap, int nbGfam, bool DeCoLTrefine, bool useWholeClass ,bool Verbose, bool SuperVerbose)
{
	// 1. creating  EquivalenceClassFamilies containiners

	vector <EquivalenceClassFamily> * ECFams = new vector <EquivalenceClassFamily>;

	//cout << "treating "<< adjacencies.size() << "adjacencies."<<endl;

	//2. going through all leaves an adding them to the corresponding EquivalenceClassFamily

	map < int, map <int , int > >::iterator it1;
	map <int , int >::iterator it2;

	for(unsigned i= 0; i < adjacencies.size(); i++)
	{
		pair <string, string> currentadj = adjacencies[i];
		//Finding which Gfams are concerned.
		int gfam1 = -1;
		int gfam2 = -1;

		int s1 = adjacenciesOrientations[i].first;
		int s2 = adjacenciesOrientations[i].second;

		map<string,int>::iterator it;

		it = LeafToGFMap.find(currentadj.first);
		if(it != LeafToGFMap.end())
			gfam1 = it->second;

		it = LeafToGFMap.find(currentadj.second);
		if(it != LeafToGFMap.end())
			gfam2 = it->second;


		if((gfam1 ==-1) || (gfam2 == -1))// at least one of the leaves was not found. Do nothing
		{
			if(Verbose)
				cout << "The adjacency " << currentadj.first << " " << currentadj.second << " had a member which did not correspond to any gene family. Ignoring that adjacency." << endl;
			continue;
		}
		//case were gfam1 and gfam2 were found
		if(SuperVerbose)
			cout << "Adj "<< i << " : " << currentadj.first << " - " << currentadj.second << " families " << gfam1 << " " << gfam2 << " | "<< s1 << " "<< s2 << endl;

		if(gfam1 > gfam2) // we will invert the values so that gfam1 is always the smallest id
		{
			int junk;
			string junkstr;
			junk = gfam1;
			gfam1 = gfam2;
			gfam2 = junk;

			junkstr = currentadj.first;
			currentadj.first = currentadj.second;
			currentadj.second = junkstr;

			junk = s1;
			s1 = s2;
			s2 = junk;

		}

		int ECindex = -1;
		for(unsigned j = 0 ; j < ECFams->size(); j++)
		{
			if( ECFams->at(j).isCompatible(gfam1, gfam2, s1, s2  ) )// found the correct
			{
				//cout << "compat "<< gfam1<<" " <<gfam2<<" "<< s1<<" "<< s2<< " <-> "<<  ECFams->at(j).getGfamily1()<<" " <<ECFams->at(j).getGfamily2()<<" "<< ECFams->at(j).getSens1()<<" "<< ECFams->at(j).getSens2() << endl;
				ECindex = j;
				break;
			}
		}


		if(ECindex == -1)
		{
			ECFams->push_back(EquivalenceClassFamily(gfam1,gfam2));
			ECFams->back().setSens1(s1);
			ECFams->back().setSens2(s2);

			ECindex = ECFams->size() -1;
		}


		bool success = ECFams->at(ECindex).CheckAddAdj(currentadj.first,currentadj.second,gfam1,gfam2, s1,s2);//adding the adj to the correct Equivalence class
		if(SuperVerbose)
			if(!success)
				cout << "Could not add adj to the Equivalence Class " << ECindex << " " << ECFams->at(ECindex).getGfamily1() << " - " << ECFams->at(ECindex).getGfamily2() << endl;
	}
	if(SuperVerbose)
		cout << "2bis"<< endl;


	if(useWholeClass)
	{
		for(unsigned gfam1 = 0; gfam1 < nbGfam; gfam1++)
		{
			for(unsigned gfam2 = gfam1 ; gfam2 < nbGfam ; gfam2++)
			{
				bool toCreate = true;
				for(unsigned j = 0 ; j < ECFams->size(); j++)
				{
					if(  (  ECFams->at(j).getGfamily1() == gfam1 ) && (ECFams->at(j).getGfamily2() == gfam2) )// found the correct
					{
						toCreate = false;
						break;
					}
				}
				
				if(toCreate)
				{
					ECFams->push_back(EquivalenceClassFamily(gfam1,gfam2));
					ECFams->back().setSens1(0);
					ECFams->back().setSens2(0);
					
				}

			}
		}

	}




	return ECFams;


}


/*
Fills  speciesC0C1 and speGeneAdjNb

Takes:
 	- adjacencies (vector< pair <string,string > >)
	- map <string,int> &LeafToSpMap
	- MySpeciesTree *speciesTree	
	- map < string, int > speciesChrNb, string: species name | int: chromosome number
	- map < int,vector <float> > speciesC0C1, int: species Id | vector <float>: value of the 6 different costs => c0/1(0,0), c0/1(1/0,0/1), c0/1(1,1)
 	- map<int,string> &species_id_name, 
	- map<int, map<string,int> > &speGeneAdjNb
	- float Break
 	- Verbose (bool)
	- SuperVerbose (bool)
	- map < string, map <string , double> > * adjacencyScores : associates somes adjacencies with a score between 0 and 1 
	- bool includeScoredAdjs : wether the adjacencies with a score != 1 should be included or not when computing the number of contigs

Returns:
	(bool): false if no problem occured; true if there a pb (less contigs than expected chromosomes.)
*/
bool computeArtDeCoMaps(vector< pair <string,string > > &adjacencies, map <string,int> &LeafToSpMap, MySpeciesTree *speciesTree,
						 map < string, int > speciesChrNb, map < int,vector <float> > &speciesC0C1, map<int,string> &species_id_name, 
						 map<int, map<string,int> > &speGeneAdjNb, float Break,bool Verbose, bool SuperVerbose, map < string, map <string , double> > * adjacencyScores, bool includeScoredAdjs)
{

	//1. initialize speGeneAdjNb by putting all genes at 0
	for( map<string,int>::iterator it = LeafToSpMap.begin() ; it != LeafToSpMap.end(); ++it) // iterate over all genes
	{
		speGeneAdjNb[ it->second ][ it->first ] = 0;
	}

	bool returnValue = false;

	//2. add the adjacencies to the mix

	for(unsigned i = 0; i < adjacencies.size(); i++)
	{
		pair <string, string> currentadj = adjacencies[i];

		//finding wether or not the adj have a score
		map < string, map <string , double> >::iterator it1;
		map <string , double>::iterator it2;
		it1 = adjacencyScores->find(currentadj.first) ;
		if(it1 != adjacencyScores->end())
		{
			it2 = it1->second.find(currentadj.second);
			if( it2 != it1->second.end() )
			{
				if( it2->second == 0)
					continue; // always ignore an adj with score 0 here
				else if(( !includeScoredAdjs )&&( it2->second < 1 ))
					continue; // we ignore the score adjs for this step
			}
		}

		//Finding which Gfams are concerned.
		int sp1 = -1;
		int sp2 = -1;

		map<string,int>::iterator it;

		it = LeafToSpMap.find(currentadj.first);
		if(it != LeafToSpMap.end())
			sp1 = it->second;

		it = LeafToSpMap.find(currentadj.second);
		if(it != LeafToSpMap.end())
			sp2 = it->second;


		if((sp1 ==-1) || (sp2 == -1))// at least one of the leaves was not found. Do nothing
		{
			if(Verbose)
				cout << "The adjacency " << currentadj.first << " " << currentadj.second << " had a member which did not correspond to any species. Ignoring that adjacency." << endl;
			continue;
		}

		//past this point sp1 and sp2 have been found
		speGeneAdjNb[ sp1 ][ currentadj.first ] += 1;
		speGeneAdjNb[ sp2 ][ currentadj.second ] += 1;

	}


    //3. Compute C0/C1 for case Ext/Ext with  stats for all species

	map < int, map <string , int > >::iterator it1;
	for(it1=speGeneAdjNb.begin();it1!=speGeneAdjNb.end();it1++)
	{
		int spe=(*it1).first;
		int nbGene0Adj=0;
		int nbGene1Adj=0;
		int nbGene2Adj=0;
		int nbGene3Adj=0;
		map <string , int >::iterator it2 = (*it1).second.begin();
		while(it2 != (*it1).second.end())
		{
			if((*it2).second == 0)
			{
				nbGene0Adj+=1;
				++it2;
			}
			else if((*it2).second == 1)
			{
				nbGene1Adj+=1;
				++it2;
			}
			else if((*it2).second == 2)
			{
				nbGene2Adj+=1;
				speGeneAdjNb[spe].erase(it2++);
			}
			else
			{
				nbGene3Adj+=1;
				speGeneAdjNb[spe].erase(it2++);
			}
		}

		string species=species_id_name[spe];
		if(Verbose)
			cout<<"\tProcessing species "<<species<<"("<<spe<<"):"<<endl;

        //if(nbGene3Adj>0)
        //	cerr<<"Scaffolding mode : genome assembly of species "<<species<<" is not linear or circular then computation of P(v1~v2) won't be exact (For more information see article on ARt-DeCo algorithm)"<<endl;

		int chr=speciesChrNb[species];
		int ctg=nbGene0Adj+(nbGene1Adj/2);

		if (Verbose)
			cout<<"\t\ta/ Compute P(v1~v2)"<<endl;


		//Compute P(v1~v2) for rho(v1~v2)=1
		float c0P_11;
		float c0P_01;
		float c0P_00;
		float c1P_11;
		float c1P_01;
		float c1P_00;

		if( ctg < chr ) // less contigs than expected chromosomes -> error!
		{
			returnValue = true;

			c0P_11 = 0;
			c0P_01 = 0;
			c0P_00 = 0;
			c1P_11 = numeric_limits<float>::max();
			c1P_01 = numeric_limits<float>::max();
			c1P_00 = numeric_limits<float>::max();

		}
		else if(ctg == 1)
		{ // fail safe if only one contig

			c0P_11 = 0;
			c0P_01 = 0;
			c0P_00 = 0;
			c1P_11 = numeric_limits<float>::max();
			c1P_01 = numeric_limits<float>::max();
			c1P_00 = numeric_limits<float>::max();


		}
		else
		{

			float Proba=float(ctg-chr)/float(2*ctg*(ctg-1));
	
			if (Verbose)
				cout<<"\t\tb/ Compute base_log"<<endl;
	
			//Compute base log to compute c0P(v1~v2) & c1P(v1~v2)      
			float base_log;
			// Case where genome assembly is complete and chr==contig
			if (Proba==0)
				base_log=2.0;
			else
				base_log=ceil(pow((1.0-Proba)/Proba,1.0/Break))*1.05;	// Computation of base_log
	
			if(Verbose)
				cout<<"\t\t\tbase_log = "<<base_log<<endl;
	
			if (Verbose)
				cout<<"\t\tc/ Compute C0 and C1 for the 3 different cases (0 vs 0, 0/1 vs 1/0 and 1 vs 1)"<<endl;
	
	
			c0P_11=-log(1.0-Proba)/log(base_log);
			if(Verbose)
				cout<<"\t\t\tc0P_11 = "<<c0P_11<<endl;
			c0P_01=-log(1.0-(2*Proba))/log(base_log);
			if(Verbose)
				cout<<"\t\t\tc0P_01 = "<<c0P_01<<endl;
			c0P_00=-log(1.0-(4*Proba))/log(base_log);
			if(Verbose)
				cout<<"\t\t\tc0P_00 = "<<c0P_00<<endl;
			c1P_11=-log(Proba)/log(base_log);
			if(Verbose)
				cout<<"\t\t\tc1P_11 = "<<c1P_11<<endl;
			c1P_01=-log(2*Proba)/log(base_log);
			if(Verbose)
				cout<<"\t\t\tc1P_01 = "<<c1P_01<<endl;
			c1P_00=-log(4*Proba)/log(base_log);
			if(Verbose)
				cout<<"\t\t\tc1P_00 = "<<c1P_00<<endl;
	
		}

		speciesC0C1[spe].push_back(c0P_11);
		speciesC0C1[spe].push_back(c0P_01);
		speciesC0C1[spe].push_back(c0P_00);
		speciesC0C1[spe].push_back(c1P_11);
		speciesC0C1[spe].push_back(c1P_01);
		speciesC0C1[spe].push_back(c1P_00);

		if(Verbose)
		{
			cout<<"\t\t=> Species "<<species<<" has "<<nbGene0Adj+nbGene1Adj+nbGene2Adj+nbGene3Adj<<" genes:"<<endl;
			cout<<"\t\t\t"<<nbGene0Adj<<" genes with 0 adjacency"<<endl;
			cout<<"\t\t\t"<<nbGene1Adj<<" genes with 1 adjacency"<<endl;
			cout<<"\t\t\t"<<nbGene2Adj<<" genes with 2 adjacencies"<<endl;
			cout<<"\t\t\t"<<nbGene3Adj<<" genes with more than 2 adjacencies"<<endl;
		}
	}

	return returnValue;
}

/*
Computes IN PLACE the EquivalenceClassFamilies's AdjMatrix

Takes:
 - vector < EquivalenceClassFamily> * ECFams : vector of equivalence class families containing refined equivalence class
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - map < string, map <string , double> > &adjacencyScores : map containing probas associated with adjacencies for ADseq
 - map<int,vector<float> > speciesC0C1 : for artdeco
 - map<int, map<string,int> > speGeneAdjNb : for artdeco
 - double Gcost : cost of an adjacency Gain
 - double Bcost  : cost of an adjacency Break
 - bool boltzmann : wether to use boltzmann computation or not
 - bool SubRecInAdj : if set to true, the weighted cost of a reconciliation event will be used to favor co-event in the adjacency matrix computation. Unavalaible for Boltzmann computation.
 - double WDupCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WLossCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WHgtCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - bool Verbose
 - bool SuperVerbose
 - double boltzmannTemp [default : 1] : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - double absencePenalty [default : -1] : if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization.
 - double adjScoreLogBase [default : 10000] : base of the log that will be used to go from adjacency confidence score to parsimony costs

*/
void ComputeEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList , map < string, map <string , double> > &adjacencyScores, map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,double WDupCost, double WLossCost, double WHgtCost, bool Verbose, bool SuperVerbose, double boltzmannTemp , double absencePenalty, double adjScoreLogBase , bool interactionMode)
{
	cout << setprecision(2);
	for(unsigned i = 0 ; i < ECFams->size(); i++)
	{
		if(Verbose)
			cout << "Adjacency matrix computation for the Equivalence class " << i << endl;
		EquivalenceClassFamily * ECF = &ECFams->at(i);

		ComputeOneEquivalenceClassFamily( ECF,  GeneFamilyList, 
										adjacencyScores,  speciesC0C1,  speGeneAdjNb,  
										Gcost,  Bcost ,  
										boltzmann,  SubRecInAdj, 
										WDupCost,  WLossCost, WHgtCost, 
										Verbose, SuperVerbose, boltzmannTemp , absencePenalty, adjScoreLogBase , interactionMode);

		if(Verbose)
			cout << "Matrix computed" << endl;
	}
	cout << setprecision(5);
}


/*
Computes IN PLACE the EquivalenceClassFamily's AdjMatrix

Takes:
 - EquivalenceClassFamily * ECF :  equivalence class family
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - map < string, map <string , double> > &adjacencyScores : map containing probas associated with adjacencies for ADseq
 - map<int,vector<float> > speciesC0C1 : for artdeco
 - map<int, map<string,int> > speGeneAdjNb : for artdeco
 - double Gcost : cost of an adjacency Gain
 - double Bcost  : cost of an adjacency Break
 - bool boltzmann : wether to use boltzmann computation or not
 - bool SubRecInAdj : if set to true, the weighted cost of a reconciliation event will be used to favor co-event in the adjacency matrix computation. Unavalaible for Boltzmann computation.
 - double WDupCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WLossCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WHgtCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - bool Verbose
 - bool SuperVerbose
 - double boltzmannTemp [default : 1] : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - double absencePenalty [default : -1] : if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization.
 - double adjScoreLogBase [default : 10000] : base of the log that will be used to go from adjacency confidence score to parsimony costs

Returns:
	vector < pair < pair<string, string> , double > >
*/
vector < pair < pair<string, string> , double > > ComputeOneEquivalenceClassFamily( EquivalenceClassFamily * ECF, vector <GeneFamily *> * GeneFamilyList, map < string, map <string , double> > &adjacencyScores, map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,double WDupCost, double WLossCost, double WHgtCost, bool Verbose, bool SuperVerbose, double boltzmannTemp , double absencePenalty, double adjScoreLogBase, bool interactionMode )
{

	ReconciledTree * Rtree1 =  GeneFamilyList->at(ECF->getGfamily1())->getRecTree();
	ReconciledTree * Rtree2 =  GeneFamilyList->at(ECF->getGfamily2())->getRecTree();

	return ComputeOneEquivalenceClassFamily( ECF,  Rtree1, Rtree2, adjacencyScores, speciesC0C1, speGeneAdjNb,  Gcost, Bcost , boltzmann, SubRecInAdj, WDupCost, WLossCost, WHgtCost, Verbose, SuperVerbose, boltzmannTemp , absencePenalty , interactionMode);
	
}

/*
Computes IN PLACE the EquivalenceClassFamily's AdjMatrix

Takes:
 - EquivalenceClassFamily * ECF :  equivalence class family
 - ReconciledTree * Rtree1 : pointer to the first reconciledtree
 - ReconciledTree * Rtree2 : pointer to the second reconciled tree
 - map < string, map <string , double> > &adjacencyScores : map containing probas associated with adjacencies for ADseq
 - map<int,vector<float> > speciesC0C1 : for artdeco
 - map<int, map<string,int> > speGeneAdjNb : for artdeco
 - double Gcost : cost of an adjacency Gain
 - double Bcost  : cost of an adjacency Break
 - bool boltzmann : wether to use boltzmann computation or not
 - bool SubRecInAdj : if set to true, the weighted cost of a reconciliation event will be used to favor co-event in the adjacency matrix computation. Unavalaible for Boltzmann computation.
 - double WDupCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WLossCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WHgtCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - bool Verbose
 - bool SuperVerbose
 - double boltzmannTemp [default : 1] : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - double absencePenalty [default : -1] : if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization.
 - double adjScoreLogBase [default : 10000] : base of the log that will be used to go from adjacency confidence score to parsimony costs

Returns:
	vector < pair < pair<string, string> , double > >
*/
vector < pair < pair<string, string> , double > > ComputeOneEquivalenceClassFamily( EquivalenceClassFamily * ECF, ReconciledTree * Rtree1, ReconciledTree * Rtree2, map < string, map <string , double> > &adjacencyScores, map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,double WDupCost, double WLossCost, double WHgtCost, bool Verbose, bool SuperVerbose, double boltzmannTemp , double absencePenalty, double adjScoreLogBase, bool interactionMode)
{
	cout << setprecision(2);

	vector < pair < pair<string, string> , double > > scoreAssociation = ECF->createAdjMatrix(adjacencyScores, speciesC0C1, speGeneAdjNb, Gcost, Bcost,  Rtree1,  Rtree2, SuperVerbose , boltzmann , boltzmannTemp , absencePenalty, adjScoreLogBase, interactionMode);
	if(!SubRecInAdj)
		ECF->computeAdjMatrix();
	else //using weighted reconciliation event costs
		ECF->computeAdjMatrix(WDupCost, WLossCost, WHgtCost);

	
	cout << setprecision(5);
	return scoreAssociation ;
}



/*
backtracks IN PLACE the EquivalenceClassFamilies's AdjMatrix

Takes:
 - vector < EquivalenceClassFamily> * ECFams : vector of equivalence class families containing refined equivalence class
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - bool boltzmann : wether to use boltzmann computation or not
 - bool Verbose
 - bool SuperVerbose
 - bool galwaysGain : there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 - double gC1Advantage : probability to choose c1 over c0 IF (and only if) they have the same score
 - bool useCounts [default = true]: whether to use count of possible backtracks or not. 

*/
void backtrackInPlaceEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double gC1Advantage, bool useCounts)
{

	for(unsigned i = 0 ; i < ECFams->size(); i++)
	{

		EquivalenceClassFamily * ECF = &ECFams->at(i);
		ReconciledTree * Rtree1 =  GeneFamilyList->at(ECF->getGfamily1())->getRecTree();
		ReconciledTree * Rtree2 =  GeneFamilyList->at(ECF->getGfamily2())->getRecTree();
		
		ECF->backtrackAdjMatrixForSelf(Rtree1, Rtree2, boltzmann, galwaysGain, gC1Advantage, useCounts);
		
		if(Verbose)
			cout << "Backtrack finished: " << ECF->getNbAdjTrees() << " trees." << endl;
		if(Verbose)
			cout << " Total of "<< ECF->getNbAdjGain() << " Adjacency Gains and " << ECF->getNbAdjBreak() << " Adjacency Breaks." << endl;
	}
}



/*
backtracks IN PLACE the EquivalenceClassFamilies's AdjMatrix

Takes:
 - vector < EquivalenceClassFamily> * ECFams : vector of equivalence class families containing refined equivalence class
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - bool boltzmann : wether to use boltzmann computation or not
 - bool Verbose
 - bool SuperVerbose
 - bool galwaysGain : there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 - double gC1Advantage : probability to choose c1 over c0 IF (and only if) they have the same score

*/
void backtrackInPlaceOneEquivalenceClassFamily(EquivalenceClassFamily * ECF, vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage)
{

	ReconciledTree * Rtree1 =  GeneFamilyList->at(ECF->getGfamily1())->getRecTree();
	ReconciledTree * Rtree2 =  GeneFamilyList->at(ECF->getGfamily2())->getRecTree();
	
	backtrackInPlaceOneEquivalenceClassFamily( ECF, Rtree1, Rtree2 , boltzmann, Verbose, SuperVerbose, galwaysGain,  gC1Advantage);

}

/*
backtracks IN PLACE the EquivalenceClassFamilies's AdjMatrix

Takes:
 - vector < EquivalenceClassFamily> * ECFams : vector of equivalence class families containing refined equivalence class
 - ReconciledTree * Rtree1 : pointer to the first reconciledtree
 - ReconciledTree * Rtree2 : pointer to the second reconciled tree
 - bool boltzmann : wether to use boltzmann computation or not
 - bool Verbose
 - bool SuperVerbose
 - bool galwaysGain : there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 - double gC1Advantage : probability to choose c1 over c0 IF (and only if) they have the same score

*/
void backtrackInPlaceOneEquivalenceClassFamily(EquivalenceClassFamily * ECF, ReconciledTree * Rtree1, ReconciledTree * Rtree2 , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage)
{
	
	ECF->backtrackAdjMatrixForSelf(Rtree1, Rtree2, boltzmann, galwaysGain, gC1Advantage);
	
	if(Verbose)
		cout << "Backtrack finished: " << ECF->getNbAdjTrees() << " trees." << endl;
	if(Verbose)
		cout << " Total of "<< ECF->getNbAdjGain() << " Adjacency Gains and " << ECF->getNbAdjBreak() << " Adjacency Breaks." << endl;

}


/*
backtracks the EquivalenceClassFamily's AdjMatrix 1 times and add it to Sample

Takes:
 - ECFsample * Sample : pointer to an ECFSample (ECFsample is the result of one pass of backtracking on and Equivalence Class Family)
 - EquivalenceClassFamily * ECF : the equivalence class family to backtrack n times
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - bool boltzmann : wether to use boltzmann computation or not
 - bool Verbose
 - bool SuperVerbose
 - bool galwaysGain : there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 - double gC1Advantage : probability to choose c1 over c0 IF (and only if) they have the same score
 - int &overflowed : nuùmber of overflowed EC

*/
void backtrackOnetimeOneEquivalenceClassFamily( ECFsample * Sample, EquivalenceClassFamily * ECF , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage, int &overflowed)
{
	ReconciledTree * Rtree1 =  GeneFamilyList->at(ECF->getGfamily1())->getRecTree();
	ReconciledTree * Rtree2 =  GeneFamilyList->at(ECF->getGfamily2())->getRecTree();


	ECF->backtrackAdjMatrix(Rtree1, Rtree2, Sample , boltzmann, overflowed, galwaysGain, gC1Advantage); // creates the adj forest 
}

/*
backtracks the EquivalenceClassFamily's AdjMatrix N times and add it to Sample

Takes:
 - NECFsample * Samples : pointer to a vector of ECFSample (ECFsample is the result of one pass of backtracking on and Equivalence Class Family)
 - EquivalenceClassFamily * ECF : the equivalence class family to backtrack n times
 - int N : number of time to backtrack
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - bool boltzmann : wether to use boltzmann computation or not
 - bool Verbose
 - bool SuperVerbose
 - bool galwaysGain : there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 - double gC1Advantage : probability to choose c1 over c0 IF (and only if) they have the same score

*/
void backtrackNtimesOneEquivalenceClassFamily(NECFsample * Samples, EquivalenceClassFamily * ECF, int N , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double  gC1Advantage)
{
	//cout << "backtrackNtimesOneEquivalenceClassFamily init" << endl;
	ReconciledTree * Rtree1 =  GeneFamilyList->at(ECF->getGfamily1())->getRecTree();
	ReconciledTree * Rtree2 =  GeneFamilyList->at(ECF->getGfamily2())->getRecTree();
	
	int overflowed = 0;

	for(unsigned i = 0 ; i < N ; i++)
	{
		//cout << "backtrackNtimesOneEquivalenceClassFamily " << i << endl;

		Samples->push_back( new ECFsample);
		ECF->backtrackAdjMatrix(Rtree1, Rtree2, Samples->back() , boltzmann, overflowed, galwaysGain, gC1Advantage); // creates the adj forest 

		if(Verbose)
			cout << "backtrack : " << i << "/" << N << "\r"; // trying something...
	}

}


/*
backtracks the EquivalenceClassFamilies's AdjMatrix N times and add it to AllSamples

Takes:
 - vector <NECFsample * > * AllSamples : pointer to a vector of pointers to vectors of ECFSample (ECFsample is the result of one pass of backtracking on and Equivalence Class Family)
 - vector < EquivalenceClassFamily> * ECFams : vector of equivalence class families containing refined equivalence class
 - int N : number of time to backtrack
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - bool boltzmann : wether to use boltzmann computation or not
 - bool Verbose
 - bool SuperVerbose
 - bool galwaysGain : there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 - double gC1Advantage : probability to choose c1 over c0 IF (and only if) they have the same score

*/
void backtrackNtimesNequivalenceClassFamilies(vector <NECFsample * > * AllSamples, vector <EquivalenceClassFamily > * ECFams, int N , vector <GeneFamily *> * GeneFamilyList , bool boltzmann, bool Verbose, bool SuperVerbose, bool galwaysGain, double gC1Advantage)
{
	for(unsigned i = 0 ; i < ECFams->size(); i++)
	{

		EquivalenceClassFamily * ECF = &ECFams->at(i);

		AllSamples->push_back(new NECFsample);
		backtrackNtimesOneEquivalenceClassFamily(AllSamples->back(), ECF, N , GeneFamilyList , boltzmann, Verbose, SuperVerbose, galwaysGain, gC1Advantage);

		if(Verbose)
			cout << "backtracked ECF " << i << " " << N << " times." << endl;

	}

}


/*
Computes and Backtrack IN PLACE the EquivalenceClassFamilies's AdjMatrix

Takes:
 - vector < EquivalenceClassFamily> * ECFams : vector of equivalence class families containing refined equivalence class
 - vector <GeneFamily *> * GeneFamilyList  : vector fo gene families
 - map < string, map <string , double> > &adjacencyScores : map containing probas associated with adjacencies for ADseq
 - map<int,vector<float> > speciesC0C1 : for artdeco
 - map<int, map<string,int> > speGeneAdjNb : for artdeco
 - double Gcost : cost of an adjacency Gain
 - double Bcost  : cost of an adjacency Break
 - bool boltzmann : wether to use boltzmann computation or not
 - bool SubRecInAdj : if set to true, the weighted cost of a reconciliation event will be used to favor co-event in the adjacency matrix computation. Unavalaible for Boltzmann computation.
 - double WDupCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WLossCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - double WHgtCost : Cost of a single reconciliation event, weighted so that this score can be added to an adjacency score.
 - bool Verbose
 - bool SuperVerbose
 - bool galwaysGain : there is always a Gain at the top of an Adjacency tree. Will add a gain to c1 at the root of the equivalence class
 - double gC1Advantage : probability to choose c1 over c0 IF (and only if) they have the same score
 - double boltzmannTemp [default : 1] : Temperature used in boltzmann computation (a temperature of 0 is not possible)
 - double absencePenalty [default : -1] : if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization.
 - double adjScoreLogBase [default : 10000] : base of the log that will be used to go from adjacency confidence score to parsimony costs

*/
void ComputeAndBacktrackEquivalenceClassFamilies(vector < EquivalenceClassFamily> * ECFams, vector <GeneFamily *> * GeneFamilyList ,map < string, map <string , double> > &adjacencyScores,  map<int,vector<float> > speciesC0C1, map<int, map<string,int> > speGeneAdjNb, double Gcost, double Bcost , bool boltzmann, bool SubRecInAdj,double WDupCost, double WLossCost, double WHgtCost, bool Verbose, bool SuperVerbose, bool galwaysGain, double gC1Advantage, double boltzmannTemp , double absencePenalty, double adjScoreLogBase, bool interactionMode )
{
	ComputeEquivalenceClassFamilies( ECFams,  GeneFamilyList , adjacencyScores, speciesC0C1, speGeneAdjNb, Gcost,  Bcost ,  boltzmann,  SubRecInAdj, WDupCost,  WLossCost,  WHgtCost,  Verbose,  SuperVerbose,  boltzmannTemp ,  absencePenalty, adjScoreLogBase , interactionMode);
	backtrackInPlaceEquivalenceClassFamilies( ECFams, GeneFamilyList , boltzmann, Verbose, SuperVerbose, galwaysGain, gC1Advantage);
}


/*
Takes:
 - AForest (vector <AdjTree * > *) : vector of adjacency trees
 - CoEventSet (vector <CoEvent> * ) : vector to put the coevents in
 - EclassId (int): id of the Equivalence Class in the EquivalenceClassFamily
 - gfam1 (int) : id of the first gene family
 - gfam2 (int) : id of the second gene family
 - Rtree1 (ReconciledTree *) : Reconciled tree first gene family
 - Rtree2 (ReconciledTree *) : Reconciled tree second gene family
 - bool ignoreTime : whether time constraints can be ignored or not (specifically for the loading of old solution where ids of non bifurcating nodes are not kept (meaning that no timeSlice can be found for their adjacencies))

*/
void PopulatesCoeventsFromAdjForest( vector <AdjTree * > * AForest, vector <CoEvent> * CoEventSet, int EclassId , int gfam1, int gfam2, ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool ignoreTime)
{
	//cout << "PopulatesCoeventsFromAdjForest "<< EclassId << ":" << gfam1 << "-" << gfam2 << endl;
	//parsing the nodes of the adj forests.
	int nbadjtrees = AForest->size();

	for(unsigned j = 0; j < nbadjtrees; j++)
	{
		//cout << "plop " << j << endl;
		AdjTree * Atree = AForest->at(j);

		vector <int > Nids = Atree->getNodesId();

		//cout << "plopBis"<< endl;

		for(unsigned n = 0; n < Nids.size(); n++)
		{
			//cout << "looking at Atree node " << Nids[n] << endl;
			//cout << "coev? " << Atree->getNodeCoEventStatus(Nids[n]) << endl;
			if( Atree->getNodeCoEventStatus(Nids[n]) ) // this is a CoEvent
			{
				//cout << "in" << Nids[n] << "testing property : " << nodeid1 << " evt :" << Atree->getNodeEvent(Nids[n]) << endl;

				//int evt = Atree->getNodeEvent(Nids[n]);
				//if( ( evt == 0 ) || ( evt == 8))
				//	continue; // adjBreak or leaf
				//1. Determining time slices...
				int UTS = -1;
				int LTS = -1;
				if(! Atree->hasNodeProperty(Nids[n], nodeid1) ) //for leaves
					continue;

				//cout << "plop" << endl;

				pair <int,int> geneids = Atree->getNodeNodeIds(Nids[n]);

				//cout << "NodeNodeIds : " << geneids.first << "-" << geneids.second << endl;

				//cout << "PopulatesCoeventsFromAdjForest" << " tree " << j << " ; Anode " << Nids[n] << " ; gnodes " << geneids.first << "," << geneids.second << " , event "<< Atree->getNodeEvent(Nids[n]) << endl;
				//cout << "ignore time " << ignoreTime << endl;
				bool tmpAccountTime = ( Rtree1->hasNode(geneids.first) ) && ( Rtree2->hasNode(geneids.second) );
				if( (!ignoreTime) && (tmpAccountTime) )
				{
					if(Rtree1->getTimeSliceStatus() > 0) //time slices matters  
					{
						UTS = Rtree1->getNodeUpperBoundaryTS(geneids.first);
						LTS = Rtree1->getNodeLowerBoundaryTS(geneids.first);
					}
					if(Rtree2->getTimeSliceStatus() > 0) //time slices matters
					{
						if(UTS == -1)//uts and lts aren't set --> set them
						{
							UTS = Rtree2->getNodeUpperBoundaryTS(geneids.second);
							LTS = Rtree2->getNodeLowerBoundaryTS(geneids.second);
						}
						else
						{
							int tmp = Rtree2->getNodeUpperBoundaryTS(geneids.second);
							if(tmp < UTS)
								UTS = tmp;//replacing by the more  restrictive upper time slice
							tmp = Rtree2->getNodeLowerBoundaryTS(geneids.second);
							if(tmp > LTS)
								LTS = tmp;//replacing by the more  restrictive lower time slice
						}
					}
				}
				//cout << "looking a pre-existing co-ev" << endl;
				//2. checking compatibility with an existing coevent
				int coeventindex = -1;
				for(unsigned c = 0; c < CoEventSet->size(); c++)
				{
					if(CoEventSet->at(c).isCompatible(Nids[n], Atree, UTS, LTS)) // found a compatible coevent
					{
						coeventindex = c;
						break;
					}
				}
				//3. adding to the found coevent or creating one
				//cout << "adding to the found coevent or creating one : " << coeventindex << endl;
				
				if(coeventindex == -1)
				{
					//creating a coevent
					CoEventSet->push_back(CoEvent());
					coeventindex = CoEventSet->size() - 1 ;
				}
				//adding to the coevent:
				if( (!ignoreTime) && (tmpAccountTime) )
					CoEventSet->at(coeventindex).addAdj(EclassId, Nids[n], Atree,  gfam1, gfam2, Rtree1, Rtree2);//also checks if the AdjNode already exists, as well as the Genes. Sets the UTS/LTS 
				else
					CoEventSet->at(coeventindex).addAdj(EclassId, Nids[n], Atree,  gfam1, gfam2);//also checks if the AdjNode already exists, as well as the Genes.

				//cout << "coevent "<< coeventindex << "now has "<<  CoEventSet->at(coeventindex).getNumberOfGene() << " genes."<<endl;
			}
		}
		//cout << "plep" << endl;
	}
	//cout << "plepEND" << endl;
}



/*
Read the species tree and do basic checks.
Does not make any subdivision. 

Takes:
 - speciesFile (string): filename of the species tree
 - dateAsBootstrap (bool) [default: false]: if true, node ordering is read directly from the nodes bootstraps

Returns;
	(MySpeciesTree *): pointer to the species tree
 */
MySpeciesTree * getSpeciesTree(string speciesFile, bool dateAsBootstrap)
{
    // read species tree
    string errString = "";
    MySpeciesTree* speciesTree = MySpeciesTree::readMySpeciesTree( 
            speciesFile.c_str(),
            errString, 
            dateAsBootstrap );
    if( errString != "" || speciesTree == NULL ) {
        cerr << "Error reading species tree: " << errString << endl;
        exit(1);
    }

   
    // species tree must be binary
    int nonBinaryCount = speciesTree->getBinaryCount();
    if( nonBinaryCount != 0 ) {
        cerr << "ERROR: Species tree is not binary." << endl;
        cerr << nonBinaryCount << " nodes are not binary." << endl;
        exit(1);
    }

    // species tree must have unique leaves
    string dupName;
    if( !speciesTree->uniqueLeaves( dupName ) ) {
        cerr << "ERROR: Species tree has duplicate leaf names: <" 
             << dupName << ">" << endl;
        exit(1);
    }

    return speciesTree;
}



/** 
 * Process species tree: trimming, costs, subdivision, date changes
 * 
 * @arg MySpeciesTree * speciesTree : the species tree
 * @arg bool dateAsBootstrap : if true, node ordering is read directly from the nodes bootstraps
 * @arg bool dated : if true the tree is considered ultrametric and will be timeSliced. 
 * @arg bool computeT : if true, transfer are computed
 * @arg bool computeTD : if true, transfer from the dead are computed; no effect if computeT is false.
 * @arg double hgtCost : cost of an horizontal gene  transfer -> used for transfer from the dead
 * @arg double lossCost : cost of a loss -> used for transfer from the dead 
 * @arg bool verbose
 *
 * NB: In DeCo framework, dated and computeT should be equal but are kept separate here for maintenance / software change purpose
 * NB: Idem with computeT and computeTD
 * 
 * @return Maximum time slice. (int)
 */
int processSpeciesTree( MySpeciesTree *speciesTree , bool dateAsBootstrap, bool dated, bool computeT, bool computeTD , double hgtCost, double lossCost ,bool verbose )
{

	//3/3/2016 : no date swap/change in DeCo for now. Uncomment and add correct arugment and missing function to enable
    //// parse data changing parameters
    //vector< pair<int, int> > dateMap;
    //string errStr = "empty error";
    //if( gStringParams.find("date.change.swap")->second != "" ) {
    //    if( !parseDateSwap( dateMap, errStr ) ) {
    //        cerr << "ERROR: " << errStr << endl;
    //        exit(1);
    //    }
    //}
    //if( gStringParams.find("date.changes")->second != "" ) {
    //    if( !parseDateChanges( changedTimeSlices, errStr ) ) {
    //        cerr << "ERROR: " << errStr << endl;
    //        exit(1);
    //    }
    //}

	//3/3/2016 : no triming in DeCo for now. Uncomment and add correct arugment and missing function to enable
    // trim tree
    //boost::unordered_map<string, int> taxaNamesGenes;
    //if( gBoolParams.find("trim.species.tree")->second ) {
    //    // create map of all taxa name in the gene trees
    //    BOOST_FOREACH( MyGeneTree *geneTree, geneTrees ) {
    //        vector<string> leaves = geneTree->getLeavesNames();
    //        BOOST_FOREACH( string leafName, leaves) {
    //            size_t pos = leafName.find( 
    //                        gStringParams.find("char.sep")->second[0] ); 
    //            string taxaName = leafName.substr(0,pos);
    //            taxaNamesGenes.insert( make_pair(taxaName,1) );
    //        }
    //    }
    //   
    //    // do undated or partially dated trimming here
    //    if( gIntParams.find("dated")->second != 2 ) {
    //        // trim species tree
    //        //if(!speciesTree->restrictTreeToLCA(taxaNamesGenes, gBoolParams.find("verbose")->second))
    //        if( !speciesTree->trimTree( taxaNamesGenes, 
    //                    gBoolParams.find("verbose")->second ) ) 
    //        {
    //            cerr << "Species tree has none of the gene taxa." << endl;
    //            exit(1);
    //        }
    //    }
    //}

    //3/3/2016 : no species specific costs in DeCo for now. Uncomment and add correct arugment and missing function to enable
    //// read variable costs
    //if( gStringParams.find("costs.file")->second != "none" ) {
    //    if( !speciesTree->assignCosts( 
    //                gStringParams.find("costs.file")->second ) )
    //        exit(1);
    //    gFixedCosts = false;
    //}

    // to add the outGroup to transfer from the dead 
    // create a sibling of the root and a new root
    if( computeTD && computeT )  
        speciesTree->addAlphaForDeadTransfer( dateAsBootstrap, hgtCost, lossCost );
   
    // must be done before subdivision
    speciesTree->compute_RealPostOrder();  

	//two parameters harshly set to a given value: species tree time slices are unchanged
	vector<int> changedTimeSlices;
	vector< pair<int, int> > dateMap;

	string errStr = "";

    // subdivide tree if dated
    if( dated ) {
        bool success = speciesTree->computeSubdivision( dateMap, 
                        dateAsBootstrap, 
                        dated, // if the species tree is dated, it is presumed to be ultrametric 
                        changedTimeSlices, errStr ); 
        if( !success ) {
            cerr << "ERROR: " << errStr << endl;
            exit(1);
        }
        if( verbose ) {
            BOOST_FOREACH(int ts, changedTimeSlices ) 
                cout << "dates changed: " << ts << endl;

            bpp::ApplicationTools::displayTime(
                    "Computing the dated subdivision done:");
            cout << "leaves and nodes: " << speciesTree->getNumberOfLeaves()
                 << " " << speciesTree->getNumberOfNodes() << endl;
        }

        //// dated tree trimming -> no trimming in DeCo for now
        //if( gBoolParams.find("trim.species.tree")->second ) {
        //    if( !gBoolParams.find("compute.TD")->second ) {
        //        cerr << "Trimming with dated species tree requires "
        //            << "transfer from the dead." << endl;
        //        // The dated tree can have transfers at each timeslice,
        //        // which is a node in the original species tree. Trimming
        //        // the species tree removes these timeslices and therefore
        //        // possiblities to transfer at these times. To compensate,
        //        // transfer from the dead can be used.
        //        exit(1);
        //    }
        //
        //    // trim species tree
        //    //if(!speciesTree->restrictTreeToLCA(taxaNamesGenes, gBoolParams.find("verbose")->second))
        //    if( !speciesTree->trimTree( taxaNamesGenes, 
        //        gBoolParams.find("verbose")->second ) ) 
        //    {
        //        cerr << "Species tree has none of the gene taxa." << endl;
        //        exit(1);
        //    }
        //}
    } else {
        speciesTree->assignNoSubdivisionTimeSlices();
    }

    // assign ids (correspondance)
    speciesTree->assignPostOrderIds();

    // create vector of nodes by time slices
    int maxTS = speciesTree->setVectorTimeSlices();

    //3/3/2016 : no species specific costs in DeCo for now. Uncomment and add correct arugment and missing function to enable
    //if( !gFixedCosts ) {
    //    // check input date costs
    //    vector<MySpeciesNode*> nodes = speciesTree->getNodes();
    //    BOOST_FOREACH( MySpeciesNode* node, nodes ) {
    //        if( node->getInfos().isAlpha ) {
    //            if( node->getInfos().hgtCost == 0
    //                || node->getInfos().lossCost == 0 ) 
    //            {
    //                cout << "BAD ALPHA COST=========" << endl;
    //                exit(1);
    //            }
    //        } else if( !node->hasFather() ) {
    //            // root doesn't matter
    //        } else if( node->getInfos().duplicationCost == 0 
    //                || node->getInfos().hgtCost == 0
    //                || node->getInfos().lossCost == 0 ) 
    //        {
    //            cout << "BAD COST=========" 
    //                 << node->getId() << endl;
    //            exit(1);
    //        }
    //    }
    //}

	//3/3/2016 : no other species tree in DeCo for now. Uncomment and add correct arugment and missing function to enable
    // process another species file
    //if( gStringParams.find("other.species.file")->second != "none" ) 
    //    processOtherSpeciesTree( speciesTree, dateMap, changedTimeSlices );

    return maxTS;
}








/*
Checks if the gene tree has enough leaves and nop duplicate leaves. exit if problem
*/
void checkGeneTree( int counter, MyGeneTree *geneTree)
{
    // needed by the algoritm!!!	
    if( geneTree->getNumberOfLeaves()<3 ) {
        cerr << "ERROR: Gene tree " << counter
             << " has only " << geneTree->getNumberOfLeaves()
             << " leaves. At least three are required." << endl;
        exit(1);
    }

    // gene trees must have unique leaves
    string dupName;
    if( !geneTree->uniqueLeaves( dupName ) ) {
        cerr << "ERROR: A gene tree " << counter
            << " has duplicate leaf names: <" 
             << dupName << ">" << endl;
        exit(1);
    }
}






/*
* @arg string ListGeneFile: filename of the file containing the filenames of the gene tree distribution files
*
* @return vector<string> : vector of the filenames of the gene tree distribution files
*/
vector<string> readGeneDistributionsFile(string ListGeneFile)
{
	//// read gene trees distributions
	ifstream fileStream(ListGeneFile.c_str());
	
	if( !fileStream.is_open() ) 
	{
	  	cout << "Could not open gene distribution files file."<< ListGeneFile <<endl;
		exit(1);
	}
	vector<string> geneFiles;
	while( !fileStream.eof() ) 
	{
		string line;
		getline( fileStream, line );
		boost::trim( line );

		if( line == "" ) 
		{
				// ignore blank lines
		}
		else
		{
			geneFiles.push_back(line);
		}
	}
	return geneFiles;
}

/*
*
* @arg vector <GeneFamily *> * GeneFamilyList : vector of GeneFamily where gene distribution will be added
* @arg MySpeciesTree *speciesTree 
* @arg vector<string> geneFiles : vector of gene distribution filenames
* @arg bool ale : true if the files are actually in the .ale format
* @arg bool reconciled : true if the files actually contains a reconciled tree in PhyloXML format (negated by ale)
* @arg char charSep : separator between species id and gene id
* @arg bool verbose 
* @arg bool superverbose
* @arg bool rooted (default = false)
*/
void readGeneDistributions(vector <GeneFamily *> * GeneFamilyList, MySpeciesTree *speciesTree ,vector<string> geneFiles, bool ale, bool reconciled, char charSep ,bool verbose, bool superverbose, bool rooted)
{
	//reading the gene trees and putting them in GeneFamilies
	bool overflow = false;
	string errString;
	for(size_t i = 0; i < geneFiles.size(); i++)
	{
		string geneFile = geneFiles[i];

		if(verbose)
		{
			cout << "Reading gene Tree Distribution :"<< geneFile << endl;
		}

		if(!ale)
		{
			if(!reconciled)
			{
				vector<MyGeneTree*> geneTrees;
			

				geneTrees = MyGeneTree::readMyGeneTrees( geneFile.c_str(), errString );//reading the tree list
				if( errString != "" ) 
				{
					cerr << "Error reading gene trees: " << errString << endl;
					exit(1);
				}
				if( verbose ) 
					cout << geneTrees.size() << " gene trees" << endl;
	
//				for(size_t j = 0 ; j < geneTrees.size();j++)
//					checkGeneTree(j,geneTrees[j]);
				if(!rooted)
					GeneFamilyList->push_back(new GeneFamily(geneTrees, charSep,verbose, superverbose));
				else // rooted case
					GeneFamilyList->push_back(new GeneFamily( *(geneTrees[0]), charSep,verbose, superverbose, rooted));
			}
			else // reconciled tree
			{
				// old way
				//ReconciledTree Rtree(geneFile , speciesTree, superverbose);
				//GeneFamilyList->push_back(new GeneFamily(Rtree, verbose, superverbose));
				//
				//if(superverbose)
				//	GeneFamilyList->back()->printRecTree();
				
				// adds as many gene family as there are reconciled tree in the file.
				
				//int previousS = GeneFamilyList->size();

				ReadRecPhyLoXMLFile(   GeneFamilyList,  geneFile , speciesTree,  verbose, superverbose); 

				//if(verbose)
				//{
				//	cout << "read " <<  GeneFamilyList->size() - previousS << " reconciled tree."<<endl;
				//}
			}
		}
		else //using ale files
		{
			GeneFamilyList->push_back(new GeneFamily(geneFile.c_str(), charSep,verbose, superverbose));
		}
		cout << "geneFile "<<geneFile<<endl;
		GeneFamilyList->at(GeneFamilyList->size() -1 )->setSourceFile(geneFile);
	}
}

/*
* @arg string AdjacencyFile : name of the file containing adjacencies (one per line, leaf name separated by a space.)
*
* @return vector< pair <string,string > > : list of adjacencies
*/
vector< pair <string,string > > readAdjacencies(string AdjacencyFile)
{

	ifstream AdjStream(AdjacencyFile.c_str());

	vector< pair <string,string > > adjacencies;
	pair <string,string> currentadj;
	
	while( !AdjStream.eof() ) 
	{
		AdjStream >> currentadj.first;

		if(AdjStream.eof()) // to avoid 
			continue;
		AdjStream >> currentadj.second;
		adjacencies.push_back(currentadj);
	}
	
	return adjacencies;
}

/*
* @arg string AdjacencyFile : name of the file containing adjacencies (one per line, leaf name separated by a space.)
* @arg map < string, map <string , double> > &adjacencyScores : map containing probas associated with adjacencies for ADseq
* @arg vector < pair <int,int> > &adjacenciesOrientations : the int indicated which extremity of the gene is actually part of the adjacency. (-1 : start ; 0 : unspecified ; 1 : stop)
*
* @return vector< pair <string,string > > : list of adjacencies
*/
vector< pair <string,string > > readAdjacencies(string AdjacencyFile , map < string, map <string , double> > &adjacencyScores, vector < pair <int,int> > &adjacenciesOrientations, bool verbose)
{
	string line;
	ifstream AdjStream(AdjacencyFile.c_str());

	vector< pair <string,string > > adjacencies;
	pair <string,string> currentadj;
	
	while( getline(AdjStream , line) ) 
	{
		boost::trim( line );
		if( line  == "" )
			continue; // ingore blanck lines

		//cout << line << endl;
		istringstream ss(line);

		vector < string > parsedLine;
		while( ! ss.eof() )
		{
			string s;
			ss >> s;
			parsedLine.push_back(s);
		}

		//for(int i = 0 ; i < parsedLine.size(); i++)
		//	cout << parsedLine[i] << " ";
		//cout << endl;

		if(parsedLine.size() < 2)
		{
			if(verbose)
				cout << "found a line with only 1 element. Ignoring."<< endl;
			continue;			
		}

		currentadj.first = parsedLine[0];
		currentadj.second = parsedLine[1];


		int sens1 = 0;
		int sens2 = 0;
		int scoreIndex = -1;
		if( parsedLine.size() == 3 ) // format : g1 g2 p
		{
			scoreIndex = 2;
		}
		else if( parsedLine.size() > 3 )// format : g1 g2 -|+ -|+ [p]
		{
			if(parsedLine.size() == 5)
			{ // presence of the score
				scoreIndex = 4;
			}

			//cout << "orientation : " << parsedLine[2] << " " << parsedLine[3] << endl;

			if( parsedLine[2] == "+" ) // --> *
				sens1 = 1; // this is the stop of the first gene
			else if( parsedLine[2] == "-" ) // <-- *
				sens1 = -1; // this is the start of the first gene

			if( parsedLine[3] == "+" )
				sens2 = -1; // this is the start of the first gene
			else if( parsedLine[3] == "-" )
				sens2 = 1; // this is the stop of the first gene
				
		}

		//cout << "  parsed : " << sens1 << " " << sens2 << " " << scoreIndex<< endl;

		if(scoreIndex != -1)
		{ // there is a score
			double score;

			if( ! ( StrToDouble( parsedLine[scoreIndex] , score ) ) )
			{
				if(verbose)
					cout << "Failed to read the adjacency ("<< currentadj.first <<"-" << currentadj.second <<") score. Using default score of 1 instead." << endl;
				score = 1;
			}

			if( (score < 0) || ( score > 1 ) )
			{
				if(verbose)
					cout << "The adjacency ("<< currentadj.first <<"-" << currentadj.second <<") score must be between 0 and 1 (read " << score <<"). Using default score of 1 instead." << endl;
				score = 1;
			}

			adjacencyScores[currentadj.first][currentadj.second] = score;
	
			//if(verbose)
			//	cout << "read adj score"<< currentadj.first << "-" << currentadj.second << " -> "<< score <<endl;

		}

		adjacenciesOrientations.push_back( pair < int , int > (sens1 , sens2) );

		adjacencies.push_back(currentadj);
	}
	
	return adjacencies;

}



/*
* @arg MySpeciesTree * speciesTree : pointer to the species tree
* @arg bool newick : if true the tree will be written in newick format. otherwise it will be in PhyloXML
* @arg string prefix : the prefix (containing path) of the filename
*
* @return : the name of the written file
*/
string WriteSpeciestree(MySpeciesTree * speciesTree, bool newick, string prefix, string char_sep)
{

	string filename = prefix;
	if(( filename[ filename.size() - 1 ] != '/')&&(filename.size() > 0))
		filename += ".";

	if(newick)
	{
		filename += "speciesTree.newick";
		
		MySpeciesTree *tree = speciesTree->getPostorderTree();

		// correspondance between leaf names and RPO
		vector <MySpeciesNode *> leaves = tree->getLeaves();

		for( unsigned i = 0 ; i < leaves.size() ; i ++)
		{
			if( leaves[i]->hasName() ) 
			{
				string name = leaves[i]->getName();
				MySpeciesNode * N = speciesTree->getNode(name);
				leaves[i]->setName( name ) ;//+ char_sep + static_cast<ostringstream*>( &(ostringstream() << N->getInfos().realPostOrder ) )->str()  );
				//cout << N->getInfos().realPostOrder << " "<< name<<endl;
			}
		}

        tree->printNewick( filename.c_str() );
        delete tree;
	}
	else
	{
		DeCoOutputManager DOM;
		ofstream ofs;

		filename += "speciesTree.phyloxml";
		ofs.open(filename.c_str(),ofstream::out);		
		DOM.WritePhyloXMLSpeTree(ofs, speciesTree);
		ofs.close();
	}
	return filename;
}


/*
* @arg vector < GeneFamily * > * GeneFamilyList 
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg string prefix : the prefix (containing path) of the filename
*
* @ret vector <string> : vector of the names of the files the trees were written in
*/
vector <string> WriteGeneFamilyReconciledTrees( vector < GeneFamily * > * GeneFamilyList, bool newick,bool hideLosses, string prefix)
{
	DeCoOutputManager DOM;

	vector <string> FILENAMES;

	for(unsigned i = 0 ; i < GeneFamilyList->size(); i++)
	{
		string filename = prefix + "geneFamily" + static_cast<ostringstream*>( &(ostringstream() << i) )->str() ;
		if(newick)
			filename += ".newick";
		else
			filename += ".phyloxml";
		ofstream ofs;
		ofs.open(filename.c_str(),ofstream::out);
		DOM.WriteRecTree(ofs, GeneFamilyList->at(i)->getRecTree(), newick, hideLosses);
		ofs.close();

		FILENAMES.push_back(filename);
	}
	return FILENAMES;
}

/*
* @arg GeneFamily *  GFam
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg string prefix : the prefix (containing path) of the filename
* @arg int index : the index of the GeneFamily 
*
* string : name of the file the tree was written in
*/
string WriteOneGeneFamilyReconciledTree( GeneFamily * GFam, bool newick,bool hideLosses, string prefix, int index) // not used anymore. 
{
	DeCoOutputManager DOM;


	string filename = prefix + "geneFamily" + static_cast<ostringstream*>( &(ostringstream() << index) )->str() ;
	if(newick)
		filename += ".newick";
	else
		filename += ".phyloxml";

	ofstream ofs;
	ofs.open(filename.c_str(),ofstream::out);

	DOM.WriteRecTree(ofs, GFam->getRecTree(), newick, hideLosses);

	ofs.close();

	return filename;
}

/*
* @arg GeneFamily *  GFam
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg string fileName : the name of the file to add the tree to
* @arg int index : the index of the GeneFamily 
* 
*/
void AddOneGeneFamilyReconciledTreeToFile( GeneFamily * GFam, bool newick,bool hideLosses, string fileName, int index)
{
	DeCoOutputManager DOM;

	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);

	DOM.WriteRecTree(ofs, GFam->getRecTree(), newick, hideLosses, index);

	ofs.close();
}



/*
Takes:
	- string filename : name of the file to open
	- bool newick: whether to write in newick or recPhyloXML
*/
void InitReconciledTreeFile(string fileName, bool newick)
{
	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out);

	if(!newick)
	{
		ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" <<endl;
		ofs << "<recPhylo xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\" xmlns=\"http://www.recg.org\" >" << endl;
	}

}

/*
Takes:
	- string filename : name of the file to open
	- bool newick: whether to write in newick or recPhyloXML
*/
void FinishReconciledTreeFile(string fileName, bool newick)
{
	if(!newick)
	{
		ofstream ofs;
		ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
		ofs << "</recPhylo>" << endl;
	}
}


/*
Takes:
 - string filenName: name of the file to wrtie in
 - bool newick : if true: newick format is used. Otherwise, some recPhyloXML derived format is used
*/
void InitAdjacencyTreeFile(string fileName, bool newick)
{
	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out);

	if(!newick)
	{
		ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" <<endl;
		ofs << "<adjPhylo xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.recg.org ./recGeneTreeXML.xsd\" xmlns=\"http://www.recg.org\" >" << endl;
	}

}


/*
* @arg string fileName : name of the file to add in
* @arg EquivalenceClassFamily * ECF : pointer the the ECF whose trees we want to write
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg double GainCost (default = 3) : cost of a single gain event
* @arg double BreakCost (default = 1) : cost of a single break event
* 
*/
void AddECFAForestToFile(string fileName, EquivalenceClassFamily * ECF, bool newick,bool hideLosses, double GainCost ,  double BreakCost)
{
	int fam1 = ECF->getGfamily1();
	int fam2 = ECF->getGfamily2();

	int sens1 = ECF->getSens1();	
	int sens2 = ECF->getSens2();


	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
	DeCoOutputManager DOM;


	if(!newick)
	{
		ofs << "  <EquivalenceClassFamily fam1=\""<<fam1;
	
		if( sens1 == 1 )
			ofs << "_stop";
		else if( sens1 == -1 )
			ofs << "_start";
	
		ofs <<"\" fam2=\""<<fam2;
	
		if( sens2 == 1 )
			ofs <<  "_stop";
		else if( sens2 == -1 )
			ofs << "_start";
	
		ofs<<"\">"<< endl;
	}
	else
	{
		//line indicating information
		ofs << ">";
	
		ofs <<  "Adjacencies between families: " << fam1;
	
		if( sens1 == 1 )
			ofs << "_stop";
		else if( sens1 == -1 )
			ofs << "_start";
	
		ofs << "-" << fam2;
	
		if( sens2 == 1 )
			ofs <<  "_stop";
		else if( sens2 == -1 )
			ofs << "_start";
	
		ofs << " number of trees: "<< ECF->getNbAdjTrees() ;
		ofs << " #Gain: " << ECF->getNbAdjGain();
		ofs << " #Break: " << ECF->getNbAdjBreak();
		ofs << " score: " << ECF->getNbAdjGain() * GainCost + ECF->getNbAdjBreak() * BreakCost ;
		ofs << endl;
	}		


	int line_indent = 2;

	for(size_t j = 0; j < ECF->getNbEqClasses() ; j++)
	{
		DOM.WriteAdjForest(ofs, ECF->getAdjForest(j),newick, hideLosses, line_indent);
	}


	if(!newick)
	{
		ofs << "  </EquivalenceClassFamily>"<< endl;
	}

	ofs.close();
	return ;
}



/*
* @arg string fileName : name of the file to add in
* @arg EquivalenceClassFamily * ECF : pointer the the ECF whose trees we want to write
* @arg NECFsample * sample : the N samples to write
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg double GainCost (default = 3) : cost of a single gain event
* @arg double BreakCost (default = 1) : cost of a single break event
* 
*/
void AddECFNsampleToFile(string fileName, EquivalenceClassFamily * ECF, NECFsample * sample, bool newick,bool hideLosses, double GainCost, double BreakCost)
{
	int fam1 = ECF->getGfamily1();
	int fam2 = ECF->getGfamily2();

	int sens1 = ECF->getSens1();	
	int sens2 = ECF->getSens2();


	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
	DeCoOutputManager DOM;


	if(!newick)
	{
		ofs << "  <EquivalenceClassFamily fam1=\""<<fam1;
	
		if( sens1 == 1 )
			ofs << "_stop";
		else if( sens1 == -1 )
			ofs << "_start";
	
		ofs <<"\" fam2=\""<<fam2;
	
		if( sens2 == 1 )
			ofs <<  "_stop";
		else if( sens2 == -1 )
			ofs << "_start";
	
		ofs<<"\">"<< endl;
	}


	int line_indent = 3;


	for(size_t sampleIndex = 0; sampleIndex < sample->size() ; sampleIndex++)
	{
		if(newick)
		{
			//line indicating
			ofs << ">";
		
			ofs <<  "Id: " << fam1;
		
			if( sens1 == 1 )
				ofs << "_stop";
			else if( sens1 == -1 )
				ofs << "_start";
		
			ofs << "-" << fam2;
		
			if( sens2 == 1 )
				ofs <<  "_stop";
			else if( sens2 == -1 )
				ofs << "_start";

			ofs << " sample: " << sampleIndex;

			int nbTrees = 0;
			int nbGain = 0;
			int nbBreak = 0;
			for(size_t j = 0; j < sample->at(sampleIndex)->size() ; j++)
			{
				nbTrees += sample->at(sampleIndex)->at(j)->size();
				for(size_t k = 0; k < sample->at(sampleIndex)->at(j)->size() ; k++)
				{
					nbBreak += sample->at(sampleIndex)->at(j)->at(k)->countNbBreak();
					nbGain += sample->at(sampleIndex)->at(j)->at(k)->countNbGain();
				}
			}

			ofs << " #Trees: "<< nbTrees ;
			ofs << " #Gain: " << nbGain;
			ofs << " #Break: " << nbBreak;
			ofs << " Score: " << nbGain * GainCost + nbBreak * BreakCost ;
			ofs << endl;

		}
		else
		{
			ofs << "    <sample number=\""<<sampleIndex<<"\">"<< endl;
		}
		

		for(size_t j = 0; j < sample->at(sampleIndex)->size() ; j++)
		{
			DOM.WriteAdjForest(ofs, sample->at(sampleIndex)->at(j),newick, hideLosses, line_indent);
		}


		if(!newick)
		{
			ofs << "    </sample>"<< endl;
		}
	}

	if(!newick)
	{
		ofs << "  </EquivalenceClassFamily>"<< endl;
	}

	ofs.close();
	return ;
}

/*
* @arg string fileName : name of the file to add in
* @arg EquivalenceClassFamily * ECF : pointer the the ECF whose trees we want to write
* @arg NECFsample * sample : the N samples to write
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg double GainCost (default = 3) : cost of a single gain event
* @arg double BreakCost (default = 1) : cost of a single break event
* @arg bool Init ( default : false) : if true, it will write some massage to signal the start of this ECF (XML only)
* @arg bool Finish ( default : false) : if true, it will write some massage to signal the end of this ECF (XML only)
* 
*/
void AddECFsampleToFile(string fileName, EquivalenceClassFamily * ECF, ECFsample * sample,
							 int sampleIndex, bool newick,bool hideLosses,  
							 double GainCost ,  double BreakCost , bool Init, bool Finish)
{
	int fam1 = ECF->getGfamily1();
	int fam2 = ECF->getGfamily2();

	int sens1 = ECF->getSens1();	
	int sens2 = ECF->getSens2();


	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
	DeCoOutputManager DOM;


	if( (!newick) && (Init) )
	{
		ofs << "  <EquivalenceClassFamily fam1=\""<<fam1;
	
		if( sens1 == 1 )
			ofs << "_stop";
		else if( sens1 == -1 )
			ofs << "_start";
	
		ofs <<"\" fam2=\""<<fam2;
	
		if( sens2 == 1 )
			ofs <<  "_stop";
		else if( sens2 == -1 )
			ofs << "_start";
	
		ofs<<"\">"<< endl;
	}


	int line_indent = 3;

	if(newick)
	{
		//line indicating
		ofs << ">";
	
		ofs <<  "Id: " << fam1;
	
		if( sens1 == 1 )
			ofs << "_stop";
		else if( sens1 == -1 )
			ofs << "_start";
	
		ofs << "-" << fam2;
	
		if( sens2 == 1 )
			ofs <<  "_stop";
		else if( sens2 == -1 )
			ofs << "_start";
		ofs << " sample: " << sampleIndex;
		int nbTrees = 0;
		int nbGain = 0;
		int nbBreak = 0;
		for(size_t j = 0; j < sample->size() ; j++)
		{
			nbTrees += sample->at(j)->size();
			for(size_t k = 0; k < sample->at(j)->size() ; k++)
			{
				nbBreak += sample->at(j)->at(k)->countNbBreak();
				nbGain += sample->at(j)->at(k)->countNbGain();
			}
		}
		ofs << " #Trees: "<< nbTrees ;
		ofs << " #Gain: " << nbGain;
		ofs << " #Break: " << nbBreak;
		ofs << " Score: " << nbGain * GainCost + nbBreak * BreakCost ;
		ofs << endl;
	}
	else
	{
		ofs << "    <sample number=\""<<sampleIndex<<"\">"<< endl;
	}
	
	for(size_t j = 0; j < sample->size() ; j++)
	{
		DOM.WriteAdjForest(ofs, sample->at(j),newick, hideLosses, line_indent);
	}

	if(!newick)
	{
		ofs << "    </sample>"<< endl;
	}




	if( (!newick) && (Finish) )
	{
		ofs << "  </EquivalenceClassFamily>"<< endl;
	}

	ofs.close();
	return ;	
}


/*
Takes:
	- string fileName
	- EquivalenceClassFamily * ECF
	- bool hideLosses
	- int fam1 
	- int fam2 
	- int ECindex
	- int sens1 = 0 
	- int sens2 = 0
	- int sample = -1
*/
void AddAForestToFileNewick(string fileName, AFOREST * forest,bool hideLosses , 
								int fam1 , int fam2 , int ECindex, int sens1 , int sens2 , int sample )
{
	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
	DeCoOutputManager DOM;

	//line indicating
	ofs << ">";

	ofs <<  "Adjacencies between families: " << fam1;

	if( sens1 == 1 )
		ofs << "_stop";
	else if( sens1 == -1 )
		ofs << "_start";

	ofs << "-" << fam2;

	if( sens2 == 1 )
		ofs <<  "_stop";
	else if( sens2 == -1 )
		ofs << "_start";


	//ofs << " EquivalenceClass: "<< ECindex;

	if( sample != -1 )
		ofs << " sample: " << sample;

	ofs << " number of trees: "<< forest->size() << endl;

	ofs << endl;

	

	ofs.close();

}

/*
Takes:
	- string fileName
	- EquivalenceClassFamily * ECF
	- bool hideLosses
	- int fam1 
	- int fam2 
	- int ECindex
	- int sens1 = 0 
	- int sens2 = 0
	- int sample = -1
*/
void AddAForestToFileXML(string fileName, AFOREST * forest,bool hideLosses ,
							 int fam1 , int fam2 , int ECindex, int sens1 , int sens2 , int sample  )
{

	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
	DeCoOutputManager DOM;

	//ofs << "    <EquivalenceClass index=\"" << ECindex << "\"";
	if(sample !=-1)
		ofs << "    <sample number=\""<<sample<<"\">"<< endl;

	int line_indent = 2;
	if(sample !=-1)
		line_indent++;

	DOM.WriteAdjForest(ofs, forest,false, hideLosses, line_indent);	

	if(sample !=-1)
		ofs << "    </sample>"<< endl;
	ofs.close();
}


/*
Takes:
 - string filenName: name of the file to wrtie in
 - bool newick : if true: newick format is used. Otherwise, some recPhyloXML derived format is used
*/
void FinishAdjacencyTreeFile(string fileName, bool newick)
{
	if(!newick)
	{
		ofstream ofs;
		ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
		ofs << "</adjPhylo>" << endl;
		ofs.close();
	}
}



/*
* @arg EquivalenceClassFamily * ECF : pointer the the ECF whose trees we want to write
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg string prefix : the prefix (containing path) of the filename
*
* @return vector of the writtent file names
*/
vector <string> WriteECFamTrees(EquivalenceClassFamily * ECF, bool newick, bool hideLosses , string prefix)
{
	DeCoOutputManager DOM;
	vector <string> fileNameList;

	int g1 = ECF->getGfamily1();
	int g2 = ECF->getGfamily2();

	string filePrefix = prefix + "EqClass_" + static_cast<ostringstream*>( &(ostringstream() << g1) )->str();

	if( ECF->getSens1() == 1 )
		filePrefix += "_stop";
	else if( ECF->getSens1() == -1 )
		filePrefix += "_start";

	filePrefix += "-" + static_cast<ostringstream*>( &(ostringstream() << g2) )->str();

	if( ECF->getSens2() == 1 )
		filePrefix += "_stop";
	else if( ECF->getSens2() == -1 )
		filePrefix += "_start";

	//filePrefix += ")";

	for(size_t j = 0; j < ECF->getNbEqClasses() ; j++)//NB: with that method, some file might not contain any tree...
	{
		string filename = filePrefix + "_" + static_cast<ostringstream*>( &(ostringstream() << j) )->str() + ".adjtrees";

		if(newick)
			filename += ".newick";
		else
			filename += ".phyloxml";

		ofstream ofs;
		ofs.open(filename.c_str(),ofstream::out);
		DOM.WriteAdjForest(ofs, ECF->getAdjForest(j),newick, hideLosses);
		ofs.close();
		fileNameList.push_back(filename);
	}
	return fileNameList;
}

/*
* @arg ECFsample * sample : pointer to the sample
* @arg int g1 : index of the first gene family. For filename purposes
* @arg int g2 : index of the first gene family. For filename purposes
* @arg int sampleIndex : index of the sample. For filename purposes
* @arg bool newick : if true the trees will be written in newick format. otherwise they will be in recPhyloXML
* @arg bool hideLosses : if true, losses and the branches leading to them wil be removed from the newick string
* @arg string prefix : the prefix (containing path) of the filename
*
* @return vector of the writtent file names
*/
vector <string> WriteECFsample(EquivalenceClassFamily * ECF,ECFsample * sample,int g1, int g2, int sampleIndex, bool newick, bool hideLosses , string prefix)
{
	DeCoOutputManager DOM;

	

	string filePrefix = prefix + "EqClass_" + static_cast<ostringstream*>( &(ostringstream() << g1) )->str();

	if( ECF->getSens1() == 1 )
		filePrefix += "_stop";
	else if( ECF->getSens1() == -1 )
		filePrefix += "_start";

	filePrefix += "-" + static_cast<ostringstream*>( &(ostringstream() << g2) )->str();

	if( ECF->getSens2() == 1 )
		filePrefix += "_stop";
	else if( ECF->getSens2() == -1 )
		filePrefix += "_start";


	filePrefix += "_sample" + static_cast<ostringstream*>( &(ostringstream() << sampleIndex) )->str() ;

	vector <string> fileNameList;

	for(size_t j = 0; j < sample->size() ; j++)//NB: with that method, some file might not contain any tree...
	{
		string filename = filePrefix + "_" + static_cast<ostringstream*>( &(ostringstream() << j) )->str() + ".adjtrees";

		if(newick)
			filename += ".newick";
		else
			filename += ".phyloxml";

		ofstream ofs;
		ofs.open(filename.c_str(),ofstream::out);
		DOM.WriteAdjForest(ofs, sample->at(j), newick, hideLosses);
		ofs.close();

		fileNameList.push_back(filename);
	}

	return fileNameList;
}


/*
* Writes adjacencies for a bunch of equivalence class families
*
* @arg vector < EquivalenceClassFamily > * ECFams
* @arg string prefix
*
*/
void WriteAdjacencies(vector < EquivalenceClassFamily > * ECFams, string prefix)
{
	DeCoOutputManager DOM;

	string filename = prefix;
	if( filename[ filename.size() - 1 ] != '/')
		filename += ".";
	filename += "adjacencies.txt";

	ofstream ofs;
	ofs.open(filename.c_str(),ofstream::out);

	for(size_t i = 0 ; i < ECFams->size();  i++)
		for(size_t j = 0; j < ECFams->at(i).getNbEqClasses() ; j++)//NB: with that method, some file might not contain any tree...
			DOM.WriteAdjForestAdjacencies(ofs, ECFams->at(i).getAdjForest(j) , ECFams->at(i).getSens1() , ECFams->at(i).getSens2() );

	ofs.close();

}


void WriteAdjacenciesOneECF( EquivalenceClassFamily * ECF, string prefix, ReconciledTree * Rtree1, ReconciledTree * Rtree2)
{
	DeCoOutputManager DOM;

	string filename = prefix;
	if( filename[ filename.size() - 1 ] != '/')
		filename += ".";
	filename += "adjacencies.txt";

	ofstream ofs;
	ofs.open(filename.c_str(),ofstream::out | ofstream::app);

	
	for(size_t j = 0; j < ECF->getNbEqClasses() ; j++)
		DOM.WriteAdjForestAdjacencies(ofs, ECF->getAdjForest(j), ECF->getSens1() , ECF->getSens2() , Rtree1, Rtree2);

	ofs.close();	
}



map < string, int > storeSpeciesChrNumber(string chrNumberFile)
{
	//// read file with chromosome number corresponding to species
	ifstream fileStream(chrNumberFile.c_str());
	
	if( !fileStream.is_open() ) 
	{
	  	cout << "Could not open chromosome number by species file."<< chrNumberFile <<endl;
		exit(1);
	}

	map < string, int > speciesChrNb;
	while (!fileStream.eof())
	{
		string species;
		unsigned int  chrNb;
		fileStream>>species;
		fileStream>>chrNb;
		speciesChrNb[species]=chrNb;
	} 
	return speciesChrNb;
}

/*
adds the name of genes at ancestral speciation and leaves to the specified file

Takes:  
	- ReconciledTree * RT : reconciled tree of the family
	- int famIndex : index of the gene family
	- string fileName : name of the file to add to
	- map<int, string> SpIdToSpName : map associating the species tree leaves RealPostOrder info to their name
*/
void WriteGeneOneGF( ReconciledTree * RT, int famIndex, string fileName, map<int, string> &SpIdToSpName)
{
	//cout << "WriteGeneOneGF "<< famIndex << " -> "<< fileName << endl;

	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);

	//1. internal 
	vector <int> nodeList = RT->getInnerNodesId();

	for(unsigned i = 0 ; i < nodeList.size(); i++)
	{
		int ev = RT->getNodeEvent(nodeList[i]);
		if( RT->isSpeciation(ev) )
		{
			ofs << RT->getNodeSpecies(nodeList[i]) << " " << famIndex << "|" << nodeList[i];

			vector <int> DescendantIds = RT->getLeafOrSpeciationDescendants(nodeList[i]);
			for(unsigned j = 0 ; j < DescendantIds.size(); j++)
			{
				if( RT->isRealLeaf(DescendantIds[j]) ) // this is a leaf
					ofs << " " << RT->getNodeName(DescendantIds[j]) ;
				else // thisis a speciation
					ofs << " " << famIndex << "|" << DescendantIds[j] ; 					
			}
			ofs << endl;
		}
	}

	nodeList.clear();

	//2. leaves
	nodeList = RT->getLeavesId();

	for(unsigned i = 0 ; i < nodeList.size(); i++)
	{
		if( RT->isRealLeaf(nodeList[i]) ) // to exclude losses
		{
			int spRPO = RT->getNodeSpecies(nodeList[i]);
			map <int,string>::iterator it = SpIdToSpName.find(spRPO);
			if( it ==  SpIdToSpName.end()) // should only occur if the map has not been constructed
			{
				ofs << spRPO;
				ofs << " " << RT->getNodeName(nodeList[i]) << endl;
			}
			else
			{
				ofs << it->second;
				ofs << " " << removePrefix(RT->getNodeName(nodeList[i]), it->second, true) << endl; // we remove the species name to avoid redundancy
			}
			
		}
			
	}


	ofs.close();
}

/*
Takes:
	- AdjTree * Atree : an adjacency tree
	- int NbSample : total number of sample done on this ECF
	- map <string, map <string, int> > & AdjIndexMap : map associating the string representation of both extremities of the adj to their index in the otehr maps
    - vector < double > & AdjInScoreList : list of inputed scores (unrepresented extant adjs and ancestral ones have an input score of 0)
    - vector < double > & AdjOutScoreList : frequency of observation of each adj in the output
    - vector < int > AdjSpeList : vector specifying the species association of an adjacency
    - ReconciledTree * Rtree1 : reconciled tree of the first family
	- ReconciledTree * Rtree2 : reconciled tree of the second family
*/
void UpdateAdjMapsFromTree( AdjTree * Atree, int NbSample, 
								map <string, map <string, int> > & AdjIndexMap,
            					vector < double > & AdjInScoreList,
            					vector < double > & AdjOutScoreList,
            					vector < int > & AdjSpeList,
            					ReconciledTree * Rtree1,
            					ReconciledTree * Rtree2)
{

	int fam1 = Atree->getGfamily1();
	int fam2 = Atree->getGfamily2();

	vector <int> nodeList = Atree->getNodesId();

	for(unsigned i = 0 ; i < nodeList.size(); i++)
	{
		int evt = Atree->getNodeEvent(nodeList[i]);
		if( ( Atree->isSpeciation(evt) ) || ( Atree->isExtant(evt) ) )
		{
			
			if(! Atree->hasNodeProperty(nodeList[i],nodeid1) ) //no nid property -> ignore
				continue;
			
			pair <int,int> gids = Atree->getNodeNodeIds(nodeList[i]);

			string name1;
			string name2;

			if( ( Atree->isExtant(evt) ) && (Rtree1 != NULL) )
				name1 = Rtree1->getNodeName( gids.first );
			else
				name1 = static_cast<ostringstream*>( &(ostringstream() << fam1) )->str() + "|" + static_cast<ostringstream*>( &(ostringstream() << gids.first) )->str();



			if( ( Atree->isExtant(evt) ) && (Rtree2 != NULL) )
				name2 = Rtree2->getNodeName( gids.second );
			else
				name2 = static_cast<ostringstream*>( &(ostringstream() << fam2) )->str() + "|" + static_cast<ostringstream*>( &(ostringstream() << gids.second) )->str();

			//cout << "found " << name1 <<" "<<name2<<endl;

			int adjIndex = -1;
			map <string, map <string, int> >::iterator it1;
			map <string, int>::iterator it2;

			it1 = AdjIndexMap.find(name1);
			if( it1 != AdjIndexMap.end() )
			{ // found the first name
				it2 = it1->second.find(name2); // search for the second name
	
				if( it2 != it1->second.end())
				{ // found the second name
					adjIndex = it2->second;
				}
			}
	
			//if(adjIndex == -1)
			//{ // reverse search: first Lnames2 then Lnames1
			//	it1 = AdjIndexMap.find(name2);
			//	if( it1 != AdjIndexMap.end() )
			//	{ // found the first name
			//		it2 = it1->second.find(name1); // search for the second name
			//
			//		if( it2 != it1->second.end())
			//		{ // found the second name
			//			adjIndex = it2->second;
			//		}
			//	}
			//}
	
			if(adjIndex == -1)
			{
				//the adj wasn't found -> add it
				adjIndex = AdjSpeList.size(); // new index

				AdjSpeList.push_back( Atree->getNodeSpecies(nodeList[i]) );
				AdjInScoreList.push_back( 0 ); // if we just discover it, it had no input score -> 0
				AdjOutScoreList.push_back( 0 );
				AdjIndexMap[name1][name2] = adjIndex;
				//cout << "created adj "<< adjIndex<<endl;
			}

			//we update the out score:
			AdjOutScoreList[adjIndex] +=  1 / (double) NbSample;

			//cout << nodeList[i] << " -> " << name1 <<" "<<name2<< " -> " << adjIndex << " -> " << AdjOutScoreList[adjIndex] << endl;

		}
	}


}

/*
Takes:
	- vector < AdjTree * > * AForest : an adjacency forest
	- int NbSample : total number of sample done on this ECF
	- map <string, map <string, int> > & AdjIndexMap : map associating the string representation of both extremities of the adj to their index in the otehr maps
    - vector < double > & AdjInScoreList : list of inputed scores (unrepresented extant adjs and ancestral ones have an input score of 0)
    - vector < double > & AdjOutScoreList : frequency of observation of each adj in the output
    - vector < int > AdjSpeList : vector specifying the species association of an adjacency
    - ReconciledTree * Rtree1 : reconciled tree of the first family
	- ReconciledTree * Rtree2 : reconciled tree of the second family
*/
void UpdateAdjMapsFromForest( vector < AdjTree * > * AForest, int NbSample, 
								map <string, map <string, int> > & AdjIndexMap,
            					vector < double > & AdjInScoreList,
            					vector < double > & AdjOutScoreList,
            					vector < int > & AdjSpeList,
            					ReconciledTree * Rtree1,
            					ReconciledTree * Rtree2
            					)
{
	for(unsigned i =0;i< AForest->size();i++)
	{
		UpdateAdjMapsFromTree( AForest->at(i),  NbSample, 
								AdjIndexMap,
            					AdjInScoreList,
            					AdjOutScoreList,
            					AdjSpeList,
            					Rtree1,
            					Rtree2);
	}
}


void AddAdjMapsToFile(string fileName , int sens1, int sens2,
								map <string, map <string, int> > & AdjIndexMap,
            					vector < double > & AdjInScoreList,
            					vector < double > & AdjOutScoreList,
            					vector < int > & AdjSpeList
								)
{

	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out | ofstream::app);

	string sens1Str = "";
	string sens2Str = "";

	if(sens1 == -1)
		sens1Str = " -";
	else if(sens1 == 1)
		sens1Str = " +";
	else if(sens2 != 0)
		sens1Str = " *";

	if(sens2 == -1)
		sens2Str = " +";
	else if(sens2 == 1)
		sens2Str = " -";
	else if(sens1 != 0)
		sens2Str = " *";

	for (map <string, map <string, int> >::iterator it=AdjIndexMap.begin(); it!=AdjIndexMap.end(); ++it)
	{
		string name1 = it->first;
		for (map <string, int>::iterator it2= it->second.begin(); it2!=it->second.end(); ++it2)
		{
			string name2 = it2->first;
			int adjIndex = it2->second;

			ofs << AdjSpeList[adjIndex] << " " << name1 << " " << name2 << sens1Str << sens2Str << " " << AdjInScoreList[adjIndex] << " " << AdjOutScoreList[adjIndex] << endl;
		}

	}



	ofs.close();	
}

/*  
Only does leaves. Doesn't check internal nodes names
*/
map <int,string> buildRPOtoSpNameMap(MySpeciesTree * Stree)
{
	map <int,string> RPOtoSpNameMap;

	vector <int> nodeList = Stree->getLeavesId();

	for(unsigned i = 0 ; i < nodeList.size();i++)
	{
		if(Stree->hasNodeName(nodeList[i]))
		{
			RPOtoSpNameMap[ Stree->getRPO( nodeList[i] ) ] = Stree->getNodeName( nodeList[i] );
			//cout << Stree->getRPO( nodeList[i] ) << " -> " << RPOtoSpNameMap[ Stree->getRPO( nodeList[i] ) ] << endl;
		}
	}
	return RPOtoSpNameMap;
}

/*
Takes:
	- string s : a string 
	- string prefix: a prefix
	- bool removeOneMore (default = false): if true, one more character will be remove a separator for instance)
Returns:
	(string) the string s without the prefix (if it is present)
*/
string removePrefix(string s, string prefix, bool removeOneMore)
{
	if( s.size() < prefix.size() ) // prefix bigger tha s -> prefix is not present
		return s;

	for(unsigned i = 0 ; i < prefix.size(); i++)//checking prefix presence
	{
		if(s[i]!= prefix[i])
			return s; // difference -> prefix absent
	}

	int removalPos = prefix.size();
	if(removeOneMore)
		removalPos++;

	return s.substr(removalPos);
}


void writeSpeciesCorrespondanceFile(string fileName, MySpeciesTree * speciesTree)
{
	ofstream ofs;
	ofs.open(fileName.c_str(),ofstream::out );

	vector <int> nodeList = speciesTree->getInnerNodesId();

	for(unsigned i = 0 ; i < nodeList.size(); i++)
	{
		if( speciesTree->getSonsId( nodeList[i] ).size() == 2 ) // 2 sons -> speciation node
		{
			vector< string > names = speciesTree->cloneSubtree(nodeList[i])->getLeavesNames(); // list of leaves under this node

			ofs << speciesTree->getRPO( nodeList[i] ) ;
			for(unsigned j = 0 ; j < names.size(); j++)
				ofs << " " << names[j];
			ofs << endl;
		}
	}

	nodeList = speciesTree->getLeavesId();
	for(unsigned i = 0 ; i < nodeList.size(); i++)
	{
		if(speciesTree->hasNodeName( nodeList[i] ))
			ofs << speciesTree->getRPO( nodeList[i] ) << " " << speciesTree->getNodeName( nodeList[i] )<<endl;
	}

}