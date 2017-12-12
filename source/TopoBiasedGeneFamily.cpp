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

This file contains a class for Gene Family as a Clade and Tripartition instance associated to a reconciled tree
that can use a biased Clade and Tripartition instance to sample tree topologies

Created the: 16-06-2016
by: Wandrille Duchemin

Last modified the: 16-06-2016
by: Wandrille Duchemin

*/

#include "TopoBiasedGeneFamily.h"

/*
For the case where ids in the adj tree does not correspond to ids in the reconciled tree

Takes:
	- AdjTree * Atree : an adjacency tree
	- int GfamIndex : wether we are interested in the ids of gfam1 or gfam2 in the adjacency tree
Returns;
	vector < vector <string> > : vector of clades represented by list of strings
*/
vector < vector <string> > TopoBiasedGeneFamily::getLeafListToKeepFromAtreeAlone(AdjTree * Atree, int GfamIndex)
{
	vector < vector <string> > LeafListList;

	Atree->getC1LeafList( GfamIndex, LeafListList );

	return LeafListList;
}

/*
For the case where ids in the adj tree does not correspond to ids in the reconciled tree

Takes:
	- AdjTree * Atree : an adjacency tree
	- int GfamIndex : wether we are interested in the ids of gfam1 or gfam2 in the adjacency tree
Returns;
	vector < vector <string> > : vector of clades represented by list of strings

*/
vector < vector <string> > TopoBiasedGeneFamily::getLeafListToKeepFromAtree(AdjTree * Atree, int GfamIndex)
{
	vector <int> nodeToKeepIds = Atree->getC1NodeIds(GfamIndex);

	map <int, vector <string> > cladeToLeaves = RecTree.cladeToLeafListAssociation();

	map <int, bool > cladeToLeavesTokeep; 
	for(unsigned i = 0 ; i < nodeToKeepIds.size() ; i++) // we filte the clades we are intrested in
	{
		int cladeId = RecTree.getNodeCladeNum(nodeToKeepIds[i]);

		if(cladeToLeavesTokeep.find(cladeId) == cladeToLeavesTokeep.end()) // new clade -> add it
		{
			if(cladeToLeaves[cladeId].size() > 1) // we aren't interested in leaf clades as they, by dfinition, always exists
			{
				cladeToLeavesTokeep[cladeId] = true;
			}
		}
	}

	vector < vector <string> > leafList;
	for (map<int,bool>::iterator it=cladeToLeavesTokeep.begin(); it!=cladeToLeavesTokeep.end(); ++it)
	{
		leafList.push_back(cladeToLeaves[it->first]);
	}

	return leafList;
}


/*
Takes:
	- vector < AdjTree * > * AdjForest : a vector of adjacency trees
	- vector<int> GfamIndex : wether the current gene family is gfam 1 (0) or gfam2 (1) in each AdjTree
	- vector <bool> isCompatible : wether the current gene family and each AdjTree are compatible (ie. ids are identical)

Returns:
	vector < vector <string> > : vector of clades represented by list of strings; WITHOUT ANY GUARANTY THAT TEH LIST IS NOT REDUNDANT
*/
vector < vector <string> > TopoBiasedGeneFamily::getLeafListToKeepFromAForest( vector < AdjTree * > * AdjForest , vector<int> GfamIndex, vector <bool> isCompatible)
{
	vector < vector <string> > leafList;

	for(unsigned i = 0 ; i < AdjForest->size() ; i++)
	{
		vector < vector <string> > tmp;
		if(isCompatible[i])
			tmp = getLeafListToKeepFromAtree(AdjForest->at(i), GfamIndex[i]);
		else
			tmp = getLeafListToKeepFromAtreeAlone(AdjForest->at(i), GfamIndex[i]);

		for(unsigned j = 0 ; j < tmp.size(); j++)
			leafList.push_back(tmp[j]);
	}

	return leafList;
}


/*
Takes:
	-

Returns;

*/
void TopoBiasedGeneFamily::setBias( vector < AdjTree * > * AdjForest , vector<int> GfamIndex, vector <bool> isCompatible)
{

	vector < vector <string> > leafListToRestrict = getLeafListToKeepFromAForest( AdjForest , GfamIndex, isCompatible);

	BiasedCCPDistrib = new MyCladesAndTripartitions(CCPDistrib);
	
	BiasedCCPDistrib->RestrictToClade( leafListToRestrict );

	BiasedCCPDistribSet = true;
}
