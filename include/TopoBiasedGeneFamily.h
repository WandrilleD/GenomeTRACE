#ifndef GENE_FAMILY_BIASED_TOPO_H_
#define GENE_FAMILY_BIASED_TOPO_H_

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

Created the: 15-06-2016
by: Wandrille Duchemin

Last modified the: 16-06-2016
by: Wandrille Duchemin

*/

#include "GeneFamily.h"
#include "AdjTree.h"

class TopoBiasedGeneFamily : public GeneFamily
{
protected:
	bool BiasedCCPDistribSet;
	MyCladesAndTripartitions * BiasedCCPDistrib;

	vector < vector <string> > getLeafListToKeepFromAtreeAlone(AdjTree * Atreeint , int GfamIndex);
	vector < vector <string> > getLeafListToKeepFromAtree(AdjTree * Atreeint , int GfamIndex);

	vector < vector <string> > getLeafListToKeepFromAForest( vector < AdjTree * > * AdjForest , vector<int> GfamIndex, vector <bool> isCompatible);

public:
		//constructors and destructors
	TopoBiasedGeneFamily(bool verbose = false, bool superverbose=false) : GeneFamily(verbose,superverbose)
	{
		BiasedCCPDistribSet = false;
	}

	~TopoBiasedGeneFamily()
	{
		if(BiasedCCPDistribSet)
			delete BiasedCCPDistrib;
	}

	//with a bunch of trees
	TopoBiasedGeneFamily(vector<MyGeneTree*> &geneTrees, char charsep,bool verbose, bool superverbose) : GeneFamily( geneTrees,  charsep, verbose, superverbose)
	{
		BiasedCCPDistribSet = false;
	}
	
	//ale file reader
	TopoBiasedGeneFamily(string aleFileName, char charsep,bool verbose, bool superverbose) : GeneFamily( aleFileName, charsep, verbose, superverbose)
	{
		BiasedCCPDistribSet = false;
	}

	//the biased CCP distrib is set from the ADJForest C1s
	MyCladesAndTripartitions * setBias( vector < AdjTree * > * AdjForest , vector<int> GfamIndex, vector <bool> isCompatible);

	void unsetBias()
	{
		if(BiasedCCPDistribSet)
		{
			delete BiasedCCPDistrib;
			BiasedCCPDistribSet=false;
		}
	}



};

#endif
