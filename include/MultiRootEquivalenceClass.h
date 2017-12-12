#ifndef MULTIROOTADJ_CLASS_H_
#define MULTIROOTADJ_CLASS_H_

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

This file contains a class for adjacency classes

Created the: 17-04-2016
by: Wandrille Duchemin

Last modified the: 21-07-2016
by: Wandrille Duchemin

*/

#include "EquivalenceClass.h"
//#include "MultiRootAdjMatrix.h"


class MultiRootEquivalenceClass : public EquivalenceClass
{
//protected:
//	MultiRootAdjMatrix * Amat;

public:
	MultiRootEquivalenceClass(int fam1, int fam2) : EquivalenceClass(fam1, fam2)
	{//cout << "MultiRootEquivalenceClass creation" << endl;
	};

	MultiRootEquivalenceClass( EquivalenceClass *EC);

	~MultiRootEquivalenceClass()
	{ //cout << "plop MREC"<<endl;
	};


/*
	vector<EquivalenceClass *> refineEqClass(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool forceLTrefining, bool verbose);


	vector<EquivalenceClass *> refineEqClassWhole(ReconciledTree * Rtree1, ReconciledTree * Rtree2, bool verbose);
*/

	//void createAdjMatrix(double Gcost, double Bcost, ReconciledTree * rtree1, ReconciledTree * rtree2, bool VERBOSE, bool boltzmann , double temp, double absencePenalty );
};

#endif