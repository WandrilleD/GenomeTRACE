#ifndef MYCLADESANDTRIPARTITIONSHELPER_H_
#define MYCLADESANDTRIPARTITIONSHELPER_H_
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


/**

@file
@author Wandrille Duchemin
@created 08-Aug-2016
@modified 09-Aug-2016


**/

/*
This file contains various functions that count number of topologies and combinations
that are involved in the computation necessary to draw trees according to a flat CCP distribution (ie. a CCP distribution that would have been built from the set of all topologies with a fixed number of leaves) or a constrained CCP distribution (same as a flat CCP distribution, but where two given set of leaves must be separated by at least a branch)
*/


#define MYINFINITY numeric_limits<double>::max()
#define MYMINFINITY -numeric_limits<double>::max()


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <iostream>

using namespace std;


//functions for flat distribution. Each have a log complement that should be prefered to avoid overflow
void getNbRootedTopo(int n, map <int , long> &NbRootedTopo);
void getNbCombinations(int n , int k, map <int , long> &NbComb);

void getLogNbRootedTopo(int n, map <int , double> &LogNbRootedTopo);
void getLogNbCombinations(int n , int k, map <int , double> &LogNbComb);

vector<long>  getWeightedSplitSizeForCladeOfSize(int n , map <int , long> &NbComb, map <int , long> &NbRootedTopo );
vector<double> getLogWeightedSplitSizeForCladeOfSize(int n , map <int , double> &LogNbComb, map <int , double> &LogNbRootedTopo );

int randomChooseWeightedSplitSizeForCladeOfSize(int n , map <int , long> &NbComb, map <int , long> &NbRootedTopo );
int randomChooseLogWeightedSplitSizeForCladeOfSize(int n , map <int , double> &LogNbComb, map <int , double> &LogNbRootedTopo );


pair < vector <string> , vector <string> > splitLeafList(vector <string> leafList, int firstSonSize);

//functions for constrained CCP distribution

double getLogNbRootedTopoWithConstraint(int x, int y , int z, map <int , double> &LogNbRootedTopo, int rootConstrained = 0);

class LogConstrainedCombinationCounter
{
	public:
		//information parameters
		int k; //the size of the set
		int N; //total number of possible element (leaves)
		int X; //Number of element labeled x  
		int Y; //Number of element labeled y  (Where x and y are arbitrary labels. A leaf can't have two labels but might have none)
		int nbUnlabeled; // number of unlabeled elements in N
		//counts
		double LognbZeroLabel; // log number of combinations where there are no labeled element in the set
		vector <double> LognbXlabel; // log number of combination where there are elements labeled x in the set (and no element labeled y). That number is equal to index + 1
		vector <double> LognbYlabel; // log number of combination where there are elements labeled x in the set (and no element labeled x). That number is equal to index + 1

		int onlyDoFor;

		LogConstrainedCombinationCounter(int setSize, int Usize, int nbX, int nbY, map< int, map <int,double> > &LogNbCombMap, int onlyDo = 0);

		double getLognbComb(int nblabeledElement, bool Xlabel);

	protected:
		double ComputeLogNbComb(int nblabeledElement, bool Xlabel, map< int, map <int,double> > &LogNbCombMap);
};





void getLogCountVectorsForCladeOfSize(int N, int k , int X, int Y,  
										map< int, map <int,double> > &LogNbCombMap, map <int , double> &LogNbRootedTopo, 
										double &logCountZero, vector<double> &vectorLogCountXlabel , vector<double> &vectorLogCountYlabel, int onlyDoFor = 0 );


void randomChooseLogWeightedContrainedSplitSize(int N, int X, int Y,int onlyDoFor, map< int, map <int,double> > &LogNbCombMap, map <int , double> &LogNbRootedTopo,
												int &k , int &x , int &y);

#endif