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


#include "MyCladesAndTripartitionsHelper.h"

/*
recursive function that fills NbRootedTopo with the number of unrooted topologies from 1 to n leaves

Takes:
	- int n : max number of leaves in the rooted topology
	- map <int , long> &NbRootedTopo : map assciating a number of leaves to a number of rooted topologies
*/
void getNbRootedTopo(int n, map <int , long> &NbRootedTopo)
{
	if( NbRootedTopo.find(n) != NbRootedTopo.end()) // n already exists -> do not compute it again
		return;

	if(n ==1) // end of the recursion
	{
		NbRootedTopo[1] = 1;
		return;
	}

	if( NbRootedTopo.find(n-1) == NbRootedTopo.end()) // n-1 has not been found -> compute it
		getNbRootedTopo(n-1 , NbRootedTopo);

	NbRootedTopo[n] = ( 2 * n - 3 ) * NbRootedTopo[ n - 1 ] ; 
	return;
}

/*
Takes:
	- int n : number of elements to choose from
	- int k : number of elements to choose
	- map <int , long> &NbComb : map associating a number k to the number of combination of k element in a set of n

*/
void getNbCombinations(int n , int k, map <int , long> &NbComb)
{
	if( NbComb.find(k) != NbComb.end()) // n already exists -> do not compute it again
		return;

	if(k ==1) // end of the recursion
	{
		NbComb[1] = (long) n;
		return;
	}

	if( NbComb.find(k-1) == NbComb.end()) // n-1 has not been found -> compute it
		getNbCombinations( n, k-1 , NbComb);

	NbComb[k] =  (long)( ( (double) ( ( n - k + 1 ) * NbComb[ k - 1 ] ) ) / ( (double) k ) ) ;
	return;
}



/*
recursive function that fills NbRootedTopo with the number of unrooted topologies from 1 to n leaves

Takes:
	- int n : max number of leaves in the rooted topology
	- map <int , double> &LogNbRootedTopo : map assciating a number of leaves to the log of a number of rooted topologies
*/
void getLogNbRootedTopo(int n, map <int , double> &LogNbRootedTopo)
{
	if( LogNbRootedTopo.find(n) != LogNbRootedTopo.end()) // n already exists -> do not compute it again
		return;

	if(n ==1) // end of the recursion
	{
		LogNbRootedTopo[1] = 0;
		return;
	}

	if( LogNbRootedTopo.find(n-1) == LogNbRootedTopo.end()) // n-1 has not been found -> compute it
		getLogNbRootedTopo(n-1 , LogNbRootedTopo);

	LogNbRootedTopo[n] = log( (double)( 2 * n - 3 ) ) + LogNbRootedTopo[ n - 1 ] ; 
	return;
}

/*
Takes:
	- int n : number of elements to choose from
	- int k : number of elements to choose
	- map <int , double> &LogNbComb : map associating a number k to the lof of a number of combination of k element in a set of n

*/
void getLogNbCombinations(int n , int k, map <int , double> &LogNbComb)
{
	if( LogNbComb.find(k) != LogNbComb.end()) // n already exists -> do not compute it again
		return;

	if(k == 0)
	{
		LogNbComb[0] = 0; // log(1)
		return; // log(1)
	}
	else if(k ==1) // end of the recursion
	{
		LogNbComb[1] = log((double) n);
		return;
	}

	if( LogNbComb.find(k-1) == LogNbComb.end()) // n-1 has not been found -> compute it
		getLogNbCombinations( n, k-1 , LogNbComb);

	LogNbComb[k] =   log( (double) ( n - k + 1 ) ) +  LogNbComb[ k - 1 ] - log( (double) k ) ;
	return;
}



/*
Takes:
- int n : size of the clade to split
- map <int , long> &NbComb : map associating a number k to the number of combination of k element in a set of n
- map <int , long> &NbRootedTopo : map assciating a number of leaves to a number of rooted topologies

Returns:
	vector <int> a vector where the element at index k correspond to the number of topologies where a clade of size n is split in two clades of size k and n-k
				By convention, the first element is always 0 except if n == 1 (in which case it is 1 and this is the only element of the vector)
*/
vector<long> getWeightedSplitSizeForCladeOfSize(int n , map <int , long> &NbComb, map <int , long> &NbRootedTopo )
{

	vector <long> WeightedSplitSizes;

	if(n == 1)
		WeightedSplitSizes.push_back(1);

	WeightedSplitSizes.push_back(0);
	//updating maps if need be
	getNbCombinations(n , n/2, NbComb);
	getNbRootedTopo(n-1, NbRootedTopo);

	for(int k = 1 ; k <= (n/2) ; k++) 
	{
		// the clade of size n is split into two clades of size k and n-k
		long nbTopoForSplitSize = NbComb[k] * NbRootedTopo[k] * NbRootedTopo[n-k] ; 
		// the number of topo for this split size correspond to all possible combination of size k times the number of rooted topologies given by each half of the split

		if( k == (n/2))
		{
			if(n%2 == 0)
			{ // if n is pair, there is a special case when k == n/2
				nbTopoForSplitSize /= 2 ; // we account for the symmetry of the combinations
			}
		}
		WeightedSplitSizes.push_back( nbTopoForSplitSize );
	}
	return WeightedSplitSizes;
}


/*
Takes:
- int n : size of the clade to split
- map <int , double> &LogNbComb : map associating a number k to the log of the number of combination of k element in a set of n
- map <int , double> &LogNbRootedTopo : map assciating a number of leaves to the log of the number of rooted topologies

Returns:
	vector <double> a vector where the element at index k correspond to the log number of topologies where a clade of size n is split in two clades of size k and n-k
				By convention, the first element is always 0 
*/
vector<double> getLogWeightedSplitSizeForCladeOfSize(int n , map <int , double> &LogNbComb, map <int , double> &LogNbRootedTopo )
{

	vector <double> WeightedSplitSizes;

	if(n == 1)
		WeightedSplitSizes.push_back(0);

	WeightedSplitSizes.push_back(0);
	//updating maps if need be
	getLogNbCombinations(n , n/2, LogNbComb);
	getLogNbRootedTopo(n-1, LogNbRootedTopo);

	for(int k = 1 ; k <= (n/2) ; k++) 
	{
		// the clade of size n is split into two clades of size k and n-k
		double nbTopoForSplitSize = LogNbComb[k] + LogNbRootedTopo[k] + LogNbRootedTopo[n-k] ; 
		// the number of topo for this split size correspond to all possible combination of size k times the number of rooted topologies given by each half of the split

		if( k == (n/2))
		{
			if(n%2 == 0)
			{ // if n is pair, there is a special case when k == n/2
				nbTopoForSplitSize -= log(2) ; // we account for the symmetry of the combinations
			}
		}
		WeightedSplitSizes.push_back( nbTopoForSplitSize );
	}
	return WeightedSplitSizes;
}


/*
Takes:
- int n : size of the clade to split
- map <int , long> &NbComb : map associating a number k to the number of combination of k element in a set of n
- map <int , long> &NbRootedTopo : map assciating a number of leaves to a number of rooted topologies

Returns:
	int : chosen size k. 
			the weight used to choose k correspond to the number of topologies where a clade of size n is split in two clades of size k and n-k

*/
int randomChooseWeightedSplitSizeForCladeOfSize(int n , map <int , long> &NbComb, map <int , long> &NbRootedTopo )
{
	if( n< 4 )
		return 1;

	//updating maps if need be
	getNbRootedTopo(n, NbRootedTopo);

	long totNb = NbRootedTopo[n]; // there is totNb possible topologies under this clade

	long chosen = (long) ((double) totNb * ((double) rand() / RAND_MAX )); //random results

	cout << "totNb :" << totNb  << " chosen : " << chosen;

	int k = 1;
	getNbCombinations(n , k, NbComb);
	long nbTopoForSplitSize = NbComb[k] * NbRootedTopo[k] * NbRootedTopo[n-k] ; 
	// the number of topo for this split size correspond to all possible combination of size k times the number of rooted topologies given by each half of the split
	if( k == (n/2))
	{
		if(n%2 == 0)
		{ // if n is pair, there is a special case when k == n/2
			nbTopoForSplitSize /= 2 ; // we account for the symmetry of the combinations
		}
	}	

	chosen -= nbTopoForSplitSize;
	while( chosen >= 0 )
	{
		k += 1;

		getNbCombinations(n , k, NbComb);
		nbTopoForSplitSize = NbComb[k] * NbRootedTopo[k] * NbRootedTopo[n-k] ; 
		// the number of topo for this split size correspond to all possible combination of size k times the number of rooted topologies given by each half of the split
		if( k == (n/2))
		{
			if(n%2 == 0)
			{ // if n is pair, there is a special case when k == n/2
				nbTopoForSplitSize /= 2 ; // we account for the symmetry of the combinations
			}
		}	

		chosen -= nbTopoForSplitSize;
	}

	cout << " k : "<<k << endl;
	return k;
}

/*
Takes:
- int n : size of the clade to split
- map <int , double> &LogNbComb : map associating a number k to the log number of combination of k element in a set of n
- map <int , double> &LogNbRootedTopo : map assciating a number of leaves to a log number of rooted topologies

Returns:
	int : chosen size k. 
			the weight used to choose k correspond to the number of topologies where a clade of size n is split in two clades of size k and n-k

*/
int randomChooseLogWeightedSplitSizeForCladeOfSize(int n , map <int , double> &LogNbComb, map <int , double> &LogNbRootedTopo )
{
	if( n< 4 )
		return 1;

	//updating maps if need be
	getLogNbRootedTopo(n, LogNbRootedTopo);

	double logtotNb = LogNbRootedTopo[n]; // there is totNb possible topologies under this clade

	double chosen = (double) rand() / RAND_MAX ; //random results

	//cout << "totNb :" << logtotNb  << " chosen : " << chosen;

	int k = 1;
	getLogNbCombinations(n , k, LogNbComb);
	double ProbaTopoForSplitSize = LogNbComb[k] + LogNbRootedTopo[k] + LogNbRootedTopo[n-k]  - logtotNb; 
	// the number of topo for this split size correspond to all possible combination of size k times the number of rooted topologies given by each half of the split
	if( k == (n/2))
	{
		if(n%2 == 0)
		{ // if n is pair, there is a special case when k == n/2
			ProbaTopoForSplitSize -= log(2) ; // we account for the symmetry of the combinations
		}
	}	
	//cout << chosen << "<->"<< exp(ProbaTopoForSplitSize) << endl;
	
	chosen -= exp(ProbaTopoForSplitSize);
	while( chosen >= 0 )
	{

		k += 1;

		getLogNbCombinations(n , k, LogNbComb);
		ProbaTopoForSplitSize = LogNbComb[k] + LogNbRootedTopo[k] + LogNbRootedTopo[n-k] - logtotNb; 
		// the number of topo for this split size correspond to all possible combination of size k times the number of rooted topologies given by each half of the split
		if( k == (n/2))
		{
			if(n%2 == 0)
			{ // if n is pair, there is a special case when k == n/2
				ProbaTopoForSplitSize -= log(2) ; // we account for the symmetry of the combinations
			}
		}	
		//cout << chosen << "<->"<< exp(ProbaTopoForSplitSize) << endl;
		chosen -= exp(ProbaTopoForSplitSize);
	}

	//cout << " k : "<<k << endl;
	return k;	
}


/*
randomly split the strings contained into leafList between two vector of size firstChildSize and ( leafList.size() - firstChildSize )

Takes:
	- vector <string> leafList
	- int firstChildSize

Returns:
	pair < vector <string> , vector <string> >
*/
pair < vector <string> , vector <string> > splitLeafList(vector <string> leafList, int firstChildSize)
{
	pair < vector <string> , vector <string> > Children;
	double n = (double) leafList.size();

	for(unsigned i = 0 ; i < leafList.size(); i++)
	{
		double r = ( (double) rand() / RAND_MAX ) * n;
		if( r <= firstChildSize)
		{
			Children.first.push_back(leafList[i]);
			firstChildSize -- ;
		}
		else
		{
			Children.second.push_back(leafList[i]);
		}
		n --;
	}
	return Children;
}



// LogConstrainedCombinationCounter functions


/*
Takes:
	- int x : number of leaf with label X
	- int y : number of leaves with label Y
	- int z : number of unlabeled leaves
	- map <int , double> &LogNbRootedTopo : map assciating a number of leaves to the log of the number of rooted topologies
	- int rootConstrained [default = 0] : 	if 0 : the root can be anywhere in the tree
											if 1 : the root must be in the X labeled subtree
											if 2 : the root must be in the y labeled subtree

Returns:
	(double) : log of the number of rooted topologies such that in their unrooted version there exists clades containing only unlabeled and all X labeled leaves (and indenticaly for Y labeled leaves)
*/
double getLogNbRootedTopoWithConstraint(int x, int y , int z, map <int , double> &LogNbRootedTopo, int rootConstrained )
{
	double result;

	if(z < 0)
	{
		cout << " !!ERROR!! : z<0 " << x<<","<<y<<","<<z <<endl;
		throw 1;
	}

	if( (x == 0) || ( y == 0 ) ) // no constraint here
	{
		getLogNbRootedTopo( x+y+z, LogNbRootedTopo);
		return LogNbRootedTopo[ x+y+z ];
	}
		

	if(z == 0) // no unlabeled leaf
	{ // the result is obtained by using:
	  //           Tx the number of rooted trees with x labeled leaves
	  //		   Ty the number of rooted trees with y labeled leaves
	  //        These two set of trees are joined by their root, giving Ty * Ty unrooted trees.
	  //        IF the root is not constrained, then:
	  //        		Each of these can be rooted in (2 * (x+y) - 3 ) ways  (ie. the nb of branches in a tree with x+y leaves )
	  //        OTHERWISE:
	  //                There is ( 2 * (a+1) - 3 ) possible root (a is equal to x if rootConstrained == 1 and a == y rootConstrained == 2 )
	  //				This is equilvalent to the number of branches + 1 in the subtree we have to root in.
		getLogNbRootedTopo( x, LogNbRootedTopo);
		getLogNbRootedTopo( y, LogNbRootedTopo);

		double nbTopo = LogNbRootedTopo[x] + LogNbRootedTopo[y];
		if(rootConstrained == 0)
			nbTopo += log( (double)( 2 * (x+y) - 3 ) ) ;
		else
		{
			int a = x;
			if(rootConstrained == 2)
				a = y;
			nbTopo += log( (double)( 2 * (a+1) - 3 ) ) ;
		}
		return  nbTopo;
	}

	// as unlabeled leaves can be placed anywhere in the tree
	// the number of trees with z unlabeled leaves is equal to ( 2*(x+y+z) - 3 ) * NbRootedTopoWithConstraint(x,y,z-1)
	return getLogNbRootedTopoWithConstraint(x,y,z-1, LogNbRootedTopo, rootConstrained) + log( (double)( 2 * (x+y+z) - 3 ) );
}


/*
Takes:
	- int setSize : size of the resulting set
	- int Usize : size of the urn
	- int nbX : number of element labeled x in the urn
	- int nbY : number of element labeled y in the urn
	- map< int, map <int,double> > &LogNbCombMap : map associtating:	a number n
																		to a map associating a number k to the log number of combination of k element in a set of n
	- int OnlyDo [default = 0] : if 0 : all counts are made
								 if 1 : the counts implying >0 y labeled leaves aren't made (and corrected for)
								 if 2 : the counts implying >0 x labeled leaves aren't made (and corrected for)


*/
LogConstrainedCombinationCounter::LogConstrainedCombinationCounter(int setSize, int Usize, int nbX, int nbY, map< int, map <int,double> > &LogNbCombMap, int onlyDo )
{
	//information parameters
	k = setSize; //the size of the set
	N = Usize; //total number of possible element (leaves)
	X = nbX; //Number of element labeled x  
	Y = nbY; //Number of element labeled y  (Where x and y are arbitrary labels. A leaf can't have two labels but might have none)
	nbUnlabeled = N - X - Y; // number of unlabeled elements in N
	onlyDoFor = onlyDo; 
	//filling counts
	

	LognbZeroLabel = MYMINFINITY;
	if( nbUnlabeled >= k )
	{
		getLogNbCombinations(nbUnlabeled , k, LogNbCombMap[k] );
		LognbZeroLabel = LogNbCombMap[nbUnlabeled][k];
	}


	bool Xlabel = true;
	int limit = min( k, X );

	if(onlyDoFor != 2)
	{
		for(int nblabeledElement = 1 ; nblabeledElement <= limit; nblabeledElement++ )
			LognbXlabel.push_back( ComputeLogNbComb(nblabeledElement , Xlabel, LogNbCombMap) );
	}

	if(onlyDoFor != 1)
	{
		Xlabel = false;
		limit = min( k, Y );
	
		for(int nblabeledElement = 1 ; nblabeledElement <= limit; nblabeledElement++ )
			LognbYlabel.push_back( ComputeLogNbComb(nblabeledElement , Xlabel, LogNbCombMap) );
	}

}

/*
Takes:
	- int nblabeledElement : number of labeled elemnt in the set of size k
	- bool Xlabel : whether the labeled elements are of label X or Y
	- map< int, map <int,double> > &LogNbCombMap : map associtating:	a number n
																		to a map associating a number k to the log number of combination of k element in a set of n

Returns:
	(double) : log of the number of combinations with nblabeledElement of label X or Y (and 0 of the other label), !!CORRECTED to avoid redundancy !!

*/
double LogConstrainedCombinationCounter::ComputeLogNbComb(int nblabeledElement, bool Xlabel, map< int, map <int,double> > &LogNbCombMap)
{


	int totLabeled = X;
	if(!Xlabel)
		totLabeled = Y;

	if( ( nblabeledElement > totLabeled ) || ( ( k - nblabeledElement ) > nbUnlabeled ) ) // not enough labeled or unlabeled element to satisfy
		return MYMINFINITY;

	// filling LogNbCombMap if necessary
	if(LogNbCombMap.find(totLabeled) == LogNbCombMap.end() )
		getLogNbCombinations(totLabeled , nblabeledElement, LogNbCombMap[totLabeled]);

	if(LogNbCombMap.find(nbUnlabeled) == LogNbCombMap.end() )
		getLogNbCombinations(nbUnlabeled , k - nblabeledElement, LogNbCombMap[nbUnlabeled]);	

	double result = LogNbCombMap[totLabeled][nblabeledElement] + LogNbCombMap[nbUnlabeled][nbUnlabeled]; 

	//if( ( k == n/2) && (n%2 == 0) )// account for the case where k == N/2 
	//{
		//-> divide the nb by 2 if the complementary set is non Mixing to avoid redundancy
	if(onlyDoFor == 0)
	{
		if( nblabeledElement == totLabeled )// the compl set is non mixing IF all labeled element of one type are on that set
			result -= log(2);
	}
	//}

	return result; 
}

/*
Takes:
	- int nblabeledElement : number of labeled elemnt in the set of size k
	- bool Xlabel : whether the labeled elements are of label X or Y
*/
double LogConstrainedCombinationCounter::getLognbComb(int nblabeledElement, bool Xlabel)
{
	if(nblabeledElement ==  0)
		return LognbZeroLabel;

	vector <double > *Vptr;

	if(Xlabel)
		Vptr = &LognbXlabel;
	else
		Vptr = &LognbYlabel;

	if(Vptr->size() < nblabeledElement)
		return MYMINFINITY;
	return Vptr->at(nblabeledElement-1);
}


/*
(non mixing clade == clade whose labeled elements all share the same label)

Takes:
	- int N :  total number of element 
	- int k :  desired number of element in the non-mixing clade 
	- int X :  total number of element labeled X
	- int Y :  total number of element labeled Y
	- map< int, map <int,double> > &LogNbCombMap : map associtating:	a number n
																		to a map associating a number k to the log number of combination of k element in a set of n
	- map <int , double> &LogNbRootedTopo : map assciating a number of leaves to the log of the number of rooted topologies

	- double &logCountZero : log of the number of topology with a non-mixing clade of size k with only unlabeled leaves
	- vector<double> &vectorLogCountXlabel : vector were the element at index i correspond to the log number of topology with a non-mixing clade of size k with i+1 leaves labeled X (and k-(i+1) unlabeled leaves )
	- vector<double> &vectorLogCountYlabel : vector were the element at index i correspond to the log number of topology with a non-mixing clade of size k with i+1 leaves labeled Y (and k-(i+1) unlabeled leaves )

	- int onlyDoFor = 0 :	if 0 : count will be computed for all possible number of x or y labeled leaves
							if 1 : count will only be computed for x labeled leaves clades (we presume that a sister clade contains x labeled leaves)
							if 2 : count will only be computed for y labeled leaves clades (we presume that a sister clade contains y labeled leaves)
*/
void getLogCountVectorsForCladeOfSize(int N, int k , int X, int Y,  
										map< int, map <int,double> > &LogNbCombMap, map <int , double> &LogNbRootedTopo, 
										double &logCountZero, vector<double> &vectorLogCountXlabel , vector<double> &vectorLogCountYlabel, int onlyDoFor )
{
	cout << "k" <<k << "N" <<N << "X" <<X << "Y" <<Y << endl;
	//1. counts all combhinations of size k in N element where at least a clade doesn't mix leaves labeled X and y (called a non mixing clade)
	LogConstrainedCombinationCounter CombinationCounter(k, N, X, Y, LogNbCombMap, onlyDoFor);

	int Z = N - X - Y;

	//2. Case where the non-mixing clade only contains unlabeled leaves
	logCountZero = CombinationCounter.LognbZeroLabel;
	
	if( logCountZero != MYMINFINITY ) // if the count is not 0 -> we multiply by the number of topology with k unlabeled leaves and the number of contrained topologies with parameters X,Y and z-k 
	{
		getLogNbRootedTopo(k,LogNbRootedTopo);
		logCountZero +=  LogNbRootedTopo[k] + getLogNbRootedTopoWithConstraint( X,Y , Z - k , LogNbRootedTopo, onlyDoFor);
	}

	// NB : the root is always constrained if the non mixing clade contains at least a labeled leaf

	if(onlyDoFor != 2)
	{
		//3. Cases where' the non-mixing clade contains a number of X-labeled leaves
		for(int i = 0 ; i < CombinationCounter.LognbXlabel.size(); i++)
		{ // i+1 X labeled leaves
			int nbLabeledXFirst = i+1;
			int nbLabeledYFirst = 0;
			int nbUnLabeledFirst = k - nbLabeledXFirst;
	
			int nbLabeledXSecond = X - nbLabeledXFirst;
			int nbLabeledYSecond = Y - nbLabeledYFirst;
			int nbUnLabeledSecond = Z - nbUnLabeledFirst;
	
			vectorLogCountXlabel.push_back( CombinationCounter.LognbXlabel[i] ) ;
			cout << "x : " << nbLabeledXFirst << " -> log combinations: " << vectorLogCountXlabel.back() << endl;
			if( vectorLogCountXlabel.back() != MYMINFINITY ) // non zero count
			{
				vectorLogCountXlabel.back() += getLogNbRootedTopoWithConstraint( nbLabeledXFirst,nbLabeledYFirst , nbUnLabeledFirst , LogNbRootedTopo, 1) +
											   getLogNbRootedTopoWithConstraint( nbLabeledXSecond,nbLabeledYSecond , nbUnLabeledSecond , LogNbRootedTopo, 1); // we constrain the root to the subtree containing x labeled leaves as its sister clade contains at least 1 x labeled leaf
			}
		}
	}

	if(onlyDoFor != 1)
	{
		//4. Cases where' the non-mixing clade contains a number of Y-labeled leaves
		for(int i = 0 ; i < CombinationCounter.LognbYlabel.size(); i++)
		{ // i+1 Y labeled leaves
			int nbLabeledXFirst = 0;
			int nbLabeledYFirst = i+1;
			int nbUnLabeledFirst = k - nbLabeledYFirst;
	
			int nbLabeledXSecond = X - nbLabeledXFirst;
			int nbLabeledYSecond = Y - nbLabeledYFirst;
			int nbUnLabeledSecond = Z - nbUnLabeledFirst;
	
			vectorLogCountYlabel.push_back( CombinationCounter.LognbYlabel[i] ) ;
			cout << "y : " << nbLabeledYFirst << " -> log combinations: " << vectorLogCountYlabel.back() << endl;
			if( vectorLogCountYlabel.back() != MYMINFINITY ) // non zero count
			{
				vectorLogCountYlabel.back() += getLogNbRootedTopoWithConstraint( nbLabeledXFirst,nbLabeledYFirst , nbUnLabeledFirst , LogNbRootedTopo, 2) +
											   getLogNbRootedTopoWithConstraint( nbLabeledXSecond,nbLabeledYSecond , nbUnLabeledSecond , LogNbRootedTopo, 2); // we constrain the root to the subtree containing y labeled leaves as its sister clade contains at least 1 y labeled leaf
			}
		}
	}
}

/*
Takes:
	- int N : total number of element 
	- int X : total number of element labeled X
	- int Y : total number of element labeled Y
	- int onlyDoFor : if 0 : count will be computed for all possible number of x or y labeled leaves
					  if 1 : count will only be computed for x labeled leaves clades (we presume that a sister clade contains x labeled leaves)
					  if 2 : count will only be computed for y labeled leaves clades (we presume that a sister clade contains y labeled leaves)

	- map< int, map <int,double> > &LogNbCombMap : map associtating:	a number n
																		to a map associating a number k to the log number of combination of k element in a set of n
	- map <int , double> &LogNbRootedTopo : map assciating a number of leaves to the log of the number of rooted topologies

	- int &k : size of the first clade of the chosen split
	- int &x : number of x labeled element in the first clade of the chosen split
	- int &y : number of y labeled element in the first clade of the chosen split

These 5 last arguments are filled during the function course.
The 3 last arguments will contain the result of the function: the information about the chosen split of a clade of size N with a constraint modelled as two set of leaves (respectively X and Y leaves).
*/
void randomChooseLogWeightedContrainedSplitSize(int N, int X, int Y, int onlyDoFor, map< int, map <int,double> > &LogNbCombMap, map <int , double> &LogNbRootedTopo,
												int &k , int &x , int &y)
{

	int Z = N - X - Y;

	double LogTotalNumberTopologies = getLogNbRootedTopoWithConstraint( X,Y , Z , LogNbRootedTopo, onlyDoFor); // getting the total

	double chosen = (double) rand() / RAND_MAX ; //random number between 0 and 1

	cout << "totNb :" << LogTotalNumberTopologies << " ( "<< exp( LogTotalNumberTopologies ) <<" ) onlyDoFor: "<< onlyDoFor <<". chosen : " << chosen << endl;

	int current_k = 1;

	double logCountZero;
	vector<double> vectorLogCountXlabel;
	vector<double> vectorLogCountYlabel;

	while(chosen > 0)
	{
		cout << "chosen "<< chosen << endl;
		cout << "getting count for k "<< current_k << endl;
		getLogCountVectorsForCladeOfSize( N, current_k, X, Y, LogNbCombMap, LogNbRootedTopo, logCountZero, vectorLogCountXlabel , vectorLogCountYlabel , onlyDoFor);	

		cout << "  all unlabeled log count : ";
		if(logCountZero == MYMINFINITY)
		{
			cout << "-inf"<<endl;
		}
		else
		{ // non zero count -> test it
			cout << logCountZero << " ( "<< exp(logCountZero) <<" ) -> " << exp( logCountZero - LogTotalNumberTopologies ) << endl;

			chosen -= exp( logCountZero - LogTotalNumberTopologies );
			if(chosen <= 0)
			{
				cout << "  --> choice ok"<<endl;
				k=current_k;
				x=0;
				y=0;
				return;
			}
		}

		for(unsigned nbx =1 ; nbx < ( vectorLogCountXlabel.size() + 1 ); nbx++)
		{

			cout << "  " << nbx << " label x log count : ";
			if(vectorLogCountXlabel[nbx -1 ] == MYMINFINITY)
			{
				cout << "-inf"<<endl;
			}
			else
			{ // non zero count -> test it
				cout << vectorLogCountXlabel[nbx -1 ] << " ( "<< exp(vectorLogCountXlabel[nbx -1 ]) <<" ) -> " << exp( vectorLogCountXlabel[nbx -1 ] - LogTotalNumberTopologies ) << endl;
	
				chosen -= exp( vectorLogCountXlabel[nbx - 1 ] - LogTotalNumberTopologies );
				if(chosen <= 0)
				{
					cout << "  --> choice ok"<<endl;
					k=current_k;
					x=nbx;
					y=0;
					return;
				}
			}
		}

		for(unsigned nby =1 ; nby < ( vectorLogCountYlabel.size() + 1 ) ; nby++)
		{

			cout << "  " << nby << " label y log count : ";
			if(vectorLogCountYlabel[nby] == MYMINFINITY)
			{
				cout << "-inf"<<endl;
			}
			else
			{ // non zero count -> test it
				cout << vectorLogCountYlabel[nby -1 ] << " ( "<< exp(vectorLogCountYlabel[nby -1 ]) <<" ) -> " << exp( vectorLogCountYlabel[nby -1 ] - LogTotalNumberTopologies ) << endl;
	
				chosen -= exp( vectorLogCountYlabel[nby - 1 ] - LogTotalNumberTopologies );
				if(chosen <= 0)
				{
					cout << "  --> choice ok"<<endl;
					k=current_k;
					x=0;
					y=nby;
					return;
				}
			}
		}

		vectorLogCountXlabel.clear();
		vectorLogCountYlabel.clear();
		current_k+=1;
	}

	cout << " not found any split !? ( chosen : " << chosen << " )" << endl;	
	k=-1;
	x=-1;
	y=-1;
	return;
}