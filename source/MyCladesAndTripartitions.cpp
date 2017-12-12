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
@created 06-Nov-2015
@modified 19-Aug-2016


**/

#include <stdlib.h>
#include <time.h>
#include <boost/foreach.hpp>

#include "MyCladesAndTripartitions.h"




/*
Here, a constraint is represented by 2 non-overlapping set of leaves (whose union does not necessarily form the set of all leaves in the CCP distribution).
A tree stisfying the constraint must contain a branch 

Takes:
	- vector <string> tested : a clade as a list of leaf labels
	- map<string, int> constraint : keys are leaf labels, values are 1 or 2 according to their appartenance to the 1st or the second leafList forming the constrain

Returns:
	(int) : -1 if tested does contain leaves from both set of the constraint
			0 if tested does not contain leaves from the sets of constraint
			1 if tested does not contain leaves from first set of constraint
			2 if tested does not contain leaves from second set of constraint
*/
int MyCladesAndTripartitions::isNotMixing(vector <string> tested, map<string, int> constraint)
{
	int n = tested.size();
	int associatedSet = 0; //0: not attributed; 1,2 -> attributed
	map<string,int>::iterator found;
	for(unsigned i = 0 ; i < n; i++)
	{
		found = constraint.find(tested[i]);
		if(found != constraint.end())
		{
			if( associatedSet == 0 )  // nothing associated yet -> associate
				associatedSet = found->second;
			else if( associatedSet != found->second) // mixing
				return -1;
		}
	}
	return associatedSet;
}


/*
Here, a constraint is represented by 2 non-overlapping set of leaves (whose union does not necessarily form the set of all leaves in the CCP distribution).
A tree stisfying the constraint must contain a branch 

Takes:
	- pair<int,int> split: a pair of clades
	- map<string, int> constraint : keys are leaf labels, values are 1 or 2 according to their appartenance to the 1st or the second leafList forming the constrain

Returns:
	(bool) : true if both clades mix leaves from the sets forming the constraint
*/
bool MyCladesAndTripartitions::isIncompatible(pair<int,int> split, map<string, int> constraint)
{
	//not the way to go I think
	bool Mix1 = (isNotMixing(mClades.getLeafListFromCladeInt( split.first ), constraint) == -1);
	bool Mix2 = (isNotMixing(mClades.getLeafListFromCladeInt( split.second ), constraint) == -1);

	return (Mix1 && Mix2) ;
}


/*
Here, a constraint is represented by 2 non-overlapping set of leaves (whose union does not necessarily form the set of all leaves in the CCP distribution).
A tree stisfying the constraint must contain a branch 

Takes:
	- vector <string> FlatLeafList : list of the leaves forming the current clade, used when drawing in the flat distribution
	- pair< vector<string> , vector<string> > constraint : two list of leaves that will be used to define aithorized and non authorized clades
												clades are authorized only if they:
																	- contain entirely one of the list in the constraint
																	- contain leaf from only one constraint
	- int onlyDoFor [default = 0]:	if 0 : all counts are made
								 	if 1 : only splits who have a clade that does not contain any leaves from the second set of the constraint will be authorized
								 	if 2 : only splits who have a clade that does not contain any leaves from the second set of the constraint will be authorized
								 	This property is used to denote that a sister clade already contains leaves from one or another set in the constraint (ie. the choice of the relatioship of the set's ancestor ( parent/child or sister ) have already been made)
	

Returns:
	pair < pair < vector<string>, vector<string> > , pair <int,int> > : 
													first: pair < vector<string>, vector<string> > : pair of list of leaves desribing the two clade of the chosen split
													second: pair <int,int> : int describing the composition of each clade of the split
																			 if -1: the clade mixes leaves from both set of the constraint
																			 if 0 : the clade does not contain any leaves from any set of the constraint
																		 	 if 1 : the clade contains at least one leaf from the first set of the constraint (and none from the second one)
								 											 if 2 : the clade contains at least one leaf from the second set of the constraint (and none from the first one)
								 	

*/
pair < pair < vector<string>, vector<string> > , pair <int,int> > MyCladesAndTripartitions::DrawInFlatConstrainedDistribution(vector <string> FlatLeafList,  pair< vector<string> , vector<string> > constraint,int onlyDoFor )
{
	pair < pair < vector<string>, vector<string> > , pair <int,int> > result;

	//1. we want to count the number of leaves that belong to either set of the constraint
	//   we also use the occasion to split leaves into three vector that we will use to create the actual split once the number of each category has been decided
	int N = 0;
	int X = 0;
	int Y = 0;

	vector <string> UnlabeledLeaves;
	vector <string> XlabeledLeaves;
	vector <string> YlabeledLeaves;

	for(unsigned i = 0 ; i < FlatLeafList.size(); i++)
	{
		N++;
		bool found = false;
		for(unsigned j = 0 ; j < constraint.first.size(); j++)
		{ // searching for the leaf in the first set
			if( FlatLeafList[i] == constraint.first[j] )
			{
				X++;
				found = true;
				XlabeledLeaves.push_back(FlatLeafList[i]);
				break;
			}
		}

		if(!found)
		{
			for(unsigned j = 0 ; j < constraint.second.size(); j++)
			{// searching for the leaf in the second set
				if( FlatLeafList[i] == constraint.second[j] )
				{
					Y++;
					found = true;
					YlabeledLeaves.push_back(FlatLeafList[i]);
					break;
				}
			}
		}

		if(!found)
			UnlabeledLeaves.push_back(FlatLeafList[i]);
	}

	// 2. we draw randomly:
	//    k -> the size the of the first child of the split
	//    x -> the number of x labeled leaves in the first child of the split
	//    y -> the number of y labeled leaves in the first child of the split
	//    NB: x*y == 0 (they can't be different from 0 at the same time) (both can be equal to 0 at the same time)

	int k,x,y;

    map< int, map <int,double> > LogNbCombMap;
    map <int , double> LogNbRootedTopo;
    
    randomChooseLogWeightedContrainedSplitSize( N,  X,  Y, onlyDoFor ,LogNbCombMap, LogNbRootedTopo, k , x , y);

    int z = k - x - y; // nb of unlabeled leaves

    // 3. from these numbers we draw the actual split
    result.first.first.reserve(k);
    result.first.second.reserve(N - k);

    //unlabeled
    pair < vector <string> , vector <string> > tmpsplit =  splitLeafList(UnlabeledLeaves, z);
    for(unsigned i = 0 ; i < tmpsplit.first.size(); i++)
    	result.first.first.push_back(tmpsplit.first[i]);
    for(unsigned i = 0 ; i < tmpsplit.second.size(); i++)
    	result.first.second.push_back(tmpsplit.second[i]);

    //X labeled
    tmpsplit =  splitLeafList(XlabeledLeaves, x);
    for(unsigned i = 0 ; i < tmpsplit.first.size(); i++)
    	result.first.first.push_back(tmpsplit.first[i]);
    for(unsigned i = 0 ; i < tmpsplit.second.size(); i++)
    	result.first.second.push_back(tmpsplit.second[i]);

    //Y labeled
    tmpsplit =  splitLeafList(YlabeledLeaves, y);
    for(unsigned i = 0 ; i < tmpsplit.first.size(); i++)
    	result.first.first.push_back(tmpsplit.first[i]);
    for(unsigned i = 0 ; i < tmpsplit.second.size(); i++)
    	result.first.second.push_back(tmpsplit.second[i]);


    //4. from the numbers we deduce the mixing property of the clades

	result.second.first = 0;
	if( x*y  > 0)
		result.second.first = -1;
    else if( x > 0 ) 
    	result.second.first = 1;
    else if( y > 0 )
    	result.second.first = 2;

    int x2 = X - x;
    int y2 = Y - y;

	result.second.second = 0;
	if( x2*y2  > 0)
		result.second.second = -1;
    else if( x2 > 0 ) 
    	result.second.second = 1;
    else if( y2 > 0 )
    	result.second.second = 2;


    return result;
}

/*
Here, a constraint is represented by 2 non-overlapping set of leaves (whose union does not necessarily form the set of all leaves in the CCP distribution).
A tree stisfying the constraint must contain a branch 

Takes:
	- int idU : current clade to split
	- pair< vector<string> , vector<string> > constraint : two list of leaves that will be used to define aithorized and non authorized clades
												clades are authorized only if they:
																	- contain entirely one of the list in the constraint
																	- contain leaf from only one constraint
	- map<int,int> &MixingMap : keys are cladeId, values indicates wether they mix leaves from the sets of the constraint or not, and if they do not, to which set they belong (if they belong to a set)
	- int onlyDoFor [default = 0]:	if 0 : all counts are made
								 	if 1 : only splits who have a clade that does not contain any leaves from the second set of the constraint will be authorized
								 	if 2 : only splits who have a clade that does not contain any leaves from the second set of the constraint will be authorized
								 	This property is used to denote that a sister clade already contains leaves from one or another set in the constraint (ie. the choice of the relatioship of the set's ancestor ( parent/child or sister ) have already been made)
	

Returns:
	(pair <int,int>) : chosen split ; -1,-1 if no observed split remained after all non-authorized were screened.
*/
pair <int,int> MyCladesAndTripartitions::getRandomCladeConstrained(int idU, pair< vector<string> , vector<string> > constraint, map<int,int> &MixingMap, int onlyDoFor)
{
	//cout << "MyCladesAndTripartitions::getRandomCladeConstrained "<< idU<< endl;
	//cout << "constraint :"<< endl;
	//for(unsigned i = 0 ; i < constraint.first.size(); i++)
	//	cout << constraint.first[i] << " - ";
	//cout << endl;

	//for(unsigned i = 0 ; i < constraint.second.size(); i++)
	//	cout << constraint.second[i] << " - ";
	//cout << endl;

	if( ( constraint.first.size() == 0 ) || ( constraint.second.size() == 0 ) )
		return getRandomClade(idU);

	map <string, int> constraintMap;
	int n = constraint.first.size();
	int set = 1;
	for(unsigned i = 0 ; i < n ; i++)
		constraintMap[constraint.first[i]] = set;

	n = constraint.second.size();
	set = 2;
	for(unsigned i = 0 ; i < n ; i++)
		constraintMap[constraint.second[i]] = set;


	
	vector <int> compatibleSplit;
	vector <double> compatibleSplitRatio;
	double compatibleRatioSum = 0;

	pair<int,int> currentSplit;
	map<int,int>::iterator mixIterator;

	int Mix1;
	int Mix2;

	int nbSplit = getSplitCount(idU);
	int counter;
	for(counter=0;counter<nbSplit;counter++)//for each split
	{

		currentSplit = getCladeSplit(idU,counter);

		if(currentSplit.first == -1)
			return currentSplit; // bypass for leaves

		mixIterator = MixingMap.find(currentSplit.first);
		if(mixIterator == MixingMap.end())
		{


			Mix1 = isNotMixing( mClades.getLeafListFromCladeInt( currentSplit.first ), constraintMap );
			MixingMap[currentSplit.first] = Mix1;

			//vector <string> leaves = mClades.getLeafListFromCladeInt( currentSplit.first );
			//cout << "mix pattern for clade "<<  currentSplit.first << " ( ";
			//for(unsigned i = 0 ; i< leaves.size(); i++)
			//	cout << leaves[i] << " ";
			//cout << ") ->" << Mix1<<endl;
		}
		else
			Mix1 = mixIterator->second;

		mixIterator = MixingMap.find(currentSplit.second);
		if(mixIterator == MixingMap.end())
		{
			Mix2 = isNotMixing( mClades.getLeafListFromCladeInt( currentSplit.second ), constraintMap );
			MixingMap[currentSplit.second] = Mix2;

			//vector <string> leaves = mClades.getLeafListFromCladeInt( currentSplit.second );
			//cout << "mix pattern for clade "<<  currentSplit.second << " ( ";
			//for(unsigned i = 0 ; i< leaves.size(); i++)
			//	cout << leaves[i] << " ";
			//cout << ") ->" << Mix2<<endl;
		}
		else
			Mix2 = mixIterator->second;


		bool accept = false;

		if( ( Mix1 == 0 ) || ( Mix2 == 0 ) ) // if at least one is neutral -> compatible
			accept = true;
		else if(onlyDoFor != 0)
		{
			if( ( Mix1 == onlyDoFor ) || ( Mix2 == onlyDoFor ) ) // given that none is neutral, we accept only if one is non mixing AND equal to onlyDoFor
				accept = true;
		}
		else
			accept = ( ( Mix1 != -1 ) || ( Mix2 != -1 ) ) ; // no constraint -> accept if at least one is non mixing

		//cout << currentSplit.first << "-" << currentSplit.second << " -> "<< Mix1 <<" "<<Mix2<< "( " << onlyDoFor << " ) -> "<< accept<<endl;

		if(accept)// compatible
		{
			compatibleSplit.push_back(counter);
			compatibleSplitRatio.push_back(getSplitRatio(idU,counter));
			compatibleRatioSum += getSplitRatio(idU,counter);
		}
	}


	if(compatibleRatioSum != 0) // some split are compatible
	{
		double r = ((double) rand()/ RAND_MAX) * compatibleRatioSum; //random draw between 0 and compatibleRatioSum
		
		for(counter = 0 ; counter < compatibleSplit.size(); counter++)
		{
			r -= compatibleSplitRatio[counter];
			if(r <= 0)
				return getCladeSplit(idU,compatibleSplit[counter]);
		}
		return getCladeSplit(idU,compatibleSplit[counter]);
	}

	// no split was compatible
	cout << "NO COMPATIBLE SPLIT FOR CLADE "<<idU<<endl;


	return pair<int,int> (-1,-1);

}


/**
Here, a constraint is represented by 2 non-overlapping set of leaves (whose union does not necessarily form the set of all leaves in the CCP distribution).
A tree stisfying the constraint must contain a branch 

Takes:
	- int &pOrd : current depth-first post-order number
	- int idU : current clade in the recursion
	- pair< vector<string> , vector<string> > constraint : two list of leaves that will be used to define aithorized and non authorized clades
												clades are authorized only if they:
																	- contain entirely one of the list in the constraint
																	- contain leaf from only one constraint
	- double probaFlat : probability to prefer a flat distribution to the one described by the CCPs
	- vector <string> FlatLeafList : list of the leaves forming the current clade, used when drawing in the flat distribution
	- map<int,int> &MixingMap : keys are cladeId, values indicates wether they mix leaves from the sets of the constraint or not, and if they do not, to which set they belong (if they belong to a set)
	- int onlyDoFor [default = 0] : if 0 : all counts are made
								 	if 1 : only splits who have a clade that does not contain any leaves from the second set of the constraint will be authorized
								 	if 2 : only splits who have a clade that does not contain any leaves from the second set of the constraint will be authorized
								 	This property is used to denote that a sister clade already contains leaves from one or another set in the constraint (ie. the choice of the relatioship of the set's ancestor ( parent/child or sister ) have already been made)


Returns:
	MyGeneNode * 
**/
MyGeneNode * MyCladesAndTripartitions::getRandomTreeConstrainedAux(int &pOrd,int idU, pair< vector<string> , vector<string> > constraint, double probaFlat, vector <string> FlatLeafList, map<int,int> &MixingMap, int onlyDoFor )
{
	//cout << "MyCladesAndTripartitions::getRandomTreeConstrainedAux "<< idU << " " << onlyDoFor << endl; 
    MyGeneNode *node = new MyGeneNode();

    bool UseConstrainedFlat = false;
    bool isLeaf = false;


	pair<int,int> cladeSplit;
	pair < vector <string> , vector <string> > ChildrenleafList;
	pair <int,int> ChildrenMixInfo;

    if(idU == -1)
    { 
    // special case where we are in an unknown clades
    //  -> we use the constrained flat distribution
    //  -> FlatLeafList contains the leaves forming the clade
		UseConstrainedFlat = true;
    }
    // we ignore probaFlat which will only be used if we revert to unconstrained distribution (when the two set of leaves have been correctly separated for instance)


	if(!UseConstrainedFlat)
	{
		cladeSplit = getRandomCladeConstrained(idU, constraint,MixingMap, onlyDoFor);

		//cout << "chosen cladeSplit " << cladeSplit.first << " " << cladeSplit.second  << endl;

		if( (cladeSplit.first == -1 ) && ( cladeSplit.second == -1 ) ) // we did not find any compatible split -> use constrained flat distribution
		{
			UseConstrainedFlat = true;
			FlatLeafList.clear();
			FlatLeafList = mClades.getLeafListFromCladeInt( idU );
	    	if(FlatLeafList.size() <= 2) // shouldn't happen here... but you never know
	    	{ // this is a leaf of a cherry -> no need to randomly choose split
	    		UseConstrainedFlat = false;
	    	}
	    }
	    else
	    {
	    	ChildrenMixInfo.first = MixingMap[ cladeSplit.first ];
	    	ChildrenMixInfo.second = MixingMap[ cladeSplit.second ];
	    }
	}

	cout << "  UseConstrainedFlat "<< UseConstrainedFlat << endl; 

	if(FlatLeafList.size() == 1)
		isLeaf = true;
	else if(UseConstrainedFlat)
	{ // we did not find any satisfying split -> rely on a constrained flat distribution (a distribution containing all the trees staisfying the constraint)
		pair < pair < vector<string>, vector<string> > , pair <int,int> > randomSplitData = DrawInFlatConstrainedDistribution(FlatLeafList,
																								constraint, onlyDoFor);
		//fill split data
		ChildrenleafList.first = randomSplitData.first.first;
		ChildrenleafList.second = randomSplitData.first.second;
		
		cladeSplit.first = mClades.getCladeIntFromLeafList(ChildrenleafList.first);
		cladeSplit.second = mClades.getCladeIntFromLeafList(ChildrenleafList.second);

		ChildrenMixInfo.first = randomSplitData.second.first;
		ChildrenMixInfo.second = randomSplitData.second.second;

	}


	if(isLeaf)
	{
        node->setName( mClades.getLeafName( idU ) );
	}
	else
	{ // we treat the split

		//1. we update onlyDoFor
		if(onlyDoFor == 0)
		{
			if( ChildrenMixInfo.first >0 ) // decided 
				onlyDoFor = ChildrenMixInfo.first;
			else if( ChildrenMixInfo.second >0 ) // decided
				onlyDoFor = ChildrenMixInfo.second;
		}
		//2. we create the sons
		//cout << ChildrenMixInfo.first<<" "<< ChildrenMixInfo.second<< "<->" << "onlyDoFor "<< onlyDoFor<<endl;

		if( ChildrenMixInfo.first != -1 ) //we continue to use the constraint only if this clade is mixing. Otherwise the consttraint is fulfilled and no supplementary computation is needed
			node->addSon( getRandomTreeAux( pOrd, cladeSplit.first  , probaFlat, ChildrenleafList.first) ); // unconstraind
		else
			node->addSon( getRandomTreeConstrainedAux(pOrd,cladeSplit.first, constraint,  probaFlat, ChildrenleafList.first, MixingMap, onlyDoFor ) );//constrained
		
		if( ChildrenMixInfo.second != -1 ) //we continue to use the constraint only if this clade is mixing. Otherwise the consttraint is fulfilled and no supplementary computation is needed
			node->addSon( getRandomTreeAux( pOrd, cladeSplit.second  , probaFlat, ChildrenleafList.second) ); // unconstrained
		else
			node->addSon( getRandomTreeConstrainedAux(pOrd,cladeSplit.second, constraint,  probaFlat, ChildrenleafList.second, MixingMap, onlyDoFor ) );// constrained
		
        node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, bpp::Number<double>(pOrd));
	}


    if(idU != -1)
    {
	    if( mBranchLengths[idU] != -1 )
    	    node->setDistanceToFather( (mBranchLengths[idU] + 1)/getCladeOccurrences(idU) );
    }

	pOrd++;
	return node;
}

/**
Takes:
	- int idU : id of a clade

Returns:
	pair<int,int> a split of that clade drawn randomly (according to its split ratio)
**/
pair<int,int> MyCladesAndTripartitions::getRandomClade(int idU)
{

	

	double r = ((double) rand()/ RAND_MAX); //random results
	double ratio;//to store the various ratios


	int nbSplit = getSplitCount(idU);
	cout << "MyCladesAndTripartitions::getRandomClade "<< idU << " nb of split :" << nbSplit << endl;
	if(nbSplit == 0)
		return pair<int,int>(-1,-1);

	int counter;
	for(counter=0;counter<nbSplit;counter++)//for each split
	{
		ratio = getSplitRatio(idU,counter);
		if( r <= ratio)
			return getCladeSplit(idU,counter);//the chosen clade is found and returned
		r = r - ratio;//not the chosen clade, we update r
	}
	return getCladeSplit(idU,nbSplit - 1);//return the last split

}

/**
Takes:
	- int idU : id of a clade

Returns:
	pair<int,int> a split of that clade that has the maximum split ratio
**/
pair<int,int> MyCladesAndTripartitions::getMaxClade(int idU)
{

	//First, we find the list of all split of clade idU with the maximum ratio
	double maxRatio = 0;
	vector<int> maxIndex;//vector to keep the indexes of the split with the max ratio

	double ratio;

	int nbSplit = getSplitCount(idU);

	//cout << idU << ":" << nbSplit << endl;
	int counter;
	for(counter=0;counter<nbSplit;counter++)//for each split
	{
		ratio = getSplitRatio(idU,counter);
		if(ratio > maxRatio)
		{
			maxRatio = ratio;
			maxIndex.clear();
			maxIndex.push_back(counter);
		}
		else if(ratio == maxRatio)
		{
			maxIndex.push_back(counter);
		}
	}


	
	if(maxIndex.size() == 1)//There is only one split with the maximum ratio: return it
	{
		return getCladeSplit(idU,maxIndex[0]);
	}

	//There is more than one split: return a random one
	int r = rand() % maxIndex.size();
	return getCladeSplit(idU,maxIndex[r]);

}


/**
Takes:
	- int &pOrd : current depth-first post-order number
	- int idU : current clade in the recursion
	- double probaFlat : probability to prefer a flat distribution to the one described by the CCPs
	- vector <string> FlatLeafList : list of the leaves forming the current clade, used when drawing in the flat distribution

Returns:
	MyGeneNode * 
**/
MyGeneNode *MyCladesAndTripartitions::getRandomTreeAux(int &pOrd,int idU, double probaFlat, vector <string> FlatLeafList )
{
    MyGeneNode *node = new MyGeneNode();

    bool UseFlat = false;
    bool isLeaf = false;

	pair<int,int> cladeSplit;
	pair < vector <string> , vector <string> > ChildrenleafList;

	cout << "MyCladesAndTripartitions::getRandomTreeAux "<< idU << " " << FlatLeafList.size()<< endl;

    if(idU == -1)
    { 
    // special case where we are in an unknown clades
    //  -> we use the flat distribution
    //  -> FlatLeafList contains the leaves forming the clade
		UseFlat = true;
    }
    else if(probaFlat > 0) // we randomly choose if we go in the flat distribution or not
    {
    	UseFlat = ( ((double) rand() / RAND_MAX ) < probaFlat );

    	if(UseFlat)
    	{
	    	FlatLeafList.clear();
	    	FlatLeafList = mClades.getLeafListFromCladeInt( idU );
	    	if(FlatLeafList.size() <= 2)
	    	{ // this is a leaf of a cherry -> no need to randomly choose split
	    		UseFlat = false;
	    	}
	    }
    }

    if(!UseFlat)
    { // we use the 

	    cladeSplit = getRandomClade( idU );

	    if( cladeSplit.first == -1 ) 
	    {
	    	FlatLeafList.clear();
	    	FlatLeafList = mClades.getLeafListFromCladeInt( idU );
	    	if(FlatLeafList.size() == 1)
	    	{ // this is a leaf of a cherry -> no need to randomly choose split
	    		isLeaf = true;
	    	}
	    	else
	    	{
	    		UseFlat = true; // some error occured so we use flat distrib to get out of it.
	    	}
	    	
	    }
	}


    if(UseFlat) //NB 
    {
    	// useflat procedure
    	cout << "USEFLAT"<< endl;

    	//1. we randomly choose the size of the son clade according to the number of possible topolgy implied by said size
    	map <int , double> LogNbComb;
    	map <int , double> LogNbRootedTopo;
    	int sizeOfChild = randomChooseLogWeightedSplitSizeForCladeOfSize(FlatLeafList.size() , LogNbComb, LogNbRootedTopo );

    	//2. split the leaf vector into two
		ChildrenleafList = splitLeafList(FlatLeafList, sizeOfChild);


		//3. check if the children correspond to known clades
		cladeSplit.first = mClades.getCladeIntFromLeafList(ChildrenleafList.first);
		cladeSplit.second = mClades.getCladeIntFromLeafList(ChildrenleafList.second);

		//cout << "USEFLAT "<< FlatLeafList.size() <<" -> " << sizeOfChild << endl;
		//cout<<"["<<FlatLeafList[0];
        //for(unsigned j = 1 ; j < FlatLeafList.size(); j++)
        //    cout << "," << FlatLeafList[j] ;
        //cout<<"]";
		//cout << " -> ";
        //cout<<"(["<<ChildrenleafList.first[0];
        //for(unsigned j = 1 ; j < ChildrenleafList.first.size(); j++)
        //    cout << "," << ChildrenleafList.first[j] ;
        //cout<<"],["<<ChildrenleafList.second[0];
        //for(unsigned j = 1 ; j < ChildrenleafList.second.size(); j++)
        //    cout << "," << ChildrenleafList.second[j] ;
        //cout << "])"<<endl;
		//cout << idU << " -> " << cladeSplit.first << "," << cladeSplit.second << endl;
    }

	cout << idU << " -> " << cladeSplit.first << "," << cladeSplit.second << endl;

	if(!isLeaf)
	{
        node->addSon( getRandomTreeAux( pOrd, cladeSplit.first  , probaFlat, ChildrenleafList.first) );
        node->addSon( getRandomTreeAux( pOrd, cladeSplit.second  , probaFlat, ChildrenleafList.second) ) ;
        node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, bpp::Number<double>(pOrd));
    } 
    else
    {
        node->setName( mClades.getLeafName( idU ) );
    }

    if(idU != -1)
    {
	    if( mBranchLengths[idU] != -1 )
    	    node->setDistanceToFather( (mBranchLengths[idU] + 1)/getCladeOccurrences(idU) );
    }


    pOrd++;

    return node;
}


/**
Takes:
	- int &pOrd : current depth-first post-order number
	- int idU : current clade in the recursion

Returns:
	MyGeneNode * 
**/
MyGeneNode * MyCladesAndTripartitions::getMaxTreeAux(int &pOrd,int idU ) 
{
    MyGeneNode *node = new MyGeneNode();
    pair<int,int> cladeSplit = getMaxClade( idU);
    //cout << idU << "->" << cladeSplit.first << "," << cladeSplit.second << endl;
    if( cladeSplit.first != -1 ) {
        node->addSon( getMaxTreeAux( pOrd, cladeSplit.first ) );
        node->addSon( getMaxTreeAux( pOrd, cladeSplit.second ) );
        node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, bpp::Number<double>(pOrd));
    } else 
        node->setName( mClades.getLeafName( idU ) );

    if( mBranchLengths[idU] != -1 )
        node->setDistanceToFather( (mBranchLengths[idU] + 1)/getCladeOccurrences(idU) );

    pOrd++;

    return node;

}

/**
RECURSIVE

Takes:
	- MyGeneNode * node: a tree to compute the likelihood from
	- double default_proba: proba to give to a split that is absent from the MyCladesAndTripartitions counts
	
Returns:
	pair<double,vector<string>>:
		double: likelihood of the node (and underlying subtree) according to the CCPs distribution.
		vector<string>: list of all the leaves in the subtree rooted at node
**/
pair<double,vector<string> > MyCladesAndTripartitions::getTreeLikelihoodAux(MyGeneNode *node, double default_proba)
{

	std::pair<double, vector<string> > toreturn;

	toreturn.first = 1;//set base proba to 1.

	if (node->isLeaf())
	{
		toreturn.second.push_back(node->getName());
		return toreturn;
	}
	
	int nbson,i, cladeId;

	nbson = node->getNumberOfSons();

	if(nbson ==1)//only one son -> this is a intermediary node of some sort 
		return getTreeLikelihoodAux(node->getSon(0),default_proba);// return the value for its only son
	else if(nbson > 2)//more than 2 sons: error
		throw bpp::Exception("MyCladesAndTripartitions::getTreeLikelihoodAux: does not work when the tree is polytomic!");

	//We are now in the caser where there is more than one son

	std::pair<double, vector<string> > tempresult;//to put the result from each son

	std::pair<int,int> split;//split

	//first son
	tempresult = getTreeLikelihoodAux(node->getSon(0),default_proba);//getting the result for this son

	toreturn.first *= tempresult.first;//updating proba

	BOOST_FOREACH( string leafName, tempresult.second)
		toreturn.second.push_back(leafName);//updating leaf list

	split.first = mClades.getCladeIntFromLeafList(tempresult.second);

	tempresult.second.clear();//clearing 

	//second son
	tempresult = getTreeLikelihoodAux(node->getSon(1),default_proba);//getting the result for this son

	toreturn.first *= tempresult.first;//updating proba
	BOOST_FOREACH(string leafName, tempresult.second)
		toreturn.second.push_back(leafName);//updating leaf list

	split.second = mClades.getCladeIntFromLeafList(tempresult.second);

	tempresult.second.clear();//clearing 

	//now getting current cladeId
	cladeId = mClades.getCladeIntFromLeafList(toreturn.second);
	
	if(cladeId == -1 )
		toreturn.first *=default_proba;
	else if(split.first != -1 || split.second != -1)//if the son clades are known

		for(i=0 ; i < getSplitCount(cladeId) ; i++)
		{
			pair <int,int> testsplit = getCladeSplit(cladeId,i);
			if(split.first == testsplit.first || split.second == testsplit.first)//found the correct split
			{

				if(cladeId != 0 || mRooted)// special condition for the root clade: should not use split ratio if the distribution is unrooted
					toreturn.first *= getSplitRatio(cladeId,i);
				else if(cladeId == 0)
					toreturn.first *= ((double) getCladeOccurrences(split.first))/ ((double) mGeneTreeCount); //  if we are the "root" clade but the distribution is unrooted, we should multiply by the number of time we saw the child clade (clade probability)
				break;
			}
		}
		if(i == getSplitCount(cladeId))//the split was not found
			toreturn.first *=default_proba;

	return toreturn;
}


double MyCladesAndTripartitions::getMaxLikelihoodAux(int idU, map <int,double> &LkhMap)
{	

	if( mClades.isLeaf(idU) ) // stopping condition
		return 1;

	map <int,double>::iterator it = LkhMap.find(idU);
	if( it != LkhMap.end() ) // found
		return it->second;

	
	double maxRatio = 0;
	double ratio;

	int nbSplit = getSplitCount(idU);

	

	int counter;
	for(counter=0;counter<nbSplit;counter++)//for each split
	{
		pair<int,int> cladeSplit = getCladeSplit(idU,counter);

		ratio = 1;
		if(idU != 0 || mRooted)// special condition for the root clade: should not use split ratio if the distribution is unrooted
			ratio *= getSplitRatio(idU,counter);
		else if(idU == 0)
			ratio *= ((double) getCladeOccurrences(cladeSplit.first))/ ((double) mGeneTreeCount); //  if we are the "root" clade but the distribution is unrooted, we should multiply by the number of time we saw the child clade (clade probability)


		ratio *= getMaxLikelihoodAux(cladeSplit.first, LkhMap) * getMaxLikelihoodAux(cladeSplit.second, LkhMap);

		if(ratio > maxRatio)
		{
			maxRatio = ratio;
		}
	}

	LkhMap[idU] = maxRatio;

	return maxRatio;

}


// Public functions



/**
 * Constructor - Compute all clades and tripartitions for the given 
 * gene trees. 
 */
MyCladesAndTripartitions::MyCladesAndTripartitions( char charSep, vector<MyGeneTree*> &geneTrees, bool verbose, bool &overflow, string &errStr, bool polytomy ) : CladesAndTripartitions( charSep, geneTrees,verbose,overflow, errStr,polytomy,NULL)
{srand (time(NULL));//init the random seed
}


/**
 * Constructor - Ale reader
 * 
 */
MyCladesAndTripartitions::MyCladesAndTripartitions( char charSep, string aleFileName, string &errStr, bool verbose ) : CladesAndTripartitions( charSep, aleFileName, errStr, verbose )
{srand (time(NULL));//init the random seed
}


/*
Takes:
	- pair< vector<string> , vector<string> > constraint : two list of leaves that will be used to define aithorized and non authorized clades
												clades are authorized only if they:
																	- contain entirely one of the list in the constraint
																	- contain leaf from only one constraint
	- double probaFlat : probability to prefer a flat distribution to the one described by the CCPs

Returns:
	(MyGeneTree *)
*/
MyGeneTree * MyCladesAndTripartitions::getRandomTreeConstrained(pair< vector<string> , vector<string> > constraint, double probaFlat)
{
	int pOrd = 0;
	int onlyDoFor = 0;
	vector <string> emptyVec;
	map<int,int> MixingMap;


	cout << "MyCladesAndTripartitions::getRandomTreeConstrained" << endl;

	MyGeneNode *node = getRandomTreeConstrainedAux(pOrd, mClades.getRootClade(), constraint, probaFlat, emptyVec, MixingMap , onlyDoFor );

	return new MyGeneTree( *node );
}

/**
Returns a tree where bipartitions have been drawn at random according to their split ratios
**/
MyGeneTree * MyCladesAndTripartitions::getRandomTree(double probaFlat)
{
    	int pOrd = 0;
    	vector <string> emptyVector;
    	MyGeneNode *node = getRandomTreeAux( pOrd, mClades.getRootClade() , probaFlat, emptyVector);
    	return new MyGeneTree( *node );
}

/**
Returns a tree where bipartition have the maximum split ratios for the encountered clades
**/
MyGeneTree * MyCladesAndTripartitions::getMaxTree()
{
	int pOrd = 0;
	MyGeneNode *node;

	if(mClades.getLeafCount() == 2)
	{
	    node = new MyGeneNode();
    
		vector<int> Scl = mClades.getSortedClades();
		for(unsigned i = 0; i < Scl.size(); i++)
		{
			if(mClades.isLeaf(Scl[i]))
			{
				MyGeneNode *sonNode = new MyGeneNode();
				sonNode->setName( mClades.getLeafName( Scl[i] ) );
				node->addSon( sonNode);
				node->setBranchProperty(bpp::TreeTools::BOOTSTRAP, bpp::Number<double>(pOrd));
			}		
		}
    }
    else
    {
		node = getMaxTreeAux( pOrd, mClades.getRootClade() );
    }
	

    return new MyGeneTree( *node );
}

/**
RECURSIVE

Takes:
	- MyGeneNode * node: a tree to compute the likelihood from
	- double default_proba: proba to give to a split that is absent from the MyCladesAndTripartitions counts

Returns:
	double: likelihood of the node (and underlying subtree) according to the CCPs distribution.

NB: requires the tree to be binary
**/
double MyCladesAndTripartitions::getTreeLikelihood(MyGeneNode *node, double default_proba)
{
	pair<double,vector<string> > result;
	result = getTreeLikelihoodAux(node,default_proba);//recursive function
	return result.first;
}

/*
Takes:
 - node (MyGeneNode *): node of a tree

Returns:
 (bool): true if the leafs of the subtree rooted at node all corresponds to existing clades in the MyCladesAndTripartitions instance
 */
bool MyCladesAndTripartitions::isCompatible(MyGeneNode * node)
{
	if(node->isLeaf())
	{
		vector <string> leafNames;
		leafNames.push_back(node->getName());

		try
		{
			mClades.getCladeIntFromLeafList(leafNames); // will throw an error if the leaf clade does not exists
		}
		catch(exception & e) // not very elegant.
		{
			return false;
		}
	}

	bool result= true;

	//recursion
	int nbson = node->getNumberOfSons();

	for(int i = 0; i < nbson; i++)
	{
		if(!isCompatible(node->getSon(i)))
		{
			result = false;
			break;
		}
	}
	return result;
}

double MyCladesAndTripartitions::getMaxLikelihood()
{
	map <int,double> LkhMap;
	return getMaxLikelihoodAux(mClades.getRootClade() , LkhMap);
}



int MyCladesAndTripartitions::RestrictToClade(vector <int> cladesToForce)
{
	//1. from clade to force to clade to delete
	// the clades to delete are the clades incompatible with the ones we want to force
	vector< vector<int> > noncompatibleLists = mClades.PUBLICgetNoncompatible();	
	map <int,bool> cladeToDeleteMAP; // we first put the clades to delete in a map to avoid doubles
	//cout << "restricting to clades" ;
	for(unsigned i = 0 ; i < cladesToForce.size(); i++)
	{
		int idU = cladesToForce[i];
		//cout << " " << idU ;
		for(unsigned j = 0 ; j < noncompatibleLists[idU].size(); j++)
		{
			cladeToDeleteMAP[ noncompatibleLists[idU][j] ] = true;
		}
	}

	for(unsigned i = 0 ; i < cladesToForce.size(); i++)// remove the clades to restrict from the clade to delete list to if they appear in the map
	{
		int idU = cladesToForce[i];
		cladeToDeleteMAP[ idU ] = false;
	}

	//cout << endl;
	vector <int> cladeToDelete;
	for (map<int,bool>::iterator it=cladeToDeleteMAP.begin(); it!=cladeToDeleteMAP.end(); ++it) 
	{
		if(it->second) // that way we don't remove a clade to restrict
	  		cladeToDelete.push_back( it->first );
  	}
	//cout << "MyCladesAndTripartitions::RestrictToClade. Deleting "<< cladeToDelete.size() << " clades"<<endl;
	//for(unsigned i = 0 ; i < cladeToDelete.size(); i++)
	//	cout << cladeToDelete[i] << " ";
	//cout << endl;
  	//2. delete the clades to delete
	return deleteClades(cladeToDelete);
}


int MyCladesAndTripartitions::RestrictToClade(vector <vector <string> > cladesToForce)
{
	vector <int> cladesToForceIDU;

	for( unsigned i = 0 ; i < cladesToForce.size() ; i++)
	{
		int idU = mClades.getCladeIntFromLeafList( cladesToForce[i] );
		if(idU >= 0 )
			cladesToForceIDU.push_back( idU );
	
		//cout << idU <<" <- ";
		//for(unsigned j = 0 ; j < cladesToForce[i].size() ; j++)
		//	cout << cladesToForce[i][j] << " ";	
		//cout << endl;
	}
	return RestrictToClade(cladesToForceIDU);
}


/*
Takes:
	- MyGeneTree * Tree : a gene Tree

Returns:
	(bool) : false if a problem occured during tree addition

*/
bool MyCladesAndTripartitions::addTreeToDistribution(MyGeneTree * Tree)
{
	bool polytomic = false;
	bool overflow = false;

	//cout << " size of mSplitsOccurrencesMap : "<< mSplitsOccurrencesMap.size() << endl;
	mSplitsOccurrencesMap.clear(); // should already be cleared, better make sure.


    computeSplits( Tree->getRootNode(), polytomic, overflow );
    if( overflow )
        return false;
	//cout << " size of mSplitsOccurrencesMap after computeSplits: "<< mSplitsOccurrencesMap.size() << endl;

	boost::unordered_map< vector <int>,double>::iterator iterSplits;
	//for(iterSplits = mSplitsOccurrencesMap.begin(); iterSplits != mSplitsOccurrencesMap.end(); iterSplits++ )
    //{
    //    vector<int> split = iterSplits->first;
    //    cout << split[0] << "->" << split[1]<<","<<split[2] ;
    //    cout << " : "<< iterSplits->second << endl;
    //}

    //adding the content of mSplitsOccurrencesMap to the counts
	for(iterSplits = mSplitsOccurrencesMap.begin(); iterSplits != mSplitsOccurrencesMap.end(); iterSplits++ )
    {
        vector<int> split = iterSplits->first;
        int cladeNum = split[0];

        pair<int,int> CladeSons(split[1], split[2]);
        int Index = -1;

		while(mCladeSplits.size() <= cladeNum)// unknown clade -> initialize its counts
		{
			mCladeSplits.push_back( vector< pair <int,int> >() );
			mSplitsRatio.push_back( vector<double>() );
			mCladeOccurrences.push_back( 0 );
		}

        double nbOccurencesOld = mCladeOccurrences[cladeNum];

        mCladeOccurrences[cladeNum] += iterSplits->second;
        
        double nbOccurencesNew = mCladeOccurrences[cladeNum];
        double NewOldRatio = nbOccurencesOld / nbOccurencesNew; // ratio to convert pre-existing split ratios into new ones


        for(unsigned i = 0 ; i < mCladeSplits[cladeNum].size(); i++)
        {
        	if( CladeSons == mCladeSplits[cladeNum][i]) // found the split
        	{
        		Index = i;
        	}

        	//updating the split ratios
        	mSplitsRatio[cladeNum][i] *= NewOldRatio;
        }
        if( Index == -1) // unknown split -> add it
        {
        	Index = mCladeSplits[cladeNum].size();
        	mCladeSplits[cladeNum].push_back( CladeSons );
        	mSplitsRatio[cladeNum].push_back( 0 );
        }
		mSplitsRatio[cladeNum][Index] += (iterSplits->second / nbOccurencesNew);
			
    }
	mSplitsOccurrencesMap.clear();

    return true;
}




/*
works if source is amalgamted trees
*/
MyCladesAndTripartitions::MyCladesAndTripartitions(MyCladesAndTripartitions * source)
{

    mRooted = false;   ///< amalgamated tree is not rooted

    mGeneTreeCount = source->getmGeneTreeCount();

    mNewickGeneTree = source->getmNewickGeneTree(); ///< save newick string for ALE output

	mClades = source->mClades;


	vector <int> sortedClades = mClades.getSortedClades();

	for(unsigned i = 0 ; i < sortedClades.size(); i++)
	{
		int cladeNum = i;// as this is how things are organized in the cladesandtripartition instance, we keep this structure

		int splitCount = source->getSplitCount( cladeNum );

		mCladeSplits.push_back( vector< pair<int,int> > () );        ///< All splits for each clade (by clade id)
		mSplitsRatio.push_back( vector <double> ()  );        ///< Occurence of each split as a vector, divided by clade occurrences.
		for(unsigned splitNum = 0 ; splitNum < splitCount; splitNum++)
		{
			mCladeSplits[cladeNum].push_back( source->getCladeSplit( cladeNum, splitNum ) ) ;
			mSplitsRatio[cladeNum].push_back( source->getSplitRatio( cladeNum, splitNum ) ) ;
		}

	    mCladeOccurrences.push_back( source->getCladeOccurrences( cladeNum ) );///< Occurence of each clade, indexed by clade number.

		mBranchLengths.push_back( source->getmBranchLengths( cladeNum ) ); ///< branch length for each clade

	}

}

