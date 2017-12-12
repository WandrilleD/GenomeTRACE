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

This file contains the functions relative to cost computation for the class AdjMatrix

Created the: 30-11-2015
by: Wandrille Duchemin

Last modified the: 09-05-2016
by: Wandrille Duchemin

*/

#include "AdjMatrix.h"


///////////////////////////////////
/////// Cost functions ////////////
///////////////////////////////////

double AdjMatrix::AggregateScore(vector <double> scores, int nbGain , int nbBreak )
{
	double s;

	if(scores.size() > 1)
	{

		// worstScore is always contaminant -> auto return worstScore (avoid overflow ?)
		if(scores[0] == worstScore)
			return worstScore;
		else if(scores[1] == worstScore)
			return worstScore;


		s = ((*this).*scoreAggregatorfunc)(scores[0], scores[1]);

		for(unsigned i = 2; i < scores.size(); i++)
		{
			if(scores[i] == worstScore)
				return worstScore;
			s = ((*this).*scoreAggregatorfunc)(s, scores[i]);
		}
	}
	else if(scores.size() == 1)
		s = scores[0];
	else
	{//giving base value?
		s = bestScore; //to adapt later for boltzmann refinement
	}

	//applying Gains and Breaks
	if(!useBoltzmann)
	{
		for(unsigned i = 0; i < nbGain; i++)
			s = ((*this).*scoreAggregatorfunc)(s, getGainCost());
	
		for(unsigned i = 0; i < nbBreak; i++)
			s = ((*this).*scoreAggregatorfunc)(s, getBreakCost());
	}
	else
		s = ((*this).*scoreAggregatorfunc)(s, getBoltzmannGainBreakCost(nbGain, nbBreak));

	return s;

}

/*
Takes:
 - Vsolution (vector <AdjSolution> ): a list of different solution

Returns
	(double): score corresponding to the best solution (non boltzmann case) or the sum of the scores of the solutions (boltzmann case)

*/
double AdjMatrix::compareScore( vector <AdjSolution>  Vsolution)
{
	vector <double> scorevector;


	for(unsigned i = 0; i < Vsolution.size(); i++)
	{
		scorevector.push_back( Vsolution.at(i).score );
	}

	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*
sets VsolutionC1 and VsolutionC0 for NodeId1 and NodeId2

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - VsolutionC1 (vector <AdjSolution> &): vector of solution that will be filled
 - VsolutionC1 (vector <AdjSolution> &): vector of solution that will be filled

*/
void AdjMatrix::computeSolution(int NodeId1, int NodeId2, vector <AdjSolution> &VsolutionC1, vector <AdjSolution> &VsolutionC0 )
{

/*	vector <int> Ids1 = Rtree1.getNodesId();

	for(unsigned i = 0; i < Ids1.size() ; i++)
	{
		cout << Ids1[i] <<",";
	}

	cout << endl;
*/
	Node * n1 = Rtree1.getNode(NodeId1);
	Node * n2 = Rtree2.getNode(NodeId2);


	// The choice is made that the TS compatibility is set by the reconciled trees rather than through the Adj algorithm.
	// This, for instance, allows a simple DeCo algorithm to perform on a full time sliced tree and respecting TS contraints

	if(!Rtree1.areTSCompatible(n1,n2))//if the node aren't TS compatible: their adjacency is impossible
	{
		if(verbose)
		{
			cout << "Nodes "<< NodeId1 << " and " << NodeId2 << " are not time compatible. Setting values accordingly." << endl;
		}

		VsolutionC1 = SolutionC1DefaultImpossibleCase();
		VsolutionC0 = SolutionC0DefaultImpossibleCaseNew(NodeId1,  NodeId2);
		//VsolutionC0 = SolutionC0DefaultImpossibleCase();

		//setRootMatrixNotPossible(NodeId1,NodeId2);//WMODIF

		return;	
	}
		// both node are TS compatible
	
	int evt1 = Rtree1.getNodeEvent(NodeId1);
	int evt2 = Rtree2.getNodeEvent(NodeId2);
	
	
	
	if(decoLTalgo)
	{//check the cases where the species can differ in DeCoLT: both Bout or Rec and different events
		if(evt1 != evt2)
		{
			bool boutwithrec = false;
			bool firstBout;
			if((Rtree1.isRec(evt1)) && (Rtree2.isBout(evt2)))
			{
				boutwithrec = true;
				firstBout = false;
			}
			else if((Rtree1.isBout(evt1)) && (Rtree2.isRec(evt2)))
			{
				boutwithrec = true;
				firstBout = true;
			}

			if(boutwithrec)
			{
				VsolutionC1 = SolutionC1BoutWithRec(NodeId1,NodeId2, firstBout);
				VsolutionC0 = SolutionC0BoutWithRec(NodeId1,NodeId2, firstBout);

				//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
				//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF
				return;
			}
		}
	}

	//now the species compatibility need to hold
	if(!Rtree1.haveSameSpecies(n1, n2))
	{
		if(verbose)
		{
			cout << "Nodes "<< NodeId1 << " and " << NodeId2 << " have different species. Setting values accordingly." << endl;
		}

		VsolutionC1 = SolutionC1DefaultImpossibleCase();
		VsolutionC0 = SolutionC0DefaultImpossibleCaseNew(NodeId1,  NodeId2);
		//VsolutionC0 = SolutionC0DefaultImpossibleCase();

		//setRootMatrixNotPossible(NodeId1,NodeId2);//WMODIF
		

		return;
	}


	//cases where species is the same

	bool firstLoss = Rtree1.isLoss(evt1);
	bool secondLoss = Rtree2.isLoss(evt2);

	if( (firstLoss) || (secondLoss) )
	{
		if( (firstLoss) && (secondLoss) )
		{
			VsolutionC1 = SolutionC1LossWithLoss(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0LossWithLoss(NodeId1,NodeId2);

			//setRootMatrixNotPossible(NodeId1,NodeId2);//WMODIF
			

			return; 
		}
		else
		{
			VsolutionC1 = SolutionC1LossWithOther(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0LossWithOther(NodeId1,NodeId2);

			//setRootMatrixNotPossible(NodeId1,NodeId2);//WMODIF
			
			return;
		}
	}


	bool firstDup = Rtree1.isDup(evt1);
	bool secondDup = Rtree2.isDup(evt2);

	if(decoLTalgo) // Sout specific part
	{

		bool firstRec = Rtree1.isRec(evt1);
		bool secondRec = Rtree2.isRec(evt2);
		bool firstSout = Rtree1.isSout(evt1);
		bool secondSout = Rtree2.isSout(evt2);
		if( (firstRec) || (secondRec) )
		{
			if( (firstRec) && (secondRec) )
			{
				VsolutionC1 = SolutionC1RecWithRec(NodeId1,NodeId2);
				VsolutionC0 = SolutionC0RecWithRec(NodeId1,NodeId2);

				//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
				//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

				return ;
			}
			else
			{


				if( (firstSout) || (secondSout) )
				{
					VsolutionC1 = SolutionC1RecWithSout(NodeId1,NodeId2,firstSout);
					VsolutionC0 = SolutionC0RecWithSout(NodeId1,NodeId2,firstSout);

				}
				else if( (firstDup) || (secondDup) )
				{
					VsolutionC1 = SolutionC1DupWithRec( NodeId1, NodeId2, firstDup );
					VsolutionC0 = SolutionC0DupWithRec( NodeId1, NodeId2, firstDup );
				}
				else
				{
					VsolutionC1 = SolutionC1RecWithOther(NodeId1,NodeId2,firstRec);
					VsolutionC0 = SolutionC0RecWithOther(NodeId1,NodeId2,firstRec);

					//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
					//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF
				}
				return ;
			}
		}
		else if( (firstSout) || (secondSout) )
		{
			if( (firstDup) || (secondDup) )
			{
				VsolutionC1 = SolutionC1DupWithSout( NodeId1, NodeId2, firstDup );
				VsolutionC0 = SolutionC0DupWithSout( NodeId1, NodeId2, firstDup );
				return;
			}

		}
		
	}




	if( (firstDup) || (secondDup) )
	{
		if( (firstDup) && (secondDup) )
		{
			VsolutionC1 = SolutionC1DupWithDup(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0DupWithDup(NodeId1,NodeId2);

			//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
			//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

			return ;
		}
		else 
		{
			VsolutionC1 = SolutionC1DupWithOther(NodeId1,NodeId2,firstDup);
			VsolutionC0 = SolutionC0DupWithOther(NodeId1,NodeId2,firstDup);

			//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
			//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

			return;
		}
	}




	if(decoLTalgo)
	{
		bool firstSout = Rtree1.isSout(evt1);
		bool secondSout = Rtree2.isSout(evt2);

		if( (firstSout) || (secondSout) )
		{
			if( (firstSout) && (secondSout) )
			{
				VsolutionC1 = SolutionC1SoutWithSout(NodeId1,NodeId2);
				VsolutionC0 = SolutionC0SoutWithSout(NodeId1,NodeId2);

				//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
				//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

				return ;
			}
			else
			{
				VsolutionC1 = SolutionC1SoutWithExtantOrSpecOrNull(NodeId1,NodeId2,firstSout);
				VsolutionC0 = SolutionC0SoutWithExtantOrSpecOrNull(NodeId1,NodeId2,firstSout);

				//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
				//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

				return ;
			}
		}

	}


	//Now: cases where both event must be equal
	if(evt1 == evt2)
	{
		if(Rtree1.isExtant(evt1))
		{
			VsolutionC1 = SolutionC1ExtantWithExtant(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0ExtantWithExtant(NodeId1,NodeId2);

			//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
			//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

			return ;	
		}
		else if(Rtree1.isSpeciation(evt1))
		{
			VsolutionC1 = SolutionC1SpecWithSpec(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0SpecWithSpec(NodeId1,NodeId2);

			//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
			//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF
			return ;	
		}
		else if(Rtree1.isNull(evt1))
		{
			VsolutionC1 = SolutionC1NullWithNull(NodeId1,NodeId2);
			VsolutionC0 = SolutionC0NullWithNull(NodeId1,NodeId2);

			//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
			//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

			return ;
		}
		else if(decoLTalgo)
		{//bout with bout and rec with rec cases
			if(Rtree1.isBout(evt1))
			{
				VsolutionC1 = SolutionC1BoutWithBout(NodeId1,NodeId2);
				VsolutionC0 = SolutionC0BoutWithBout(NodeId1,NodeId2);

				//setRootMatrixPossible(NodeId1,NodeId2);//WMODIF
				//scanSolutionForRooting(NodeId1,NodeId2);//WMODIF

				return ;		
			}
		}
	}


	//we are out of the case covered -> return some default value
	if(verbose)
	{
		cout << "Nodes "<< NodeId1 << " (" << evt1 << ")" << " and " << NodeId2 << " (" << evt2 << ")" << " do not fall in any cases. Setting values accordingly." << endl;
	}
	VsolutionC1 = SolutionC1DefaultImpossibleCase();//impossible adjacency
	VsolutionC0 = SolutionC0DefaultImpossibleCaseNew(NodeId1,  NodeId2);
	//VsolutionC0 = SolutionC0DefaultImpossibleCase();
	//setRootMatrixNotPossible(NodeId1,NodeId2);//WMODIF
	return;


}



// simple DeCo cases
/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector <AdjSolution> AdjMatrix::SolutionC1ExtantWithExtant(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1ExtantWithExtant"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,true) );

	Vsolution.back().score = getC1(NodeId1,NodeId2); // the score of that solution

	return Vsolution;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector <AdjSolution> AdjMatrix::SolutionC0ExtantWithExtant(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0ExtantWithExtant"<< endl;


	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = getC0(NodeId1,NodeId2); // the score of that solution

	return Vsolution;
}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
*/
vector <AdjSolution> AdjMatrix::SolutionC1LossWithOther(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1LossWithOther"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = bestScore;

	return Vsolution;

}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector <AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
*/
vector <AdjSolution> AdjMatrix::SolutionC0LossWithOther(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0LossWithOther"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = bestScore;

	return Vsolution;

}




/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: add something later to account for possible co-events
*/
vector<AdjSolution> AdjMatrix::SolutionC1LossWithLoss(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1LossWithLoss"<< endl;
	
	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,true) );

	Vsolution.back().score = WeightedLossCost; // Always the WeightedLossCost as this is either bestScore (default) or a value to add/multiply to bestscore = that value

	return Vsolution;

}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

NB: add something later to account for possible co-events
*/
vector<AdjSolution> AdjMatrix::SolutionC0LossWithLoss(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0LossWithLoss"<< endl;

	vector<AdjSolution> Vsolution;

	//only one possible solution
	Vsolution.push_back( AdjSolution(0,0,false) );

	Vsolution.back().score = bestScore;

	return Vsolution;	
}

/*
Takes:
 - NodeIdDup (int): node id in RtreeDup
 - NodeIdOther (int): node id in RtreeOther
 - RtreeDup (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the dupplication is in the tree1 (true) or not (false)

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeIdDup and NodeIdOther
*/
vector<AdjSolution> AdjMatrix::SolutionD1(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert )
{

	vector<AdjSolution> Vsolution;
	bool iscoevent = false;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId = RtreeDup->getSonsId(NodeIdDup);

	//cout << "PLOP " << NodeIdDup << " " << NodeIdOther << " " << invert << " : " << DupSonsId[0] << " " << DupSonsId[1] <<  endl;

	//There is 4 possible cases
	vector <double> caseScore;

	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab1.id1 = NodeIdOther;
		Sc1ab1.id2 = DupSonsId[0];
	}
	else
	{
		Sc1ab1.id1 = DupSonsId[0];
		Sc1ab1.id2 = NodeIdOther;
	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);



	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab2.id1 = NodeIdOther;
		Sc1ab2.id2 = DupSonsId[1];
	}
	else
	{
		Sc1ab2.id1 = DupSonsId[1];
		Sc1ab2.id2 = NodeIdOther;
	}
	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);

	

	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab1.id1 = NodeIdOther;
		Sc0ab1.id2 = DupSonsId[0];
	}
	else
	{
		Sc0ab1.id1 = DupSonsId[0];
		Sc0ab1.id2 = NodeIdOther;
	}
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab2.id1 = NodeIdOther;
		Sc0ab2.id2 = DupSonsId[1];
	}
	else
	{
		Sc0ab2.id1 = DupSonsId[1];
		Sc0ab2.id2 = NodeIdOther;
	}
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);


	//case 1 -> child 1 of the duplication conserved an adjacency but not child 2
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	Vsolution.push_back(*currentSolution);

	
	delete currentSolution;


	caseScore.clear();
	//case 2 -> child 2 of the duplication conserved an adjacency but not child 1
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab2 );
	currentSolution->components.push_back( Sc0ab1 );
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 1 of the duplication conserved an adjacency and child 2 too -> need a gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency -> need a break
	currentSolution = new AdjSolution(0,1,iscoevent) ;  // 0 gain 1 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;

	return Vsolution;
}


/*
Takes:
 - NodeIdDup (int): node id in RtreeDup
 - NodeIdOther (int): node id in RtreeOther
 - RtreeDup (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the dupplication is in the tree1 (true) or not (false)

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeIdDup and NodeIdOther
*/
vector<AdjSolution> AdjMatrix::SolutionD0(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert )
{
	vector<AdjSolution> Vsolution;
	bool iscoevent = false;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId = RtreeDup->getSonsId(NodeIdDup);

	//There is 4 possible cases
	vector <double> caseScore;

	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab1.id1 = NodeIdOther;
		Sc1ab1.id2 = DupSonsId[0];
	}
	else
	{
		Sc1ab1.id1 = DupSonsId[0];
		Sc1ab1.id2 = NodeIdOther;
	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);



	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab2.id1 = NodeIdOther;
		Sc1ab2.id2 = DupSonsId[1];
	}
	else
	{
		Sc1ab2.id1 = DupSonsId[1];
		Sc1ab2.id2 = NodeIdOther;
	}
	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);

	

	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab1.id1 = NodeIdOther;
		Sc0ab1.id2 = DupSonsId[0];
	}
	else
	{
		Sc0ab1.id1 = DupSonsId[0];
		Sc0ab1.id2 = NodeIdOther;
	}
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab2.id1 = NodeIdOther;
		Sc0ab2.id2 = DupSonsId[1];
	}
	else
	{
		Sc0ab2.id1 = DupSonsId[1];
		Sc0ab2.id2 = NodeIdOther;
	}
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);



	//case 1 -> child 1 of the duplication conserved an adjacency but not child 2 -> need a gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 2 -> child 2 of the duplication conserved an adjacency but not child 1 -> need a gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab2 );
	currentSolution->components.push_back( Sc0ab1 );
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 1 of the duplication conserved an adjacency and child 2 too -> need 2 gain
	currentSolution = new AdjSolution(2,0,iscoevent) ;  // 2 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency 
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;

	return Vsolution;
}


/*
Both nodes are duplications that are presumed simultaneous

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionD12(int NodeId1, int NodeId2)
{
	vector<AdjSolution> Vsolution;
	bool iscoevent = true;

	AdjSolution * currentSolution;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId1 = Rtree1.getSonsId(NodeId1);
	vector <int> DupSonsId2 = Rtree2.getSonsId(NodeId2);


	//There is 16 possible cases representing all c1 and c0 combination of the sons of 1 an 2
	vector <double> caseScore;

	for(unsigned c1a1b1 = 0; c1a1b1 < 2; c1a1b1 ++)
			for(unsigned c1a1b2 = 0; c1a1b2 < 2; c1a1b2 ++)
					for(unsigned c1a2b1 = 0; c1a2b1 < 2; c1a2b1 ++)
							for(unsigned c1a2b2 = 0; c1a2b2 < 2; c1a2b2 ++)
							{
								currentSolution = new AdjSolution();
								currentSolution->coevent = iscoevent;

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[0],DupSonsId2[0]));
								if(!c1a1b1)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[0]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[0]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[0],DupSonsId2[1]));
								if(!c1a1b2)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[1]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[1]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[1],DupSonsId2[0]));
								if(!c1a2b1)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[0]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[0]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[1],DupSonsId2[1]));
								if(!c1a2b2)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[1]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[1]));

								//determining the number of gains and breaks
								int nbadj = c1a1b1 + c1a1b2 + c1a2b1 + c1a2b2;
								switch( nbadj )
								{

									case 0: // none linked -> 2 break

										currentSolution->NbGain = 0;
										currentSolution->NbBreak = 2;
										break;

									case 1: // only one linked -> 1 break
										currentSolution->NbGain = 0;
										currentSolution->NbBreak = 1;	
										break;
									case 3: // three linked -> 1 gain
										currentSolution->NbGain = 1;
										currentSolution->NbBreak = 0;	
										break;
									case 4: // all linked -> 2 gain
										currentSolution->NbGain = 2;
										currentSolution->NbBreak = 0;	
										break;
									default: //exactly 2 adjs
										currentSolution->NbGain = 0;
										currentSolution->NbBreak = 0;
										if(( c1a1b1 == c1a1b2 ) || (c1a1b1 == c1a2b1 )) // scenarios where 1 child has two adjs -> one gain and one break
										{
											currentSolution->NbGain = 1;
											currentSolution->NbBreak = 1;
										}
								}

								if(WeightedDupCost != bestScore) //subtracting a Duplication to the score if that option was set
									caseScore.push_back(WeightedDupCost);

								//getting the new score
								currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);
								
								if(Vsolution.empty())
									Vsolution.push_back(*currentSolution);
								else
								{
									if(currentSolution->score <= Vsolution.front().score)
										Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
									else
										Vsolution.push_back(*currentSolution);
								}

								//clearing
								delete currentSolution;
								caseScore.clear();
							}


	return Vsolution;
}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
vector<AdjSolution> AdjMatrix::SolutionC1DupWithOther(int NodeId1, int NodeId2, bool firstDup )
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithOther"<< endl;	

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;
	
	vector<AdjSolution> Vsolution;
	if(!firstDup) // NodeId2 is the duplication
		Vsolution = SolutionD1(NodeId2, NodeId1, rtree2, rtree1, !firstDup);
	else // nodeId1 is the duplication
		Vsolution = SolutionD1(NodeId1, NodeId2, rtree1, rtree2, !firstDup);

	return Vsolution;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
vector<AdjSolution> AdjMatrix::SolutionC0DupWithOther(int NodeId1, int NodeId2, bool firstDup )
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0DupWithOther"<< endl;	

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;
	
	vector<AdjSolution> Vsolution;
	if(!firstDup) // NodeId2 is the duplication
		Vsolution = SolutionD0(NodeId2, NodeId1, rtree2, rtree1, !firstDup);
	else // nodeId1 is the duplication
		Vsolution = SolutionD0(NodeId1, NodeId2, rtree1, rtree2, !firstDup);

	return Vsolution;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC1SpecWithSpec(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SpecWithSpec"<< endl;	

	vector<AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	bool iscoevent = true;
	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);

	double sa1,sa2,sb1,sb2;

	if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[0]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[1]) ) ))//case where son1 of node 1 corresponds to son1 of node 2 // using the OR here does not change a thing in the case of the speciation, but will allow to treat them in the same way as synchronous speciation out
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[0];
		sb2 = sonsId2[1];
	}
	else if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[1]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[0]) ) )) //sanity check
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[1];
		sb2 = sonsId2[0];	
	}
	else
	{
		throw Exception("AdjMatrix::C1SpecWithSpec : no children corresponds to each other!");
	}



	double c1a1b1 = getC1(sa1, sb1);
	AdjScore Sc1a1b1(true ,sa1,sb1);
	double c0a1b1 = getC0(sa1, sb1);
	AdjScore Sc0a1b1(false,sa1,sb1);
	double c1a2b2 = getC1(sa2, sb2);
	AdjScore Sc1a2b2(true ,sa2,sb2);
	double c0a2b2 = getC0(sa2, sb2);
	AdjScore Sc0a2b2(false,sa2,sb2);

	vector <double> scorevector;
	vector <double> caseScore;

	//case 1 -> both adjacency have been kept
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc1a2b2 );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 2 -> adjacency have been kept in child 1 but not child 2 -> 1 break
	currentSolution = new AdjSolution(0,1,iscoevent) ;  // 0 gain 1 break 
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc0a2b2 );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 3 -> adjacency have been kept in child 2 but not child 1 -> 1 break
	currentSolution = new AdjSolution(0,1,iscoevent) ;  // 0 gain 1 break 
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc1a2b2 );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 4 -> adjacency have not been kept in child 2 and child 1 -> 2 break
	currentSolution = new AdjSolution(0,2,iscoevent) ;  // 0 gain 2 break 
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc0a2b2 );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacencies between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC0SpecWithSpec(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SpecWithSpec"<< endl;	

	vector<AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	bool iscoevent = false;
	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);

	double sa1,sa2,sb1,sb2;

	if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[0]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[1]) ) ))//case where son1 of node 1 corresponds to son1 of node 2 // using the OR here does not change a thing in the case of the speciation, but will allow to treat them in the same way as synchronous speciation out
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[0];
		sb2 = sonsId2[1];
	}
	else if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[1]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[0]) ) )) //sanity check
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[1];
		sb2 = sonsId2[0];	
	}
	else
	{
		throw Exception("AdjMatrix::C1SpecWithSpec : no children corresponds to each other!");
	}



	double c1a1b1 = getC1(sa1, sb1);
	AdjScore Sc1a1b1(true ,sa1,sb1);
	double c0a1b1 = getC0(sa1, sb1);
	AdjScore Sc0a1b1(false,sa1,sb1);
	double c1a2b2 = getC1(sa2, sb2);
	AdjScore Sc1a2b2(true ,sa2,sb2);
	double c0a2b2 = getC0(sa2, sb2);
	AdjScore Sc0a2b2(false,sa2,sb2);

	vector <double> scorevector;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept -> 2 gain
	currentSolution = new AdjSolution(2,0,iscoevent) ;  // 2 gain 0 break 
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc1a2b2 );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 2 -> adjacency have been kept in child 1 but not child 2 -> 1 Gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc0a2b2 );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 3 -> adjacency have been kept in child 2 but not child 1 -> 1 gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc1a2b2 );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 4 -> adjacency have not been kept in child 2 and child 1
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc0a2b2 );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0a2b2 );

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();


	return Vsolution;

}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: add something later to account for co-event
*/
vector<AdjSolution> AdjMatrix::SolutionC1DupWithDup(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithDup"<< endl;	

	vector<AdjSolution> Vsolution,tmp;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	//D1 -> the duplication of NodeId1 was before the duplication of NodeId2
	Vsolution = SolutionD1(NodeId1, NodeId2, rtree1, rtree2, false);

	//D2 -> the duplication of NodeId2 was before the duplication of NodeId1
	tmp = SolutionD1(NodeId2, NodeId1, rtree2, rtree1, true);

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	//D12 -> both duplications were simultaneous
	tmp = SolutionD12(NodeId1,NodeId2);

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacencies between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC0DupWithDup(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0DupWithDup"<< endl;	

	vector<AdjSolution> Vsolution,tmp;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	//D0 -> the duplication of NodeId1 was before the duplication of NodeId2
	Vsolution = SolutionD0(NodeId1, NodeId2, rtree1, rtree2, false);

	//D0 -> the duplication of NodeId2 was before the duplication of NodeId1
	tmp = SolutionD0(NodeId2, NodeId1, rtree2, rtree1, true);

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	return Vsolution;
}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC1NullWithNull(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1NullWithNull"<< endl;

	vector<AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	//NodeId1 and NodeId2 only have 1 child
	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);


	//2possible cases
	vector <double> caseScore;

	//case 1 -> an adjacency exists between their children
	currentSolution = new AdjSolution(0,0,false);
	currentSolution->components.push_back(AdjScore(true, sonsId1[0], sonsId2[0]));
	caseScore.push_back( getC1(sonsId1[0], sonsId2[0]) );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // 0 gain 0 break

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 1 -> no adjacency exists between their children -> we need a break
	currentSolution = new AdjSolution(0,1,false);
	currentSolution->components.push_back(AdjScore(false, sonsId1[0], sonsId2[0]));
	caseScore.push_back( getC0(sonsId1[0], sonsId2[0]) );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // 0 gain 1 break

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacencies between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC0NullWithNull(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0NullWithNull"<< endl;


	vector<AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	//NodeId1 and NodeId2 only have 1 child
	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);


	//2possible cases
	vector <double> caseScore;

	//case 1 -> an adjacency exists between their children -> need 1 gain
	currentSolution = new AdjSolution(1,0,false);
	currentSolution->components.push_back(AdjScore(true, sonsId1[0], sonsId2[0]));
	caseScore.push_back( getC1(sonsId1[0], sonsId2[0]) );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // 0 gain 0 break

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 1 -> no adjacency exists between their children
	currentSolution = new AdjSolution(0,0,false);
	currentSolution->components.push_back(AdjScore(false, sonsId1[0], sonsId2[0]));
	caseScore.push_back( getC0(sonsId1[0], sonsId2[0]) );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // 0 gain 1 break
	
	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
vector<AdjSolution> AdjMatrix::SolutionC1DupWithRec(int NodeId1, int NodeId2, bool firstDup )
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithRec"<< endl;	
	
	vector<AdjSolution> Vsolution;

	Vsolution = SolutionC1DupWithOther(NodeId1, NodeId2, firstDup ); // cases where the dup happened first

	vector<AdjSolution> tmp =  SolutionC1RecWithOther(NodeId1, NodeId2, !firstDup);

	cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithRec --> tmp: "<< tmp.size()<< " " ;
	cout << tmp.back().components.back().id1 << "-" << tmp.back().components.back().id2 << endl;

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}	


	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
vector<AdjSolution> AdjMatrix::SolutionC0DupWithRec(int NodeId1, int NodeId2, bool firstDup )
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithRec"<< endl;	
	
	vector<AdjSolution> Vsolution;

	Vsolution = SolutionC0DupWithOther(NodeId1, NodeId2, firstDup ); // cases where the dup happened first

	vector<AdjSolution> tmp =  SolutionC0RecWithOther(NodeId1, NodeId2, !firstDup);

	cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithRec --> tmp: "<< tmp.size()<< " " ;
	cout << tmp.back().components.back().id1 << "-" << tmp.back().components.back().id2 << endl;

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}	


	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
vector<AdjSolution> AdjMatrix::SolutionC1DupWithSout(int NodeId1, int NodeId2, bool firstDup )
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithSout"<< endl;	
	
	vector<AdjSolution> Vsolution;

	Vsolution = SolutionC1DupWithOther(NodeId1, NodeId2, firstDup ); // cases where the dup happened first

	vector<AdjSolution> tmp = SolutionC1SoutWithExtantOrSpecOrNull(NodeId1, NodeId2, !firstDup); // cases where the dup happened second


	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}	

	//cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithSout --> tmp: "<< Vsolution.size()<< " " ;
	//for(unsigned i = 0; i < Vsolution.size();i++)
	//{
	//	for(unsigned j = 0; j < Vsolution[i].components.size();j++)
	//	{
	//		cout << Vsolution[i].components[j].id1 << "-" << Vsolution[i].components[j].id2 << " ";
	//	}
	//}
	//cout << endl;



	return Vsolution;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
vector<AdjSolution> AdjMatrix::SolutionC0DupWithSout(int NodeId1, int NodeId2, bool firstDup )
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithSout"<< endl;	
	
	vector<AdjSolution> Vsolution;

	Vsolution = SolutionC0DupWithOther(NodeId1, NodeId2, firstDup );

	vector<AdjSolution> tmp = SolutionC0SoutWithExtantOrSpecOrNull(NodeId1, NodeId2, !firstDup);


	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}	
	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstSout (bool): wether the Sout is NodeId1 or NodeId2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacencies between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC1SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SoutWithExtantOrSpecOrNull"<< endl;

	pair <int,int> orderedIds;
	orderedIds = CostSoutWithExtantOrSpecOrNullAux( NodeId1,  NodeId2, firstSout);

	//now that we have the correct son, we can apply the formula
	//only one case, with no gain or break -> simple transmission

	vector<AdjSolution> Vsolution;

	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().components.push_back(AdjScore(true, orderedIds.first, orderedIds.second ));
	Vsolution.back().score = getC1(orderedIds.first, orderedIds.second); // 0 gain 0 break

	return Vsolution;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstSout (bool): wether the Sout is NodeId1 or NodeId2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacencies between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC0SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SoutWithExtantOrSpecOrNull"<< endl;

	pair <int,int> orderedIds;
	orderedIds = CostSoutWithExtantOrSpecOrNullAux( NodeId1,  NodeId2, firstSout);

	//now that we have the correct son, we can apply the formula
	//only one case, with no gain or break -> simple transmission

	vector<AdjSolution> Vsolution;

	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().components.push_back(AdjScore(false, orderedIds.first, orderedIds.second));
	Vsolution.back().score = getC0(orderedIds.first, orderedIds.second); // 0 gain 0 break

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacencies between NodeId1 and NodeId2
Presumes that both Speciation out event occurs at the same time
*/
vector<AdjSolution> AdjMatrix::SolutionC1SoutWithSoutSynchronous(int NodeId1, int NodeId2) 
{
	//presumes that both Speciation out event occurs at the same time
	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> SonsIds1 = Rtree1.getSonsId(NodeId1);
	vector < int> SonsIds2 = Rtree2.getSonsId(NodeId2);

	double sa1,sat,sb1,sbt; //sat and sbt are the transferred children

	bool iscoevent = true;


	int CorrectSonInd;
	int sp;
	unsigned i;
	//we want to get the Son of the Sout that stayed in the species tree for the node 1
	CorrectSonInd = -1;
	sp =  Rtree1.getNodeSpecies(NodeId1);
	for(i = 0; i < SonsIds1.size(); i++)
	{
		if(Rtree1.getNodeSpecies(SonsIds1[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::SolutionC1SoutWithSoutSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sat = SonsIds1[1];
		sa1 = SonsIds1[0];
	}
	else // the non transferred son as index 1
	{
		sat = SonsIds1[0];
		sa1 = SonsIds1[1];	
	}

	//we want to get the Son of the Sout that stayed in the species tree for the node 2
	CorrectSonInd = -1;
	sp =  Rtree2.getNodeSpecies(NodeId2);
	for(i = 0; i < SonsIds2.size(); i++)
	{
		if(Rtree2.getNodeSpecies(SonsIds2[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::SolutionC1SoutWithSoutSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sbt = SonsIds2[1];
		sb1 = SonsIds2[0];
	}
	else // the non transferred son as index 1
	{
		sbt = SonsIds2[0];
		sb1 = SonsIds2[1];
	}


	double c1a1b1 = getC1(sa1, sb1);
	AdjScore Sc1a1b1(true,sa1, sb1);
	double c0a1b1 = getC0(sa1, sb1);
	AdjScore Sc0a1b1(false,sa1,sb1);
	double c1atbt = getC1(sat, sbt);
	AdjScore Sc1atbt(true,sat, sbt);
	double c0atbt = getC0(sat, sbt);
	AdjScore Sc0atbt(false,sat,sbt);
	//each occurence of c0a1b1 costs a break; c0atbt are free because of the possibly asynchronous transfer events

	vector <AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc1atbt );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // 0 gain 0 break

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 2 -> adjacency have been kept in child 1 but not child t
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc0atbt );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // 0 gain 0 break

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 3 -> adjacency have been kept in child t but not child 1 -> 1 break
	currentSolution = new AdjSolution(0,1,iscoevent);
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc1atbt );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();

	//case 4 -> adjacency have not been kept in child t and child 1 -> 1 break
	currentSolution = new AdjSolution(0,1,iscoevent);
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc0atbt );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacencies between NodeId1 and NodeId2
Presumes that both Speciation out events do not occur at the same time
*/
vector<AdjSolution> AdjMatrix::SolutionC1SoutWithSoutaSynchronous(int NodeId1, int NodeId2) 
{
	//presumes that both Speciation out events don't occurs at the same time
	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> SonsIds1 = Rtree1.getSonsId(NodeId1);
	vector < int> SonsIds2 = Rtree2.getSonsId(NodeId2);

	double sa1,sat,sb1,sbt; //sat and sbt are the transferred children

	bool iscoevent = false;

	int CorrectSonInd;
	int sp;
	unsigned i;
	//we want to get the Son of the Sout that stayed in the species tree for the node 1
	CorrectSonInd = -1;
	sp =  Rtree1.getNodeSpecies(NodeId1);
	for(i = 0; i < SonsIds1.size(); i++)
	{
		if(Rtree1.getNodeSpecies(SonsIds1[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::SolutionC1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sat = SonsIds1[1];
		sa1 = SonsIds1[0];
	}
	else // the non transferred son as index 1
	{
		sat = SonsIds1[0];
		sa1 = SonsIds1[1];	
	}

	//we want to get the Son of the Sout that stayed in the species tree for the node 2
	CorrectSonInd = -1;
	sp =  Rtree2.getNodeSpecies(NodeId2);
	for(i = 0; i < SonsIds2.size(); i++)
	{
		if(Rtree2.getNodeSpecies(SonsIds2[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::SolutionC1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sbt = SonsIds2[1];
		sb1 = SonsIds2[0];
	}
	else // the non transferred son as index 1
	{
		sbt = SonsIds2[0];
		sb1 = SonsIds2[1];
	}

	// the two Sout are asynchronous. This means that one arrived before the other.

	/* NEWWMODIF
	//adjacencies between transfered sons
	double c1atbt = getC1(sat, sbt);
	AdjScore Sc1atbt(true,sat, sbt);
	double c0atbt = getC0(sat, sbt);
	AdjScore Sc0atbt(false,sat,sbt);
	*/

	// here we compare the transfered sons and not one transfered son with a Sout node (contrary to what we do in the double duplication case)
	// this is due to the fact that the transfered son and the Sout node re not comparable
	// but also because doing so would block backtrack: the equivalence classes would depend on the choices done during the backtrack.

	//CASETYPE 1 : NodeId1 was the first sout
	double c1a1bp = getC1(sa1, NodeId2);
	AdjScore Sc1a1bp(true,sa1, NodeId2);
	double c0a1bp = getC0(sa1, NodeId2);
	AdjScore Sc0a1bp(true,sa1, NodeId2);


	double c1atbp = getC1(sat, NodeId2);
	AdjScore Sc1atbp(true,sat, NodeId2);
	double c0atbp = getC0(sat, NodeId2);
	AdjScore Sc0atbp(true,sat, NodeId2);



	vector <AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept -> 1 gain in the transfered sons
	currentSolution = new AdjSolution(1,0,iscoevent);
	currentSolution->components.push_back( Sc1a1bp );
	//currentSolution->components.push_back( Sc1atbt );
	currentSolution->components.push_back( Sc1atbp );
	caseScore.push_back( c1a1bp );
	//caseScore.push_back( c1atbt );
	caseScore.push_back( c1atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 2 -> only the adjacencies in the same species has been kept -> 0 gain nor loss
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back( Sc1a1bp );
	//currentSolution->components.push_back( Sc0atbt );
	currentSolution->components.push_back( Sc0atbp );
	caseScore.push_back( c1a1bp );
	//caseScore.push_back( c0atbt );
	caseScore.push_back( c0atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 3 -> only the adjacency in transfered sons has been kept -> 1 gain and 1 break 
	currentSolution = new AdjSolution(1,1,iscoevent);
	currentSolution->components.push_back( Sc0a1bp );
	currentSolution->components.push_back( Sc1atbp );
	caseScore.push_back( c0a1bp );
	caseScore.push_back( c1atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();


	//case 4 -> no adjacencies have been kept -> 0 gain , 1 break
	currentSolution = new AdjSolution(0,1,iscoevent);
	currentSolution->components.push_back( Sc0a1bp ); // 1 break
	currentSolution->components.push_back( Sc0atbp );
	caseScore.push_back( c0a1bp );
	caseScore.push_back( c0atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//CASETYPE 2 : NodeId2 was the first sout
	double c1apb1 = getC1(NodeId1, sb1);
	AdjScore Sc1apb1(true,NodeId1, sb1);
	double c0apb1 = getC0(NodeId1, sb1);
	AdjScore Sc0apb1(true,NodeId1, sb1);
	double c1apbt = getC1(NodeId1, sbt);
	AdjScore Sc1apbt(true,NodeId1, sbt);
	double c0apbt = getC0(NodeId1, sbt);
	AdjScore Sc0apbt(true,NodeId1, sbt);



	//case 1 -> both adjacencies have been kept -> 1 gain in the transfered sons
	currentSolution = new AdjSolution(1,0,iscoevent);
	currentSolution->components.push_back( Sc1apb1 );
	currentSolution->components.push_back( Sc1apbt ); //1 gain
	caseScore.push_back( c1apb1 );
	caseScore.push_back( c1apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 2 -> only the adjacencies in the same species has been kept -> 0 gain nor loss
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back( Sc1apb1 ); 
	currentSolution->components.push_back( Sc0apbt );
	caseScore.push_back( c1apb1 );
	caseScore.push_back( c0apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 3 -> only the adjacency in transfered sons has been kept -> 1 gain and 1 break 
	currentSolution = new AdjSolution(1,1,iscoevent);
	currentSolution->components.push_back( Sc0apb1 ); // 1 break
	currentSolution->components.push_back( Sc1apbt ); // 1 gain
	caseScore.push_back( c0apb1 );
	caseScore.push_back( c1apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();


	//case 4 -> no adjacencies have been kept -> 0 gain , 1 break
	currentSolution = new AdjSolution(0,1,iscoevent);
	currentSolution->components.push_back( Sc0apb1 ); // 1 break
	currentSolution->components.push_back( Sc0apbt );
	caseScore.push_back( c0a1bp );
	caseScore.push_back( c0apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();



	// These formula are erroneous

	/*
	double c1a1b1 = getC1(sa1, sb1);
	AdjScore Sc1a1b1(true,sa1, sb1);
	double c0a1b1 = getC0(sa1, sb1);
	AdjScore Sc0a1b1(false,sa1,sb1);
	double c1atbt = getC1(sat, sbt);
	AdjScore Sc1atbt(true,sat, sbt);
	double c0atbt = getC0(sat, sbt);
	AdjScore Sc0atbt(false,sat,sbt);
	
	//each occurence of c0a1b1 costs a break, each occurence of c1atbt costs a gain.
	

	vector <AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept -> 1 gain
	currentSolution = new AdjSolution(1,0,iscoevent);
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc1atbt );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	Vsolution.push_back(*currentSolution);

	caseScore.clear();

	//case 2 -> adjacency have been kept in child 1 but not child t
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back( Sc1a1b1 );
	currentSolution->components.push_back( Sc0atbt );
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	caseScore.clear();

	//case 3 -> adjacency have been kept in child t but not child 1 -> 1 gain and 1 break
	currentSolution = new AdjSolution(1,1,iscoevent);
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc1atbt );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	//case 4 -> adjacency have not been kept in child t and child 1 -> 1 break
	currentSolution = new AdjSolution(0,1,iscoevent);
	currentSolution->components.push_back( Sc0a1b1 );
	currentSolution->components.push_back( Sc0atbt );
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0atbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	caseScore.clear();*/

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1SoutWithSout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SoutWithSout"<< endl;

	vector<AdjSolution> Vsolution,tmp;

	Vsolution = SolutionC1SoutWithSoutSynchronous( NodeId1,NodeId2); // co event case
	tmp = SolutionC1SoutWithSoutaSynchronous(NodeId1,NodeId2); // non co event case

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}	


	//cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SoutWithSout --> Vsolution: "<< Vsolution.size()<< " " ;
	//for(unsigned i = 0; i < Vsolution.size();i++)
	//{
	//	for(unsigned j = 0; j < Vsolution[i].components.size();j++)
	//	{
	//		cout << Vsolution[i].components[j].id1 << "-" << Vsolution[i].components[j].id2 << " ";
	//	}
	//}
	//cout << endl;



	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacencies between NodeId1 and NodeId2

NB: same formula wether the Souts are simultaneous or not 
*/
vector<AdjSolution> AdjMatrix::SolutionC0SoutWithSout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0SoutWithSout"<< endl;

	/////////// VVV Not sure we can.
	//we can treat that event with the speciation formulae (which only need one of the two children to have the same species)
	//return SolutionC0SpecWithSpec(NodeId1, NodeId2);


	vector < int> SonsIds1 = Rtree1.getSonsId(NodeId1);
	vector < int> SonsIds2 = Rtree2.getSonsId(NodeId2);

	double sa1,sat,sb1,sbt; //sat and sbt are the transferred children

	bool iscoevent = false;

	int CorrectSonInd;
	int sp;
	unsigned i;
	//we want to get the Son of the Sout that stayed in the species tree for the node 1
	CorrectSonInd = -1;
	sp =  Rtree1.getNodeSpecies(NodeId1);
	for(i = 0; i < SonsIds1.size(); i++)
	{
		if(Rtree1.getNodeSpecies(SonsIds1[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::SolutionC1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sat = SonsIds1[1];
		sa1 = SonsIds1[0];
	}
	else // the non transferred son as index 1
	{
		sat = SonsIds1[0];
		sa1 = SonsIds1[1];	
	}

	//we want to get the Son of the Sout that stayed in the species tree for the node 2
	CorrectSonInd = -1;
	sp =  Rtree2.getNodeSpecies(NodeId2);
	for(i = 0; i < SonsIds2.size(); i++)
	{
		if(Rtree2.getNodeSpecies(SonsIds2[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::SolutionC1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sbt = SonsIds2[1];
		sb1 = SonsIds2[0];
	}
	else // the non transferred son as index 1
	{
		sbt = SonsIds2[0];
		sb1 = SonsIds2[1];
	}



	//CASETYPE 1 : NodeId1 was the first sout
	double c1a1bp = getC1(sa1, NodeId2);
	AdjScore Sc1a1bp(true,sa1, NodeId2);
	double c0a1bp = getC0(sa1, NodeId2);
	AdjScore Sc0a1bp(true,sa1, NodeId2);


	double c1atbp = getC1(sat, NodeId2);
	AdjScore Sc1atbp(true,sat, NodeId2);
	double c0atbp = getC0(sat, NodeId2);
	AdjScore Sc0atbp(true,sat, NodeId2);



	vector <AdjSolution> Vsolution;
	AdjSolution *currentSolution;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept -> 2 gains
	currentSolution = new AdjSolution(2,0,iscoevent);
	currentSolution->components.push_back( Sc1a1bp ); // 1 gain
	//currentSolution->components.push_back( Sc1atbt );
	currentSolution->components.push_back( Sc1atbp ); // 1 gain
	caseScore.push_back( c1a1bp );
	//caseScore.push_back( c1atbt );
	caseScore.push_back( c1atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 2 -> only the adjacencies in the same species has been kept -> 1 gain
	currentSolution = new AdjSolution(1,0,iscoevent);
	currentSolution->components.push_back( Sc1a1bp ); //1 gain
	//currentSolution->components.push_back( Sc0atbt );
	currentSolution->components.push_back( Sc0atbp );
	caseScore.push_back( c1a1bp );
	//caseScore.push_back( c0atbt );
	caseScore.push_back( c0atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 3 -> only the adjacency in transfered sons has been kept -> 1 gain and 0 break 
	currentSolution = new AdjSolution(1,0,iscoevent);
	currentSolution->components.push_back( Sc0a1bp );
	currentSolution->components.push_back( Sc1atbp ); // 1 gain
	caseScore.push_back( c0a1bp );
	caseScore.push_back( c1atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();


	//case 4 -> no adjacencies have been kept -> 0 gain , 0 break
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back( Sc0a1bp ); 
	currentSolution->components.push_back( Sc0atbp );
	caseScore.push_back( c0a1bp );
	caseScore.push_back( c0atbp );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//CASETYPE 2 : NodeId2 was the first sout
	double c1apb1 = getC1(NodeId1, sb1);
	AdjScore Sc1apb1(true,NodeId1, sb1);
	double c0apb1 = getC0(NodeId1, sb1);
	AdjScore Sc0apb1(true,NodeId1, sb1);
	double c1apbt = getC1(NodeId1, sbt);
	AdjScore Sc1apbt(true,NodeId1, sbt);
	double c0apbt = getC0(NodeId1, sbt);
	AdjScore Sc0apbt(true,NodeId1, sbt);



	//case 1 -> both adjacencies have been kept -> 2 gains
	currentSolution = new AdjSolution(2,0,iscoevent);
	currentSolution->components.push_back( Sc1apb1 ); //1 gain
	currentSolution->components.push_back( Sc1apbt ); //1 gain
	caseScore.push_back( c1apb1 );
	caseScore.push_back( c1apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 2 -> only the adjacencies in the same species has been kept -> 0 gain nor loss
	currentSolution = new AdjSolution(1,0,iscoevent);
	currentSolution->components.push_back( Sc1apb1 ); //1 gain
	currentSolution->components.push_back( Sc0apbt );
	caseScore.push_back( c1apb1 );
	caseScore.push_back( c0apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	//case 3 -> only the adjacency in transfered sons has been kept -> 1 gain 
	currentSolution = new AdjSolution(1,0,iscoevent);
	currentSolution->components.push_back( Sc0apb1 ); 
	currentSolution->components.push_back( Sc1apbt ); // 1 gain
	caseScore.push_back( c0apb1 );
	caseScore.push_back( c1apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();


	//case 4 -> no adjacencies have been kept -> 0 gain
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back( Sc0apb1 ); 
	currentSolution->components.push_back( Sc0apbt );
	caseScore.push_back( c0a1bp );
	caseScore.push_back( c0apbt );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); 

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	return Vsolution;

}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: add something later for co-event checking (a priori this formulae are good for both as long as a gain is added externally (ie. after the adj tree is built))
*/
vector<AdjSolution> AdjMatrix::SolutionC1RecWithSout(int NodeId1, int NodeId2, bool firstSout)
{
	vector<AdjSolution> Vsolution = SolutionC1RecWithOther(NodeId1, NodeId2, !firstSout);

	vector<AdjSolution> tmp = SolutionC1SoutWithExtantOrSpecOrNull(NodeId1, NodeId2, firstSout);	

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	return Vsolution;
} 


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

NB: add something later for co-event checking (a priori this formulae are good for both as long as a gain is added externally (ie. after the adj tree is built))
*/
vector<AdjSolution> AdjMatrix::SolutionC0RecWithSout(int NodeId1, int NodeId2, bool firstSout)
{
	vector<AdjSolution> Vsolution = SolutionC0RecWithOther(NodeId1, NodeId2, !firstSout);

	vector<AdjSolution> tmp = SolutionC0SoutWithExtantOrSpecOrNull(NodeId1, NodeId2, firstSout);	

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	return Vsolution;
} 


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

NB: add something later for co-event checking (a priori this formulae are good for both as long as a gain is added externally (ie. after the adj tree is built))
*/
vector<AdjSolution> AdjMatrix::SolutionC1RecWithRec(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1RecWithRec"<< endl;

	vector<AdjSolution> Vsolution;
	AdjSolution *currentSolution;

	bool iscoevent = true;

	//both reception node only a have  child
	int sonId1 = Rtree1.getSonsId(NodeId1)[0];
	int sonId2 = Rtree2.getSonsId(NodeId2)[0];

	double c1 = getC1(sonId1 , sonId2);
	AdjScore Sc1(true ,sonId1, sonId2);
	double c0 = getC0(sonId1 , sonId2);
	AdjScore Sc0(false,sonId1, sonId2);

	vector <double > caseScore;

	//two cases
	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back(Sc1);
	caseScore.push_back(c1);
	if(WeightedHgtCost != bestScore)
		caseScore.push_back(WeightedHgtCost);

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	currentSolution = new AdjSolution(0,1,iscoevent);// C0 under -> one break
	currentSolution->components.push_back(Sc0);
	caseScore.push_back(c0);
	if(WeightedHgtCost != bestScore)
		caseScore.push_back(WeightedHgtCost);

	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

NB: the synchronous case recquires the node to be linked
*/
vector<AdjSolution> AdjMatrix::SolutionC0RecWithRec(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0RecWithRec"<< endl;



	vector<AdjSolution> Vsolution;
	AdjSolution *currentSolution;

	bool iscoevent = false;

	//both reception node only a have  child
	int sonId1 = Rtree1.getSonsId(NodeId1)[0];
	int sonId2 = Rtree2.getSonsId(NodeId2)[0];

	double c1 = getC1(sonId1 , sonId2);
	AdjScore Sc1(true ,sonId1, sonId2);
	double c0 = getC0(sonId1 , sonId2);
	AdjScore Sc0(false,sonId1, sonId2);

	vector <double > caseScore;

	//two cases
	currentSolution = new AdjSolution(1,0,iscoevent); //C1 under -> 1 gain
	currentSolution->components.push_back(Sc1);
	caseScore.push_back(c1);
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	currentSolution = new AdjSolution(0,0,iscoevent);
	currentSolution->components.push_back(Sc0);
	caseScore.push_back(c0);
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();

	// asynchronous cases: 
	vector<AdjSolution> tmp = SolutionC0RecWithOther( NodeId1, NodeId2, true);
	*currentSolution = tmp[0];

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);	



	tmp = SolutionC0RecWithOther( NodeId1, NodeId2, false);

	*currentSolution = tmp[0];

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	//cout << "AdjMatrix::SolutionC0RecWithRec " << NodeId1 << "-" << NodeId2 << "  " << Vsolution.size() << endl;

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstRec (bool): wether the Rec is NodeId1 or NodeId2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacencies between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC1RecWithOther(int NodeId1, int NodeId2, bool firstRec)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1RecWithOther"<< endl;

	pair <int,int> orderedIds;

	if(firstRec)//use the son of the reception for the transmission
	{
		orderedIds.first = Rtree1.getSonsId(NodeId1)[0];
		orderedIds.second = NodeId2;
	}
	else
	{
		orderedIds.first = NodeId1;
		orderedIds.second = Rtree2.getSonsId(NodeId2)[0];
	}

	//now that we have the correct son, we can apply the formula
	//only one case, with no gain or break -> simple transmission

	vector<AdjSolution> Vsolution;

	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().components.push_back(AdjScore(true, orderedIds.first, orderedIds.second));
	Vsolution.back().score = getC1(orderedIds.first, orderedIds.second);

	return Vsolution;

}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstRec (bool): wether the Rec is NodeId1 or NodeId2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacencies between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionC0RecWithOther(int NodeId1, int NodeId2, bool firstRec)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0RecWithOther"<< endl;

	pair <int,int> orderedIds;

	if(firstRec)//use the son of the reception for the transmission
	{
		orderedIds.first = Rtree1.getSonsId(NodeId1)[0];
		orderedIds.second = NodeId2;
	}
	else
	{
		orderedIds.first = NodeId1;
		orderedIds.second = Rtree2.getSonsId(NodeId2)[0];
	}

	//now that we have the correct son, we can apply the formula
	//only one case, with no gain or break -> simple transmission

	vector<AdjSolution> Vsolution;

	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().components.push_back(AdjScore(false, orderedIds.first, orderedIds.second));
	Vsolution.back().score = getC0(orderedIds.first, orderedIds.second);

	return Vsolution;

}

/*

Takes:
 - NodeIdBout (int): node id in RtreeBout
 - NodeIdOther (int): node id in RtreeOther
 - RtreeBout (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the Bifurcation Out is in the tree1 (true) or not (false)

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeIdBout and NodeIdOther
*/
vector<AdjSolution> AdjMatrix::SolutionB1(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert)
{
	vector<AdjSolution> Vsolution;
	bool iscoevent = false;

	//we first get the ids of the children of NodeIdDup
	vector <int> BoutSonsId = RtreeBout->getSonsId(NodeIdBout);

	//There is 4 possible cases
	vector <double> caseScore;

	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab1.id1 = NodeIdOther;
		Sc1ab1.id2 = BoutSonsId[0];
	}
	else
	{
		Sc1ab1.id1 = BoutSonsId[0];
		Sc1ab1.id2 = NodeIdOther;
	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);



	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab2.id1 = NodeIdOther;
		Sc1ab2.id2 = BoutSonsId[1];
	}
	else
	{
		Sc1ab2.id1 = BoutSonsId[1];
		Sc1ab2.id2 = NodeIdOther;
	}
	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);

	

	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab1.id1 = NodeIdOther;
		Sc0ab1.id2 = BoutSonsId[0];
	}
	else
	{
		Sc0ab1.id1 = BoutSonsId[0];
		Sc0ab1.id2 = NodeIdOther;
	}
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab2.id1 = NodeIdOther;
		Sc0ab2.id2 = BoutSonsId[1];
	}
	else
	{
		Sc0ab2.id1 = BoutSonsId[1];
		Sc0ab2.id2 = NodeIdOther;
	}
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);


	//case 1 -> child 1 of the bout conserved an adjacency but not child 2
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 2 -> child 2 of the bout conserved an adjacency but not child 1
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab2 );
	currentSolution->components.push_back( Sc0ab1 );
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 1 of the bout conserved an adjacency and child 2 too -> need a gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency -> no break because we are in the dead
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();

	return Vsolution;
}

/*

Takes:
 - NodeIdBout (int): node id in RtreeBout
 - NodeIdOther (int): node id in RtreeOther
 - RtreeBout (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the Bifurcation Out is in the tree1 (true) or not (false)

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeIdBout and NodeIdOther
*/
vector<AdjSolution> AdjMatrix::SolutionB0(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert)
{
	vector<AdjSolution> Vsolution;
	bool iscoevent = false;

	//we first get the ids of the children of NodeIdDup
	vector <int> BoutSonsId = RtreeBout->getSonsId(NodeIdBout);

	//There is 4 possible cases
	vector <double> caseScore;

	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab1.id1 = NodeIdOther;
		Sc1ab1.id2 = BoutSonsId[0];
	}
	else
	{
		Sc1ab1.id1 = BoutSonsId[0];
		Sc1ab1.id2 = NodeIdOther;
	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);



	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc1ab2.id1 = NodeIdOther;
		Sc1ab2.id2 = BoutSonsId[1];
	}
	else
	{
		Sc1ab2.id1 = BoutSonsId[1];
		Sc1ab2.id2 = NodeIdOther;
	}
	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);

	

	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab1.id1 = NodeIdOther;
		Sc0ab1.id2 = BoutSonsId[0];
	}
	else
	{
		Sc0ab1.id1 = BoutSonsId[0];
		Sc0ab1.id2 = NodeIdOther;
	}
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);

	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	if(invert)
	{
		Sc0ab2.id1 = NodeIdOther;
		Sc0ab2.id2 = BoutSonsId[1];
	}
	else
	{
		Sc0ab2.id1 = BoutSonsId[1];
		Sc0ab2.id2 = NodeIdOther;
	}
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);

//	cout << "Bout case0 " << NodeIdOther << " " << BoutSonsId[0] << " " << BoutSonsId[1] << " " << invert <<endl;
//	cout << "Bout case0 " << c1ab1 << " " << c1ab2 << " " << c0ab1 << " " << c0ab2 << endl;

	//case 1 -> child 1 of the bout conserved an adjacency but not child 2 -> 1 gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 2 -> child 2 of the bout conserved an adjacency but not child 1 -> 1 gain
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab2 );
	currentSolution->components.push_back( Sc0ab1 );
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 1 of the bout conserved an adjacency and child 2 too -> need 2 gain
	currentSolution = new AdjSolution(2,0,iscoevent) ;  // 2 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();

	return Vsolution;
}


/*
Both nodes are Bifurcation Out that are presumed simultaneous

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2
*/
vector<AdjSolution> AdjMatrix::SolutionB12(int NodeId1, int NodeId2)
{
	vector<AdjSolution> Vsolution;
	bool iscoevent = true;

	AdjSolution * currentSolution;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId1 = Rtree1.getSonsId(NodeId1);
	vector <int> DupSonsId2 = Rtree2.getSonsId(NodeId2);


	//There is 16 possible cases representing all c1 and c0 combination of the sons of 1 an 2
	vector <double> caseScore;

	for(unsigned c1a1b1 = 0; c1a1b1 < 2; c1a1b1 ++)
			for(unsigned c1a1b2 = 0; c1a1b2 < 2; c1a1b2 ++)
					for(unsigned c1a2b1 = 0; c1a2b1 < 2; c1a2b1 ++)
							for(unsigned c1a2b2 = 0; c1a2b2 < 2; c1a2b2 ++)
							{
								currentSolution = new AdjSolution();
								currentSolution->coevent = iscoevent;

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[0],DupSonsId2[0]));
								if(!c1a1b1)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[0]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[0]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[0],DupSonsId2[1]));
								if(!c1a1b2)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[1]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[1]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[1],DupSonsId2[0]));
								if(!c1a2b1)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[0]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[0]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[1],DupSonsId2[1]));
								if(!c1a2b2)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[1]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[1]));

								//determining the number of gains and breaks
								int nbadj = c1a1b1 + c1a1b2 + c1a2b1 + c1a2b2;
								switch( nbadj )
								{

									case 0: // none linked -> 0 break because we are in the dead

										currentSolution->NbGain = 0;
										currentSolution->NbBreak = 0;
										break;

									case 1:// only one linked -> 0 break because we are in the dead
										currentSolution->NbGain = 0;
										currentSolution->NbBreak = 0;
										break;
									case 3: // three linked -> 1 gain
										currentSolution->NbGain = 1;
										currentSolution->NbBreak = 0;	
										break;
									case 4: // all linked -> 2 gain
										currentSolution->NbGain = 2;
										currentSolution->NbBreak = 0;	
										break;
									default: //exactly 2 adjs
										currentSolution->NbGain = 0;
										currentSolution->NbBreak = 0;
										if(( c1a1b1 == c1a1b2 ) || (c1a1b1 == c1a2b1 )) // scenarii where 1 child has two adjs -> one gain and no break because we are in the dead
										{
											currentSolution->NbGain = 1;
											currentSolution->NbBreak = 0;
										}
								}
								//getting the new score
								currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);
								
								if(Vsolution.empty())
									Vsolution.push_back(*currentSolution);
								else
								{
									if(currentSolution->score <= Vsolution.front().score)
										Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
									else
										Vsolution.push_back(*currentSolution);
								}

								//clearing
								delete currentSolution;
								caseScore.clear();
							}

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1BoutWithBout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1BoutWithBout"<< endl;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	vector<AdjSolution> Vsolution,tmp;
	//B1 -> the Bout of NodeId1 was before the Bout of NodeId2
	Vsolution = SolutionB1(NodeId1, NodeId2, rtree1, rtree2, false);

	//B2 -> the Bout of NodeId2 was before the Bout of NodeId1
	tmp = SolutionB1(NodeId2, NodeId1, rtree2, rtree1, true);

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	//B12 -> both Bout were simultaneous
	tmp = SolutionB12(NodeId1,NodeId2);

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	return Vsolution;
}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC0BoutWithBout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0BoutWithBout"<< endl;


	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	vector<AdjSolution> Vsolution,tmp;
	//B1 -> the Bout of NodeId1 was before the Bout of NodeId2
	Vsolution = SolutionB0(NodeId1, NodeId2, rtree1, rtree2, false);

	//B2 -> the Bout of NodeId2 was before the Bout of NodeId1
	tmp = SolutionB0(NodeId2, NodeId1, rtree2, rtree1, true);

	for(unsigned i = 0; i < tmp.size();i++)
	{
		if(tmp[i].score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),tmp[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(tmp[i]);

	}

	return Vsolution;
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstBout (bool): whether nodeId1 is the Bout  (true) or not (false)

Returns:
	(vector<AdjSolution>): the different solutions where there is an adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC1BoutWithRec(int NodeId1, int NodeId2, bool firstBout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1BoutWithRec"<< endl;


	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;



	if(firstBout)	// the Bout is NodeId1
		return SolutionB1(NodeId1, NodeId2, rtree1, rtree2, !firstBout);

	// the Bout is NodeId2
	return SolutionB1(NodeId2, NodeId1, rtree2, rtree1, !firstBout);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstBout (bool): whether nodeId1 is the Bout  (true) or not (false)

Returns:
	(vector<AdjSolution>): the different solutions where there is no adjacency between NodeId1 and NodeId2

*/
vector<AdjSolution> AdjMatrix::SolutionC0BoutWithRec(int NodeId1, int NodeId2, bool firstBout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0BoutWithRec"<< endl;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	if(firstBout)	// the Bout is NodeId1
		return SolutionB0(NodeId1, NodeId2, rtree1, rtree2, !firstBout);

	// the Bout is NodeId2
	return SolutionB0(NodeId2, NodeId1, rtree2, rtree1, !firstBout);
}

/*
Returns:
	(vector<AdjSolution>): a vector containing an unique solution with worst score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC1DefaultImpossibleCase()
{
	vector <AdjSolution> Vsolution;
	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().score = worstScore;

	return Vsolution;
}



/*
Returns:
	(vector<AdjSolution>): a vector containing an unique solution with best score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC0DefaultImpossibleCase()
{
	vector <AdjSolution> Vsolution;

	Vsolution.push_back( AdjSolution(0,0,false) );
	Vsolution.back().score = bestScore;

	return Vsolution;
}


/*
Function specific to the case where the nodes have only one child
Two snodes can only be synchronous if they have the same event, which also means that they have the same number of children

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool inTheDead: whether 

Returns:
	(vector<AdjSolution>): a vector containing an unique solution with best score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC1synchronousOneChild(int NodeId1, int NodeId2, bool inTheDead)
{

	//1. determining coevent status
	bool iscoevent = true; // synchronous implies that this is a co-event


	vector <AdjSolution> Vsolution;

	vector <int> SonsIds1 = Rtree1.getSonsId(NodeId1); //should be size 1
	vector <int> SonsIds2 = Rtree2.getSonsId(NodeId2); //should be size 1


	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1a1b1;
	Sc1a1b1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0a1b1;
	Sc0a1b1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)

	Sc1a1b1.id1 = SonsIds1[0];
	Sc1a1b1.id2 = SonsIds2[0];

	Sc0a1b1.id1 = SonsIds1[0];
	Sc0a1b1.id2 = SonsIds2[0];

	double c1a1b1 = getC1(Sc1a1b1.id1,Sc1a1b1.id2);
	double c0a1b1 = getC0(Sc0a1b1.id1,Sc0a1b1.id2);


	int nbBreak = 0;

	// 2 possible solutions
	//case 1 -> the adj is not present >> break
	if(!inTheDead)
		nbBreak = 1;


	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 1 break if in the dead
	currentSolution->components.push_back( Sc0a1b1 );
	caseScore.push_back( c0a1b1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	Vsolution.push_back(*currentSolution);
	delete currentSolution;
	caseScore.clear();

	//case 2 -> the adj is present present << 0 gain
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1a1b1 );
	caseScore.push_back( c1a1b1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	return Vsolution;
}



/*
Function specific to the case where the nodes have two children
Two snodes can only be synchronous if they have the same event, which also means that they have the same number of children

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool inTheDead: whether 

Returns:
	(vector<AdjSolution>): a vector containing an unique solution with best score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC1synchronousTwoChildren(int NodeId1, int NodeId2, bool inTheDead)
{
	vector<AdjSolution> Vsolution;
	bool iscoevent = true;

	AdjSolution * currentSolution;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId1 = Rtree1.getSonsId(NodeId1);
	vector <int> DupSonsId2 = Rtree2.getSonsId(NodeId2);


	//There is 16 possible cases representing all c1 and c0 combination of the sons of 1 an 2
	vector <double> caseScore;

	for(unsigned c1a1b1 = 0; c1a1b1 < 2; c1a1b1 ++)
			for(unsigned c1a1b2 = 0; c1a1b2 < 2; c1a1b2 ++)
					for(unsigned c1a2b1 = 0; c1a2b1 < 2; c1a2b1 ++)
							for(unsigned c1a2b2 = 0; c1a2b2 < 2; c1a2b2 ++)
							{
								currentSolution = new AdjSolution();
								currentSolution->coevent = iscoevent;

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[0],DupSonsId2[0]));
								if(!c1a1b1)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[0]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[0]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[0],DupSonsId2[1]));
								if(!c1a1b2)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[1]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[1]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[1],DupSonsId2[0]));
								if(!c1a2b1)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[0]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[0]));

								currentSolution->components.push_back(AdjScore(true,DupSonsId1[1],DupSonsId2[1]));
								if(!c1a2b2)
								{
									currentSolution->components.back().C1 = false;
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[1]));
								}
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[1]));

								//determining the number of gains and breaks
								int nbadj = c1a1b1 + c1a1b2 + c1a2b1 + c1a2b2;
								switch( nbadj )
								{

									case 0: // none linked -> 2 break unless we are in the dead

										currentSolution->NbGain = 0;
										if(inTheDead)
											currentSolution->NbBreak = 0;
										else
											currentSolution->NbBreak = 2;
										break;

									case 1:// only one linked -> 1 break unles we are in the dead
										currentSolution->NbGain = 0;
										if(inTheDead)
											currentSolution->NbBreak = 0;
										else
											currentSolution->NbBreak = 1;
										break;
									case 3: // three linked -> 1 gain
										currentSolution->NbGain = 1;
										currentSolution->NbBreak = 0;	
										break;
									case 4: // all linked -> 2 gain
										currentSolution->NbGain = 2;
										currentSolution->NbBreak = 0;	
										break;
									default: //exactly 2 adjs
										currentSolution->NbGain = 0;
										currentSolution->NbBreak = 0;
										if(( c1a1b1 == c1a1b2 ) || (c1a1b1 == c1a2b1 )) // scenarii where 1 child has two adjs -> one gain and 1 break unless we are in the dead
										{
											currentSolution->NbGain = 1;
											if(inTheDead)
												currentSolution->NbBreak = 0;
											else
												currentSolution->NbBreak = 1;
										}
								}
								//getting the new score
								currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak);
								
								if(Vsolution.empty())
									Vsolution.push_back(*currentSolution);
								else
								{
									if(currentSolution->score <= Vsolution.front().score)
										Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
									else
										Vsolution.push_back(*currentSolution);
								}

								//clearing
								delete currentSolution;
								caseScore.clear();
							}

	return Vsolution;
}

/*
Function specific to the case where the more ancient node has only one child

Takes:
 - int NodeId1: an id in the first tree
 - int NodeId2: an id in the second tree
 - bool NodeId2First: true if nodeid2 is before nodeid1
 - bool inTheDead: whether 

Returns:
	(vector<AdjSolution>): a vector containing an unique solution with best score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC1asynchronousOneChild(int NodeId1, int NodeId2, bool NodeId2First, bool inTheDead)
{
	vector <AdjSolution> Vsolution;

	vector <int> SonsIds;

	if(!NodeId2First) // node 1 before
	 	SonsIds = Rtree1.getSonsId(NodeId1);
	else // node 2 before
		SonsIds = Rtree2.getSonsId(NodeId2);

	bool iscoevent = false; // never a co-event

	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	if(NodeId2First) // node 2 before
	{
		Sc1ab1.id1 = NodeId1;
		Sc1ab1.id2 = SonsIds[0];

		Sc0ab1.id1 = NodeId1;
		Sc0ab1.id2 = SonsIds[0];

	}
	else // node 1 before
	{
		Sc1ab1.id1 = SonsIds[0];
		Sc1ab1.id2 = NodeId2;

		Sc0ab1.id1 = SonsIds[0];
		Sc0ab1.id2 = NodeId2;


	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);


	int nbBreak = 0;

	// 2 possible solutions
	//case 1 -> the adj is not present << nothing happens
	if(!inTheDead)
		nbBreak = 1;


	currentSolution = new AdjSolution(0,nbBreak,iscoevent) ;  // 0 gain 1 break if in the dead
	currentSolution->components.push_back( Sc0ab1 );
	caseScore.push_back( c0ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	Vsolution.push_back(*currentSolution);
	delete currentSolution;
	caseScore.clear();

	//case 2 -> the adj is present present << 0 gain
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	caseScore.push_back( c1ab1 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution
	
	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);


	delete currentSolution;
	caseScore.clear();
	return Vsolution;

}

/*
Returns:
	(vector<AdjSolution>): a vector containing an unique solution with best score and no components

*/
vector<AdjSolution> AdjMatrix::SolutionC0DefaultImpossibleCaseNew(int NodeId1, int NodeId2)
{
	vector <AdjSolution> Vsolution;

	//NodeId1 before
	vector <AdjSolution> Vsolution1 = SolutionC0DefaultImpossibleCaseAux( NodeId1, NodeId2, false);

	for(unsigned i = 0; i < Vsolution1.size();i++)
	{
		Vsolution.push_back(Vsolution1[i]);
		//cout << "AdjMatrix::SolutionC0DefaultImpossibleCaseNew " << NodeId1 << "-" << NodeId2 << " n1 before " << Vsolution1[i].score << endl;
	}




	//NodeId2 before
	vector <AdjSolution> Vsolution2 = SolutionC0DefaultImpossibleCaseAux( NodeId1, NodeId2, true );

	for(unsigned i = 0; i < Vsolution2.size();i++)
	{
		if(Vsolution.size() != 0)
		{
			if(Vsolution2[i].score <= Vsolution.front().score)
			{
				Vsolution.insert(Vsolution.begin(),Vsolution2[i]);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
			}
		}
		else
		{
			Vsolution.push_back(Vsolution2[i]);
		}
		//cout << "AdjMatrix::SolutionC0DefaultImpossibleCaseNew " << NodeId1 << "-" << NodeId2 << " n2 before " << Vsolution2[i].score << endl;
	}

	if(Vsolution.size() == 0) // leaf|loss  with  leaf|loss
		Vsolution = SolutionC0DefaultImpossibleCase();

	return Vsolution;
}



vector<AdjSolution> AdjMatrix::SolutionC0DefaultImpossibleCaseAux(int NodeId1, int NodeId2, bool invert)
{
	vector <AdjSolution> Vsolution;
	
	vector <int> SonsIds;

	if(!invert) // node 1 before
	 	SonsIds = Rtree1.getSonsId(NodeId1);
	else // node 2 before
		SonsIds = Rtree2.getSonsId(NodeId2);


	if(SonsIds.size() == 0) // cases of losses and leaves
		return Vsolution;



	bool iscoevent = false; // never a co-event

	vector <double> caseScore;
	AdjSolution * currentSolution;
	
	//getting the different scores that we need
	AdjScore Sc1ab1;
	Sc1ab1.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab1;
	Sc0ab1.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	AdjScore Sc1ab2;
	Sc1ab2.C1 = true;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)
	AdjScore Sc0ab2;
	Sc0ab2.C1 = false;//true if the score is in C1 (adjacency presence); false if it is in C0 (adjacency absence)


	if(invert) // node 2 before
	{
		Sc1ab1.id1 = NodeId1;
		Sc1ab1.id2 = SonsIds[0];

		Sc0ab1.id1 = NodeId1;
		Sc0ab1.id2 = SonsIds[0];

	}
	else // node 1 before
	{
		Sc1ab1.id1 = SonsIds[0];
		Sc1ab1.id2 = NodeId2;

		Sc0ab1.id1 = SonsIds[0];
		Sc0ab1.id2 = NodeId2;


	}
	double c1ab1 = getC1(Sc1ab1.id1,Sc1ab1.id2);
	double c0ab1 = getC0(Sc0ab1.id1,Sc0ab1.id2);


	if(SonsIds.size() == 1) // ONLY one son 
	{
		// 2 possible solutions

		//case 1 -> the adj is not present << nothing happens
		currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
		currentSolution->components.push_back( Sc0ab1 );
		caseScore.push_back( c0ab1 );
		currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

		Vsolution.push_back(*currentSolution);

		delete currentSolution;
		caseScore.clear();

		//case 2 -> the adj is present present << 1 gain
		currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
		currentSolution->components.push_back( Sc1ab1 );
		caseScore.push_back( c1ab1 );
		currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

		if(currentSolution->score <= Vsolution.front().score)
			Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
		else
			Vsolution.push_back(*currentSolution);


		delete currentSolution;
		caseScore.clear();

		return Vsolution;
	}


	/// two sons (more than that is not covered yet)

	if(invert) // node 2 before
	{
		Sc1ab2.id1 = NodeId1;
		Sc1ab2.id2 = SonsIds[1];

		Sc0ab2.id1 = NodeId1;
		Sc0ab2.id2 = SonsIds[1];
	}
	else // ndoe 1 before
	{
		Sc1ab2.id1 = SonsIds[1];
		Sc1ab2.id2 = NodeId2;

		Sc0ab2.id1 = SonsIds[1];
		Sc0ab2.id2 = NodeId2;
	}

	double c1ab2 = getC1(Sc1ab2.id1,Sc1ab2.id2);
	double c0ab2 = getC0(Sc0ab2.id1,Sc0ab2.id2);


	///4 cases
	//case 1 -> neither child 1 nor child 2 kept an adjacency 
	currentSolution = new AdjSolution(0,0,iscoevent) ;  // 0 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 2 -> child 1  kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc0ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution


	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 3 -> child 2 kept an adjacency 
	currentSolution = new AdjSolution(1,0,iscoevent) ;  // 1 gain 0 break 
	currentSolution->components.push_back( Sc0ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	//case 4 -> child 1 and child 2 kept an adjacency
	currentSolution = new AdjSolution(2,0,iscoevent) ;  // 2 gain 0 break 
	currentSolution->components.push_back( Sc1ab1 );
	currentSolution->components.push_back( Sc1ab2 );
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	currentSolution->score = AggregateScore(caseScore,currentSolution->NbGain,currentSolution->NbBreak); // the score of that solution

	if(currentSolution->score <= Vsolution.front().score)
		Vsolution.insert(Vsolution.begin(),*currentSolution);//the new score is lower or equal to the smallest score --> insert in first position so that the first element of Vsolution is the one with lowest score
	else
		Vsolution.push_back(*currentSolution);

	delete currentSolution;
	caseScore.clear();
	
	return Vsolution;
}



/////////////////////////////////////////////////////////////////////////////////////////
/////////// OLD VERSION OF COST FUNCTION THAT ONLY RETURNS DOUBLE ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////



/*
Set C1 and C0 for NodeId1 , NodeId2

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

*/
pair <double, double> AdjMatrix::computeScore(int NodeId1, int NodeId2)
{


	Node * n1 = Rtree1.getNode(NodeId1);
	Node * n2 = Rtree2.getNode(NodeId2);



	pair <double, double> scores;

	// The choice is made that the TS compatibility is set by the reconciled trees rather than through the Adj algorithm.
	// This, for instance, allows a simple DeCo algorithm to perform on a full time sliced tree and respecting TS contraints

	if(!Rtree1.areTSCompatible(n1,n2))//if the node aren't TS compatible: their adjacency is impossible
	{
		if(verbose)
		{
			cout << "Nodes "<< NodeId1 << " and " << NodeId1 << " are not time compatible. Setting values accordingly." << endl;
		}

		scores.first = worstScore;
		scores.second = bestScore;
		return scores;
	}
	// both node are TS compatible

	int evt1 = Rtree1.getNodeEvent(NodeId1);
	int evt2 = Rtree2.getNodeEvent(NodeId2);



	if(decoLTalgo)
	{//check the cases where the species can differ in DeCoLT: both Bout or Rec and different events
		if(evt1 != evt2)
		{
			bool boutwithrec = false;
			bool firstBout;
			if((Rtree1.isRec(evt1)) && (Rtree2.isBout(evt2)))
			{
				boutwithrec = true;
				firstBout = false;
			}
			else if((Rtree1.isBout(evt1)) && (Rtree2.isRec(evt2)))
			{
				boutwithrec = true;
				firstBout = true;
			}

			if(boutwithrec)
			{
				scores.first =  C1BoutWithRec(NodeId1,NodeId2, firstBout);
				scores.second = C0BoutWithRec(NodeId1,NodeId2, firstBout);
				return scores ;
			}
		}
	}


	//now the species compatibility need to hold
	if(!Rtree1.haveSameSpecies(n1, n2))
	{
		if(verbose)
		{
			cout << "Nodes "<< NodeId1 << " and " << NodeId2 << " have different species. Setting values accordingly." << endl;
		}

		scores.first = worstScore;
		scores.second = bestScore;
		return scores;
	}


	//cases where species is the same

	bool firstLoss = Rtree1.isLoss(evt1);
	bool secondLoss = Rtree2.isLoss(evt2);

	if( (firstLoss) || (secondLoss) )
	{
		if( (firstLoss) && (secondLoss) )
		{
			scores.first =  C1LossWithLoss(NodeId1,NodeId2);
			scores.second = C0LossWithLoss(NodeId1,NodeId2);
			return scores; 
		}
		else
		{
			scores.first =  C1LossWithOther(NodeId1,NodeId2);
			scores.second = C0LossWithOther(NodeId1,NodeId2);
			return scores;
		}
	}



	bool firstDup = Rtree1.isDup(evt1);
	bool secondDup = Rtree2.isDup(evt2);

	if( (firstDup) || (secondDup) )
	{
		if( (firstDup) && (secondDup) )
		{
			scores.first =  C1DupWithDup(NodeId1,NodeId2);
			scores.second = C0DupWithDup(NodeId1,NodeId2);
			return scores ;
		}
		else
		{
			scores.first =  C1DupWithOther(NodeId1,NodeId2,firstDup);
			scores.second = C0DupWithOther(NodeId1,NodeId2,firstDup);
			return scores;
		}
	}


	if(decoLTalgo) // Sout specific part
	{
		bool firstSout = Rtree1.isSout(evt1);
		bool secondSout = Rtree2.isSout(evt2);

		if( (firstSout) || (secondSout) )
		{
			if( (firstSout) && (secondSout) )
			{
				scores.first =  C1SoutWithSout(NodeId1,NodeId2);
				scores.second = C0SoutWithSout(NodeId1,NodeId2);
				return scores;
			}
			else
			{
				scores.first =  C1SoutWithExtantOrSpecOrNull(NodeId1,NodeId2,firstSout);
				scores.second = C0SoutWithExtantOrSpecOrNull(NodeId1,NodeId2,firstSout);
				return scores;
			}
		}
	}


	//Now: cases where both event must be equal
	if(evt1 == evt2)
	{
		if(Rtree1.isExtant(evt1))
		{
			scores.first =  C1ExtantWithExtant(NodeId1,NodeId2);
			scores.second = C0ExtantWithExtant(NodeId1,NodeId2);
			return scores;	
		}
		else if(Rtree1.isSpeciation(evt1))
		{
			scores.first =  C1SpecWithSpec(NodeId1,NodeId2);
			scores.second = C0SpecWithSpec(NodeId1,NodeId2);
			return scores;	
		}
		else if(Rtree1.isNull(evt1))
		{
			scores.first =  C1NullWithNull(NodeId1,NodeId2);
			scores.second = C0NullWithNull(NodeId1,NodeId2);
			return scores;
		}
		else if(decoLTalgo)
		{//bout with bout and rec with rec cases
			if(Rtree1.isBout(evt1))
			{
				scores.first =  C1BoutWithBout(NodeId1,NodeId2);
				scores.second = C0BoutWithBout(NodeId1,NodeId2);
				return scores;		
			}
			else if(Rtree1.isRec(evt1))
			{
				scores.first =  C1RecWithRec(NodeId1,NodeId2);
				scores.second = C0RecWithRec(NodeId1,NodeId2);
				return scores;
			}
		}
	}


	//we are out of the case covered -> return some default value
	if(verbose)
	{
		cout << "Nodes "<< NodeId1 << " (" << evt1 << ")" << " and " << NodeId2 << " (" << evt2 << ")" << " do not fall in any cases. Setting values accordingly." << endl;
	}
	scores.first = worstScore;//impossible adjacency
	scores.second = bestScore;
	return scores;
}





/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

*/
double AdjMatrix::C1ExtantWithExtant(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1ExtantWithExtant"<< endl;
	return getC1(NodeId1,NodeId2);
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2

*/
double AdjMatrix::C0ExtantWithExtant(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0ExtantWithExtant"<< endl;

	return getC0(NodeId1,NodeId2);
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
*/
double AdjMatrix::C1LossWithOther(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1LossWithOther"<< endl;
	return bestScore;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the loss; Other may be anything but a loss
*/
double AdjMatrix::C0LossWithOther(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0LossWithOther"<< endl;
	return bestScore;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

NB: add something later to account for possible co-events
*/
double AdjMatrix::C1LossWithLoss(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1LossWithLoss"<< endl;
	return bestScore;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2

NB: add something later to account for possible co-events
*/
double AdjMatrix::C0LossWithLoss(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0LossWithLoss"<< endl;
	return bestScore;
}

/*
Takes:
 - NodeIdDup (int): node id in RtreeDup
 - NodeIdOther (int): node id in RtreeOther
 - RtreeDup (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the dupplication is in the tree1 (true) or not (false)

Returns:
	(double): the cost of having an adjacency between NodeIdDup and NodeIdOther
*/
double AdjMatrix::D1(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert )
{
	vector <double> scorevector;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId = RtreeDup->getSonsId(NodeIdDup);

	//There is 4 possible cases
	vector <double> caseScore;
	
	double c1ab1 = getC1(DupSonsId[0],NodeIdOther,invert);
	double c1ab2 = getC1(DupSonsId[1],NodeIdOther,invert);
	double c0ab1 = getC0(DupSonsId[0],NodeIdOther,invert);
	double c0ab2 = getC0(DupSonsId[1],NodeIdOther,invert);

	//case 1 -> child 1 of the duplication conserved an adjacency but not child 2
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();
	//case 2 -> child 2 of the duplication conserved an adjacency but not child 1
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();
	//case 3 -> child 1 of the duplication conserved an adjacency and child 2 too -> need a gain
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency -> need a break
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 1 break


	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*
Takes:
 - NodeIdDup (int): node id in RtreeDup
 - NodeIdOther (int): node id in RtreeOther
 - RtreeDup (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the dupplication is in the tree1 (true) or not (false)

Returns:
	(double): the cost of having no adjacency between NodeIdDup and NodeIdOther
*/
double AdjMatrix::D0(int NodeIdDup, int NodeIdOther, ReconciledTree * RtreeDup, ReconciledTree * Rtreeother,bool invert )
{
	vector <double> scorevector;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId = RtreeDup->getSonsId(NodeIdDup);

	//There is 4 possible cases
	vector <double> caseScore;
	
	double c1ab1 = getC1(DupSonsId[0],NodeIdOther,invert);
	double c1ab2 = getC1(DupSonsId[1],NodeIdOther,invert);
	double c0ab1 = getC0(DupSonsId[0],NodeIdOther,invert);
	double c0ab2 = getC0(DupSonsId[1],NodeIdOther,invert);

	//case 1 -> child 1 of the duplication conserved an adjacency but not child 2 -> need a gain
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();
	//case 2 -> child 2 of the duplication conserved an adjacency but not child 1 -> need a gain
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();
	//case 3 -> child 1 of the duplication conserved an adjacency and child 2 too -> need 2 gain
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	scorevector.push_back(AggregateScore(caseScore,2,0)); // 2 gain 0 break

	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break


	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*
Both nodes are duplications that are presumed simultaneous

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2
*/
double AdjMatrix::D12(int NodeId1, int NodeId2)
{
	vector <double> scorevector;

	//we first get the ids of the children of NodeIdDup
	vector <int> DupSonsId1 = Rtree1.getSonsId(NodeId1);
	vector <int> DupSonsId2 = Rtree2.getSonsId(NodeId2);


	//There is 16 possible cases representing all c1 and c0 combination of the sons of 1 an 2
	vector <double> caseScore;
	int nbGain;
	int nbBreak;

	for(unsigned c1a1b1 = 0; c1a1b1 < 2; c1a1b1 ++)
			for(unsigned c1a1b2 = 0; c1a1b2 < 2; c1a1b2 ++)
					for(unsigned c1a2b1 = 0; c1a2b1 < 2; c1a2b1 ++)
							for(unsigned c1a2b2 = 0; c1a2b2 < 2; c1a2b2 ++)
							{

								if(!c1a1b1)
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[0]));
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[0]));

								if(!c1a1b2)
									caseScore.push_back(getC0(DupSonsId1[0],DupSonsId2[1]));
								else
									caseScore.push_back(getC1(DupSonsId1[0],DupSonsId2[1]));

								if(!c1a2b1)
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[0]));
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[0]));

								if(!c1a2b2)
									caseScore.push_back(getC0(DupSonsId1[1],DupSonsId2[1]));
								else
									caseScore.push_back(getC1(DupSonsId1[1],DupSonsId2[1]));

								//determining the number of gains and breaks
								int nbadj = c1a1b1 + c1a1b2 + c1a2b1 + c1a2b2;
								switch( nbadj )
								{

									case 0: // none linked -> 2 break

										nbGain = 0;
										nbBreak = 2;
										break;

									case 1: // only one linked -> 1 break
										nbGain = 0;
										nbBreak = 1;	
										break;
									case 3: // three linked -> 1 gain
										nbGain = 1;
										nbBreak = 0;	
										break;
									case 4: // all linked -> 2 gain
										nbGain = 2;
										nbBreak = 0;	
										break;
									default: //exactly 2 adjs
										nbGain = 0;
										nbBreak = 0;
										if(( c1a1b1 == c1a1b2 ) || (c1a1b1 == c1a2b1 )) // scenarios where 1 child has two adjs -> one gain and one break
										{
											nbGain = 1;
											nbBreak = 1;
										}
								}
								
								if(WeightedDupCost != bestScore) //subtracting a Duplication to the score if that option was set
									caseScore.push_back(WeightedDupCost);

								//getting the new score
								scorevector.push_back(AggregateScore(caseScore,nbGain,nbBreak));

								//clearing
								caseScore.clear();
							}


	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
double AdjMatrix::C1DupWithOther(int NodeId1, int NodeId2, bool firstDup )
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithOther"<< endl;	

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	double score;
	if(!firstDup) // NodeId2 is the duplication
		score = D1(NodeId2, NodeId1, rtree2, rtree1, !firstDup);
	else // nodeId1 is the duplication
		score = D1(NodeId1, NodeId2, rtree1, rtree2, !firstDup);

	return score;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstDup (bool): true if NodeId1 is the duplication

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2

NB: Do not presume wether NodeId1 or NodeId2 is the duplication; Other may be Extant, Speciation, Null or SpeciationOut
*/
double AdjMatrix::C0DupWithOther(int NodeId1, int NodeId2, bool firstDup)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0DupWithOther"<< endl;	

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	double score;
	if(!firstDup) // NodeId2 is the duplication
		score = D0(NodeId2, NodeId1, rtree2, rtree1, !firstDup);
	else // nodeId1 is the duplication
		score = D0(NodeId1, NodeId2, rtree1, rtree2, !firstDup);

	return score;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2
*/
double AdjMatrix::C1SpecWithSpec(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SpecWithSpec"<< endl;	


	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);

	double sa1,sa2,sb1,sb2;

	if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[0]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[1]) ) ))//case where son1 of node 1 corresponds to son1 of node 2 // using the OR here does not change a thing in the case of the speciation, but will allow to treat them in the same way as synchronous speciation out
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[0];
		sb2 = sonsId2[1];
	}
	else if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[1]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[0]) ) )) //sanity check
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[1];
		sb2 = sonsId2[0];	
	}
	else
	{
		throw Exception("AdjMatrix::C1SpecWithSpec : no children corresponds to each other!");
	}


	double c1a1b1 = getC1(sa1,sb1);
	double c0a1b1 = getC0(sa1,sb1);
	double c1a2b2 = getC1(sa2,sb2);
	double c0a2b2 = getC0(sa2,sb2);

	vector <double> scorevector;
	vector <double> caseScore;

	//case 1 -> both adjacency have been kept
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1a2b2 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();

	//case 2 -> adjacency have been kept in child 1 but not child 2 -> 1 break
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0a2b2 );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 1 break

	caseScore.clear();

	//case 3 -> adjacency have been kept in child 2 but not child 1 -> 1 break
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1a2b2 );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 1 break

	//case 4 -> adjacency have not been kept in child 2 and child 1 -> 2 break
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0a2b2 );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 2 break


	caseScore.clear();

	return ((*this).*scoreComparatorfunc)(scorevector);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2
*/
double AdjMatrix::C0SpecWithSpec(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0SpecWithSpec"<< endl;	
	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);

	double sa1,sa2,sb1,sb2;

	if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[0]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[1]) ) ))//case where son1 of node 1 corresponds to son1 of node 2 // using the OR here does not change a thing in the case of the speciation, but will allow to treat them in the same way as synchronous speciation out
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[0];
		sb2 = sonsId2[1];
	}
	else if(( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[0]) , Rtree2.getNode(sonsId2[1]) ) ) || ( Rtree1.haveSameSpecies( Rtree1.getNode(sonsId1[1]) , Rtree2.getNode(sonsId2[0]) ) )) //sanity check
	{
		sa1 = sonsId1[0];
		sa2 = sonsId1[1];
		sb1 = sonsId2[1];
		sb2 = sonsId2[0];	
	}
	else
	{
		throw Exception("AdjMatrix::C1SpecWithSpec : no children corresponds to each other!");
	}


	double c1a1b1 = getC1(sa1,sb1);
	double c0a1b1 = getC0(sa1,sb1);
	double c1a2b2 = getC1(sa2,sb2);
	double c0a2b2 = getC0(sa2,sb2);

	vector <double> scorevector;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept -> 2 gain
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1a2b2 );
	scorevector.push_back(AggregateScore(caseScore,2,0)); // 2 gain 0 break

	caseScore.clear();

	//case 2 -> adjacency have been kept in child 1 but not child 2 -> 1 Gain
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0a2b2 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();

	//case 3 -> adjacency have been kept in child 2 but not child 1 -> 1 Gain
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1a2b2 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	//case 4 -> adjacency have not been kept in child 2 and child 1
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0a2b2 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 2 break


	caseScore.clear();

	return ((*this).*scoreComparatorfunc)(scorevector);
}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2

NB: add something later to account for co-event
*/
double AdjMatrix::C1DupWithDup(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1DupWithDup"<< endl;	

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	vector<double> scorevector;
	//D1 -> the duplication of NodeId1 was before the duplication of NodeId2
	scorevector.push_back( D1(NodeId1, NodeId2, rtree1, rtree2, false));

	//D2 -> the duplication of NodeId2 was before the duplication of NodeId1
	scorevector.push_back( D1(NodeId2, NodeId1, rtree2, rtree1, true));

	//D12 -> both duplications were simultaneous
	scorevector.push_back( D12(NodeId1,NodeId2) ) ;

	return ((*this).*scoreComparatorfunc)(scorevector);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2
*/
double AdjMatrix::C0DupWithDup(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0DupWithDup"<< endl;	

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	vector<double> scorevector;
	//D0 -> the duplication of NodeId1 was before the duplication of NodeId2
	scorevector.push_back( D0(NodeId1, NodeId2, rtree1, rtree2, false));

	//D0 -> the duplication of NodeId2 was before the duplication of NodeId1
	scorevector.push_back( D0(NodeId2, NodeId1, rtree2, rtree1, true));

	return ((*this).*scoreComparatorfunc)(scorevector);
}



//DeCoLT specific cases

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2
*/
double AdjMatrix::C1NullWithNull(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1NullWithNull"<< endl;

	vector<double> scorevector;

	//NodeId1 and NodeId2 only have 1 child
	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);


	//2possible cases
	vector <double> caseScore;

	//case 1 -> an adjacency exists between their children
	caseScore.push_back( getC1(sonsId1[0], sonsId2[0]) );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();

	//case 1 -> no adjacency exists between their children -> we need a break
	caseScore.push_back( getC0(sonsId1[0], sonsId2[0]) );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 1 break

	caseScore.clear();

	return ((*this).*scoreComparatorfunc)(scorevector);	
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2
*/
double AdjMatrix::C0NullWithNull(int NodeId1, int NodeId2)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0NullWithNull"<< endl;

	vector<double> scorevector;

	//NodeId1 and NodeId2 only have 1 child
	vector < int> sonsId1 = Rtree1.getSonsId(NodeId1);
	vector < int> sonsId2 = Rtree2.getSonsId(NodeId2);


	//2possible cases
	vector <double> caseScore;

	//case 1 -> an adjacency exists between their children -> need 1 gain
	caseScore.push_back( getC1(sonsId1[0], sonsId2[0]) );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();

	//case 1 -> no adjacency exists between their children
	caseScore.push_back( getC0(sonsId1[0], sonsId2[0]) );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();

	return ((*this).*scoreComparatorfunc)(scorevector);	
}

/*
Find the correct pair of id to compare in the case where one is a Sout

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstSout (bool): wether the Sout is NodeId1 or NodeId2

Returns:
	(pair <int,int>): pair of ids to get a score from
*/
pair <int,int> AdjMatrix::CostSoutWithExtantOrSpecOrNullAux(int NodeId1, int NodeId2, bool firstSout)
{
	int SoutNodeId;
	ReconciledTree * RtreeSout;
	vector<int> SonsIds;
	int CorrectSonId = -1;

	//determining which node is the Sout
	if(firstSout)
	{
		SoutNodeId = NodeId1;
		*RtreeSout = Rtree1;
	}
	else
	{
		SoutNodeId = NodeId2;
		*RtreeSout = Rtree2;	
	}

	int sp = RtreeSout->getNodeSpecies(SoutNodeId);

	SonsIds = RtreeSout->getSonsId(SoutNodeId);

	//we want to get the Son of the Sout that stayed in the species tree
	for(unsigned i = 0; i < SonsIds.size(); i++)
	{
		if(RtreeSout->getNodeSpecies(SonsIds[i]) == sp)
		{
			CorrectSonId = SonsIds[i];
			break;
		}
	}

	//cout << NodeId1 << "," << NodeId2 << " " << firstSout << "->" << sp << "-" << CorrectSonId << endl;
	//sanity check
	if(CorrectSonId == -1)
		throw Exception("AdjMatrix::CostSoutWithExtantOrSpecOrNullAux : found no son of the Sout node that stayed in the same species.");

	if(firstSout)
		NodeId1 = CorrectSonId;
	else
		NodeId2 = CorrectSonId;

	pair <int,int> IdsToReturn (NodeId1,NodeId2);

	return IdsToReturn;
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstSout (bool): wether the Sout is NodeId1 or NodeId2

Returns:
	(double): the cost of having an adjacencies between NodeId1 and NodeId2
*/
double AdjMatrix::C1SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SoutWithExtantOrSpecOrNull"<< endl;

	pair <int,int> orderedIds;
	orderedIds = CostSoutWithExtantOrSpecOrNullAux( NodeId1,  NodeId2, firstSout);

	//now that we have the correct son, we can apply the formula
	//only one case, with no gain or break -> simple transmission

	return getC1(orderedIds.first, orderedIds.second);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstSout (bool): wether the Sout is NodeId1 or NodeId2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2
*/
double AdjMatrix::C0SoutWithExtantOrSpecOrNull(int NodeId1, int NodeId2, bool firstSout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0SoutWithExtantOrSpecOrNull"<< endl;

	pair <int,int> orderedIds;
	orderedIds = CostSoutWithExtantOrSpecOrNullAux(NodeId1, NodeId2,  firstSout);

	//now that we have the correct son, we can apply the formula
	//only one case, with no gain or break -> simple transmission

	return getC0(orderedIds.first, orderedIds.second);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2
Presumes that both Speciation out event occurs at the same time -> treat as a speciation
*/
double AdjMatrix::C1SoutWithSoutSynchronous(int NodeId1, int NodeId2) 
{
	//presumes that both Speciation out event occurs at the same time
	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> SonsIds1 = Rtree1.getSonsId(NodeId1);
	vector < int> SonsIds2 = Rtree2.getSonsId(NodeId2);

	double sa1,sat,sb1,sbt; //sat and sbt are the transferred children



	int CorrectSonInd;
	int sp;
	unsigned i;
	//we want to get the Son of the Sout that stayed in the species tree for the node 1
	CorrectSonInd = -1;
	sp =  Rtree1.getNodeSpecies(NodeId1);
	for(i = 0; i < SonsIds1.size(); i++)
	{
		if(Rtree1.getNodeSpecies(SonsIds1[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::C1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sat = SonsIds1[1];
		sa1 = SonsIds1[0];
	}
	else // the non transferred son as index 1
	{
		sat = SonsIds1[0];
		sa1 = SonsIds1[1];	
	}

	//we want to get the Son of the Sout that stayed in the species tree for the node 2
	CorrectSonInd = -1;
	sp =  Rtree2.getNodeSpecies(NodeId2);
	for(i = 0; i < SonsIds2.size(); i++)
	{
		if(Rtree2.getNodeSpecies(SonsIds2[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::C1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sbt = SonsIds2[1];
		sb1 = SonsIds2[0];
	}
	else // the non transferred son as index 1
	{
		sbt = SonsIds2[0];
		sb1 = SonsIds2[1];
	}


	double c1a1b1 = getC1(sa1,sb1);
	double c0a1b1 = getC0(sa1,sb1);
	double c1atbt = getC1(sat,sbt);
	double c0atbt = getC0(sat,sbt);
	//each occurence of c0a1b1 costs a break; c0atbt are free because of the possibly asynchronous transfer events

	vector <double> scorevector;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1atbt );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();

	//case 2 -> adjacency have been kept in child 1 but not child t
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0atbt );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();

	//case 3 -> adjacency have been kept in child t but not child 1 -> 1 break
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1atbt );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 1 break

	//case 4 -> adjacency have not been kept in child t and child 1 -> 1 break
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0atbt );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 1 break

	caseScore.clear();

	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2

Presumes that both Speciation out event are two separate events -> the transfered children aren't expected to be neighbors
*/
double AdjMatrix::C1SoutWithSoutaSynchronous(int NodeId1, int NodeId2)
{

	//Four possible cases once we've matched the children of NodeId1 and NodeId2 by compatibility

	vector < int> SonsIds1 = Rtree1.getSonsId(NodeId1);
	vector < int> SonsIds2 = Rtree2.getSonsId(NodeId2);

	double sa1,sat,sb1,sbt; //sat and sbt are the transferred children



	int CorrectSonInd;
	int sp;
	unsigned i;
	//we want to get the Son of the Sout that stayed in the species tree for the node 1
	CorrectSonInd = -1;
	sp =  Rtree1.getNodeSpecies(NodeId1);
	for(i = 0; i < SonsIds1.size(); i++)
	{
		if(Rtree1.getNodeSpecies(SonsIds1[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::C1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sat = SonsIds1[1];
		sa1 = SonsIds1[0];
	}
	else // the non transferred son as index 1
	{
		sat = SonsIds1[0];
		sa1 = SonsIds1[1];	
	}

	//we want to get the Son of the Sout that stayed in the species tree for the node 2
	CorrectSonInd = -1;
	sp =  Rtree2.getNodeSpecies(NodeId2);
	for(i = 0; i < SonsIds2.size(); i++)
	{
		if(Rtree2.getNodeSpecies(SonsIds2[i]) == sp)
		{
			CorrectSonInd = i;
			break;
		}
	}
	//sanity check
	if(CorrectSonInd == -1)
		throw Exception("AdjMatrix::C1SoutWithSoutaSynchronous : found no son of the Sout node that stayed in the same species.");
	else if (CorrectSonInd == 0) // the non transferred son as index 0
	{
		sbt = SonsIds2[1];
		sb1 = SonsIds2[0];
	}
	else // the non transferred son as index 1
	{
		sbt = SonsIds2[0];
		sb1 = SonsIds2[1];
	}


	double c1a1b1 = getC1(sa1,sb1);
	double c0a1b1 = getC0(sa1,sb1);
	double c1atbt = getC1(sat,sbt);
	double c0atbt = getC0(sat,sbt);
	//each occurence of c0a1b1 costs a break, each occurence of c1atbt costs a gain.

	vector <double> scorevector;
	vector <double> caseScore;

	//case 1 -> both adjacencies have been kept -> 1 gain
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c1atbt );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();

	//case 2 -> adjacency have been kept in child 1 but not child t
	caseScore.push_back( c1a1b1 );
	caseScore.push_back( c0atbt );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();

	//case 3 -> adjacency have been kept in child 2 but not child 1 -> 1 break + 1 gain
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c1atbt );
	scorevector.push_back(AggregateScore(caseScore,1,1)); // 1 gain 1 break

	//case 4 -> adjacency have not been kept in child 2 and child 1 -> 1 break
	caseScore.push_back( c0a1b1 );
	caseScore.push_back( c0atbt );
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 1 break

	caseScore.clear();

	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

*/
double AdjMatrix::C1SoutWithSout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1SoutWithSout"<< endl;

	vector <double > scorevector;

	scorevector.push_back(C1SoutWithSoutSynchronous( NodeId1,NodeId2)); // co event case
	scorevector.push_back(C1SoutWithSoutaSynchronous(NodeId1,NodeId2)); // non co event case

	return ((*this).*scoreComparatorfunc)(scorevector);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacencies between NodeId1 and NodeId2

NB: same formula wether the Souts are simultaneous or not 
*/
double AdjMatrix::C0SoutWithSout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0SoutWithSout"<< endl;
	//we can treat that event with the speciation formulae (which only need one of the two children to have the same species)
	return C0SpecWithSpec(NodeId1, NodeId2);
}



/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

NB: add something later for co-event checking (a priori this formulae are good for both as long as a gain is added externally (ie. after the adj tree is built))
*/
double AdjMatrix::C1RecWithRec(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1RecWithRec"<< endl;
	//both reception node only a have  child
	int sonId1 = Rtree1.getSonsId(NodeId1)[0];
	int sonId2 = Rtree2.getSonsId(NodeId2)[0];

	double c1 = getC1(sonId1, sonId2);
	double c0 = getC0(sonId1, sonId2);

	vector <double > scorevector;
	vector <double > caseScore;

	//two cases
	
	caseScore.push_back(c1);
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();

	caseScore.push_back(c0);
	scorevector.push_back(AggregateScore(caseScore,0,1)); // 0 gain 1 break

	return ((*this).*scoreComparatorfunc)(scorevector);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2

NB: the synchronous case recquires the node to be linked
*/
double AdjMatrix::C0RecWithRec(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0RecWithRec"<< endl;

	//both reception node only a have  child
	int sonId1 = Rtree1.getSonsId(NodeId1)[0];
	int sonId2 = Rtree2.getSonsId(NodeId2)[0];

	double c1 = getC1(sonId1, sonId2);
	double c0 = getC0(sonId1, sonId2);

	vector <double > scorevector;
	vector <double > caseScore;

	//two cases
	
	caseScore.push_back(c1);
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();

	caseScore.push_back(c0);
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	return ((*this).*scoreComparatorfunc)(scorevector);
}



/*

Takes:
 - NodeIdBout (int): node id in RtreeBout
 - NodeIdOther (int): node id in RtreeOther
 - RtreeBout (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the Bifurcation Out is in the tree1 (true) or not (false)

Returns:
	(double): the cost of having an adjacency between NodeIdBout and NodeIdOther
*/
double AdjMatrix::B1(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert)
{
	vector <double> scorevector;

	//we first get the ids of the children of NodeIdDup
	vector <int> BoutSonsId = RtreeBout->getSonsId(NodeIdBout);

	//There is 4 possible cases
	vector <double> caseScore;
	
	double c1ab1 = getC1(BoutSonsId[0],NodeIdOther,invert);
	double c1ab2 = getC1(BoutSonsId[1],NodeIdOther,invert);
	double c0ab1 = getC0(BoutSonsId[0],NodeIdOther,invert);
	double c0ab2 = getC0(BoutSonsId[1],NodeIdOther,invert);

	//case 1 -> child 1 of the bout conserved an adjacency but not child 2
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();
	//case 2 -> child 2 of the bout conserved an adjacency but not child 1
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break

	caseScore.clear();
	//case 3 -> child 1 of the bout conserved an adjacency and child 2 too -> need a gain
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency -> no break because we are in the dead
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break


	return ((*this).*scoreComparatorfunc)(scorevector);
}



/*
Both nodes are Bifurcation Out that are presumed simultaneous

Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2
*/
double AdjMatrix::B12(int NodeId1, int NodeId2)
{
	vector <double> scorevector;

	//we first get the ids of the children of NodeIdDup
	vector <int> BoutSonsId1 = Rtree1.getSonsId(NodeId1);
	vector <int> BoutSonsId2 = Rtree2.getSonsId(NodeId2);

	//There is 16 possible cases representing all c1 and c0 combination of the sons of 1 an 2
	vector <double> caseScore;
	int nbGain;
	int nbBreak = 0; // 0 break because we are in the dead

	for(unsigned c1a1b1 = 0; c1a1b1 < 2; c1a1b1 ++)
			for(unsigned c1a1b2 = 0; c1a1b2 < 2; c1a1b2 ++)
					for(unsigned c1a2b1 = 0; c1a2b1 < 2; c1a2b1 ++)
							for(unsigned c1a2b2 = 0; c1a2b2 < 2; c1a2b2 ++)
							{

								if(!c1a1b1)
									caseScore.push_back(getC0(BoutSonsId1[0],BoutSonsId2[0]));
								else
									caseScore.push_back(getC1(BoutSonsId1[0],BoutSonsId2[0]));

								if(!c1a1b2)
									caseScore.push_back(getC0(BoutSonsId1[0],BoutSonsId2[1]));
								else
									caseScore.push_back(getC1(BoutSonsId1[0],BoutSonsId2[1]));

								if(!c1a2b1)
									caseScore.push_back(getC0(BoutSonsId1[1],BoutSonsId2[0]));
								else
									caseScore.push_back(getC1(BoutSonsId1[1],BoutSonsId2[0]));

								if(!c1a2b2)
									caseScore.push_back(getC0(BoutSonsId1[1],BoutSonsId2[1]));
								else
									caseScore.push_back(getC1(BoutSonsId1[1],BoutSonsId2[1]));

								//determining the number of gains and breaks
								int nbadj = c1a1b1 + c1a1b2 + c1a2b1 + c1a2b2;
								switch( nbadj )
								{

									case 0: // none linked -> 0 break because we are in the dead
										nbGain = 0;
										break;
									case 1: // only one linked -> 0 break because we are in the dead
										nbGain = 0;	
										break;
									case 3: // three linked -> 1 gain
										nbGain = 1;	
										break;
									case 4: // all linked -> 2 gain
										nbGain = 2;	
										break;
									default: //exactly 2 adjs
										nbGain = 0;
										if(( c1a1b1 == c1a1b2 ) || ( c1a1b1 == c1a2b1 )) // scenarios where 1 child has two adjs -> one gain and no break because we are in the dead
										{
											nbGain = 1;
										}
								}
								
								//getting the new score
								scorevector.push_back(AggregateScore(caseScore,nbGain,nbBreak));

								//clearing
								caseScore.clear();
							}


	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

*/
double AdjMatrix::C1BoutWithBout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1BoutWithBout"<< endl;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	vector<double> scorevector;
	//B1 -> the Bout of NodeId1 was before the Bout of NodeId2
	scorevector.push_back( B1(NodeId1, NodeId2, rtree1, rtree2, false));

	//B2 -> the Bout of NodeId2 was before the Bout of NodeId1
	scorevector.push_back( B1(NodeId2, NodeId1, rtree2, rtree1, true));

	//B12 -> both Bout were simultaneous
	scorevector.push_back( B12(NodeId1,NodeId2));

	return ((*this).*scoreComparatorfunc)(scorevector);
}

/*

Takes:
 - NodeIdBout (int): node id in RtreeBout
 - NodeIdOther (int): node id in RtreeOther
 - RtreeBout (ReconciledTree *): a reocnciled tree
 - Rtreeother (ReconciledTree *): a reconciled tree
 - invert (bool): wether the Bifurcation Out is in the tree1 (true) or not (false)

Returns:
	(double): the cost of having no adjacency between NodeIdBout and NodeIdOther
*/
double AdjMatrix::B0(int NodeIdBout, int NodeIdOther, ReconciledTree * RtreeBout, ReconciledTree * Rtreeother, bool invert)
{
	vector <double> scorevector;

	//we first get the ids of the children of NodeIdDup
	vector <int> BoutSonsId = RtreeBout->getSonsId(NodeIdBout);

	//There is 4 possible cases
	vector <double> caseScore;
	
	double c1ab1 = getC1(BoutSonsId[0],NodeIdOther,invert);
	double c1ab2 = getC1(BoutSonsId[1],NodeIdOther,invert);
	double c0ab1 = getC0(BoutSonsId[0],NodeIdOther,invert);
	double c0ab2 = getC0(BoutSonsId[1],NodeIdOther,invert);

	//each c1 cost a gain

	//case 1 -> child 1 of the bout conserved an adjacency but not child 2 -> 1 gain
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();
	//case 2 -> child 2 of the bout conserved an adjacency but not child 1 -> 1 gain
	caseScore.push_back( c1ab2 );
	caseScore.push_back( c0ab1 );
	scorevector.push_back(AggregateScore(caseScore,1,0)); // 1 gain 0 break

	caseScore.clear();
	//case 3 -> child 1 of the bout conserved an adjacency and child 2 too -> need 2 gain
	caseScore.push_back( c1ab1 );
	caseScore.push_back( c1ab2 );
	scorevector.push_back(AggregateScore(caseScore,2,0)); // 2 gain 0 break

	caseScore.clear();
	//case 4 -> neither child 1 nor child 2 kept an adjacency
	caseScore.push_back( c0ab1 );
	caseScore.push_back( c0ab2 );
	scorevector.push_back(AggregateScore(caseScore,0,0)); // 0 gain 0 break


	return ((*this).*scoreComparatorfunc)(scorevector);	
}

/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2

*/
double AdjMatrix::C0BoutWithBout(int NodeId1, int NodeId2) 
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0BoutWithBout"<< endl;

	vector<double> scorevector;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;



	//B0 -> the Bout of NodeId1 was before the Bout of NodeId2
	scorevector.push_back( B0(NodeId1, NodeId2, rtree1, rtree2, false));

	//B0 -> the Bout of NodeId2 was before the Bout of NodeId1
	scorevector.push_back( B0(NodeId2, NodeId1, rtree2, rtree1, true));

	return ((*this).*scoreComparatorfunc)(scorevector);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstBout (bool): whether nodeId1 is the Bout  (true) or not (false)

Returns:
	(double): the cost of having an adjacency between NodeId1 and NodeId2

*/
double AdjMatrix::C1BoutWithRec(int NodeId1, int NodeId2, bool firstBout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C1BoutWithRec"<< endl;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;



	if(firstBout)	// the Bout is NodeId1
		return B1(NodeId1, NodeId2, rtree1, rtree2, !firstBout);

	// the Bout is NodeId2
	return B1(NodeId2, NodeId1, rtree2, rtree1, !firstBout);
}


/*
Takes:
 - NodeId1 (int): node id in Rtree1
 - NodeId2 (int): node id in Rtree2
 - firstBout (bool): whether nodeId1 is the Bout  (true) or not (false)

Returns:
	(double): the cost of having no adjacency between NodeId1 and NodeId2

*/
double AdjMatrix::C0BoutWithRec(int NodeId1, int NodeId2, bool firstBout)
{
	if(verbose)
		cout << NodeId1 << " , " << NodeId2 << " -> " << "C0BoutWithRec"<< endl;

	ReconciledTree * rtree1;
	ReconciledTree * rtree2;
	rtree1 = &Rtree1;
	rtree2 = &Rtree2;


	if(firstBout)	// the Bout is NodeId1
		return B0(NodeId1, NodeId2, rtree1, rtree2, !firstBout);

	// the Bout is NodeId2
	return B0(NodeId2, NodeId1, rtree2, rtree1, !firstBout);
}


