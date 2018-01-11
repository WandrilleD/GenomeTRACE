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

This file contains utils specific to Wmodif stuff like redrawing reconciliations , ...

Created the: 24-03-2016
by: Wandrille Duchemin

Last modified the: 11-01-2018
by: Wandrille Duchemin

*/

#include "WmodifUtils.h"





/*
Reset the ECFs whose gfam 1 or 2 is GeneFamilyId

Takes:
 - GeneFamilyId (int) : id of the gene family that will be reset
 - ECFams  (vector < EquivalenceClassFamily > *) : vector of Equivalence class families 

*/
void prepRecRedrawEqClass(int GeneFamilyId, vector < EquivalenceClassFamily > * ECFams)
{
    //cout << "prepping EqClass" << endl;
	for(size_t i = 0; i < ECFams->size(); i++)
	{

		if((ECFams->at(i).getGfamily1() == GeneFamilyId) || (ECFams->at(i).getGfamily2() == GeneFamilyId))
        {
            //WMODIF
            pair <int,int> gfams ;
            gfams.first = ECFams->at(i).getGfamily1();
            gfams.second = ECFams->at(i).getGfamily2();


            vector <pair <string, string> > adjacencies =  ECFams->at(i).getAdjs();

            //cout << "before assign:" << ECFams->at(i).getNbAdj();

            //cout << "assign" << endl;

            ECFams->at(i) = EquivalenceClassFamily(gfams.first,gfams.second);
            //cout << "after assign:" << ECFams->at(i).getNbAdjTrees()<<endl;

            //cout << "before"<<ECFams->at(i).getNbAdj() << endl;
            ECFams->at(i).CheckAddAdjList(adjacencies , gfams); // adds several pair of leaf names
            //cout << "after"<<ECFams->at(i).getNbAdj() << endl;
            //WMODIF



			//ECFams->at(i).reset();
        }
	}
    //cout << "end prepping EqClass" << endl;
}

/*
remove from the Coevents the adjs whose gfam 1 or 2 is GeneFamilyId

Takes:
 - GeneFamilyId (int) : id of the gene family that will be reset
 - CoEventSet (vector < CoEvent > *) : set of co-events to update


*/
void prepRecRedrawCoEvent(int GeneFamilyId, vector <CoEvent> * CoEventSet, bool verbose)
{
	for(size_t i = 0; i < CoEventSet->size(); i++)
	{
		bool changed = CoEventSet->at(i).removeGeneFam(GeneFamilyId);
		if(verbose)
			if(changed)
				cout << "changed CoEvent " << i << endl;
	}
}


/*
remove from the Coevents the adjs whose 

Takes:
 - ECFidToReset (int) : id of the ECF that will be redrawn
 - CoEventSet (vector < CoEvent > *) : set of co-events to update


*/
void prepECFRedrawCoEvent(int ECFidToReset, vector <CoEvent> * CoEventSet, bool verbose)
{
    for(size_t i = 0; i < CoEventSet->size(); i++)
    {
        bool changed = CoEventSet->at(i).removeECF(ECFidToReset);
        if(verbose)
            if(changed)
                cout << "changed CoEvent " << i << endl;
    }
}

/*
Takes:
 - EquivalenceClassFamily * ECF : the requivalence class family to refine, compute and backtrack
 - ReconciledTree * Rtree1 : reconciled tree of the gfam1 of ECF
 - ReconciledTree * Rtree2 : reconciled tree of the gfam2 of ECF
 - double Again : cost of an adjacency gain
 - double Bgain : cost of an adjacency break
 - double absencePenalty : penalty to be given if an adjacency is absent (if -1, then this option is ignored)
 - bool doAllPair : forces refinement to consider equaivalence classes without any adjacency
 - bool withTransfer : if transfers are allowed or not
 - bool substractRecoToAdj : if true, the weighted cost of an evolutionnary event will be subtracted to the adjacency scenario cost whenever a co-event is found
 - double weightedDupCost : cost of a gene duplication weighted to be compared with adjacency histories costs
 - double weightedLossCost : cost of a gene loss weighted to be compared with adjacency histories costs
 - double weightedHGTCost : cost of a gene transfer weighted to be compared with adjacency histories costs
 - bool alwaysGainAtTop : if true, there will always be a gain at the top of the tree
 - double C1Advantage : define the probability to choose C1 over C0 if they have the same score at the root of the adjacency tree (0.5 means equal chance)
 - bool superverbose : prints LOTS of info.

*/
void refiningAndComputingECFamily(EquivalenceClassFamily * ECF, 
									ReconciledTree * Rtree1, ReconciledTree * Rtree2,
									double Again, double Bgain, 
									double absencePenalty, bool doAllPair, bool withTransfer, bool substractRecoToAdj, 
									double weightedDupCost, double weightedLossCost, double weightedHGTCost, 
									bool alwaysGainAtTop, double C1Advantage,bool superverbose,
                                    map<int,vector<float> > speciesC0C1,
                                    map<int, map<string,int> > speGeneAdjNb,
                                    bool interactionMode
                                    )
{

    if(superverbose)
        cout << "refiningAndComputingECFamily. "<< ECF->getNbAdj() << " adjs." << endl;

    ECF->refine(Rtree1, Rtree2, doAllPair , withTransfer , superverbose); 


    //cout << "refiningAndComputingECFamily " << ECF->getNbAdjTrees() << " adj trees before matrix creation and stuff" << endl;    
    
    
    if(superverbose)
    {
        cout << ECF->getNbAdj() << " adjs in " << ECF->getNbEqClasses() << " Eclass." << endl;
        for(unsigned j = 0; j <  ECF->getNbEqClasses() ; j++)
            cout << ECF->getEClass(j)->getNbAdj() << "-" ;
        cout << endl;
    }

    
    
    map < string, map <string , double> > adjacencyScores; // will stay empty, just here to ensure compat with deco's current version
    //cout << "refiningAndComputingECFamily " << absencePenalty << endl;
    ECF->createAdjMatrix( adjacencyScores, speciesC0C1, speGeneAdjNb,// dummy args because of artdeco
                          Again, Bgain ,
                          Rtree1,  Rtree2,
                          superverbose ,
                          false , 1 , // boltzmann stuff
                          absencePenalty,10000,
                          interactionMode);


    if(!substractRecoToAdj)
        ECF->computeAdjMatrix();
    else //using weighted reconciliation event costs
        ECF->computeAdjMatrix( weightedDupCost, weightedLossCost, weightedHGTCost );

    bool useCount = false;
    ECF->backtrackAdjMatrixForSelf(Rtree1, Rtree2, 
                                    false, // boltzmann
                                    alwaysGainAtTop , C1Advantage, useCount);

    if(superverbose)
        cout << "refiningAndComputingECFamily " << ECF->getNbAdjTrees() << " adj trees" << endl;


    //cout << "dumping for ECF : " << ECF->getGfamily1() << "-" << ECF->getGfamily2()  << endl;
    ECF->dumpStuff();//accepted memleak
    //cout << "dumping done" << endl;

}

/*
Takes:
 - vector <EquivalenceClassFamily> * ECFams : the requivalence class family list
 - vector < GeneFamily * > * GeneFamilyList : List of Gene Families
 - int GfamIdToReset : id of the resetted Gfam: ECFs with that id as gfam 1 or 2 will be computed
 - double Again : cost of an adjacency gain
 - double Bgain : cost of an adjacency break
 - double absencePenalty : penalty to be given if an adjacency is absent (if -1, then this option is ignored)
 - bool doAllPair : forces refinement to consider equaivalence classes without any adjacency
 - bool withTransfer : if transfers are allowed or not
 - bool substractRecoToAdj : if true, the weighted cost of an evolutionnary event will be subtracted to the adjacency scenario cost whenever a co-event is found
 - double weightedDupCost : cost of a gene duplication weighted to be compared with adjacency histories costs
 - double weightedLossCost : cost of a gene loss weighted to be compared with adjacency histories costs
 - double weightedHGTCost : cost of a gene transfer weighted to be compared with adjacency histories costs
 - bool alwaysGainAtTop : if true, there will always be a gain at the top of the tree
 - double C1Advantage : define the probability to choose C1 over C0 if they have the same score at the root of the adjacency tree (0.5 means equal chance)
 - bool verbose : prints some info.
 - bool superverbose : prints LOTS of info.

*/
void refiningAndComputingECFamilyList(vector <EquivalenceClassFamily> * ECFams, vector < GeneFamily * > * GeneFamilyList,
										int GfamIdToReset,
										double Again, double Bgain, 
										double absencePenalty, bool doAllPair, bool withTransfer, bool substractRecoToAdj, 
										double weightedDupCost, double weightedLossCost, double weightedHGTCost, 
										bool alwaysGainAtTop, double C1Advantage, bool verbose, bool superverbose,
                                        map<int,vector<float> > speciesC0C1,
                                        map<int, map<string,int> > speGeneAdjNb,
                                        bool interactionMode 
                                        )
{
    
    for(unsigned i = 0; i <  ECFams->size(); i++)
    {
        EquivalenceClassFamily * ECF = &ECFams->at(i);
        int gfam1 = ECF->getGfamily1();
        int gfam2 = ECF->getGfamily2();
 
        

        if( (gfam1 == GfamIdToReset) || (gfam2 == GfamIdToReset) )
        {
           ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
           ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();

           
           
            if(verbose)
                cout << "refining equivalence class " << i << "(" << gfam1 << "," << gfam2 << "). " << ECF->getNbAdj() << " adjacencies."<< endl; 

            refiningAndComputingECFamily(ECF, Rtree1, Rtree2,
      									Again, Bgain, 
										absencePenalty, doAllPair, withTransfer, substractRecoToAdj, 
										weightedDupCost, weightedLossCost, weightedHGTCost, 
										alwaysGainAtTop, C1Advantage, superverbose,
                                        speciesC0C1,
                                        speGeneAdjNb,
                                        interactionMode 
                                        );


            if(verbose)
                cout << "-> " << ECF->getNbEqClasses() << " equivalence classes." << endl;
            if(verbose)
                cout << "Backtrack finished: " << ECF->getNbAdjTrees() << " trees." << endl;
            if(verbose)
                cout << " Total of "<< ECF->getNbAdjGain() << " Adjacency Gains and " << ECF->getNbAdjBreak() << " Adjacency Breaks." << endl;
         }
	}
}


void GettingCoEvents(vector <CoEvent> * CoEventSet, vector <EquivalenceClassFamily> * ECFams, vector < GeneFamily * > * GeneFamilyList)
{
 
    for(unsigned i = 0 ; i < ECFams->size(); i++)
    {
        int gfam1 = ECFams->at(i).getGfamily1();
        int gfam2 = ECFams->at(i).getGfamily2();

        ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
        ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();

        //cout << "getting co-ev between "<< gfam1 << "-" << gfam2 << endl;

        for(unsigned j = 0; j < ECFams->at(i).getNbEqClasses() ; j++)
        {
            PopulatesCoeventsFromAdjForest( ECFams->at(i).getAdjForest(j), CoEventSet, i ,  gfam1,  gfam2,  Rtree1,  Rtree2);
        }

    }

    //cout << "created Co-events" << endl;

    //cout << "before refine : " << CoEventSet->size() << " coevents" << endl;
    //refineCoEvents( CoEventSet);
    //cout << "after refine : " << CoEventSet->size() << " coevents" << endl;

    //vector <bool> hasAncestorSonPairs = checkingCoEventsAncestorSonConstraint(  CoEventSet ,  GeneFamilyList);

    //cout << "checked Co-events" << endl;    
    
    //for(unsigned i = 0; i < hasAncestorSonPairs.size(); i++)
    //{
    //    if(hasAncestorSonPairs[i])
    //    {
    //        cout << "PB : co-event "<<i;
    //        cout << " ( event: " << CoEventSet->at(i).getEvent() << " , species: " << CoEventSet->at(i).getEvent();
    //        cout << " ) has some ancestor-son pair." << endl;
    //    }
    //}

}

/*
Splits the coevents by connex component
*/
void refineCoEvents(vector <CoEvent> * CoEventSet)
{
    vector <CoEvent> newCoEventSet;

    for(unsigned i = 0 ; i < CoEventSet->size() ; i++)
    {
        vector <CoEvent> coevs =  CoEventSet->at(i).splitByConnexComponent();
        for( unsigned j = 0 ; j < coevs.size() ; j++)
        {
            newCoEventSet.push_back(coevs[j]);
            //cout << "added coevent of size "<< newCoEventSet.back().getNumberOfGene() <<endl;
        }
    }
    //replacing the old coevent by the new
    CoEventSet->clear();
    for(unsigned i = 0 ; i < newCoEventSet.size() ; i++)
        CoEventSet->push_back(newCoEventSet[i]);



}

/*
Takes:
    - vector <CoEvent> * CoEventSet : set f coevents
    - vector < GeneFamily * > * GeneFamilyList : set of GeneFamilies

Returns:
    (vector <bool>) : vector which element are true if the correesponding co-"event contains at least one pait of ancestor-son

*/
vector <bool> checkingCoEventsAncestorSonConstraint( vector <CoEvent> * CoEventSet , vector < GeneFamily * > * GeneFamilyList)
{
    // 1. making a vector of reconciled trees
    vector <ReconciledTree *> RecTrees;

    for(unsigned i = 0 ; i < GeneFamilyList->size(); i++) 
    {
        RecTrees.push_back( GeneFamilyList->at(i)->getRecTree() );
    }

    //2. checking each co-event
    vector <bool> hasAncestorSonPair ; 

    for(unsigned i = 0 ; i < CoEventSet->size(); i++) 
    {
        hasAncestorSonPair.push_back(false);

        //2.1 getting the ancestor son pairs
        vector < vector < pair <int,int> > > AncestorSonPairs = CoEventSet->at(i).getAncestorSonPairs( RecTrees );

        //2.2 checking that the list are empty
        for(unsigned j = 0 ; j < AncestorSonPairs.size(); j++)
        {
            if(AncestorSonPairs[j].size() > 0)
            { // found non empty pair list -> pb
                hasAncestorSonPair[i] = true;
                break;
            }
        }

    }
    return hasAncestorSonPair;
}



/*
Draw *IN PLACE* a ne reconciliation for a given Gene Family (speciafied by its id).

*/
void ReDrawGeneFamilyProper(int GfamIdToReset, int redrawAlgo, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams,
                        MySpeciesTree * speciesTree, bool datedSpeciesTree, bool boundedTS, bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        bool verbose, bool superverbose ,
                        bool trySeveralECFSolutions, double RecWeight,
                        double probaFlat,
                        double probaC1Bias , MyCladesAndTripartitions * biasedCCP,
                        string gPathPrefix, int nbTry, int maxNbTry ,
                        bool noCoev , double DLRecCoevTemp,
                        string DTLRecCoevExecPath
                         )
{


    if(redrawAlgo == 0)
        GeneFamilyList->at(GfamIdToReset)->changeReconciliation( speciesTree, datedSpeciesTree, boundedTS); // just re-draw a tree in the already computed DTLMatrix. Won't work if no DTLMatrix have been computed
    else if(redrawAlgo == 1)
    { // if the topology is reset, then the reconciliation must be too

        double r = ( (double) rand() / RAND_MAX );
        if( r <= probaC1Bias)
        {
             GeneFamilyList->at(GfamIdToReset)->setUnrootedTree( *(biasedCCP->getRandomTree(probaFlat)) );
        }
        else
        {
            GeneFamilyList->at(GfamIdToReset)->makeUnrootedTree(true, probaFlat); // drawing a random tree in the CCP distribution
        }


        GeneFamilyList->at(GfamIdToReset)->makeReconciliation(speciesTree,
                                            withTransfer, withTransfer,  // transfer and transferLoss
                                            DupliCost, HGTCost, LossCost , // costs
                                            maxTS, TopoWeight, datedSpeciesTree, //here, dated.species.tree means that species tree is subdivided
                                            false, //gBoolParams.find("try.all.amalgamation")->second); // -> we don't want to amalgamate otherwise the drawn tree wuldn't be kept
                                            true); // to get random reconciliations
    }
    else if(redrawAlgo == 2)
    {
        cout << "using D(T)LRecCoev"<<endl;
        string DLRecCoevFolder = gPathPrefix;
        if( DLRecCoevFolder[ DLRecCoevFolder.size() - 1 ] != '/')
                DLRecCoevFolder += "/";
        DLRecCoevFolder += "DTLRecCoev/";
        DTLRecCoevWrapper wrapper( DLRecCoevFolder , DTLRecCoevExecPath);

        wrapper.setGFamToRedraw(GeneFamilyList->at(GfamIdToReset));
        wrapper.setRelativeTopoWeight(TopoWeight / RecWeight);
        wrapper.setDupCost(DupliCost);
        wrapper.setLossCost(LossCost);
        wrapper.setTemperature(DLRecCoevTemp);
        wrapper.setnbSample(maxNbTry);

        wrapper.setTransfer(withTransfer);
        wrapper.setDated(datedSpeciesTree);
        wrapper.setHGTCost(HGTCost);
        
        if(noCoev)
            wrapper.setNoCoev();

        wrapper.setSpTree(speciesTree);
    
        if(nbTry == 0)
        {
            bool allowSelf = false; // don't allow to be a self guide
            int guideGF = chooseOneLinkedfamily(GfamIdToReset,  ECFams, allowSelf );

            wrapper.MakePath();

            if(guideGF != -1)
                wrapper.writeCoevFile(GeneFamilyList->at(guideGF));
            else
                wrapper.writeCoevFile();
    
            //cout << "check 3.2"<<endl;
            wrapper.writeSpeciesFile();
    
            //cout << "check 4"<<endl;
            
            wrapper.setCommand();
            wrapper.launchCommand();
        }


        wrapper.readTree(nbTry);


    }
    else
    {
        cerr << "unknown redrawing algorithm selected..."<< endl;
        exit(1);
    }

    if(boundedTS)
        GeneFamilyList->at(GfamIdToReset)->setReconciledTreeTimeSlicesToBTS(speciesTree); // set the tree to a bounded time slice one (BTS)
    else if(datedSpeciesTree)
        GeneFamilyList->at(GfamIdToReset)->setReconciledTreeTimeSlicesToTS(speciesTree); // set the tree to a complete time sliced one (TS)
    else
        GeneFamilyList->at(GfamIdToReset)->setReconciledTreeTimeSlicesToNoTS();//no TS

    if(superverbose)
        GeneFamilyList->at(GfamIdToReset)->printRecTree();


}


/*
Draw *IN PLACE* a ne reconciliation for a given Gene Family (speciafied by its id).


Takes:
 - int GfamIdToReset : id of the gene family whose rconciliation to redraw
 - int redrawAlgo : wether a random topology should be dranw in the CCP distribution of the GeneFamily ()if false, the topology will remain the same)
 - vector < GeneFamily * > * GeneFamilyList : List of Gene Families
 - vector < EquivalenceClassFamily > * ECFams
 - vector < CoEvent > * CoEventSet : set of co-events to update
 - MySpeciesTree * speciesTree 
 - bool datedSpeciesTree
 - bool boundedTS
 - bool withTransfer
 - double DupliCost
 - double HGTCost
 - double LossCost
 - int maxTS
 - double TopoWeight
 - double Again : cost of an adjacency gain
 - double Abreak : cost of an adjacency break
 - double absencePenalty : penalty to be given if an adjacency is absent (if -1, then this option is ignored)
 - bool doAllPair : forces refinement to consider equaivalence classes without any adjacency
 - bool substractRecoToAdj : if true, the weighted cost of an evolutionnary event will be subtracted to the adjacency scenario cost whenever a co-event is found
 - double weightedDupCost : cost of a gene duplication weighted to be compared with adjacency histories costs
 - double weightedLossCost : cost of a gene loss weighted to be compared with adjacency histories costs
 - double weightedHGTCost : cost of a gene transfer weighted to be compared with adjacency histories costs
 - bool alwaysGainAtTop : if true, there will always be a gain at the top of the tree
 - double C1Advantage : define the probability to choose C1 over C0 if they have the same score at the root of the adjacency tree (0.5 means equal chance)
 - bool superverbose : prints LOTS of info.
 - bool trySeveralECFSolutions
 - double Temperature
 - double RecWeight 
 - bool deleteCoev
...

 - double probaFlat : probability to draw clades in a flat CCP distribution rather than the one given as input
 - double probaC1Bias = 0 
 - MyCladesAndTripartitions * biasedCCP = NULL
 - string gPathPrefix : path to ouput in
 - int nbTry [default = 0] : index of the current try 
 - int maxNbTry [default = 1] : total number of try for this family

*/
void ReDrawGeneFamily(int GfamIdToReset, int redrawAlgo, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                        MySpeciesTree * speciesTree, bool datedSpeciesTree, bool boundedTS, bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        double Again, double Abreak, double absencePenalty, bool doAllPair, bool substractRecoToAdj,
                        double weightedDupCost, double weightedLossCost, double weightedHGTCost, bool alwaysGainAtTop, double C1Advantage,
                        bool verbose, bool superverbose ,
                        bool trySeveralECFSolutions, double Temperature, double RecWeight, bool deleteCoev,
                        map<int,vector<float> > speciesC0C1,
                        map<int, map<string,int> > speGeneAdjNb, 
                        double probaFlat,
                        double probaC1Bias , MyCladesAndTripartitions * biasedCCP,
                        string gPathPrefix, int nbTry, int maxNbTry,
                        bool noCoev , double DLRecCoevTemp,
                        string DTLRecCoevExecPath,
                        bool interactionMode
                         )
{

    //1. prepping
    
    prepRecRedrawEqClass(GfamIdToReset, ECFams);

    if(deleteCoev)
        prepRecRedrawCoEvent(GfamIdToReset, CoEventSet, superverbose);

    //2. re-drawing a reconciliation
    ReDrawGeneFamilyProper( GfamIdToReset,  redrawAlgo, 
                        GeneFamilyList, ECFams,
                        speciesTree, datedSpeciesTree, boundedTS, withTransfer,
                        DupliCost, HGTCost, LossCost, maxTS, TopoWeight,
                        verbose, superverbose ,
                        trySeveralECFSolutions, RecWeight,
                        probaFlat,
                        probaC1Bias, biasedCCP,
                        gPathPrefix, nbTry, maxNbTry,
                        noCoev, DLRecCoevTemp,
                        DTLRecCoevExecPath
                         );


    //3. refining and computing equivalence classes
    
    refiningAndComputingECFamilyList( ECFams, GeneFamilyList,
                                                GfamIdToReset,
                                                Again, Abreak,
                                                absencePenalty, doAllPair , 
                                                withTransfer , substractRecoToAdj, 
                                                weightedDupCost, weightedHGTCost, weightedLossCost,
                                                alwaysGainAtTop, C1Advantage,
                                                verbose,superverbose,
                                                speciesC0C1,
                                                speGeneAdjNb ,
                                                interactionMode
                                                );



    //4. Co-Event Building
    if(deleteCoev)
        delete CoEventSet;//simpler for now to destroy and then rebuild from scratch...

    //CoEventSet = new vector <CoEvent>;
    GettingCoEvents( CoEventSet, ECFams, GeneFamilyList );
    //cout << "nb of coevents" << CoEventSet->size() << endl;
    //for(unsigned i = 0 ; i < CoEventSet->size(); i++)
    //{
    //    cout << CoEventSet->at(i).getEvent() << " ";
    //}
    //cout << endl;
    //double CoEventScore =  computeCoEventScore(*CoEventSet, DupliCost, LossCost, HGTCost, RecWeight);
    //cout << "coev score " << CoEventScore << endl;



    if(trySeveralECFSolutions)
    {    //5. trying to get a better solution on ECFs
        double CoEventScore =  computeCoEventScore(*CoEventSet, DupliCost, LossCost, HGTCost, RecWeight);

        for(size_t i = 0; i < ECFams->size(); i++) // this method is fairly crude and only tries to update once. Alternative could include : doing the same thing n time; doing this as long as something is updated (with max iteration)
        {
            if((ECFams->at(i).getGfamily1() == GfamIdToReset) || (ECFams->at(i).getGfamily2() == GfamIdToReset)) //this ECF is concerned by 
            {
                if(tryECFreDraw(i , Temperature, CoEventScore, GeneFamilyList, ECFams, CoEventSet, alwaysGainAtTop, C1Advantage, superverbose, DupliCost, LossCost, HGTCost, RecWeight) )
                { // we have chosen the new solution for this ECF -> update coev score
                    CoEventScore =  computeCoEventScore(*CoEventSet, DupliCost, LossCost, HGTCost, RecWeight);
                }
            }
        }
    }

    if(superverbose)
    {
        cout << "nb of coevents" << CoEventSet->size() << endl;
        for(unsigned i = 0 ; i < CoEventSet->size(); i++)
        {
            cout << CoEventSet->at(i).getEvent() << " ";
        }
        cout << endl;
        double CoEventScore =  computeCoEventScore(*CoEventSet, DupliCost, LossCost, HGTCost, RecWeight);
        cout << "coev score " << CoEventScore << endl;
    }
}


/*
Replace the reconciled tree by the given one and recompute adjacency histories and co-events

Takes:
- int GfamIdToReset : id of the gene family whose rconciliation to redraw
- ReconciledTree * newRTree : reconciled tree to put in place of the old one
- double treeLkh : likelihood of that tree
- vector < GeneFamily * > * GeneFamilyList : List of Gene Families
- vector < EquivalenceClassFamily > * ECFams : list of equivalence class families
- vector < CoEvent > * CoEventSet : set of co-events to update
- bool withTransfer : whether there is transfer or not
- double DupliCost : cost of a single duplication
- double HGTCost : cost of a single lateral gene transfer
- double LossCost : cost of a single gene loss
- int maxTS : maximum time slice in the species tree
- double TopoWeight : weight of the topology in the system score
- double Again : cost of a single adjacency gain
- double Abreak : cost of a single adjacency break
- double absencePenalty : penalty for an absence of adjacency
- bool doAllPair : whether to do all possible equivalence class, even when they don't have any adjacency
- bool substractRecoToAdj : whether to substract the weighted cost of a single adjacency to the score of a co-event in the adjacency matrix
- double weightedDupCost : cost of a weighted single gene duplication
- double weightedLossCost : cost of a weighted single gene loss
- double weightedHGTCost : cost of a weighted single lateral gene transfer
- bool alwaysGainAtTop : whether to always add the gain of an adjacency gain at the top of the adjacency trees
- double C1Advantage : probability to choose C1 over C0  IF both have thee same cost at the "root" of the equivalence class matrix
- bool verbose : prints some info
- bool superverbose : prints tons of info (debugging only)

*/
void ReplaceGeneFamilyRecTree(int GfamIdToReset, ReconciledTree * newRTree, double treeLkh, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                        bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        double Again, double Abreak, double absencePenalty, bool doAllPair, bool substractRecoToAdj,
                        double weightedDupCost, double weightedLossCost, double weightedHGTCost, bool alwaysGainAtTop, double C1Advantage,
                        bool verbose , bool superverbose,
                        map<int,vector<float> > speciesC0C1,
                        map<int, map<string,int> > speGeneAdjNb,
                        bool interactionMode 
                         )
{
    //1. prepping
    //clock_t begin;
    //clock_t end;
    //double elapsed_secs;

    //if(verbose)
    //{
    //    begin = clock();
    //}

    prepRecRedrawEqClass(GfamIdToReset, ECFams);
    prepRecRedrawCoEvent(GfamIdToReset, CoEventSet, superverbose);
    

    //if(verbose)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "time elapsed for preppin:" << elapsed_secs << endl;
    //    begin = clock();
    //}

    
    //2. replacing the reconciliation
    GeneFamilyList->at(GfamIdToReset)->setReconciliation(*newRTree);
    GeneFamilyList->at(GfamIdToReset)->setRecScore(DupliCost, HGTCost, LossCost);

    GeneFamilyList->at(GfamIdToReset)->setTreeLikelihood(treeLkh);

    //if(verbose)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "time elapsed for reconciliation:" << elapsed_secs << endl;
    //    begin = clock();
    //}


    //3. refining and computing equivalence classes
    
    refiningAndComputingECFamilyList( ECFams, GeneFamilyList,
                                                GfamIdToReset,
                                                Again, Abreak,
                                                absencePenalty, doAllPair , 
                                                withTransfer , substractRecoToAdj, 
                                                weightedDupCost, weightedHGTCost, weightedLossCost,
                                                alwaysGainAtTop, C1Advantage,
                                                superverbose,false,
                                                speciesC0C1,
                                                speGeneAdjNb,
                                                interactionMode 
                                                );


    //if(verbose)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "time elapsed for equivalence class:" << elapsed_secs << endl;
    //    begin = clock();
    //}

    //4. Co-Event Building
    delete CoEventSet;//simpler for now to destroy and then rebuild from scratch...
    CoEventSet = new vector <CoEvent>;
    GettingCoEvents( CoEventSet, ECFams, GeneFamilyList );   

    //if(verbose)
    //{
    //    end = clock();
    //    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //    cout << "time elapsed for co-event:" << elapsed_secs << endl;
    // 
    //} 
}

/*
re-do the backtrack of an equivalence class families and generate approrpiate co-event

Takes:
 - int ECFidToReset : id of the ECF to reset
 - vector < GeneFamily * > * GeneFamilyList : list of gene families
 - vector < EquivalenceClassFamily > * ECFams : list of Equivalence Class Families
 - vector < CoEvent > * CoEventSet : list of co-events !!!!!!!!!!!!!!!!!!!!!!!!!! SHOULD BE AN EMPTY VECTOR !!!!!!!!!!!!!!!
 - bool alwaysGainAtTop :  whether to always add the gain of an adjacency gain at the top of the adjacency trees
 - double C1Advantage : probability to choose C1 over C0  IF both have thee same cost at the "root" of the equivalence class matrix
 - bool superverbose : prints tons of info (debugging only)

*/
void ReDrawECFamily(int ECFidToReset , 
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool alwaysGainAtTop, double C1Advantage,bool superverbose)
{

    //1.prep
    //prepECFRedrawCoEvent(ECFidToReset, CoEventSet, superverbose);

    //2.backtrack
    EquivalenceClassFamily * ECF = &ECFams->at(ECFidToReset);
    int gfam1 = ECF->getGfamily1();
    int gfam2 = ECF->getGfamily2();
 
    ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
    ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();

    bool useCount = false;
    ECF->backtrackAdjMatrixForSelf(Rtree1, Rtree2, 
                                    false, // boltzmann
                                    alwaysGainAtTop , C1Advantage, useCount);


    //3. Co-Event Building

    GettingCoEvents( CoEventSet, ECFams, GeneFamilyList );
    //cout <<  "CoEventSet" << CoEventSet->size() << endl;
}


/*
replace the adjacency forests of a given equivalence class family along with a set of co-events

Takes:
 - int ECFidToReset : id of the ECF to reset
 - vector < vector< AdjTree *> * > * NewAdjacencyTreesVector : adjacency forests to put in place of the current one of the ECF to reset
 - vector < CoEvent > * NewCoEventSet : set of co-vent to put in place of the existing one
 - vector < GeneFamily * > * GeneFamilyList : list of Equivalence Class Families
 - vector < EquivalenceClassFamily > * ECFams : list of co-events !!!!!!!!!!!!!!!!!!!!!!!!!! SHOULD BE AN EMPTY VECTOR !!!!!!!!!!!!!!!
 - vector < CoEvent > * CoEventSet :  whether to always add the gain of an adjacency gain at the top of the adjacency trees
 - bool superverbose : prints tons of info (debugging only)

*/
void ReplaceECFamily(int ECFidToReset , vector < vector< AdjTree *> * > * NewAdjacencyTreesVector, vector < CoEvent > * NewCoEventSet,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool superverbose)
{
    for(unsigned i = 0 ; i <  NewAdjacencyTreesVector->size(); i++)
    {
        ECFams->at(ECFidToReset).setAdjForest(i, NewAdjacencyTreesVector->at(i));
    } 

    delete CoEventSet;
    CoEventSet = NewCoEventSet;
    if(CoEventSet->size() == 0)
        GettingCoEvents( CoEventSet, ECFams, GeneFamilyList );
}


/*
Chooses whether or not to accept a solution over an old one base on their score respective score difference

acceptation probability = exp( (scoreInit - scoreProposed) / Temperature )
NB:  if scoreInit > scoreProposed, the new solution is automatically accepted

Takes:
- double scoreInit : score of the initial solution
- double scoreProposed : score of the ne solution
- double Temperature : high temp -> worst solution have a higher chance of being selected

Returns:
    (bool): treu if the the solution is accepted
*/
bool acceptProposition(double scoreInit, double scoreProposed, double Temperature)
{

    
    if(scoreInit > scoreProposed)
        return true;

    if(Temperature <= 0)
        return false;

    double Paccept = exp( (scoreInit - scoreProposed) / Temperature ) ;
    double r = ((double) rand()/ RAND_MAX); //random results

    if( r < Paccept)
        return true;
    return false;
}


/*
Tries to re-do the backtrack of an equivalence class families and generate appropriate co-event
and chooses whether to keep the new solution or not

Takes:
 - int ECFidToReset : id of the ECF to reset
 - double Temperature :  high temp -> worst solution have a higher chance of being selected
 - double oldCoEventScore : co-event score of the previous (current) solution; used to choose whether to accept or no the new solution
 - vector < GeneFamily * > * GeneFamilyList : list of gene families
 - vector < EquivalenceClassFamily > * ECFams : list of equivalence class families
 - vector < CoEvent > * CoEventSet : list of co-events
 - bool alwaysGainAtTop : whether to always add the gain of an adjacency gain at the top of the adjacency trees
 - double C1Advantage : probability to choose C1 over C0  IF both have thee same cost at the "root" of the equivalence class matrix
 - bool superverbose : prints tons of info (debugging only)
 - double DupliCost : cost of a single duplication
 - double LossCost : cost of a single gene loss
 - double HGTCost : cost of a single lateral gene transfer
 - double RecWeight : weight of the reconciliation in the global system score

Returns:
    (bool): treu if the the solution is accepted
*/
bool tryECFreDraw(int ECFidToReset , double Temperature, double oldCoEventScore,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool alwaysGainAtTop, double C1Advantage,bool superverbose,
                    double DupliCost, double LossCost, double HGTCost, double RecWeight)
{
    if(superverbose)
    {
        cout << "nb of coevents" << CoEventSet->size() << endl;
        for(unsigned i = 0 ; i < CoEventSet->size(); i++)
        {
            cout << CoEventSet->at(i).getEvent() << " ";
        }
        cout << endl;
    }

    //1. keeping the initial solution

    vector < CoEvent > * INITCoEventSet = CoEventSet;
    
    bool deleteCoev = false;
    
    int nbEclass = ECFams->at(ECFidToReset).getNbEqClasses();

    vector < vector< AdjTree *> * > * AdjacencyTreesVectorINIT = new vector < vector< AdjTree *> * >;
    for(unsigned i = 0; i < nbEclass; i++)
        AdjacencyTreesVectorINIT->push_back( ECFams->at(ECFidToReset).getAdjForest(i) );
    
    //2. redrawing
    if(deleteCoev)
        delete CoEventSet;//simpler for now to destroy and then rebuild from scratch...

    CoEventSet = new vector <CoEvent>;


    ReDrawECFamily(ECFidToReset , 
                     GeneFamilyList,  ECFams,  CoEventSet,
                     alwaysGainAtTop,  C1Advantage, superverbose);

    if(superverbose)
        cout <<  "CoEventSet" << CoEventSet->size() << endl;

    //3. getting new score. Only the co-event score changes...
    double CoEventScore =  computeCoEventScore(*CoEventSet, DupliCost, LossCost, HGTCost, RecWeight);
    
    if(superverbose)
        cout << "new coevent score for ECF " << ECFidToReset << " : " << CoEventScore << " <-> " << oldCoEventScore << endl; 
        

    //4. accepting or not.
    bool accepted = acceptProposition(oldCoEventScore, CoEventScore, Temperature);

    if(!accepted)
    {
        //going back
        ReplaceECFamily( ECFidToReset , AdjacencyTreesVectorINIT, INITCoEventSet,
                     GeneFamilyList,  ECFams,  CoEventSet,
                     superverbose);

    }
    else if(superverbose)
        cout << "accepted new solution for ECFam " << ECFidToReset << endl;

    return accepted;
}


GfamSave createHistoric(int GfamIdToReset , 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * * CoEventSetPointer,
                        double Again, double Abreak,
                        double DupliCost, double HGTCost, double LossCost,
                        double RecWeight
                        )
{
    //1. keeping the initial solution
    GfamSave historic;

    

    //ReconciledTree OldRecTree = GeneFamilyList->at(GfamIdToReset)->getRecTreeCopy();
    //historic.recTree = &OldRecTree;


    historic.recTree = new ReconciledTree();
    *(historic.recTree) = GeneFamilyList->at(GfamIdToReset)->getRecTreeCopy() ;

    historic.treeLkh = GeneFamilyList->at(GfamIdToReset)->getTreeLikelihood();
    
    

    //select the ECF that will be re-drawn
    
    //cout << "GfamIdToReset " << GfamIdToReset << endl;

    historic.AdjScore = 0;//only for the changed ecfs

    for(size_t i = 0; i < ECFams->size(); i++)
    {
        if((ECFams->at(i).getGfamily1() == GfamIdToReset) || (ECFams->at(i).getGfamily2() == GfamIdToReset))
        {
            //cout << " "<< i << "->" << ECFams->at(i).getGfamily1() << '-' << ECFams->at(i).getGfamily2() << endl;

            historic.ECFamsMap[i] = ECFams->at(i);
            historic.AdjScore += ECFams->at(i).getNbAdjGain() * Again;
            historic.AdjScore += ECFams->at(i).getNbAdjBreak() * Abreak;
            //cout << "getting score for ecf " << i << endl;
            //cout << " nbtrees " << historic.ECFamsMap[i].getNbAdjTrees() << endl;
            //cout << historic.ECFamsMap[i].getNbAdjGain() <<"*" <<Again << endl;
            //cout << historic.ECFamsMap[i].getNbAdjBreak() <<"*"<<Abreak << endl;
        }
    }

    //cout << "ECF in historic: " << historic.ECFamsMap.size() <<" / " << ECFams->size()<< endl;


    historic.CoEventSet = *CoEventSetPointer;

    

    historic.RecScore = GeneFamilyList->at(GfamIdToReset)->getRecScore();//only for the changed tree
    
    historic.CoEventScore = computeCoEventScore(*(*CoEventSetPointer), DupliCost, LossCost, HGTCost, RecWeight);

    return historic;
}


/*

Takes:
 - int GfamIdToReset : index of  the gene family to try to change
 - int redrawAlgo : 0:change reconiliation only; 1:sample in ccp; 2:DLRecCoev
 - bool trySeveralECFSolutions : if several backtracking solutions should be tried 
 - double Temperature : temperature to choose new solution (a higher the temperature means that solution with a higher score will have a higher probability to be chosen)
 - vector < GeneFamily * > * GeneFamilyList : list of gene family
 - vector < EquivalenceClassFamily > * ECFams : list of equivalence class families
 - vector < CoEvent > * * CoEventSetPointer : set of co-events
 - MySpeciesTree * speciesTree : species tree
 - int maxTS : maximum time slice of the species tree
 - bool datedSpeciesTree : false if the species tree is undated, true if it is dated
 - bool boundedTS : if true, bounded time slice will be used
 - bool withTransfer : wether there is transfers or not
 - double DupliCost : cost of a single gene duplication
 - double HGTCost : cost of a single gene transfer
 - double LossCost : cost of a single loss
 - double TopoWeight : weight of the topology score in the global system score
 - double RecWeight : weight of the reconciliation score in the global system score
 - double AdjWeight : weight of the adjacency score in the global system score
 - double Again : cost of a single adjacency gain
 - double Abreak : cost of a single adjacency breakage
 - double absencePenalty : penality to apply in the case of an absence of adjacency
 - bool doAllPair : if true, all possible pairs of gene family are considered as equivalence class families, even if no adjacencies support them
 - bool substractRecoToAdj : if true, the weighted cost of an event is substracted to the cost of a potential co-event in the adjacency matrices
 - bool alwaysGainAtTop : wether to aautomatically put a gain at the top of any adjacency tree or not
 - double C1Advantage : probability of choosing presence over absence of an adjacency at the root of an adjacency matrix if both have the same probability
 - bool verbose : prints some info
 - bool superverbose  : prints a lot of info
 - double probaC1Bias : [default =0] probability to sample the tree topology in a CCP biased by the adjacency trees obtained before
 - MyCladesAndTripartitions * biasedCCP : [default = NULL] CCP distribution biased by the adjacency trees obtained before


*/
bool tryGfamRedraw(int GfamIdToReset, int redrawAlgo, bool trySeveralECFSolutions, double Temperature, 
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * * CoEventSetPointer, MySpeciesTree * speciesTree,
                    int maxTS, bool datedSpeciesTree, bool boundedTS, bool withTransfer, 
                    double DupliCost, double HGTCost, double LossCost,
                    double TopoWeight, double RecWeight, double AdjWeight,
                    double Again, double Abreak,
                    double absencePenalty, bool doAllPair, bool substractRecoToAdj, 
                    bool alwaysGainAtTop, double C1Advantage,
                    bool verbose, bool superverbose,
                    map<int,vector<float> > speciesC0C1,
                    map<int, map<string,int> > speGeneAdjNb,
                    double probaC1Bias , MyCladesAndTripartitions * biasedCCP ,
                    string gPathPrefix,
                    bool interactionMode
                    )
{
    
    //1. keeping the initial solution
    double normTreeLkhINIT = - log( GeneFamilyList->at(GfamIdToReset)->getNormalizedTreeLikelihood());

    GfamSave historic =  createHistoric(GfamIdToReset , 
                                        GeneFamilyList, ECFams, CoEventSetPointer,
                                        Again, Abreak,
                                        DupliCost, HGTCost, LossCost,
                                        RecWeight
                                        );

    double historicScore = normTreeLkhINIT * TopoWeight + historic.RecScore * RecWeight + historic.AdjScore * AdjWeight + historic.CoEventScore ;

    //cout << "historic.ECFamsMap.size() " << historic.ECFamsMap.size() << endl;

    *CoEventSetPointer = new vector< CoEvent >;
    //cout << "historic.CoEventSet " << historic.CoEventSet->size() << "<>"<< historic.CoEventScore <<  "-" << CoEventSet->size() <<endl;

    bool probaFlat = 0; //dummy arg

    //2. make the proposition
    ReDrawGeneFamily( GfamIdToReset, redrawAlgo, 
                        GeneFamilyList, ECFams, *CoEventSetPointer,
                        speciesTree, datedSpeciesTree, boundedTS, withTransfer,
                        DupliCost, HGTCost, LossCost, maxTS, TopoWeight,
                        Again, Abreak, absencePenalty, doAllPair, substractRecoToAdj,
                        DupliCost * RecWeight / AdjWeight, LossCost * RecWeight / AdjWeight, HGTCost * RecWeight / AdjWeight,
                        alwaysGainAtTop, C1Advantage,
                        verbose, superverbose ,
                        trySeveralECFSolutions, Temperature, RecWeight, false,
                        speciesC0C1,
                        speGeneAdjNb, 
                        probaFlat,
                        probaC1Bias , biasedCCP ,
                        gPathPrefix,
                        interactionMode
                        );

    //for(unsigned i = 0 ; i < (*CoEventSetPointer)->size(); i++)
    //{
    //    cout << (*CoEventSetPointer)->at(i).getEvent() << " ";
    //}
    //cout << endl;
    //double CoEventScore =  computeCoEventScore(*(*CoEventSetPointer), DupliCost, LossCost, HGTCost, RecWeight);
    //cout << "coev score " << CoEventScore << endl;


    //3. compute the new scores.

    double topoScore = (- log( GeneFamilyList->at(GfamIdToReset)->getNormalizedTreeLikelihood())) * TopoWeight;
    double recScore = RecWeight * GeneFamilyList->at(GfamIdToReset)->getRecScore();//only for the changed tree
    //cout << "recScore" << RecWeight << "*" << GeneFamilyList->at(GfamIdToReset)->getRecScore() << endl;//only for the changed tree


    double adjScore = 0;

    map<int, EquivalenceClassFamily>::iterator it;
    for(it = historic.ECFamsMap.begin() ; it != historic.ECFamsMap.end() ; ++it)
    {
        int ecId = it->first;

        adjScore += ((double) ECFams->at(ecId).getNbAdjGain() * Again);
        adjScore += ((double) ECFams->at(ecId).getNbAdjBreak() * Abreak);
    }

    adjScore *= AdjWeight;

    double coEventScore = computeCoEventScore(*(*CoEventSetPointer), DupliCost, LossCost, HGTCost, RecWeight);

    if(verbose)
        cout << "comparing scores for gfam " << GfamIdToReset << " : " << topoScore << " + " << recScore << " + " << adjScore << " + " <<  coEventScore << " = " << topoScore + recScore + adjScore + coEventScore << " <-> " << historicScore << endl; 

    //4. accepting or no the new proposition
    bool accepted = acceptProposition( historicScore ,
                                    topoScore + recScore + adjScore + coEventScore,
                                    Temperature);



    if(!accepted)
    {
        //5. return to the old solution
        GeneFamilyList->at(GfamIdToReset)->setReconciliation(*(historic.recTree));
        GeneFamilyList->at(GfamIdToReset)->setRecScore(DupliCost, HGTCost, LossCost);

        GeneFamilyList->at(GfamIdToReset)->setTreeLikelihood(historic.treeLkh);

        //cout << "rec nodes in historic" << historic.recTree->getNumberOfNodes() << endl;
        delete historic.recTree;
        //cout << "rec nodes in solution" << GeneFamilyList->at(GfamIdToReset)->getRecTree()->getNumberOfNodes() << endl;

        //ECF
        //cout << "copying map of size "<<historic.ECFamsMap.size() << endl;
        for(it = historic.ECFamsMap.begin() ; it != historic.ECFamsMap.end() ; ++it)
        {
            //cout << it->first << " "<< it->second.getNbAdjTrees() << endl;

            int ecId = it->first;

            //cout << "assign"<<endl;
            //ECFams->at(ecId) = it->second; // CLONEMOD
            ECFams->at(ecId) = EquivalenceClassFamily();
            //cout << "cloning"<<endl;
            ECFams->at(ecId).clone( &(it->second) );

            //cout << "replaced ECF "<<  ecId << " "<< ECFams->at(ecId).getNbAdjTrees() << endl;
        }

        //cout << "clearing map of size " << historic.ECFamsMap.size() <<endl;
        historic.ECFamsMap.clear();

        (*CoEventSetPointer)->clear();

        delete *CoEventSetPointer;

        *CoEventSetPointer = historic.CoEventSet; // MEMWORK


        if(verbose)
            cout << "Not accepted new solution for Gfam " << GfamIdToReset << endl;

    }
    else 
    {
        if(verbose)
            cout << "accepted new solution for Gfam " << GfamIdToReset << endl;


        //cout << "rec nodes in historic" << historic.recTree->getNumberOfNodes() << endl;
        historic.clearP();
        delete historic.recTree;
        //historic.CoEventSet->clear();
        //delete historic.CoEventSet;


        //cout << historic.CoEventSet->size() << endl;
        //cout << "<<<<<" << endl;

    }

    return accepted;
}

/*

Takes:
 - int GfamIdToReset : index of  the gene family to try to change
 - int maxNbTry : maximum number of try to attempt
 - int redrawAlgo : 0:change reconiliation only; 1:sample in ccp; 2:DLRecCoev
 - bool trySeveralECFSolutions : if several backtracking solutions should be tried 
 - double Temperature : temperature to choose new solution (a higher the temperature means that solution with a higher score will have a higher probability to be chosen)
 - vector < GeneFamily * > * GeneFamilyList : list of gene family
 - vector < EquivalenceClassFamily > * ECFams : list of equivalence class families
 - vector < CoEvent > ** CoEventSetPointer : set of co-events
 - MySpeciesTree * speciesTree : species tree
 - int maxTS : maximum time slice of the species tree
 - bool datedSpeciesTree : false if the species tree is undated, true if it is dated
 - bool boundedTS : if true, bounded time slice will be used
 - bool withTransfer : wether there is transfers or not
 - double DupliCost : cost of a single gene duplication
 - double HGTCost : cost of a single gene transfer
 - double LossCost : cost of a single loss
 - double TopoWeight : weight of the topology score in the global system score
 - double RecWeight : weight of the reconciliation score in the global system score
 - double AdjWeight : weight of the adjacency score in the global system score
 - double Again : cost of a single adjacency gain
 - double Abreak : cost of a single adjacency breakage
 - double absencePenalty : penality to apply in the case of an absence of adjacency
 - bool doAllPair : if true, all possible pairs of gene family are considered as equivalence class families, even if no adjacencies support them
 - bool substractRecoToAdj : if true, the weighted cost of an event is substracted to the cost of a potential co-event in the adjacency matrices
 - bool alwaysGainAtTop : wether to aautomatically put a gain at the top of any adjacency tree or not
 - double C1Advantage : probability of choosing presence over absence of an adjacency at the root of an adjacency matrix if both have the same probability
 - bool verbose : prints some info
 - bool superverbose  : prints a lot of info
 - double probaC1Bias : [default =0] probability to sample the tree topology in a CCP biased by the adjacency trees obtained before
 - MyCladesAndTripartitions * biasedCCP : [default = NULL] CCP distribution biased by the adjacency trees obtained before

*//*
bool updateGfam(int GfamIdToReset, int maxNbTry, int redrawAlgo, bool trySeveralECFSolutions, double Temperature,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * * CoEventSetPointer, MySpeciesTree * speciesTree,
                    int maxTS, bool datedSpeciesTree, bool boundedTS, bool withTransfer, 
                    double DupliCost, double HGTCost, double LossCost,
                    double TopoWeight, double RecWeight, double AdjWeight,
                    double Again, double Abreak,
                    double absencePenalty, bool doAllPair, bool substractRecoToAdj, 
                    bool alwaysGainAtTop, double C1Advantage,
                    bool verbose, bool superverbose,
                    map<int,vector<float> > speciesC0C1,
                    map<int, map<string,int> > speGeneAdjNb,
                    double probaC1Bias, MyCladesAndTripartitions * biasedCCP,
                    string gPathPrefix
                    )
{
    

    int nbtry = 0;
    bool accept = false;

    clock_t begin;
    clock_t end;
    double elapsed_secs;
    if(verbose)
    {
        begin = clock();
    }


    while(nbtry < maxNbTry)
    {
        accept = tryGfamRedraw( GfamIdToReset, redrawAlgo, trySeveralECFSolutions, Temperature,
                                GeneFamilyList, ECFams, CoEventSetPointer, speciesTree,
                                maxTS, datedSpeciesTree, boundedTS, withTransfer, 
                                DupliCost , HGTCost, LossCost,
                                TopoWeight, RecWeight, AdjWeight,
                                Again, Abreak,
                                absencePenalty, doAllPair, substractRecoToAdj,
                                alwaysGainAtTop , C1Advantage,
                                verbose, superverbose,
                                speciesC0C1,
                                speGeneAdjNb,
                                probaC1Bias, biasedCCP,
                                gPathPrefix
                                 );
        if(accept)
            break;

        //for(unsigned i = 0 ; i != ECFams->size() ; i++)
        //{
        //    cout << "ATree ECF "<<  i << " "<< ECFams->at(i).getNbEqClasses()<< " "<< ECFams->at(i).getNbAdj() << " "<< ECFams->at(i).getNbAdjTrees() << endl;
        //}

        if(superverbose)
        {
            //5. score update
            double SystemScore = 0;
        
            double TopoScore = computeTopoScore( *GeneFamilyList, TopoWeight);
            double ReconScore = computeReconciliationScore( *GeneFamilyList, RecWeight );
            double AdjScore =  computeAdjacenciesScore( ECFams,Again, Abreak , AdjWeight);
            double CoEventScore =  computeCoEventScore(*(*CoEventSetPointer), DupliCost , LossCost , HGTCost , RecWeight );
        
            SystemScore = TopoScore + ReconScore + AdjScore + CoEventScore ;
        
    
            cout << "current score : " << TopoScore << " + " << ReconScore << " + " << AdjScore << " + " << CoEventScore << " = " << SystemScore  << endl;
        }
        nbtry++;
    }

    if(verbose)
    {
        if(nbtry == maxNbTry)
            cout << "did not accept after " << nbtry +1 << "try."<< endl;
        else
            cout << "accepted after try " << nbtry +1 << endl;
    }


    //cout << "nb of coevents" << (*CoEventSetPointer)->size() << endl;
    //for(unsigned i = 0 ; i < (*CoEventSetPointer)->size(); i++)
    //{
    //    cout << (*CoEventSetPointer)->at(i).getEvent() << " ";
    //}
    //cout << endl;
    //double CoEventScore =  computeCoEventScore(*(*CoEventSetPointer), DupliCost, LossCost, HGTCost, RecWeight);
    //cout << "coev score " << CoEventScore << endl;

    //// not elegant...
    //*CoEventSetPointer = new vector<CoEvent>;
    //GettingCoEvents( *CoEventSetPointer, ECFams, GeneFamilyList );


    //cout << "nb of coevents" << (*CoEventSetPointer)->size() << endl;
    //for(unsigned i = 0 ; i < (*CoEventSetPointer)->size(); i++)
    //{
    //    cout << (*CoEventSetPointer)->at(i).getEvent() << " ";
    //}
    //cout << endl;
    //CoEventScore =  computeCoEventScore(*(*CoEventSetPointer), DupliCost, LossCost, HGTCost, RecWeight);
    //cout << "coev score " << CoEventScore << endl;


    if(verbose)
    {
        //5. score update
        double SystemScore = 0;
    
        double TopoScore = computeTopoScore( *GeneFamilyList, TopoWeight);
        double ReconScore = computeReconciliationScore( *GeneFamilyList, RecWeight );
        double AdjScore =  computeAdjacenciesScore( ECFams,Again, Abreak , AdjWeight);
        double CoEventScore =  computeCoEventScore(*(*CoEventSetPointer), DupliCost , LossCost , HGTCost , RecWeight );
    
        SystemScore = TopoScore + ReconScore + AdjScore + CoEventScore ;
    

        cout << "current score : " << TopoScore << " + " << ReconScore << " + " << AdjScore << " + " << CoEventScore << " = " << SystemScore  << endl;
    }


    if(verbose)
    {
        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout << "time elapsed for sampling:" << elapsed_secs << "s. That's "<< elapsed_secs / (nbtry+1) << "s / try" << endl;
    }


    return  accept;
}*/

/*
Takes:
 - string filename: name of the file containing the trees
 - int gfam1 : id of the first gfam
 - int gfam2 : id of the second gfam
 - bool gainAtRoot : wetehr there is a gain at the root
 - bool VERBOSE
*/
vector <AdjTree * > * readAdjForest(string filename, int gfam1, int gfam2, bool gainAtRoot, bool VERBOSE)
{
 
    vector <AdjTree * > * Aforest = new vector <AdjTree * > ;


    ifstream fileStream(filename.c_str());
    
    if( !fileStream.is_open() ) 
    {
        cout << "Could not open adjacency tree file : "<< filename <<endl;
        return Aforest;
        //exit(1);
    }

    string line;
    getline( fileStream, line );

    while( !fileStream.eof() ) 
    {

        
        if( line == "" ) 
        {
                // ignore blank lines
        }
        else
        {

            map <string, string> * properties = new map <string,string>;
            string * value = new string(""); 
            string Tname = InterpretLineXML(line,properties,value);

            if(Tname.compare("clade") == 0) //--> line signing the beginning of a tree
            {
                Aforest->push_back( new AdjTree(fileStream, gainAtRoot, gfam1, gfam2, VERBOSE) );
                if(VERBOSE)
                    cout << "ADDED A TREE"<<endl;
            }

            delete properties;
            delete value;
        }
        getline( fileStream, line );
        //cout << line << endl;
    }

    fileStream.close();

    return Aforest;
}



/*
Takes:
 - string filename: name of the file containing the trees
 - vector< EquivalenceClassFamily > * ECFams
 - bool gainAtRoot : whether there is a gain at the root
 - bool VERBOSE

returns:
 (bool) : true if no pb ; false otherwise
*/
bool readECFamTreesOneFile(string filename, vector< EquivalenceClassFamily > * ECFams, bool gainAtRoot, bool VERBOSE)
{

    ifstream fileStream(filename.c_str());
    
    if( !fileStream.is_open() ) 
    {
        cerr << "Could not open adjacency tree file : "<< filename  <<endl;

        return false;
    }

    string line;
    getline( fileStream, line );

    vector <AdjTree * > * Aforest = NULL;

    int current = -1;

    int fam1 = -1 ;
    int fam2 = -1 ;

    while( !fileStream.eof() ) 
    {

        
        if( line == "" ) 
        {
                // ignore blank lines
        }
        else
        {

            map <string, string> * properties = new map <string,string>;
            string * value = new string(""); 
            string Tname = InterpretLineXML(line,properties,value);


            if(Tname.compare("/EquivalenceClassFamily") == 0)
            {
                if(current != -1)
                {
                    ECFams->at(current).setAdjForest( 0 , Aforest );
                    if(VERBOSE)
                        cout << "read adj forest for ecf between " << ECFams->at(current).getGfamily1() << "-" << ECFams->at(current).getGfamily2() << endl;
                }
                fam1 = -1 ;
                fam2 = -1 ;
                current = -1;
                Aforest = NULL;
            }
            else if(Tname.compare("EquivalenceClassFamily") == 0)
            {
                Aforest = new vector <AdjTree * > ;
                char * pEnd;
                fam1 =  (int) strtol(properties->find("fam1")->second.c_str(), &pEnd, 10) ;
                fam2 =  (int) strtol(properties->find("fam2")->second.c_str(), &pEnd, 10) ;
                if(VERBOSE)
                    cout << "reading ecf " << fam1  << "-" << fam2 << endl;

                for(unsigned i = 0 ; i < ECFams->size() ; i++ )
                {
                    int f1 = ECFams->at(i).getGfamily1();
                    int f2 = ECFams->at(i).getGfamily2();
                    
                    if( ( f1 == fam1 )&&( f2 == fam2 ) )
                    {
                        current = i;
                        break;
                    }
                    else if( ( f1 == fam2 )&&( f2 == fam1 ) )
                    {
                        current = i;
                        break;
                    }
                }
                //cout <<  " --> index " << current << endl;

            }
            else if(Tname.compare("clade") == 0) //--> line signing the beginning of a tree
            {
                Aforest->push_back( new AdjTree(fileStream, gainAtRoot, fam1, fam2, VERBOSE) );
                if(VERBOSE)
                    cout << "ADDED A TREE"<<endl;
            }

            delete properties;
            delete value;
        }
        getline( fileStream, line );
        //cout << line << endl;
    }

    fileStream.close();
    return true;
}


/*
* @arg vector < GeneFamily * > * GeneFamilyList 
* @arg string prefix : the prefix (containing path) of the filename
* @arg MySpeciesTree * Stree : pointer to the species tree
* @arg bool VERBOSE : write loading details
*/
void LoadGeneFamilyReconciledTrees( vector < GeneFamily * > * GeneFamilyList, string prefix, MySpeciesTree * Stree, bool VERBOSE , bool OneFile )
{


    if(!OneFile)
    {
        for(unsigned i = 0 ; i < GeneFamilyList->size(); i++)
        {
            string filename = prefix + "geneFamily" + static_cast<ostringstream*>( &(ostringstream() << i) )->str() ;
            filename += ".phyloxml";
    
            if(VERBOSE)
                cout << "Reading in file " << filename << endl;
    
            ReconciledTree Rtree(filename , Stree, VERBOSE);
    
            GeneFamilyList->at(i)->setReconciliation(Rtree);
    
        }
    }
    else
    {
        string fileName = prefix ;
        if( fileName.length() >0 )
        {
            if(fileName[ fileName.length() -1 ] != '/')
                fileName += ".";
        }
        fileName += "reconciliations.xml";
        cout << "reading in "<< fileName << endl;

        ifstream fileStream(fileName.c_str());
        if( !fileStream.is_open() ) 
        {
            if(VERBOSE)
                cout << "LoadGeneFamilyReconciledTrees : could not open reconciled tree file : " << fileName << " will try on separated files."<< endl;

            LoadGeneFamilyReconciledTrees( GeneFamilyList, prefix, Stree,  VERBOSE , false );

            return;
        }

        string RecGeneTreeTag = "recGeneTree";
        string startTag = "clade";

        int current = 0;

        while( goToNextTag(fileStream,RecGeneTreeTag) )  //goes to the next RecGeneTreeTag -> next reconciled gene tree
        {
            if( goToNextTag(fileStream,startTag) ) // go to clade -> the root of the reconciled gene tree
            {
                ReconciledTree Rtree(fileStream, Stree, VERBOSE);
                GeneFamilyList->at(current)->setReconciliation(Rtree);
                current++;
            }
            else
            { // no clade ? pb but ignore...
                if(VERBOSE)
                    cout << "found a "<< RecGeneTreeTag << " without any " << startTag << " in the file " << fileName << " -> ignoring that tree."<< endl;
            }
        }

        if( current < GeneFamilyList->size() )
        {
            cerr << "!!ERROR!! while loading the reconciled gene trees xml file. not enough recGeneTree (" << current << " found ; "<< GeneFamilyList->size() << " expected)" << endl;
        }
    }
}




/*
* @arg EquivalenceClassFamily * ECF : pointer the the ECF whose trees we want to write
* @arg string prefix : the prefix (containing path) of the filename
* @arg bool gainAtRoot : wether there is automatically a gain at the root of the trees or not
* @arg bool VERBOSE
*/
void LoadECFamTrees(EquivalenceClassFamily * ECF, string prefix,bool gainAtRoot, bool VERBOSE)
{



    int g1 = ECF->getGfamily1();
    int g2 = ECF->getGfamily2();

    string filePrefix = prefix + "EqClass_" + static_cast<ostringstream*>( &(ostringstream() << g1) )->str() + "-" + static_cast<ostringstream*>( &(ostringstream() << g2) )->str() ;



    for(size_t j = 0; j < ECF->getNbEqClasses() ; j++)//NB: with that method, some file might not contain any tree...
    {
        string filename = filePrefix + "_" + static_cast<ostringstream*>( &(ostringstream() << j) )->str() + ".adjtrees";

        filename += ".phyloxml";

        if(VERBOSE)
            cout << "Reading in file " << filename << endl;

        vector <AdjTree * > * newAForest = readAdjForest(filename, g1, g2, gainAtRoot, VERBOSE);

        if(VERBOSE)
            cout << "Finished reading in file " << filename << endl;


        ECF->setAdjForest( j, newAForest);

    }
}

/*
get the list of gene families index 
that share at least an adjacency with the given gene family
and corresponding Equivalence class families

Takes:
 - int GeneFamilyId : a gene family index
 - vector < EquivalenceClassFamily > * ECFams : pointer to the vector of equivalence class families
 - vector <int> &LinkedGFs  : address of the vector that will contains the indexes of linfes GFs
 - vector <int> &LinkedECFs : address of the vector that will contains the indexes of linfes ECFs
 - bool allowSelf [default = false] : if true, allow the given GeneFamilyId to be part of the result (in case of an adjacency between two genes of the same family)
*/
void getLinkedFamilies(int GeneFamilyId, vector < EquivalenceClassFamily > * ECFams, vector <int> &LinkedGFs , vector <int> &LinkedECFs, bool allowSelf )
{

    for(unsigned i = 0 ; i < ECFams->size();i++)
    {
        bool add = false;
        int LinkedGfamIndex = -1;
        if( ECFams->at(i).getGfamily1() == GeneFamilyId)
        {
            if(!allowSelf)
            {
                if(ECFams->at(i).getGfamily2() == GeneFamilyId)// self link -> ignore
                {
                    continue; 
                }
            }
            LinkedGfamIndex = ECFams->at(i).getGfamily2();
            add = true;
        }
        else if( ECFams->at(i).getGfamily2() == GeneFamilyId )
        {
            LinkedGfamIndex = ECFams->at(i).getGfamily1();
            add = true;
        }

        if(add)
        {
            LinkedECFs.push_back( i );
            LinkedGFs.push_back( LinkedGfamIndex );
        }
    }
}


/*
get the list of gene families index 
that share at least an adjacency with the given gene family
and corresponding Equivalence class families

Takes:
 - int GeneFamilyId : a gene family index
 - vector < EquivalenceClassFamily > * ECFams : pointer to the vector of equivalence class families
 - bool allowSelf [default = false] : if true, allow the given GeneFamilyId to be part of the result (in case of an adjacency between two genes of the same family)

Returns:
    int : index of the chosen linked family ; -1 if None was possible
*/
int chooseOneLinkedfamily(int GeneFamilyId, vector < EquivalenceClassFamily > * ECFams, bool allowSelf )
{
    vector <int> LinkedGFs;
    vector <int> LinkedECFs;

    getLinkedFamilies( GeneFamilyId,  ECFams, LinkedGFs , LinkedECFs,  allowSelf );

    if(LinkedGFs.size() == 0)
        return -1;

    int total = 0;
    vector<int> Ps;
    Ps.reserve(LinkedGFs.size());

    for(unsigned i = 0; i < LinkedGFs.size();i++ )
    {
        //cout << LinkedGFs[i] << "-> (" << ECFams->at(LinkedECFs[i]).getGfamily1() <<","<<ECFams->at(LinkedECFs[i]).getGfamily2()<<") : " << ECFams->at(LinkedECFs[i]).getNbAdj()<<endl;

        int nbAdjs = ECFams->at(LinkedECFs[i]).getNbAdj();
        total += nbAdjs;
        Ps.push_back(nbAdjs);
    }

    int r = rand() % total;
    //cout << r <<" << "<<total << endl;
    
    int i = 0;
    r -= Ps[i];
    while( r >= 0 )
    {
        i++;
        r -= Ps[i];
    }
    
    //cout << "chosen "<< i << "->"<<LinkedGFs[i]<<endl;
    return LinkedGFs[i];
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
void AddECFToFile(string fileName, EquivalenceClassFamily * ECF, 
                             bool newick,bool hideLosses,  
                             double GainCost ,  double BreakCost )
{
    int fam1 = ECF->getGfamily1();
    int fam2 = ECF->getGfamily2();

    int sens1 = ECF->getSens1();    
    int sens2 = ECF->getSens2();


    ofstream ofs;
    ofs.open(fileName.c_str(),ofstream::out | ofstream::app);
    DeCoOutputManager DOM;


    if( !newick )
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
        ofs << " sample: " << "1";
        int nbTrees = ECF->getNbAdjTrees();
        int nbGain = ECF->getNbAdjGain();
        int nbBreak = ECF->getNbAdjBreak();

        ofs << " #Trees: "<< nbTrees ;
        ofs << " #Gain: " << nbGain;
        ofs << " #Break: " << nbBreak;
        ofs << " Score: " << nbGain * GainCost + nbBreak * BreakCost ;
        ofs << endl;
    }
    else
    {
        ofs << "    <sample number=\""<<"1"<<"\">"<< endl;
    }


    int nbEqClass = ECF->getNbEqClasses();    
    for(size_t j = 0; j < nbEqClass ; j++)
    {
        DOM.WriteAdjForest(ofs, ECF->getAdjForest(j),newick, hideLosses, line_indent);
    }

    if(!newick)
    {
        ofs << "    </sample>"<< endl;
    }




    if( !newick )
    {
        ofs << "  </EquivalenceClassFamily>"<< endl;
    }

    ofs.close();
    return ;    
}

/*
Takes:
    - adjacencies (vector< pair <string,string > >)
    - GeneFamilyList (vector <GeneFamily *> *)
    - Verbose (bool)
    - SuperVerbose (bool)

Returns:
 - (vector < EquivalenceClassFamily > * ) : a pointer to a vector of refined Equivalence class families
*/
vector < EquivalenceClassFamily > * BasicCreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies, 
                                                                    vector <GeneFamily *> * GeneFamilyList, 
                                                                    bool Verbose, bool SuperVerbose)
{
    // 1. creating  EquivalenceClassFamilies containiners

    vector <EquivalenceClassFamily> * ECFams = new vector <EquivalenceClassFamily>;
    map < int, map <int , int > > Gfam1ToGfam2ToIndex; // to easily find the already created EquivalenceClasses

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
            ECindex = Gfam1ToGfam2ToIndex[gfam1][gfam2];
        }

        bool success = ECFams->at(ECindex).CheckAddAdj(currentadj.first,currentadj.second,gfam1,gfam2);//adding the adj to the correct Equivalence class
        if(SuperVerbose)
            if(!success)
                cout << "Could not add adj to the Equivalence Class " << ECindex << " " << ECFams->at(ECindex).getGfamily1() << " - " << ECFams->at(ECindex).getGfamily2() << endl;

    }


    return ECFams;
}