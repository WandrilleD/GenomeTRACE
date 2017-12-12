#ifndef WMODIF_UTIL_H_
#define WMODIF_UTIL_H_

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

Last modified the: 12-12-2017
by: Wandrille Duchemin

*/


#include "GeneFamily.h"
//#include "TopoBiasedGeneFamily.h"
#include "EquivalenceClassFamily.h"
#include "CoEvent.h"
#include "DeCoUtils.h"
#include "XMLUtils.h"
#include "DTLRecCoevWrapper.h"


class GfamSave
{
public:

    ReconciledTree * recTree;
    double treeLkh;
    double RecScore;
    map <int , EquivalenceClassFamily > ECFamsMap;
    double AdjScore;
    
    vector <CoEvent> * CoEventSet;
    double CoEventScore;

    GfamSave(){}
    ~GfamSave(){}

    void clearP()
    {
        //delete recTree;
        CoEventSet->clear();
        delete CoEventSet;
    }

};




void prepRecRedrawEqClass(int GeneFamilyId, vector < EquivalenceClassFamily > * ECFams);
void prepRecRedrawCoEvent(int GeneFamilyId, vector <CoEvent> * CoEventSet, bool verbose);

void prepECFRedrawCoEvent(int ECFidToReset, vector <CoEvent> * CoEventSet, bool verbose);

void refiningAndComputingECFamily(EquivalenceClassFamily * ECF, 
									ReconciledTree * Rtree1, ReconciledTree * Rtree2,
									double Again, double Bgain, 
									double absencePenalty, bool doAllPair, bool withTransfer, bool substractRecoToAdj, 
									double weightedDupCost, double weightedLossCost, double weightedHGTCost, 
									bool alwaysGainAtTop, double C1Advantage,bool superverbose,
                                    map<int,vector<float> > speciesC0C1,
                                    map<int, map<string,int> > speGeneAdjNb );


void refiningAndComputingECFamilyList(vector <EquivalenceClassFamily> * ECFams, vector < GeneFamily * > * GeneFamilyList,
										int GfamIdToReset,
										double Again, double Bgain, 
										double absencePenalty, bool doAllPair, bool withTransfer, bool substractRecoToAdj, 
										double weightedDupCost, double weightedLossCost, double weightedHGTCost, 
										bool alwaysGainAtTop, double C1Advantage, bool verbose, bool superverbose,
                                        map<int,vector<float> > speciesC0C1,
                                        map<int, map<string,int> > speGeneAdjNb 
                                        );


void GettingCoEvents(vector <CoEvent> * CoEventSet, vector <EquivalenceClassFamily> * ECFams, vector < GeneFamily * > * GeneFamilyList);

void refineCoEvents(vector <CoEvent> * CoEventSet);

vector <bool> checkingCoEventsAncestorSonConstraint( vector <CoEvent> * CoEventSet , vector < GeneFamily * > * GeneFamilyList);





void ReDrawGeneFamilyProper(int GfamIdToReset, int redrawAlgo, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams,
                        MySpeciesTree * speciesTree, bool datedSpeciesTree, bool boundedTS, bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        bool verbose, bool superverbose ,
                        bool trySeveralECFSolutions, double RecWeight,
                        double probaFlat,
                        double probaC1Bias , MyCladesAndTripartitions * biasedCCP,
                        string gPathPrefix, int nbTry = 0, int maxNbTry = 1,
                        bool noCoev=false , double DLRecCoevTemp=0.5,
                        string DTLRecCoevExecPath=""
                         );

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
                        double probaC1Bias = 0 , MyCladesAndTripartitions * biasedCCP = NULL,
                        string gPathPrefix = "", int nbTry = 0, int maxNbTry = 1,
                        bool noCoev =false, double DLRecCoevTemp=0.5,
                        string DTLRecCoevExecPath=""
                        );


void ReplaceGeneFamilyRecTree(int GfamIdToReset, ReconciledTree * newRTree, double treeLkh, 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                        bool withTransfer,
                        double DupliCost, double HGTCost, double LossCost, int maxTS, double TopoWeight,
                        double Again, double Abreak, double absencePenalty, bool doAllPair, bool substractRecoToAdj,
                        double weightedDupCost, double weightedLossCost, double weightedHGTCost, bool alwaysGainAtTop, double C1Advantage,
                        bool verbose, bool superverbose,
                        map<int,vector<float> > speciesC0C1,
                        map<int, map<string,int> > speGeneAdjNb 
                         );


/*void ReDrawECFamily(int ECFidToReset , 
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool alwaysGainAtTop, double C1Advantage,bool superverbose);*/

void ReplaceECFamily(int ECFidToReset , vector < vector< AdjTree *> * > * NewAdjacencyTreesVector, vector < CoEvent > * NewCoEventSet,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool superverbose);


bool acceptProposition(double scoreInit, double scoreProposed, double Temperature=1);


bool tryECFreDraw(int ECFidToReset , double Temperature, double oldCoEventScore,
                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * CoEventSet,
                    bool alwaysGainAtTop, double C1Advantage,bool superverbose,
                    double DupliCost, double LossCost, double HGTCost, double RecWeight);

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
                    double probaC1Bias = 0 , MyCladesAndTripartitions * biasedCCP = NULL,
                    string gPathPrefix = ""
                    );


//bool updateGfam(int GfamIdToReset, int maxNbTry, int redrawAlgo, bool trySeveralECFSolutions, double Temperature,
//                    vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * * CoEventSetPointer, MySpeciesTree * speciesTree,
//                    int maxTS, bool datedSpeciesTree, bool boundedTS, bool withTransfer, 
//                    double DupliCost, double HGTCost, double LossCost,
//                    double TopoWeight, double RecWeight, double AdjWeight,
//                    double Again, double Abreak,
//                    double absencePenalty, bool doAllPair, bool substractRecoToAdj, 
//                    bool alwaysGainAtTop, double C1Advantage,
//                    bool verbose, bool superverbose,
//                    map<int,vector<float> > speciesC0C1,
//                    map<int, map<string,int> > speGeneAdjNb,
//                    double probaC1Bias = 0 , MyCladesAndTripartitions * biasedCCP = NULL,
//                    string gPathPrefix = ""
//                    );




vector <AdjTree * > * readAdjForest(string filename, int gfam1, int gfam2, bool gainAtRoot, bool VERBOSE);

void LoadGeneFamilyReconciledTrees( vector < GeneFamily * > * GeneFamilyList, string prefix, MySpeciesTree * Stree, bool VERBOSE , bool OneFile = false);
void LoadECFamTrees(EquivalenceClassFamily * ECF, string prefix,bool gainAtRoot, bool VERBOSE);


bool readECFamTreesOneFile(string filename, vector< EquivalenceClassFamily > * ECFams, bool gainAtRoot, bool VERBOSE);

void getLinkedFamilies(int GeneFamilyId, vector < EquivalenceClassFamily > * ECFams, vector <int> &LinkedGFs , vector <int> &LinkedECFs, bool allowSelf =false);

int chooseOneLinkedfamily(int GeneFamilyId, vector < EquivalenceClassFamily > * ECFams,  bool allowSelf );

void AddECFToFile(string fileName, EquivalenceClassFamily * ECF, 
                             bool newick,bool hideLosses,  
                             double GainCost ,  double BreakCost );


GfamSave createHistoric(int GfamIdToReset , 
                        vector < GeneFamily * > * GeneFamilyList, vector < EquivalenceClassFamily > * ECFams, vector < CoEvent > * * CoEventSetPointer,
                        double Again, double Abreak,
                        double DupliCost, double HGTCost, double LossCost,
                        double RecWeight
                        );


vector < EquivalenceClassFamily > * BasicCreateEquivalenceClassFamilies(vector< pair <string,string > > adjacencies, 
                                                                    vector <GeneFamily *> * GeneFamilyList, 
                                                                    bool Verbose, bool SuperVerbose);

#endif