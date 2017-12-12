/*

@file
@author Celine Scornavacca
@author Edwin Jacox

@section LICENCE
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

@section DESCRIPTION
*/

#include <ctime>
#include <boost/foreach.hpp>

#include "IDTLMatrix.h"


//#define DEBUG

struct time_sort {
inline bool operator() (const MySpeciesNode * a, const MySpeciesNode * b ) {
    return b->getInfos().timeSlice > a->getInfos().timeSlice;
}
};


double IDTLMatrix::computeTransferCost(
    int idUsource,  ///< gene source clade
    int idUtarget,  ///< gene target clade
    int timeSlice,  ///< current time slice
    int idX,         ///< speciesClade
    int &targetSpecies )
{
    targetSpecies = mBestReceiver[timeSlice][idUtarget];  
    if( targetSpecies == idX ) 
        targetSpecies = mSecondBestReceiver[timeSlice][idUtarget]; 
    
    if( targetSpecies == -1 ) // no best receiver 
        return std::numeric_limits<double>::max();

    double cost = mHGTCost
                + mMatrix.getValueSure( idUsource, idX) 
                + mMatrix.getValueSure( idUtarget, targetSpecies);
//cout << "xfer source: " << idUsource << "," << idX << "=" << mMatrix.getValueSure( idUsource, idX ) << endl;
//cout << "xfer target: " << idUtarget << "," << targetSpecies << "=" << mMatrix.getValueSure( idUtarget, targetSpecies) << endl;

    return cost;
}

void IDTLMatrix::setBacktrack(
        int u1,
        int x1,
        int u2,
        int x2 ) 
{
    mU1 = u1;
    mX1 = x1;
    mU2 = u2;
    mX2 = x2;
}

/**
 * Compute the optima for the matrix cell specified by
 * gene and species clade ids, which is the optimal costs 
 * without considering the transfer loss cost.
 *
 * If best splits are needed, the split and event associated
 * with the best cost is saved.
 *
 * @return optimal cost
 */
string IDTLMatrix::computeOptimaForCell(
        int idU,    ///< gene clade
        int timeSlice,  ///< current time slice
        int idX,        ///< species clade
        bool backtrack )
{
    double optCost = std::numeric_limits<double>::max();
#ifdef DEBUG
cout << "  ===computeOptimaForCell idX=" << idX << endl;
#endif

    string event = "U";

    if( timeSlice == 0 
        && mSpeciesTree->getNodeById( idX ) != NULL
        && mCladesTrips->mClades.isLeaf( idU ) 
        && mSpeciesTree->getNodeById( idX )->isLeaf()
        && mCladesTrips->mClades.getSpeciesName(idU) 
            == mSpeciesTree->getNodeById(idX)->getName() ) 
    {
#ifdef DEBUG
cout << "  LEAF: " << mSpeciesTree->getNodeById(idX)->getName() << " idX=" << idX << endl;
#endif
        mMatrix.setValue( idU, idX, 0 );
        if( backtrack ) setBacktrack( -1, -1, -1, -1 );
//cout << idU << "," << idX << "," << timeSlice << ":" << 0 << endl;
        return "C";
    }

    // Normal speciation + lost
    for( size_t i=0; i<mSpeciesTree->getSplits(idX).size(); i++ ) {
        pair<int,int> spSplit = mSpeciesTree->getSplits(idX)[i];
#ifdef DEBUG
cout << "    species " << spSplit.first << "/" << spSplit.second << endl;
#endif
        if( spSplit.second == -1 ) {
            double nullCost = mMatrix.getValueSure( idU, spSplit.first );
            if( optCost > nullCost ) {
                optCost = nullCost;
                event = 'N';
                if( backtrack ) setBacktrack( idU, spSplit.first, -1, -1 );
            }
#ifdef DEBUG
cout << "   NULL: " << nullCost << endl;
#endif
        } else {
            double spLossCost = mMatrix.getValueSure( idU, spSplit.first )
                          + mLossCost;
#ifdef DEBUG
cout << "   spLossCost1 " << spLossCost << " from " << spSplit.second << endl;
#endif
            if( optCost > spLossCost ) {
                optCost = spLossCost;
                event = "SL";
                if( backtrack ) setBacktrack( idU, spSplit.first, -1, -1 );
            }
            spLossCost = mMatrix.getValueSure( idU, spSplit.second ) 
                       + mLossCost;
#ifdef DEBUG
cout << "   spLossCost2 " << spLossCost << " from " << spSplit.first  << endl;
#endif
            if( optCost > spLossCost ) {
                optCost = spLossCost;
                event = "SL";
                if( backtrack ) setBacktrack( idU, spSplit.second, -1, -1 );
            }
        }
    }

    int splitCount = mCladesTrips->getSplitCount( idU );
    for( int splitIdx=0; splitIdx<splitCount; splitIdx++ ) {
        pair<int,int> geneSplit = mCladesTrips->getCladeSplit( idU, splitIdx );
#ifdef DEBUG
cout << "    gene " << geneSplit.first << "/" << geneSplit.second << endl;
#endif
        if( geneSplit.first == -1 )
            continue; // leaf

        // duplication
        double dupliCost = mDupliCost
                         + mMatrix.getValueSure( geneSplit.first, idX )
                         + mMatrix.getValueSure( geneSplit.second, idX );
#ifdef DEBUG
cout << "   DUPLI: " << dupliCost << " split:" << geneSplit.first << "/" << geneSplit.second << endl;
#endif
        if( optCost > dupliCost ) {
            optCost = dupliCost;
            event = "D";
            if( backtrack ) 
                setBacktrack( geneSplit.first, idX, geneSplit.second, idX );
        }

        // transfer
        if( mComputeT ) {
            int targetSpecies = -1;
            double transCost = computeTransferCost( geneSplit.first,
                              geneSplit.second, timeSlice, idX, targetSpecies );
            if( optCost > transCost ) {
                optCost = transCost;
                event = "T";
                if( backtrack ) setBacktrack( geneSplit.first, idX, 
                                        geneSplit.second, targetSpecies);
            }
#ifdef DEBUG
cout << "   TRANS1: " << transCost << endl;
#endif
            transCost = computeTransferCost( geneSplit.second,
                              geneSplit.first, timeSlice, idX, targetSpecies );
            if( optCost > transCost ) {
                optCost = transCost;
                event = "T";
                if( backtrack ) setBacktrack( geneSplit.second, idX, 
                                        geneSplit.first, targetSpecies);
            }
#ifdef DEBUG
cout << "   TRANS2: " << transCost << endl;
#endif
        }

        // speciation
        for( size_t i=0; i<mSpeciesTree->getSplits(idX).size(); i++ ) {
            pair<int,int> spSplit = mSpeciesTree->getSplits(idX)[i];
            if( spSplit.second == -1 )
                continue;
            double spCost;
            spCost = mMatrix.getValueSure( geneSplit.first, spSplit.first )
                   + mMatrix.getValueSure( geneSplit.second, spSplit.second );
#ifdef DEBUG
cout << "   SP1: " << spCost << endl;
#endif
            if( optCost > spCost ) {
                optCost = spCost;
                event = "S";
                if( backtrack ) setBacktrack( geneSplit.first, spSplit.first,
                                             geneSplit.second, spSplit.second );
            }
            spCost = mMatrix.getValueSure( geneSplit.first, spSplit.second )
                   + mMatrix.getValueSure( geneSplit.second, spSplit.first );
#ifdef DEBUG
cout << "   SP2: " << spCost << endl;
#endif
            if( optCost > spCost ) {
                optCost = spCost;
                event = "S";
                if( backtrack ) setBacktrack( geneSplit.first, spSplit.second,
                                             geneSplit.second, spSplit.first );
            }
        }
    
        // ILS 
        for( size_t i=0; i<mSpeciesTree->getIlsSplits(idX).size(); i++ ) {
            pair<int,int> ilsSplit = mSpeciesTree->getIlsSplits(idX)[i];
            if( ilsSplit.second == -1 )
// Is this possible
                continue;
            double ilsCost;
            ilsCost = mMatrix.getValueSure( geneSplit.first, ilsSplit.first )
                    + mMatrix.getValueSure( geneSplit.second, ilsSplit.second )
                    + mILScost;
#ifdef DEBUG
cout << "   ILS1: " << ilsCost << endl;
#endif
            if( optCost > ilsCost ) {
                optCost = ilsCost;
                event = "I";
                if( backtrack ) setBacktrack( geneSplit.first, ilsSplit.first,
                                             geneSplit.second, ilsSplit.second);
            }
            ilsCost = mMatrix.getValueSure( geneSplit.first, ilsSplit.second )
                    + mMatrix.getValueSure( geneSplit.second, ilsSplit.first )
                    + mILScost;
#ifdef DEBUG
cout << "   ILS2: " << ilsCost << endl;
#endif
            if( optCost > ilsCost ) {
                optCost = ilsCost;
                event = "I";
                if( backtrack ) setBacktrack( geneSplit.first, ilsSplit.second,
                                             geneSplit.second, ilsSplit.first );
            }
        }
    }

    // ILS + lost
    for( size_t i=0; i<mSpeciesTree->getIlsSplits(idX).size(); i++ ) {
        pair<int,int> ilsSplit = mSpeciesTree->getIlsSplits(idX)[i];
#ifdef DEBUG
cout << "    ils " << ilsSplit.first << "/" << ilsSplit.second << endl;
#endif
        if( ilsSplit.second == -1 ) {
// Is this possible?
            double nullCost = mMatrix.getValueSure( idU, ilsSplit.first );
#ifdef DEBUG
cout << "   NULL22222222222222222222222222: " << nullCost << endl;
#endif
            if( optCost > nullCost ) {
                optCost = nullCost;
                event = "N";
                if( backtrack ) setBacktrack( idU, ilsSplit.first, -1, -1 );
            }
        } else {
            double ilsLossCost = mMatrix.getValueSure( idU, ilsSplit.first )
                          + mLossCost + mILScost;
#ifdef DEBUG
cout << "   ilsLossCost1 " << ilsLossCost << endl;
cout << "           " << mMatrix.getValueSure( idU, ilsSplit.first ) << endl;
#endif
            if( optCost > ilsLossCost ) {
                optCost = ilsLossCost;
                event = "IL";
                if( backtrack ) setBacktrack( idU, ilsSplit.first, -1, -1 );
            }
            ilsLossCost = mMatrix.getValueSure( idU, ilsSplit.second )
                          + mLossCost + mILScost;
#ifdef DEBUG
cout << "   ilsLossCost2 " << ilsLossCost << endl;
cout << "           " << mMatrix.getValueSure( idU, ilsSplit.second ) << endl;
#endif
            if( optCost > ilsLossCost ) {
                optCost = ilsLossCost;
                event = "IL";
                if( backtrack ) setBacktrack( idU, ilsSplit.second, -1, -1 );
            }
        }
    }

#ifdef DEBUG
cout << "  computeOptimaForCell idX=" << idX << "  optCost=" << optCost << " event=" << event << endl;
#endif
    
    mMatrix.setValue( idU, idX, optCost );
//cout << idU << "," << idX << "," << timeSlice << ":" << optCost << endl;

    return event;
}







/**
 * Compute and save bestReceiver and second best recievers.
 *
 */
void IDTLMatrix::computeBestReceivers( 
    int timeSlice,      ///< time slice to save
    int idU,            ///< gene id
    vector<int> &speciesInTS) ///< nodes to considier
{

    int bestReceiver = -1;
    int secondBestReceiver = -1;
    double bestCost = std::numeric_limits<double>::max();
    double secondBestCost = std::numeric_limits<double>::max();

    // Save other bests for graph construction.
    // if best==secondBest, save all nodes equal to best cost
    // if best!=secondBest, save all nodes equal to second best cost
    bool first = true;
    for( size_t i=0; i<speciesInTS.size(); i++ ) {
        int idX = speciesInTS[i];

        if( mSpeciesTree->getNodeById(idX) == NULL ) 
            continue;

        double spNodeOpt = mMatrix.getValueSure( idU, idX );
        if( first ) {
            first = false;
            bestReceiver = idX;
            bestCost = spNodeOpt;  
        } else if( COST_GREATER( bestCost, spNodeOpt ) ) {

            // move best to second best
            if( timeSlice > 0 &&
                !COST_EQUAL( bestCost, numeric_limits<double>::max() ) )  
            {
                // put old best as second best id and save old 
                // if not timeSlice=0 || not MAX
                secondBestReceiver = bestReceiver;
                secondBestCost = bestCost;
            }

            // new best
            bestReceiver = idX;
            bestCost = spNodeOpt;  

        } else if( timeSlice == 0 
                && COST_EQUAL( spNodeOpt, numeric_limits<double>::max() ) )
        {
            //  don't add max cost to second best in time slice 0 
        } else if( COST_GREATER( secondBestCost, spNodeOpt ) ) {
            // save best receiver for best recevier
            secondBestReceiver = idX;
            secondBestCost = spNodeOpt;
        } 

    }

    if( first )
        throw bpp::Exception ("DTLMatrix::computeBestReceivers:"
                          "best receiver not set" );

    if( bestReceiver != -1 ) {
        mBestReceiver[timeSlice][idU] = bestReceiver; 
        mBestReceiverCost[timeSlice][idU] = bestCost;
    } 

    if( secondBestReceiver != -1 ) {
        mSecondBestReceiver[timeSlice][idU] = secondBestReceiver; 
        mSecondBestReceiverCost[timeSlice][idU] = secondBestCost;
    }
}


IDTLMatrix::IDTLMatrix( MySpeciesTree *speciesTree, CladesAndTripartitions *cat,
               bool computeT, bool computeTL,
               double dupliCost, double hgtCost, double lossCost, 
               double ilsCost ) 
{
    mSpeciesTree = speciesTree;
    mCladesTrips = cat;
    mComputeT = computeT;
    mComputeTL = computeTL;
    mDupliCost = dupliCost;
    mHGTCost = hgtCost;
    mLossCost = lossCost;
    mILScost = ilsCost;

    mMaxTS = mSpeciesTree->getRootNode()->getInfos().timeSlice;

//    clock_t start = clock();
	mSTnodes = mSpeciesTree->getNumberOfIds();
//    double duration = clock()-start;
//cout << "orig ticks: " << duration << endl;

    init();
}

// Basically the same as DTLMatrix without:
//  subopt
//  iterations/updateCosts
//  recalc (tsStart,tsEnd)
void IDTLMatrix::calculateMatrix( 
    bool verbose,   ///< print extra info
    bool unDated, ///< species tree subdivided
    bool partialDates ) ///< partially dated species tree (if not fully dated)
{
    vector< vector<int> > cladesBySize = mCladesTrips->getCladesBySize();
    BOOST_FOREACH( vector<int> &sameSizeClades, cladesBySize ) {
        for( size_t p=0; p<sameSizeClades.size(); p++ ) {
            int idU = sameSizeClades[p];
            if( unDated || partialDates )
                //calculateMatrixNoSub( idU );
                throw bpp::Exception( "IDTLMatrix::calculateMatrix:"
                        " undated not implemented" );
            else
                // loop over all time slices
                for( int ts=0; ts<=mMaxTS; ts++ ) 
                    calculateMatrixTS( idU, ts );
        }
    }

}
double IDTLMatrix::computeTL(
    int idU,    ///< clade id
    int timeSlice, ///< current time slice
    int idX )    ///< species id
{
    double bestRecCost = mBestReceiverCost[timeSlice][idU]; 
    int bestReceiverId = mBestReceiver[timeSlice][idU];
    if( bestReceiverId == idX ) {
        //if the bestReceiver is x, we take the second one ..
        bestRecCost = mSecondBestReceiverCost[timeSlice][idU]; 
        bestReceiverId = mSecondBestReceiver[timeSlice][idU];
    }

    // Transfer from optimal node. Add loss cost if this is not alpha.
    if( bestReceiverId == -1 ) 
        return -1;
    double cost = bestRecCost + mHGTCost + mLossCost;		
#ifdef DEBUG
cout << "   tlCost " << cost << " from " << idX << "<-" << bestReceiverId << endl;
#endif
    return cost; 
}

// Differs from DTLMatrix:
//  subopts
//  triplets
//  speciesNodesTS are from clades
//  No alpha
void IDTLMatrix::calculateMatrixTS( 
    int idU,    ///< clade id
    int timeSlice ) ///< current time slice
{
#ifdef DEBUG
cout << "calculatematrixTS idU=" << idU << " ts=" << timeSlice << endl;
#endif

    // sorting by id ensures that children are seen before parent clades
    vector<int> tsClades = mSpeciesTree->getVectorWithTS( timeSlice );
//    sort( tsClades.begin(), tsClades.end() );

    BOOST_FOREACH( int idX, tsClades ) 
        computeOptimaForCell( idU, timeSlice, idX );

    // calculate and set best receivers (first and second least cost nodes)
    if( tsClades.size() > 1 ) 
        computeBestReceivers( timeSlice, idU, tsClades );

    // Calculate transfer loss costs and set matrix cell.
    if( mComputeTL ) {
        for( size_t i=0; i<tsClades.size(); i++ ) {
            int idX = tsClades[i];
            double opt = mMatrix.getValueSure( idU, idX );

            // check for better opt by transfer loss
            if( opt != 0 && tsClades.size() > 1 ) {
                double cost = computeTL( idU, timeSlice, idX );
                if( cost != -1 && cost < opt ) 
                    mMatrix.setValue( idU, idX, cost ); 
            }
        }
    }
}


void IDTLMatrix::getReconciliation( int u, int x ) {


    if( u==-1 ) {
        cout << "=======RECONCILIATION================================" << endl;
        getBestCost( u, x );
    }
cout << "=== " << u << "," << x << "  " << mMatrix.getValueSure( u, x ) << endl;
    int ts = mSpeciesTree->getNodeById(x)->getInfos().timeSlice;
    double tlCost = computeTL( u, ts, x );
    string event = computeOptimaForCell( u, ts, x, true );
cout << "timeSlice=" << ts << endl;
cout << "tlCost=" << tlCost << endl;
cout << "event=" << event << "  " << mU1 << "," << mX1 << " and " << mU2 << "," << mX2 << endl;
    int u2 = mU2;
    int x2 = mX2;
    string eventStr = bpp::TextTools::toString(x)  + "," + event + "," + bpp::TextTools::toString(mX1) + "," + bpp::TextTools::toString(mX2) + "@1";
//print reconciliation and species tree (using what ids?)
    if( mU1 != -1 ) 
        getReconciliation( mU1, mX1 );
    if( u2 != -1 ) 
        getReconciliation( u2, x2 );

// 1. **Use compute split - need 4 numbers (2 matrix values) (or a code?)
// 2. Duplicate logic
}
