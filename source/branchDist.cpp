
//#define VERBOSE


//
// File: branchDist.cpp
// Created by: Edwin Jacox
// Created on: Dec 10 2013
//

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
Calculation of branch score distance (generalization of the Robinson-Foulds 
distance) between two trees. This measure is based on the clusters/biparitions 
(bips) of the trees. At it's simplest, the score is the sum of squares, 
over all bips, of the difference in branch length. If a bip exists in one tree,
but not the other, the branch length of the non-existent bip is zero. All 
branch lengths are set to 1 to derive the Robinson-Foulds distance.

The simplest algorithm finds all of the bips for both trees, then identifies
matching bips. The bottleneck is finding matching bips. 
The optimized algorithm (for one-to-one) is based on the Day (1985) 
algorithm. bips are hashed/mapped, providing efficient access. Furthemore, 
leaf names are converted to integers for more efficient hashing. Leaf 
name/number associations are stored in a map. The Day algorithm maps just 
one tree. The leaf numbers are based on the one tree and assigned such 
that all leafs in a bips are contiguous integers. In this way, a bip 
can be represented by the min and max leaf numbers.

For all-to-all (or networks), it could be more efficient to create bip 
representations once for each tree using a bitmap (a bit is set for 
each leaf present in bip). These are hashed for efficient access.


The input is two or more phylogenetic trees.
 */

// From the STL:
#include <iostream>
#include <sstream>
//#include <unordered_map>

using namespace std;

// From PhylLib:
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/BipartitionTools.h>


// From Utils:
#include <Bpp/Io/FileTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Utils/AttributesTools.h>

using namespace bpp;


// set to map or unordered_map (requires stl11)
#define MY_MAP map

void help()
{
  ApplicationTools::displayMessage("__________________________________________________________________________" );
  ApplicationTools::displayMessage("This program takes as input phylogenetic trees in newick format");
  ApplicationTools::displayMessage("(in one or two files) and outputs their RF distance and branch scores.");
  ApplicationTools::displayMessage("first.tree.file            | [path] toward the first tree (newick)     " );
  ApplicationTools::displayMessage("second.tree.file           | [path] toward the second tree (newick)     " );
  ApplicationTools::displayMessage("common.names           | [bool] Restricts the trees to the common set of labels of all trees" );
  ApplicationTools::displayMessage("one.to.many           | Computes the score for the first tree against all of the other input trees" );
  ApplicationTools::displayMessage("many.to.many           | Computes the score for every pair of input trees." );
  ApplicationTools::displayMessage("___________________________|___________________________________________" );
}




/** TO DO
 *  - MORE TESTS: 5 sizes, same and diff leaf names, tree with no distances?
 *  - excpetion handling
 *  - check for memory leaks
 *
 *  - convert to bio++ trees and read trees - NOT NECESSARY
 *  - (possible) add an option to return clusters
 *  - unwind traversal for speed? - NO
    - (possible) Create a class/function to do the hashing of the bip.
     - use same map if contiguous
     - else a different map 
        - CHECK CLUSTER AND CHOOSEMAP, use same numbering so one distance vector!
    - (possible) clean up 
       create classes and shorten functions: bipMap and nameMap become members
*/


/** TESTS
 * g1 g2 0.00619433
 * g1 g2a 0.00731645
 * g1 g2b 0.00508027
 * g1-g2b g2a 0.00571145, 0.00731645
 * gen50-1 gen50-2 score=98, RF=92 time-1000loops = 2 seconds
 * gen100 score=200, RF=194, time1000 = 3 seconds
 * gen300 score=600, RF=594, time1000 = 9 seconds
 * gen1000 score=2000, RF=1994, time1000 = 34 seconds
 *
 */




/*
 * @brief Traverse tree, adding or deleting names. A tree traversal is done
 * explicitly rather than a function call like getLeaves() to ensure the 
 * ordering matches other function in branch distance, which is important
 * for the Day algorithm implementation because bipartitions must have a contiguous
 * numbering.
 * @param node - root of the tree or subtree 
 * @param nameMap - mapping of leaf names to integers 
 * @param leafCnt - number of entries for assign leafNum
 */
void assignNumericLeafNames( const Node *node, 
                            MY_MAP<string,size_t> &nameMap, 
                             size_t &leafCnt ) {
   
    int numSons = node->getNumberOfSons();
    for(int i=0; i<numSons; i++)
        assignNumericLeafNames( node->getSon(i), nameMap, leafCnt );

    if( !numSons ) {
        // leaf 
        nameMap.insert( make_pair(node->getName(), leafCnt++) );
    } 
}

/*
 * @brief Set hasIds elements to zero if leaf id is not in the tree.
 * @param tree 
 * @param nameMap - mapping of leaf names to integers
 * @param hasId - a vector indicating common leaf ids
 */
void findCommonNames( TreeTemplate<Node> *tree, 
                      MY_MAP<string,size_t> &nameMap, 
                      vector<unsigned char> &hasId) {

    //create array, corresponding to hasId, of leaf names found
    vector<unsigned char> foundIt( hasId.size() );
    vector <Node *> leaves = tree->getLeaves();
	for (size_t i=0; i<leaves.size(); i++) {
        MY_MAP<string,size_t>::const_iterator pos = nameMap.find( leaves[i]->getName() );
        if (pos != nameMap.end()) {
            foundIt[pos->second] = 1;
        }
    }

    // mark as zero those ids not found in the tree
    for( size_t i=0; i<foundIt.size(); i++ ) 
        if( !foundIt[i] ) 
            hasId[i] = 0;
}


/*
 * @brief Create nameMap of with common names. First tree determines ordering. 
 * Leaves not in the other trees are removed the map. This creates
 * non-contigous leaf numbers, which requires that the numbers be redone.
 * @param tree 
 * @param nameMap - mapping of leaf names to integers
 */
void createCommonNames( const vector<TreeTemplate<Node> *>& trees, 
                        MY_MAP<string,size_t> &nameMap ) {

    size_t leafCnt = 0;
    assignNumericLeafNames( trees[0]->getRootNode(), nameMap, leafCnt );

    vector<unsigned char> hasIds(nameMap.size()+1, 1);
    for( size_t i=1; i<trees.size(); i++) {
        findCommonNames( trees[i], nameMap, hasIds ); 
    }

    // remove those not in array from nameMap
    MY_MAP<string,size_t>::iterator iter;
    for (iter = nameMap.begin(); iter!=nameMap.end(); ) {
        if( !hasIds[iter->second] ) {
            nameMap.erase(iter++);
        } else {
            ++iter;
        }
    }
}



struct ChildInfo {
    int minLeafNum;
    int maxLeafNum;
    int bipSize;
    int bipNum;
};

/*
 * @brief Function do to the traversal for mapBipartitionsDay.
 *
 * Read tree and create clusters/bipartiions (bip) for each node. 
 * Each bip is stored as minimum leaf number, maximum leaf number and size. 
 * The bip are given ids (bipNum) in bipMap, which is
 * also used to find matching bip from other trees. The id is used
 * as an index into the distances vector.
 *
 * The algorithm is recursive. Each call is at a node and creates a bip,
 * except for the root and root's second son.
 *
 *    Assign all bipartition numbers for index into distances, 
 *    but we don't need to insert into bipMap for the last tree.
 *
 * @param node - node of a tree 
 * @param bipMap - map of bipartitions (min/size) to a bipNum
 * @param nameMap - maps leaf names to numbers
 * @param createNameMap - indicates to generate name map
 * @param leafCnt - counter for generating nameMap
 * @param firstTree - indicates to reunumber name map (if not creating it),
 *                    and to add to bipMap
 * @param bipCnt - for indexing distances array
 * @param distances - a vector of distances, indexed by bipNum
 * @return struct of ChildInfo
 */
ChildInfo mapBipartitionsDayAux( const Node *node, 
                      MY_MAP<string,size_t> &bipMap, 
                      MY_MAP<string,size_t> &nameMap, 
                      const bool createNameMap,
                      size_t &leafCnt,
                      const bool firstTree,
                      size_t &bipCnt,
                      vector<double> &distances ) {

    int minLeafNum = -1;
    int maxLeafNum = -1;
    size_t bipSize = 0; // leaves beneath this node

    int numValidSons = 0; // number of sons from nameMap

    int firstValidChildInfo_bipNum = -1;
    ChildInfo secondValidChildInfo;
    secondValidChildInfo.bipSize = -1; 
    if( node->isLeaf() ) {
        // leaf 
        if( createNameMap ) {
            // give leaf a new number id
            MY_MAP<string,size_t>::iterator pos = nameMap.find(node->getName());
            if (pos != nameMap.end()) {
                minLeafNum = maxLeafNum = leafCnt++;
                nameMap.insert( make_pair(node->getName(), minLeafNum) );
                bipSize = 1;
#ifdef VERBOSE
                cout << "nameMap " << node->getName() << " -> " 
                     << minLeafNum << endl;
#endif
            } 
            // else   not part of common set of leaves

        } else {
            // find leaf number id
            MY_MAP<string,size_t>::iterator pos = nameMap.find(node->getName());
            if (pos != nameMap.end()) {
                bipSize = 1;
                if( firstTree ) {
                    minLeafNum = maxLeafNum = leafCnt++;
                    pos->second = minLeafNum;
#ifdef VERBOSE
                    cout << "nameMap " << node->getName() << " -> " 
                        << minLeafNum << endl;
#endif
                } else {
                    minLeafNum = maxLeafNum = pos->second;
                }
            } 
            // else   not part of common set of leaves
        }
    } else {
        for( size_t i=0; i<node->getNumberOfSons(); i++) {
            // recurse
            ChildInfo childInfo = mapBipartitionsDayAux( 
                    node->getSon(i), bipMap, nameMap, createNameMap, 
                    leafCnt, firstTree, bipCnt, distances ); 

            // accumulate child info into current bipartition 
            if( childInfo.bipSize != 0 ) {
                numValidSons++;

                if( numValidSons==1 || childInfo.minLeafNum < minLeafNum ) 
                    minLeafNum = childInfo.minLeafNum;
                if( numValidSons==1 || childInfo.maxLeafNum > maxLeafNum ) 
                    maxLeafNum = childInfo.maxLeafNum;
                bipSize += childInfo.bipSize;

                // save info for collapsing and merging rooted tree
                if( numValidSons == 1 ) {
                    firstValidChildInfo_bipNum = childInfo.bipNum;
                } else if( numValidSons == 2 ) {
                    secondValidChildInfo = childInfo;
                }
            }
        }
    }


    // get distance for this node/bipartition if it is not the root
    int bipNum = -1;  // root unassigned
    if( bipSize == nameMap.size() ) {
        // this is the real root 
        if( numValidSons == 2 
            && firstValidChildInfo_bipNum != secondValidChildInfo.bipNum ) 
        {
            // rooted so merge distances
            distances[firstValidChildInfo_bipNum] += 
                distances[secondValidChildInfo.bipNum];

            // leave second child unused
            distances[secondValidChildInfo.bipNum] = 0;

            // recreate IdStr to remap bipNum
            stringstream secondChildBipIdStr;
            secondChildBipIdStr << secondValidChildInfo.minLeafNum 
                << "-" << secondValidChildInfo.maxLeafNum 
                << "-" << secondValidChildInfo.bipSize;
            bipMap[secondChildBipIdStr.str()] = firstValidChildInfo_bipNum;
        }
    } else if( bipSize > 0 && node->hasFather() ) {

// Are there trees with no distances for RF? If so, set distance to 1 if no
// distance. CASES! TEST THEM! Or always set to 1 if no distance.
        double dist = 1;
        if( node->hasDistanceToFather() ) 
            dist = node->getDistanceToFather();

        stringstream bipIdStr;
        bipIdStr << minLeafNum << "-" << maxLeafNum << "-" << bipSize;


        // set info
        if( numValidSons == 1 ) {
            // This bipartition exists because there was only one child.
            // Collapse the branch by adding the distances.
            bipNum = firstValidChildInfo_bipNum;
            distances[bipNum] += dist;
        } else {
            // check if bip has been assigned already
            MY_MAP<string,size_t>::iterator iter;
            if( !firstTree )  
                iter = bipMap.find( bipIdStr.str() );

            if( firstTree || iter == bipMap.end( )) {
                // bip not in map, add it 
                bipNum = bipCnt++;
                distances.push_back(dist);
#ifdef VERBOSE
                cout << "NEW CLUSTER " << bipNum << endl;
#endif
                if( firstTree ) 
                    bipMap.insert( make_pair(bipIdStr.str(), bipNum) );
                    // else second tree and don't care because there will be 
                    // no look up of this value
            } else {
                // modify existing
                bipNum = iter->second;
                // normally, this shouldn't be set already, but
                // can occur if this is a root son that has
                // been adjusted before
                distances[bipNum] += dist; 
            }
        }



#ifdef VERBOSE
        string name = "internal";
        if( node->hasName() ) {
            name = node->getName();
        }
        cout << "node " << name << "(" << numValidSons << ")"
            << " bip " << bipIdStr.str() << " -> " << bipNum 
            << " = " << distances[bipNum] << endl;
#endif


    } 


    ChildInfo info;
    info.minLeafNum = minLeafNum;
    info.maxLeafNum = maxLeafNum;
    info.bipSize = bipSize;
    info.bipNum = bipNum;

    return info;
}


/*
 * @brief Prepare each tree for a call to mapBipartitionsDayAux, which
 * does the work to get the distances for each cluster.
 *
 * @param trees
 * @param firstTreeIdx - use only trees with or after idx
 * @param nameMap - maps leaf names to numbers
 * @param createNameMap
 * @param allDistances - distance vectors
 * @return number of bipartitions in all trees (double counts root sons)
 */
size_t mapBipartitionsDay( const vector<TreeTemplate<Node> *>& trees, 
                        size_t firstTreeIdx,
                         MY_MAP<string,size_t> &nameMap,
                         bool createNameMap,
                         vector< vector<double> > &allDistances )
{

    // associates bipartition with numbers
    MY_MAP<string,size_t> bipMap;

    size_t bipCnt = 0; // number of bipartitions found so far
    bool firstTree = true; // renumber in first pass to get contiguous
    for( size_t i=firstTreeIdx; i<trees.size(); i++ ) {

#ifdef VERBOSE
        cout << "======== tree " << i << " ================" << endl;
#endif
        Node *root = trees[i]->getRootNode();


        vector<double> distances(bipCnt);
        size_t leafCnt = 0;
        mapBipartitionsDayAux( root, bipMap, nameMap,
                                createNameMap, leafCnt, firstTree, 
                                bipCnt, distances );
        allDistances.push_back( distances );

        createNameMap = false; // only need first time, if at all
        firstTree = false;
    }

    return bipCnt;

}

struct Scores {
    double long branchScore;
    double long normalizedScore;
    size_t RFscore;
};


/*
 * @brief Score all trees against first tree.
 *
 * @param trees
 * @params nameMap
 * @param firstTreeIdx - use only trees with or after idx
 * @param createNameMap - Ignore labels not shared by all trees.
 * @return vector of structures with various scores
 */
vector<Scores> branchDist( const vector<TreeTemplate<Node> *>& trees, 
                    size_t firstTreeIdx,
                    MY_MAP<string,size_t> nameMap,
                    bool createNameMap ) {


    vector< vector<double> > allDistances; 
    size_t bipCount= mapBipartitionsDay( trees, firstTreeIdx,
                                        nameMap, createNameMap, allDistances );

    // Calculate distances
    vector<Scores> allScores;
    vector<double> firstDists = allDistances.front();
    vector< vector<double> >::iterator secondDistIter = allDistances.begin();
    secondDistIter++; // skip to second
    for (; secondDistIter != allDistances.end(); secondDistIter++) {

        Scores scores;
        scores.branchScore = 0;
        scores.RFscore = 0;
#ifdef VERBOSE
        cout << "SCORE:" << endl;
#endif
        for( size_t i=0; i<bipCount; i++ ) {
            double long sum = 0;
            int RFsum = 0;
            double dist1 = 0;
            if( i<firstDists.size() )
                dist1 = firstDists[i];
            sum += dist1;
            if( dist1 ) 
                RFsum += 1;

            double dist2 = 0;
            if( i<secondDistIter->size() )
                dist2 = (*secondDistIter)[i];
            sum -= dist2;
            if( dist2 ) 
                RFsum -= 1;

#ifdef VERBOSE
            cout << i << " " << dist1 << " " << dist2 << " = " << sum << endl;
#endif
            scores.branchScore += pow((long double) sum, 2 );
            scores.RFscore += pow( RFsum, (long double) 2 );
        }
        scores.normalizedScore = scores.branchScore /(2*nameMap.size()-3);

        allScores.push_back( scores );
    }


    return allScores;






}



/* this function reads a list of trees written in a file (path) in a newich format, separed by semicolons and returns a list of TreeTemplate<Node>*/

vector < TreeTemplate<Node> *>  readTrees(const string & path) throw (bpp::Exception) {
    // Checking the existence of specified file
    
    ifstream file(path.c_str(), ios::in);
    if (! file) { throw bpp::IOException ("\nError reading file.\n"); }
    
    vector < TreeTemplate<Node> *> trees;
    string temp, description;// Initialization
    // Main loop : for all file lines
    
    while (! file.eof()) {
        temp=FileTools::getNextLine(file);
        if(temp.size()!=0){
            string::size_type index = temp.find(";");
            if(index== string::npos) 
                throw bpp::Exception("readTrees(). Bad format: no semi-colon found.");
            if(index >= 0 && index < temp.size()) {
                description += temp.substr(0, index + 1);	
                TreeTemplate<Node> * tree = 
                    TreeTemplateTools::parenthesisToTree(description,true);    
                trees.push_back(tree);
                description = temp.substr(index);	
            } 
            else description += temp;
        }
    }
    file.close();
    return trees;	
};




int main(int args, char ** argv){

	#ifdef VERBOSE 
	cout << "******************************************************************" << endl;
  	cout << "*               Bio++ branch score distance, version 0.2.0 *" << endl;
  	cout << "* Author: C. Scornavacca                    Created     21/10/13 *" << endl;
  	cout << "*    and  E. Jacox                          Last Modif. 17/12/13 *" << endl;
  	cout << "******************************************************************" << endl;
  	cout << endl;
	#endif 
	
	if(args == 1)
  	{
    	help();
    	exit(0);
  	}
  
	try {
	  
		

		#ifdef VERBOSE
		cout << "Parsing options:" << endl;
		#endif 
		
		// Get the parameters from command line:
		map<string, string> cmdParams = 
            AttributesTools::getAttributesMap(
                    AttributesTools::getVector(args, argv), "=");
		
		// Look for a specified file with parameters:
		map<string, string> params;
		if(cmdParams.find("param") != cmdParams.end())
		{
			string file = cmdParams["param"];
			if(!FileTools::fileExists(file))
			{
				cerr << "Parameter file not found." << endl;
			    exit(-1);
			}
			else
			{
				params = AttributesTools::getAttributesMapFromFile(file, "=");
		  	    // Actualize attributes with ones passed to command line:
		  	  AttributesTools::actualizeAttributesMap(params, cmdParams);
			}
		}
		else
		{
			params = cmdParams;
		}
		
		    
		string listPath1 = ApplicationTools::getAFilePath("first.tree.file", params);
		#ifdef VERBOSE
		ApplicationTools::displayResult("Input first file", listPath1);
		#endif 
		if(listPath1 == "none") 
            throw bpp::Exception("You must provide the first tree file.");
	
	  	vector<TreeTemplate<Node> *> trees = readTrees(listPath1);
    
		if(cmdParams.find("second.tree.file") != cmdParams.end())  {
		    string listPath2 = ApplicationTools::getAFilePath(
                    "second.tree.file", params);
		    #ifdef VERBOSE
		    ApplicationTools::displayResult("Input second file", listPath2);
		    #endif 

// Bpp tree reader OK now. Switch to that.
// How should I handle multiple tree inputs?
            vector<TreeTemplate<Node> *> treesTemp = readTrees(listPath2);

            for(unsigned int y=0;y< treesTemp.size();y++){
                trees.push_back(treesTemp[y]);
            }
        }

// options: 1n,nn, common, error if more than 2 and not first two
		if(trees.size()<2) {
            std::stringstream ss;
			ss <<  "This program expects at least two trees as input, but you gave " 
               << trees.size() << " trees";
            throw bpp::Exception( ss.str() );
        } else if(trees.size()>2
		    && cmdParams.find("one.to.many") == cmdParams.end()
		    && cmdParams.find("many.to.many") == cmdParams.end())  {
            std::stringstream ss;
			ss <<  "This program expects exactly two trees as input unless"
                << " the one.to.one or one.to.many options is specified."
                " You gave " 
               << trees.size() << " trees";
            throw bpp::Exception( ss.str() );
        }
    
        // find common set of taxa if true
        bool checkCommonNames = ApplicationTools::getBooleanParameter(
                            "common.names", cmdParams, true, "", true, false );
        MY_MAP<string,size_t> nameMap; // associates leaf names with numbers
        /* Find common names if necessary */
        if( checkCommonNames ) {
            createCommonNames( trees, nameMap );
        } // else let branchDist create nameMap

        size_t lastIdx = 0;
		if( cmdParams.find("many.to.many") != cmdParams.end()) 
            lastIdx = trees.size()-1;

        cout << "firstTree-secondTree: Robinson-Foulds, branch distance, "
             << "normalized branch distance" << endl;
        for( size_t firstTreeIdx = 0; firstTreeIdx<=lastIdx; firstTreeIdx++ ) {
                                                    
            vector<Scores> scores = branchDist( trees, firstTreeIdx, 
                                                nameMap, !checkCommonNames );
            vector<Scores>::iterator iter;
            size_t cnt = firstTreeIdx+2;
            for ( iter = scores.begin(); iter!=scores.end(); ++iter ) {
                cout << firstTreeIdx+1 << "-" << cnt 
                    << ": " << iter->RFscore 
                    << ", " << iter->branchScore
                    << ", " <<  iter->normalizedScore << endl;
                cnt++;
            }

            // remove first element
        }
	
// Where should this clean up go? Needs to be handled during exceptions.
        for(unsigned int i = 0; i < trees.size(); i++) delete trees[i];

	} catch (std::exception&  e) {
    	cout << e.what() << endl;
    	exit(-1);
	}


	return (0);
};	
	
		
		

