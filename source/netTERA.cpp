/**

@file netTERA.cpp
@author Celine Scornavacca
@author Edwin Jacox
@version 1.2 
@date 01/07/2015

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

**/

#include <boost/foreach.hpp>
#include <Bpp/Utils/AttributesTools.h>

#include "NetAlg.h"

string version = "1.0";


//name, type, default, description
std::map<string,bool> gBoolParams;
std::map<string,int> gIntParams;
std::map<string,double> gDoubleParams;
std::map<string,string> gStringParams; // string, char, and path 


// Ordered for the help message.
// Parameters with an empty description are not displayed in the help.
const int gParameterCount = 11;
string const gParameters[gParameterCount][4] = {
 {"species.file", "path", "required", "species tree file (newick)" },
 {"gene.file", "path", "required", "gene trees file (newick)"},
 {"verbose", "bool", "false", "show progress and timing"},
 {"force.binary", "bool", "false", "make gene trees binary if they are not"},
  {"transfer.cost", "double", "3", "cost of a transfer"}, 
 {"dupli.cost", "double", "2", "cost of a duplication"},
 {"loss.cost", "double", "1", "cost of a loss"},
 {"char.sep", "char", "_", "char separating gene names in gene tree files"},
 {"best.switch", "bool", "true", "run minimum switching algorithm"},
 {"min.recon", "bool", "false", "run minimum network reconciliation algorithm"},
 {"run.brute", "bool", "false", "run brute minimum switching algorithm"}
};

void help() {
    cout << "The program outputs the cost of the most parsimonious"
            " reconcilation between the species tree and the given"
            " gene networks." << endl;
    cout << "version: " << version << endl;
    cout << "__________________________________________________________"
            "_____________________" << endl;
    cout << "___________________________________INPUT__________________"
            "_____________________" << endl;
    cout << "format: parameter.name | [type,default] <description>" << endl;

    for( size_t i=0; i<gParameterCount; i++ ) {
        if( gParameters[i][3] != "" )
            cout << gParameters[i][0] << " | [" << gParameters[i][1] << ","
                 << gParameters[i][2] << "] " << gParameters[i][3] << endl;
    }

}


/**
 * read input parameters
 */
void readParameters( 
        map<string,string> params ) ///< input paramsters
{
    map<string,string>::iterator iterStr;
    bool problem = false;

    // check for unknown parameters
    map<string,bool> allParams;
    for( size_t i=0; i<gParameterCount; i++ ) 
        allParams[gParameters[i][0]] = true;
    for( iterStr = params.begin(); iterStr != params.end(); ++iterStr ) {
        map<string,bool>::iterator iterBool = allParams.find( iterStr->first );
        if( iterBool == allParams.end() ) {
            problem = true;
            cerr << "Unknown parameter: " << iterStr->first << endl;
        }
    }
        

    for( size_t i=0; i<gParameterCount; i++ ) {
        string name = gParameters[i][0];
        string type = gParameters[i][1];
        string value;

        // get value 
        iterStr = params.find( name );	
        if( iterStr == params.end() ) {
            if( gParameters[i][1] == "required" ) {
                problem = true;
                cerr << name << " must be specified" << endl;
            }
            value = gParameters[i][2]; // default
        } else {
            value = iterStr->second;
        }


        if( type == "bool" ) {
            if ((value == "true") 
                || (value == "TRUE")
                || (value == "t")
                || (value == "T")
                || (value == "yes")
                || (value == "YES")
                || (value == "y")
                || (value == "Y")
                || (value == "1") ) 
            {
                gBoolParams[name] = true;
            } 
            else if ((value == "false") 
                || (value == "FALSE")
                || (value == "f")
                || (value == "F")
                || (value == "no")
                || (value == "NO")
                || (value == "n")
                || (value == "N")
                || (value == "0") ) 
            {
                gBoolParams[name] = false;
            }
            else {
                cerr << "Invalid boolean value (" << value << ") for " 
                     << name << endl;
                problem = true;
            }
        } else if( type == "string" ) {
            gStringParams[name] = value;
        } else if( type == "char" ) {
            if( value.size() != 1 ) {
                cerr << "Invalid char value (" << value << ") for " 
                     << name << endl;
                problem = true;
            } else {
                gStringParams[name] = value;
            }

        } else if( type == "int" ) {
            gIntParams[name] = TextTools::toInt( value );
        } else if( type == "double" ) {
            gDoubleParams[name] = TextTools::toDouble( value );
        } else if( type == "path" ) {
            if( value == "required" ) {
                cerr << name << " is a required parameter." << endl;
                problem = true;
            } else if( value != "none" ) {
                if( !FileTools::fileExists( value ) ) {
                    cerr << value << " does not exist." << endl;
                    problem = true;
                } else {
                    gStringParams[name] = value;
                }
            } else 
                gStringParams[name] = "none";
        } else {
            cerr << "Found type " << type << endl;
            throw bpp::Exception( "readParameters: unknown type" );
        }
    }

    if( problem )
        exit(1);
}



/////////////////////////////////////////////////
// Main 
/////////////////////////////////////////////////

int main(int args, char ** argv)
{
	if(args == 1)
	{
		help();
		exit(0);
	}
	
	try {
        // fill global parameter variables
        map<string, string> params = AttributesTools::parseOptions(args, argv);
        readParameters( params );

        if( gBoolParams.find("verbose")->second ) {
	cout << "******************************************************************"
         << endl;
	cout << "*                 netTERA, version " << version << "          *"
        << endl;
	cout << "* Authors: C. Scornavacca, E. Jacox         Created     16/07/15 *"
        << endl;
	cout << "*                                           Last Modif. 30/03/16 *"
        << endl;
	cout << "******************************************************************"
        << endl;
	cout << endl;
        }


        if( gBoolParams.find("verbose")->second ) 
            ApplicationTools::startTimer();


        ////////////////////////////////////////
        // Get trees and check them and other input 
        ////////////////////////////////////////


        // read species tree
        string errString = "";
        MyNetwork* speciesNetwork = MyNetwork::readMyNetwork( 
                gStringParams.find("species.file")->second.c_str(),
                                     errString );
        if( errString != "" || speciesNetwork == NULL ) {
            cerr << "Error reading species tree: " << errString << endl;
            exit(1);
        }

        speciesNetwork->assignNetworkPostOrderIds();
        

        vector<MySpeciesNode*> allNodes = speciesNetwork->getNodes();
        
        string errStr;
        
        BOOST_FOREACH( MySpeciesNode *node, allNodes ) {
			int sonCount = node->getNumberOfSons();
			if( sonCount > 2 ) {
				errStr = "A network node has more than two children";
				return false;
			}
			for( int i=0; i<sonCount; i++ ) {
				MySpeciesNode *son = node->getSon( i );
				if(son->getInfos().primaryFather !=NULL){
					if( node->getId() != son->getInfos().primaryFather->getId() ) {
						son->getInfos().secondaryFather = node;
					}
				}
			}

    	}
    	
        // species tree must be binary, except for hybrid nodes (2 parents)
        if( !speciesNetwork->checkNetwork( errStr ) ) {
            cerr << "ERROR: Invalid species network: " << errStr << endl;
            exit(1);
        }

        if( !speciesNetwork->checkNetwork( errStr ) ) {
            cerr << "ERROR: second check failed: " << errStr << endl;
            exit(1);
        }
        
                

    

        // map of taxaNames required for the function restrictTreeToASetOfTaxa
        boost::unordered_map<string, int> taxaNamesSpecies;
        vector<string> leafNames = speciesNetwork->getLeavesNames();
        BOOST_FOREACH( string leafName, leafNames ) {
            size_t pos = leafName.find( "#" ); 
            leafName = leafName.substr(0,pos);
            taxaNamesSpecies.insert( make_pair(leafName,1) );
        }


		//reading the gene trees 
        vector<MyGeneTree*> geneTrees = MyGeneTree::readMyGeneTrees(
                        gStringParams.find("gene.file")->second.c_str(),
                        errString );
        if( errString != "" ) {
            cerr << "Error reading gene trees: " << errString << endl;
            exit(1);
        }
        if( geneTrees.size() < 1 ) {
            cerr << "No gene trees found" << endl;
            exit(1);
        }
        if( geneTrees.size() > 1 ) 
            cout << "More than one gene tree, using the first one." << endl;
        MyGeneTree *geneTree = geneTrees[0];




        //////////////////////////////////////////////////////////
		// Gene trees
        //////////////////////////////////////////////////////////
        
        if( !geneTree->restrictTreeToASetOfTaxa( 
            taxaNamesSpecies, gStringParams.find("char.sep")->second[0],
            gBoolParams.find("verbose")->second )  )
        {
            cerr << "ERROR: Gene tree has no valid leaves." << endl;
            exit(1);
        }

// Is there a limit?
        // needed by the algoritm!!!	
        if( geneTree->getNumberOfLeaves()<3 ) {
            cerr << "ERROR: Gene tree has only " 
                 << geneTree->getNumberOfLeaves()
                 << " leaves. At least three are required." << endl;
            exit(1);
        }

        if( !geneTree->isBinary() ) {
            if( gBoolParams.find("force.binary")->second )
                geneTree->makeBinary();
            else {
                cerr << "ERROR: Gene tree is not binary." << endl;

        		vector<MyGeneNode*> nodes = geneTree->getNodes();		
        		for( size_t i=0; i<nodes.size(); i++ ) {
            		if( !nodes[i]->isLeaf() && nodes[i]->getNumberOfSons() != 2 ) { 
                		cout << nodes[i]->getNumberOfSons() << endl;
           	 		}    
       			}

                exit(1);
            }
        }

        // gene trees must have unique leaves
        string dupName;
        if( !geneTree->uniqueLeaves( dupName ) ) {
            cerr << "ERROR: A gene tree has duplicate leaf names: <" 
                 << dupName << ">" << endl;
            exit(1);
        }

        NetAlg netAlg( speciesNetwork, geneTree,
                       gStringParams.find("char.sep")->second[0],
                       gDoubleParams.find("dupli.cost")->second,
                       gDoubleParams.find("loss.cost")->second,
                       gDoubleParams.find("transfer.cost")->second );
                       

        
        
           
        if( gBoolParams.find("run.brute")->second ) {
            cout << "======== BRUTE ALG====== " << endl;
            int numLosses = 0;
            int numDupli = 0;
            int numTransfers= 0;

            double cost = netAlg.runBrute( numLosses, numDupli, numTransfers);        
                                
            cout << numLosses << " losses, " << numDupli
                 << " duplications and " << numTransfers << " transfers" << endl;
            cout << "cost = " << cost << endl;
        }

        if( gBoolParams.find("best.switch")->second ) {
            cout << "===============SWITCH ALG============" << endl;

            int numLosses = 0;
            int numDupli = 0;
            int numTransfers= 0;
             vector<std::pair <int,int> >  edgesBestSwitchings;
            double costMinSwitch = netAlg.runMinSwitch( numLosses, numDupli , numTransfers, edgesBestSwitchings );
            cout << numLosses << " losses, " << numDupli
                 << " duplications and " << numTransfers << " transfers" << endl;
            cout << "cost = " << costMinSwitch << endl;
        }
        
        if( gBoolParams.find("min.recon")->second ) {
            cout << "=============== MIN RECON ALG============" << endl;
            double costMinRec = netAlg.runMinRecon();
            cout << "cost = " << costMinRec << endl;
        }

        BOOST_FOREACH( MyGeneTree *tree, geneTrees ) 
            delete tree;
        delete speciesNetwork;

        if( gBoolParams.find("verbose")->second )
		    ApplicationTools::displayTime("Done:");

	}  catch(exception & e) {
		cerr << e.what() << endl;
		exit(-1);
	}
	
	return 0;
}
