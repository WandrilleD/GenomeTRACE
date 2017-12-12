// Created by: Celine Scornavacca


/**

@file ilsTERA.cpp
@author Celine Scornavacca
@author Edwin Jacox
@version 0.1
@date 22/10/2015

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

#include <sstream>
#include <sys/stat.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Utils/AttributesTools.h>

#include "CladesAndTripartitions.h"
#include "IDTLMatrix.h"
#include <Bpp/Phyl/Io/Newick.h>
			
string version = "0.1";

//name, type, default, description
std::map<string,bool> gBoolParams;
std::map<string,int> gIntParams;
std::map<string,double> gDoubleParams;
std::map<string,string> gStringParams; // string, char, and path 


// Ordered for the help message.
// Parameters with an empty description are not displayed in the help.
const int gParameterCount = 35;
string const gParameters[gParameterCount][4] = {
 {"species.file", "path", "required", "species tree file (newick)" },
 {"gene.file", "path", "required", "gene trees file (newick)"},
 {"verbose", "bool", "false", "show progress and timing"},
 {"force.binary", "bool", "false", "make gene trees binary if they are not"},
 {"dupli.cost", "double", "2", "cost of a duplication"},
 {"HGT.cost", "double", "3", "cost of an HGT"}, 
 {"loss.cost", "double", "1", "cost of a loss"},
 {"ils.cost", "double", "1", "cost of a incomplete lineage sorting"},
 {"ils.cutoff", "double", "0", 
     "branch length cutoff for incomplete lineage sorting"},
 {"compute.T", "bool", "true", "compute transfers"},
 {"compute.TL", "bool", "true", "use transfer losses"},
 {"char.sep", "char", "_", "char separating gene names in gene tree files"},
 {"gene.mapping.file", "path", "none", 
     "file with mapping between gene names and species names"}, 
// {"transfer.dead", "bool", "true", "allow transfers from the dead"}, 
 {"dates.as.bootstraps", "bool", "false", 
     "read the species node ordering directly from the bootstrap values"},
 {"ultrametric.only", "bool", "true", 
    "return an error if a dated species tree is not ultrametric"},
 {"print.matrix.file", "string", "none", "file for output matrix"}, 
 {"print.newick.file", "string", "none", 
    "base file name for output of gene/species tree"}, 
 {"print.newick.dir", "string", "none", 
     "directory to output modified gene/species tree"}, 
 {"amalgamate", "bool", "false", "try all amalgamations"}, 
 {"force.rooting", "bool", "false", 
     "unroot gene, but don't do an amalalgamation"}, 
 {"weight.amalgamation", "double", "0", 
     "weight multipler for clade costs (amalgamation)"}, 
 {"fully.dated", "bool", "true", "tree is fully dated, subdivide the tree"}, 
 {"partially.dated", "bool", "false", 
     "consider the species tree as partially dated"}, 
 {"trim.species.tree", "bool", "false", 
        "Remove elements of the species tree above LCA of shared taxa"
        " with the gene trees."}, 
 {"allow.polytomy", "bool", "false", 
     "create combination clades from polytomies (amalgamations only)"},
 {"collapse.tree", "bool", "false", 
     "collapse nodes below a threshold to create polytomic trees"},
 {"collapse.threshold", "double", "0.5", 
     "collapse trees if nodes are below threshold"},
 {"collapse.mode", "int", "1", "0=distances, 1=bootstrap values"},
 {"tree.limit", "double", "0", 
     "maximum number of possible polytomic trees (0=no limit)"},
 {"degree.limit", "int", "0", 
     "maximum out degree of collapsed polytomic trees (0=no limit)"},
 {"collapse.decrement", "double", "0.1", 
     "amount to decrement collapse.threshold of tree.limit is surpassed"},
 {"reroot.file", "string", "none", 
    "output a rerooted input gene tree based on the amalgamation to this file"},
 {"reroot.proportion", "int", "2", 
     "ratio of distances of root sons in rerooted tree"},
 {"use.bootstrap.weighting", "bool", "false", 
     "weight split ratios with bootstrap values"},


 // THESE DON'T MAKE SENSE - unique dated idX ids are all collapsed
 {"internal.graph.ids", "bool", "false", ""}
        // "otherwise they match the reconciliation"},
        
};



void help() {
    cout << "The program outputs the cost of the most parsimonious"
            " reconcilation between the species tree and the ALE"
            " corresponding to the given gene trees." << endl;
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




// change to boost::filesystem?
static bool do_mkdir( const char *path ) {
    struct stat st;
    if( stat(path, &st) != 0 ) {
        if( mkdir( path, 0777 ) != 0 && errno != EEXIST)
            return false;
    } else if( !S_ISDIR(st.st_mode) ) { // check if it is a directory
        return false;
    }

    return true;
}
/**
 ** mkpath - ensure all directories in path exist
 ** Algorithm takes the pessimistic view and works top-down to ensure
 ** each directory in path exists, rather than optimistically creating
 ** the last element and working backwards.
 **/
int mkpath( string path ) // mode_t mode )
{
    char *copypath = strdup(path.c_str());
    char *sp;
    int status = 0;
    char *pp = copypath;
    while( status == 0 && (sp = strchr(pp, '/')) != 0 ) {
        if (sp != pp) { 
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath ); //, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    //if (status == 0)
        status = do_mkdir(path.c_str()); //, mode);
    free(copypath);
    return status;
}

struct MatchPathSeparator {
        bool operator()( char ch ) const { return ch == '/'; }
};

// change to boost::filesystem
static string basename( string const& pathname ) {
    return string( find_if( pathname.rbegin(), pathname.rend(),
                   MatchPathSeparator() ).base(), pathname.end() );
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
            gIntParams[name] = bpp::TextTools::toInt( value );
        } else if( type == "double" ) {
            gDoubleParams[name] = bpp::TextTools::toDouble( value );
        } else if( type == "path" ) {
            if( value == "required" ) {
                cerr << name << " is a required parameter." << endl;
                problem = true;
            } else if( value != "none" ) {
                if( !bpp::FileTools::fileExists( value ) ) {
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



    // PARAMETER CHECKS
    
    string newickDir = gStringParams.find("print.newick.dir")->second;
    if( newickDir != "none" ) {
        if( !mkpath( newickDir ) ) {
            cout << "ERROR: " << newickDir << " is not a directory." 
                 << endl;
            exit(1);
        }
    }


    if( gBoolParams.find("collapse.tree")->second )
        gBoolParams["allow.polytomy"] = true;

    if( gBoolParams.find("partially.dated")->second ) {
        gBoolParams["fully.dated"] = false;
        gBoolParams["dates.as.bootstraps"] = true;
    }



}

/**
 * Return the memory used in MB.
 */
void printMemory( string msg ){

    int result = -1;

#ifdef __linux__
    // This gets the value in KB.
    // VmPeak = peak virtual memory (all memory requested)
    // VmRSS = resident memory (currently in RAM)
    FILE* file = fopen("/proc/self/status", "r");
    char line[128];
    while( fgets(line, 128, file) != NULL ) {
        if (strncmp(line, "VmPeak:", 6) == 0)
            break;
    }
    fclose(file);

    // parse line
    int len = strlen(line);
    result = 0;
    for( int i=0; i<len; i++ ) {
        if( line[i] >= '0' && line[i] <= '9' ) {
            result = result*10 + line[i]-'0';
        }
    }
    result = round( result/1000 );
   
    if( gBoolParams.find("verbose")->second )
        cout << msg << " memory usage (MB): " << result << endl;

#elif __APPLE__
/*
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    if( KERN_SUCCESS != task_info(mach_task_self(),
        TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count) )
    {
            return -1;
    }
    // resident size is in t_info.resident_size;
    // virtual size is in t_info.virtual_size;
*/
#endif
}


void readGeneMappingFile(
    const char *mappingFileName, ///< path to file with mapping
    vector<string> &speciesNames, ///< array to fill
    vector<string> &geneNames,    ///< array to fill
    string &errStr )    ///< returns any errors
{
    ifstream fileStream( mappingFileName, ios::in );  
    if( !fileStream) { 
        errStr = "Failed to open file.";
        return;
    }

    string description;
	while( !fileStream.eof() ) {
        string temp;
        // Copy current line to temporary string
        getline(fileStream, temp, '\n');  

        boost::char_separator<char> sep(" ");
        Tokenizer tok( temp, sep );
        Tokenizer::iterator iter=tok.begin();
        if( *iter == "" ) // blank line
            continue;
        geneNames.push_back( *iter );
        iter++;
        if( iter==tok.end() ) {
            errStr = "Invalid line, missing species name"; 
            return;
        }
        speciesNames.push_back( *iter );
    }
}



/**
 * Run MPR and output results
 * 
 * @return a backtrack gene tree.
 */
MyGeneTree *run( 
    bool amalgamation, ///< do amalgamation if true
    MySpeciesTree *speciesTree, ///< species tree
    vector<MyGeneTree*> geneTrees, ///< gene trees
    int maxTS,  ///< maximum time slice
    int counter, ///< which gene tree this is in the gene tree file
    bool moreThanOneTree, ///< true if there is more than one gene tree
    bool returnTree,  ///< return backtrack gene tree
    bool polytomic,  ///< true if tree is non-binary
    bool noCost=false, ///< print cost
    MyGeneTree *bootstrapGeneTree=NULL, ///< gene tree with bootstrap values
    MyGeneTree *polytomyTree=NULL ) ///< polytomy tree to expand
{
    //compute the clades and tripartitions
    CladesAndTripartitions *cladesAndTripartitions;
    if( !polytomic && !amalgamation ) {
        cladesAndTripartitions = new CladesAndTripartitions( 
                   gStringParams.find("char.sep")->second[0], *(geneTrees[0]) );
    } else {
        bool overflow = false;
        string errStr = "";
        cladesAndTripartitions = 
            new CladesAndTripartitions( 
                    gStringParams.find("char.sep")->second[0], 
                    geneTrees, gBoolParams.find("verbose")->second, 
                    overflow, errStr, polytomic, bootstrapGeneTree ); 
        if( overflow ) {
            delete cladesAndTripartitions;
            cerr << "Too many possible gene trees due to polytomies" << endl;
            exit( 1 );
        }
        if( errStr != "" ) {
            delete cladesAndTripartitions;
            cerr << errStr << endl;
            exit( 1 );
        }
        if( gBoolParams.find("verbose")->second ) {
            bpp::ApplicationTools::displayTime("Computing ALE done:");
            printMemory( "Computing ALE" );
        }
    }

    if( gStringParams.find("gene.mapping.file")->second != "none" ) {
        vector<string> speciesNames;
        vector<string> geneNames;
        string errStr = "";
        readGeneMappingFile( 
                gStringParams.find("gene.mapping.file")->second.c_str(),
                speciesNames, geneNames, errStr );
        if( errStr == "" )
            cladesAndTripartitions->mClades.mapSpeciesNames( 
                                    speciesNames, geneNames, errStr );
        if( errStr != "" ) {
            delete cladesAndTripartitions;
            cerr << "Error in file " 
                << gStringParams.find("gene.mapping.file")->second
                 << ": " << errStr << endl;
            exit( 1 );
        }
    }


    //////////////////////////////////////////////////////////
    // Create the matrix
    //////////////////////////////////////////////////////////


    // compute MPR cost
    IDTLMatrix *idtlMatrix 
            = new IDTLMatrix( speciesTree, 
            cladesAndTripartitions, 
            gBoolParams.find("compute.T")->second, 
            gBoolParams.find("compute.TL")->second, 
            gDoubleParams.find("dupli.cost")->second, 
            gDoubleParams.find("HGT.cost")->second, 
            gDoubleParams.find("loss.cost")->second, 
            gDoubleParams.find("ils.cost")->second
            );


    idtlMatrix->calculateMatrix( 
            gBoolParams.find("verbose")->second, 
            !gBoolParams.find("fully.dated")->second,
            gBoolParams.find("partially.dated")->second);


    double bestCost = idtlMatrix->getBestCost();
    if( !noCost ) {
        if( bestCost == std::numeric_limits<double>::max() ) 
            cout << "Cost of a most parsimonious reconciliation: OVERFLOW" 
                 << endl;
        else
            cout << "Cost of a most parsimonious reconciliation: " 
                 << bestCost << endl;
    }
//idtlMatrix->getReconciliation();
//idtlMatrix->printRealMatrixValues();
    
    printMemory( "Matrix" );



    /////////////////////////////////////////////////
    //  output
    /////////////////////////////////////////////////

    string ext = "";
    if( !gBoolParams.find("amalgamate")->second && moreThanOneTree )
        ext = "_" + bpp::TextTools::toString( counter );


    // matrix - only print once
    string matrixFileStr = gStringParams.find("print.matrix.file")->second;
    if( matrixFileStr != "none" ) {
        idtlMatrix->printMatrixCSV( 
                (matrixFileStr+ext+".csv").c_str(),
                gStringParams.find("species.file")->second.c_str(),
                gStringParams.find("gene.file")->second.c_str() );
    }


    ///////////////////////////////////////////////////////////
    ////////// Resulting gene tree 
    ///////////////////////////////////////////////////////////
    
    bool append = false; // append newick or rerooted if not the first
    if( counter > 1 )
        append = true;

    
    // newick gene tree
    MyGeneTree *treeToReturn = NULL;
    if( returnTree 
        || gStringParams.find("print.newick.file")->second != "none" 
        || gStringParams.find("print.newick.dir")->second != "none" ) 
    {
        string pathName = gStringParams.find("print.newick.file")->second;
        if( gStringParams.find("print.newick.file")->second == "none" ) 
            pathName = gStringParams.find("print.newick.dir")->second 
                   + "/" + basename( gStringParams.find("gene.file")->second ); 

        MyGeneTree *tree = NULL;
        MyGeneNode *node = NULL;
        if( polytomic || amalgamation ) {
            node = idtlMatrix->backtrack( false );
            // create a tree to delete node structre easily
            tree = new MyGeneTree(*node); 
        } else {
            tree = cladesAndTripartitions->getTree();
            node = tree->getRootNode();
        }

        if( gStringParams.find("print.newick.file")->second != "none" 
            || gStringParams.find("print.newick.dir")->second != "none" ) 
        {
            tree->printNewick( pathName, append );
        }

        if( !returnTree )
            delete tree;
        else
            treeToReturn = tree;
    }

    if( gStringParams.find("reroot.file")->second != "none" 
            && amalgamation && geneTrees.size() == 1 ) 
    {
        MyGeneNode *newRoot = idtlMatrix->backtrack( false );
        vector<string> sonLeaves = 
                bpp::TreeTemplateTools::getLeavesNames( *(newRoot->getSon(0)) );
        MyGeneTree *geneTreeCopy = new MyGeneTree( *(geneTrees[0]) );
        bool success = geneTreeCopy->reroot( sonLeaves, 
                            gIntParams.find("reroot.proportion")->second );
        if( success )
            geneTreeCopy->printNewick( 
                    gStringParams.find("reroot.file")->second, append );
        else 
            cout << "Failed to reroot tree. =======================" << endl;
    }



    printMemory( "Final" );

    delete cladesAndTripartitions;
    delete idtlMatrix;


    return treeToReturn;
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
        map<string, string> params 
            = bpp::AttributesTools::parseOptions(args, argv);
        readParameters( params );

        if( gBoolParams.find("verbose")->second ) {
	        cout << "**********************************************"
                    "********************" << endl;
	        cout << "*                 ilsTERA, version " 
                 << version << "          *" << endl;
	        cout << "* Authors: C. Scornavacca, E. Jacox         "
                    "Created     20/10/15 *" << endl;
	        cout << "*                                           "
                    "Last Modif. 20/10/15 *" << endl;
	        cout << "**********************************************"
                    "********************" << endl;
	        cout << endl;
        }


        if( gBoolParams.find("verbose")->second ) 
            bpp::ApplicationTools::startTimer();


        ////////////////////////////////////////
        // Get trees and check them and other input 
        ////////////////////////////////////////

        // read species tree
        string errString = "";
        MySpeciesTree* speciesTree = MySpeciesTree::readMySpeciesTree( 
                gStringParams.find("species.file")->second.c_str(),
                errString, 
                gBoolParams.find("dates.as.bootstraps")->second );
        if( errString != "" || speciesTree == NULL ) {
            cerr << "Error reading species tree: " << errString << endl;
            exit(1);
        }

        //to have random dates
        //TreeTools::computeBranchLengthsGrafen( *speciesTree, 1, true);
		//TreeTools::convertToClockTree( *speciesTree, 
        //                          speciesTree->getRootNode()->getId());
		//Newick * print = new Newick(false,false);
		//print->write( *speciesTree, cout);		

		
        // species tree must be binary
        int nonBinaryCount = speciesTree->getBinaryCount();
        if( nonBinaryCount != 0 ) {
            cerr << "ERROR: Species tree is not binary." << endl;
            cerr << nonBinaryCount << " nodes are not binary." << endl;
            exit(1);
        }

        // species tree must have unique leaves
        string dupName;
        if( !speciesTree->uniqueLeaves( dupName ) ) {
            cerr << "ERROR: Species tree has duplicate leaf names: <" 
                 << dupName << ">" << endl;
            exit(1);
        }

        // map of taxaNames required for the function restrictTreeToASetOfTaxa
        boost::unordered_map<string, int> taxaNamesSpecies;
        vector<string> leaves = speciesTree->getLeavesNames();
        BOOST_FOREACH( string leaf, leaves) 
            taxaNamesSpecies.insert(make_pair(leaf,1));

		//reading the gene trees
        vector<MyGeneTree*> geneTrees;
        geneTrees = MyGeneTree::readMyGeneTrees( 
                gStringParams.find("gene.file")->second.c_str(),
                errString );
        if( errString != "" ) {
            cerr << "Error reading gene trees: " << errString << endl;
            exit(1);
        }
        if( gBoolParams.find("verbose")->second ) 
            cout << geneTrees.size() << " gene trees" << endl;

        boost::unordered_map<string, int> taxaNamesGenes;
        if( gBoolParams.find("trim.species.tree")->second ) {
            // create map of all taxa name in the gene trees
            BOOST_FOREACH( MyGeneTree *geneTree, geneTrees ) {
                vector<string> leaves = geneTree->getLeavesNames();
                BOOST_FOREACH( string leafName, leaves) {
                    size_t pos = leafName.find( 
                                gStringParams.find("char.sep")->second[0] ); 
                    string taxaName = leafName.substr(0,pos);
                    taxaNamesGenes.insert( make_pair(taxaName,1) );
                }
            }
            
            if( !gBoolParams.find("fully.dated")->second ) {
                // trim species tree
                if( !speciesTree->trimTree( taxaNamesGenes, 
                            gBoolParams.find("verbose")->second ) ) 
                {
                    cerr << "Species tree has none of the gene taxa." << endl;
                    exit(1);
                }
            }
        }


        ////////////////////////////////////////
        // Process species tree
        ////////////////////////////////////////

/*
        // to add the outGroup to transfer from the dead 
        // create a sibling of the root and a new root
        if( gBoolParams.find("transfer.dead")->second )  
            speciesTree->addAlphaForDeadTransfer( 
                    gBoolParams.find("dates.as.bootstraps")->second,
                    gDoubleParams.find("HGT.cost")->second, 
                    gDoubleParams.find("loss.cost")->second );
*/       

        // must be done before subdivision
        speciesTree->compute_RealPostOrder();  

        // subdivision
        if( gBoolParams.find("fully.dated")->second ) {
            string errStr = "empty error";
            vector<int> changedTimeSlices;
            vector< pair<int,int> > dateMap;
            bool success = speciesTree->computeSubdivision( dateMap, 
                            gBoolParams.find("dates.as.bootstraps")->second, 
                            gBoolParams.find("ultrametric.only")->second,
                            changedTimeSlices, errStr ); 
            if( !success ) {
                cerr << "ERROR: " << errStr << endl;
                exit(1);
            }

            if( gBoolParams.find("trim.species.tree")->second ) {
                /*
                if( !gBoolParams.find("transfer.dead")->second ) {
                    cerr << "Trimming with dated species tree requires "
                        << "transfer from the dead." << endl;
                    // The dated tree can have transfers at each timeslice,
                    // which is a node in the original species tree. Trimming
                    // the species tree removes these timeslices and therefore
                    // possiblities to transfer at these times. To compensate,
                    // transfer from the dead can be used.
                    exit(1);
                }
                */

                // trim species tree
                //if(!speciesTree->restrictTreeToLCA(taxaNamesGenes, gBoolParams.find("verbose")->second))
                if( !speciesTree->trimTree( taxaNamesGenes, 
                    gBoolParams.find("verbose")->second ) ) 
                {
                    cerr << "Species tree has none of the gene taxa." << endl;
                    exit(1);
                }
            }
        } else {
            speciesTree->assignNoSubdivisionTimeSlices();
        }

        
		int maxTS = speciesTree->setVectorTimeSlices();
		
        // assign ids (correspondance)
        speciesTree->assignPostOrderIds();

        speciesTree->computeSpeciesCladesAndSplits( 
                gDoubleParams.find("ils.cutoff")->second,
                15, true );


        //////////////////////////////////////////////////////////
		// Gene trees
        //////////////////////////////////////////////////////////
        int counter = 0;
        BOOST_FOREACH( MyGeneTree *geneTree, geneTrees ) {
            counter++;

            if( !geneTree->restrictTreeToASetOfTaxa( 
                        taxaNamesSpecies, 
                        gStringParams.find("char.sep")->second[0], 
                        gBoolParams.find("verbose")->second ) ) 
            {
                cerr << "ERROR: Gene tree " << counter
                     << " has no valid leaves." << endl;
                exit(1);
            }

            // needed by the algoritm!!!	
            if( geneTree->getNumberOfLeaves()<3 ) {
                cerr << "ERROR: Gene tree " << counter
                     << " has only " << geneTree->getNumberOfLeaves()
                     << " leaves. At least three are required." << endl;
                //exit(1);
                continue;
            }

            // gene trees must have unique leaves
            string dupName;
            if( !geneTree->uniqueLeaves( dupName ) ) {
                cerr << "ERROR: A gene tree " << counter
                    << " has duplicate leaf names: <" 
                     << dupName << ">" << endl;
                exit(1);
            }


            bool polytomic = false; // non-binary = false
            bool amalgamation = false;
            MyGeneTree *bootstrapGeneTree = NULL;
            if( gBoolParams.find("allow.polytomy")->second ) {
                amalgamation = false;
                polytomic = true;


                if( gBoolParams.find("use.bootstrap.weighting")->second ) 
                    bootstrapGeneTree = new MyGeneTree( *geneTree );
                if( gBoolParams.find("collapse.tree")->second ) {
                    if( gIntParams.find("collapse.mode")->second != 0 
                        && gIntParams.find("collapse.mode")->second != 1 ) 
                    { 
                        cerr << "ERROR: collapse.mode not 0 or 1" << endl;
                        exit(1);
                    }
                    double threshold = 
                        gDoubleParams.find("collapse.threshold")->second;
                    int maxDegree = 0;
                    bool overflow;
                    bool firstLoop = true;
                    bool continueLoop = true;
                    unsigned long long numberOfTrees = 0;
                    MyGeneTree *geneTreeCopy = NULL;
                    while( firstLoop || continueLoop ) {
                        if( !firstLoop ) {
                            if( gDoubleParams.find("tree.limit")->second!=0 ) 
                                cout << "number of trees (" << numberOfTrees
                                    << ") greater than limit, trying threshold "
                                    << threshold << endl;
                            if( gIntParams.find("degree.limit")->second!=0 ) 
                                cout << "max out degeree (" << maxDegree 
                                    << ") greater than limit, trying threshold "
                                    << threshold << endl;
                            // revert to copy
                            delete geneTree;
                            geneTree = geneTreeCopy;
                            geneTreeCopy = NULL;
                        }

                        if( gDoubleParams.find("tree.limit")->second != 0 
                            || gIntParams.find("degree.limit")->second != 0 ) 
                        {
                            // copy original tree in case 
                            // gDoubleParams.find("tree.limit")->second exceeded
                             geneTreeCopy = new MyGeneTree( *geneTree );
                        }
                        maxDegree = geneTree->collapseOnTree( 
                                threshold, 
                                gIntParams.find("collapse.mode")->second );
                        numberOfTrees
                            = CladesAndTripartitions::findPolytomies( 
                                        geneTree->getRootNode(), overflow );
                        continueLoop = false;
                        if( gIntParams.find("degree.limit")->second != 0 
                            || gDoubleParams.find("tree.limit")->second != 0 ) 
                        {
                            if( overflow 
                               || (gDoubleParams.find("tree.limit")->second!=0 
                                   && numberOfTrees>gDoubleParams.find(
                                        "tree.limit")->second)
                               || (gIntParams.find("degree.limit")->second!=0 
                                   && maxDegree>gIntParams.find(
                                       "degree.limit")->second) )
                            {
                                continueLoop = true;
                            }
                        } else if( overflow ) {
                            cerr << "Overflowed number of possible polytomic"
                                 << " trees - exiting." << endl;
                            exit(1);
                        }
                        threshold -= 
                            gDoubleParams.find("collapse.decrement")->second; 
                        firstLoop = false;
                    }
                    if( geneTreeCopy != NULL )
                        delete geneTreeCopy;
//geneTree->printNewick( "collapsedTree" );
                }
            } else {
                // Rearrange tree if it is unrooted (3 root sons)
                // and use clades by setting amalgamation = true.
                amalgamation = geneTree->rootTree();
                if( gBoolParams.find("force.rooting")->second )
                    amalgamation = false;
                
                if( !geneTree->isBinary() ) {
                    if( gBoolParams.find("force.binary")->second )
                        geneTree->makeBinary();
                    else {
                        cerr << "ERROR: Gene tree " << counter
                             << " is not binary." << endl;
                        exit(1);
                    }
                }
            }


           
            // do each tree individually if this isn't an amalgamation
            if( !gBoolParams.find("amalgamate")->second ) {
                vector<MyGeneTree*> singleGeneTree;
                singleGeneTree.push_back( geneTree );
                string rStr = "";
                if( amalgamation ) {
                    rStr = " unrooted"; 
                    /*if( gBoolParams.find("construct.graph")->second ) {
                        cerr << "ERROR: Cannot construct a graph for unrooted"
                                " gene trees (tree " << counter
                                << ")" << endl;
                        exit(1);
                    }*/
                }
                if( gBoolParams.find("verbose")->second )
                    cout << "Gene tree " << counter << rStr << endl;
                bool returnTree = false;
                MyGeneTree *tree = run( amalgamation, speciesTree, 
                        singleGeneTree, 
                        maxTS, 
                        //changedTimeSlices, 
                        counter, geneTrees.size() > 1,
                        //constructGraph, 
                        returnTree, polytomic, false,
                        bootstrapGeneTree );
                if( returnTree ) {
                    // construct graph from backtracked tree
                    vector<MyGeneTree*> backtrackGeneTree;
                    backtrackGeneTree.push_back( tree );
                    run( false, speciesTree, backtrackGeneTree, maxTS, 
                            //changedTimeSlices, 
                            0, false, 
                            //true, 
                            false, false,
                            true, NULL, geneTree );
                    delete tree;
                }
                if( gBoolParams.find("verbose")->second )
                    cout << endl;
            }

            if( bootstrapGeneTree != NULL )
                delete bootstrapGeneTree;
		}
        if( geneTrees.size() == 0 ) {
            cout << "No gene trees found." << endl;
            exit(1);
        }

        if( gBoolParams.find("verbose")->second ) {
            bpp::ApplicationTools::displayTime("Reading the gene trees done:");
        }


        if( gBoolParams.find("amalgamate")->second ) {
            run( true, speciesTree, geneTrees, maxTS, 
                   //changedTimeSlices, 
                   0, false, false, 
                   //gBoolParams.find("construct.graph")->second, 
                   false );
        }

        BOOST_FOREACH( MyGeneTree *tree, geneTrees ) 
            delete tree;

        // print a species file with post order ids if requested for graphs
        if( gStringParams.find("print.newick.dir")->second != "none" )
        {
            string pathName = 
                gStringParams.find("print.newick.file")->second + "_species";
            if( gStringParams.find("print.newick.file")->second == "none" ) 
                pathName = gStringParams.find("print.newick.dir")->second 
                     + "/" + basename( 
                    gStringParams.find("species.file")->second.c_str() );

            MySpeciesTree *tree = speciesTree->getPostorderTree();
            tree->printNewick( pathName );
            delete tree;
        }
	 
        delete speciesTree;

        if( gBoolParams.find("verbose")->second )
		    bpp::ApplicationTools::displayTime("Done:");

	}  catch(exception & e) {
		cerr << e.what() << endl;
		exit(-1);
	}
	
	return 0;
}








