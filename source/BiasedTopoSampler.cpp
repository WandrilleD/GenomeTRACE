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

This file contains a small executable to play with biased topology sampling in a CCP distribution

Created the: 15-06-2016
by: Wandrille Duchemin

Last modified the: 15-06-2016
by: Wandrille Duchemin

*/


#include "MyCladesAndTripartitions.h"
#include "MyGeneTree.h"

#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>


#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Utils/AttributesTools.h>


const int gParameterCount = 1;
string const gParameters[gParameterCount][4] = {
// files
 {"gene.trees", "path", "required", "a file with some gene trees" }
};


//name, type, default, description
map<string,bool> gBoolParams;
map<string,int> gIntParams;
map<string,double> gDoubleParams;
map<string,string> gStringParams; // string, char, and path 



/**
 * read input parameters
 */
void readParameters( 
        map<string,string> &params ) ///< input paramsters
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
            cerr << "Found type <" << type << "> for " << name << endl;
            throw bpp::Exception( "readParameters: unknown type" );
        }
    }

    if( problem )
        exit(1);

    // PARAMETER CHECKS -> forcing some dependent parameters

}

/////////////////////////////////////////////////
// Main 
/////////////////////////////////////////////////

int main(int args, char ** argv)
{

	if(args == 1)
	{
		exit(0);
	}
	
	try 
	{
        // fill global parameter variables
        map<string, string> params = 
            bpp::AttributesTools::parseOptions(args, argv);
        readParameters( params );

		char charSep = '_';
		bool overflow = false;
		string errString = "";
		bool verbose = true;
		///////////////////////////////////////
		/// Reading data
		///////////////////////////////////////

		vector<MyGeneTree*> geneTrees;
	
		geneTrees = MyGeneTree::readMyGeneTrees( gStringParams.find("gene.trees")->second.c_str(), errString );//reading the tree list
		if( errString != "" ) 
		{
			cerr << "Error reading gene trees: " << errString << endl;
			exit(1);
		}
		if( verbose ) 
			cout << geneTrees.size() << " gene trees" << endl;


		MyCladesAndTripartitions CCPDistrib( charSep, geneTrees, verbose, overflow, errString, false );

		CCPDistrib.printMe();

		if( errString != "" ) 
		{
			cerr << "Error creating CCP distrib: " << errString << endl;
			exit(1);
		}
		if(overflow)
		{
			cout << "overflow!" << endl;
			exit(1);
		}


		//CCPDistrib.printNonCompatible();

		vector < vector <string> > cladesToForce;
		cladesToForce.push_back( vector <string>() );

		cladesToForce[0].push_back("A_a");
		cladesToForce[0].push_back("B_b");
		cladesToForce[0].push_back("E_e");

		CCPDistrib.RestrictToClade( cladesToForce );


		CCPDistrib.printMe();

		
		MyGeneTree * Gtree;
		for(unsigned i = 0; i < 10 ; i++)
		{
			Gtree = CCPDistrib.getRandomTree();
			Gtree->printNewick("res.txt",true);
		}


    }
	catch(exception & e)
	{
		cerr << e.what() << endl;
		exit(-1);
	}
	
	return 0;


}
