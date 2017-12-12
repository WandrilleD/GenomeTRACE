// Created by: Wandrille Duchemin

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
@author Celine Scornavacca
@author Edwin Jacox
@version 1
@date 10-02-2017

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


#include "DeCoUtils.h"
#include "WmodifUtils.h"

#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>


#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Utils/AttributesTools.h>



string version = "1";

//name, type, default, description
map<string,bool> gBoolParams;
map<string,int> gIntParams;
map<string,double> gDoubleParams;
map<string,string> gStringParams; // string, char, and path 



const int gParameterCount = 45;
string const gParameters[gParameterCount][4] = {
// files
 {"parameter.file", "path", "none", "a file with input parameters" },
 {"species.file", "path", "required", "species tree file (newick). NB: default behavior wants it to be ultrametric if transfers are used (use dated.species.tree=0 to circumvent)" },
 {"gene.distribution.file", "path", "required", "gene distribution files file (one file name per line)"},
 {"adjacencies.file","path","required","adjacencies file (one adjacency per line; leafnames separated by a space)"},

// basic options
 {"with.transfer","bool","true","allow transfers the reconciliation and adjacency histories reconstruction."},
 {"dated.species.tree","bool","true","the species tree is ultrametric and dites will be used to subdivide the trees in time slices."},

// format of input
 {"char.sep", "char", "_", "char separating gene names in gene tree files. One character only."},
 {"ale","bool","false","gene tree distribution are ALE files"},
 {"already.reconciled","bool","false","gene tree distribution are reconciled gene trees in recPhyloXML format. Will skip the reconciliation phase"},
 {"verbose","bool","false","show progress and timing"},//10
 {"superverbose","bool","false","show lots of text about progress and timing"},

//reconciliation options
 {"dupli.cost", "double", "2", "cost of a duplication"},
 {"HGT.cost", "double", "3", "cost of an HGT"}, 
 {"loss.cost", "double", "1", "cost of a loss"},
 {"try.all.amalgamation","bool","true","try all possible amalgamation when reconciling gene trees. Otherwise the best possible tree is used"},

//DeCo option
 {"AGain.cost","double","2","cost of an adjacency gain"},
 {"ABreak.cost","double","1","cost of an adjacency break"},
 {"C1.Advantage","double","0.5","between 0 and 1. Probability to choose C1 (presence of adjacency) over C0 (absence of adjacency) in case of a score tie at the root of an equivalence class"},

//Boltzmann
// {"use.boltzmann","bool","false","use Boltzmann sampling for the adjacencies history computation"},
// {"boltzmann.temperature","double","1","Temperature to use in the Boltzmann sampling (if used)"},
// {"nb.sample","int","1","number of samples to get from the adjacency matrix in the case of Boltzmann sampling"},

//output options
 {"write.newick","bool","false","use newick format rather than phyloXML-like format"},
 {"hide.losses.newick","bool","false","if true, losses and the branches leading to them wil be removed from the newick string"},//20
 {"write.adjacencies","bool","true","write the adjacencies inferred in ancestral species"},
 {"output.dir", "string", "none", "directory for printed files"}, 
 {"output.prefix", "string", "none", "A prefix to prepend to all output files."},


// scaffolding
  {"scaffolding.mode","bool","false","use scaffolding algorithm to improve extant genomes scaffolding/assembly"},
  {"chromosome.file","path","none","used with the scaffolding.mode option. A file containing the number of chromosome in each each species (one species per line, each line comprised name of the species followed by the number of chromosome, separated by a tabulation)"},

//advanced
 {"all.pair.equivalence.class","bool","false","compute adjacency histories for all pair of gene families (even if they share no adjacencies)."},
 {"bounded.TS","bool","false","use bounded time slices in adjacency history computations"},
 {"always.AGain","bool","true","always put an Ajdacency Gain at the top of an equivalence class tree"},
 {"absence.penalty","double","-1","if set to -1 (the default), nothing changes. Otherwise, specify the cost of having an adjacency at a pair of leaves which are not part of the list of adjacencies given at initialization"},
 {"substract.reco.to.adj","bool","false","if set to 1, the weighted cost of a reconciliation event will be used to favor co-event in the adjacency matrix computation. Unavalaible for Boltzmann computation."},
 
 {"Topology.weight","double","1","weight of the topology in the global score" },
 {"Reconciliation.weight","double","1","weight of the reconciliation in the global score" },//30
 {"Adjacency.weight","double","1","weight of the adjacency in the global score" },

 {"sampling.temp","double","1","temperature to accept a new sampled solution" },
 {"sampling.nbMax","int","0","maximum number of sample to draw" },


//loading previous instance
 {"load.save" , "bool" , "false" ,"wether to load a previous instance or not"},
 {"load.folder" , "string" , "none","name of folder where a previous instance can be found (used with load.save=1)."},
 {"load.prefix" , "string" , "none","prefix of file of a previous instance (used with load.save=1)."},

 {"proba.flat","double","0.01","probability to draw bipartitions in a flat tree distribution rather than the CCP distribution when redrawing gene topologies" },
 {"proba.bias","double","0","probability to bias the CCP distribution to force it to keep the clade where an ancestral adjacency was computed (no effect is 0)" },
 {"continue.sampling.while.accepted","bool","false","If true, the system will continue to run if at least 1 new gene family was succesfully updated during this round" },
 {"redraw.algo","int","0","0 : only change reconciliation (doesn't work right now) ; 1 : ccp only redraw ; 2 : externally calls DTLRecCoev"},

 {"DTLRecCoevExecPath","string","", "used only when redraw.algo=2. specifies the path where the external DTLRecCoev executable can be found"},

 {"DTLRecCoev.NoCoev","bool","false","doesn't take potential co-event in account when externally calling DTLRecCoev"},
 {"DTLRecCoev.Temp","double","0.5","temperature of DTLRecCoev"}

};


void help() {
    cout << "This progam computes and outputs the most parsimonious"
            " history of a set of given adjacencies between extant genes,"
            " according to their gene trees and their reconciliation with"
            " a species tree." << endl;
    cout << "version: " << version << endl;
    cout << "__________________________________________________________"
            "_____________________" << endl;
    cout << "___________________________________INPUT__________________"
            "_____________________" << endl;
    cout << "format: parameter.name | [type,default] <description>" << endl;

    for( size_t i=0; i<gParameterCount; i++ ) {
        if( gParameters[i][3] != "" )
        {
            cout << gParameters[i][0] ;
            for (size_t j = 0; j <  50 - gParameters[i][0].size() ; j++)
            	cout << " ";
        	cout << "| [" << gParameters[i][1] << ","
                 << gParameters[i][2] << "] " << gParameters[i][3] << endl;
		}
    }

}

string gPathPrefix = ""; // path and prefix for all output files
string outputDir ="";


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

/*struct MatchPathSeparator {
        bool operator()( char ch ) const { return ch == '/'; }
};*/


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

    // check for a parameter file
    iterStr = params.find( "parameter.file" );	
    if( iterStr != params.end() ) {
        map<string,string> fileParams;
        bpp::AttributesTools::getAttributesMapFromFile( iterStr->second,
                                                       fileParams, "=" );
        // add file parameters if not already there
        map<string,string>::iterator fIter;
        for( fIter=fileParams.begin(); fIter!=fileParams.end(); ++fIter) 
        {
            iterStr = params.find( fIter->first );
            if( iterStr == params.end() ) 
                params[fIter->first] = fIter->second;
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

    if(gBoolParams.find("superverbose")->second && !gBoolParams.find("verbose")->second)
    {
    	cout << "superverbose : forcing verbose=true" << endl;
    	gBoolParams["verbose"]=true;
    }

    if( !gBoolParams.find("with.transfer")->second) //no transfer
    {
    	if(gBoolParams.find("bounded.TS")->second)
    	{
    		cout << "No transfer : forcing bounded.TS=false" << endl;
    		gBoolParams["bounded.TS"] = false;
    	}
    	if(gBoolParams.find("dated.species.tree")->second)
    	{
    		cout << "No transfer : forcing dated.species.tree=false" << endl;
    		gBoolParams["dated.species.tree"] = false;	
    	}
    }
    else // with transfers.
    {
    	if(!gBoolParams.find("dated.species.tree")->second)
    		if(gBoolParams.find("bounded.TS")->second)
    		{
    			cout << "Non dated species tree : forcing bounded.TS=false" << endl;
    			gBoolParams["bounded.TS"] = false;
    		}
    }



	if( gBoolParams.find("ale")->second && gBoolParams.find("already.reconciled")->second)
	{
		cout << "already.reconciled : forcing ale=false" << endl;
		gBoolParams["ale"] = false;	
	}


	if( gBoolParams.find("substract.reco.to.adj")->second  && gDoubleParams.find("Adjacency.weight")->second == 0) //substract reco score but adj weight = 0 -> division by 0!!
	{
		cout << "substract.reco.to.adj=true and Adjacency.weight=0. Will set Adjacency.weight to 0.000001 in Adjacency matrix computation to avoid division by 0." << endl;
		// modification on the fly
	}

	//ensuring 0 <= C1.advantage <= 1
	if(gDoubleParams.find("C1.Advantage")->second < 0 )
	{
		cout << "Negative value of C1.Advantage: forcing C1.Advantage=0" << endl;
		gDoubleParams["C1.Advantage"] = 0;
	}
	else if (gDoubleParams.find("C1.Advantage")->second > 1 )
	{
		cout << "value of C1.Advantage > 1: forcing C1.Advantage=1" << endl;
		gDoubleParams["C1.Advantage"] = 1;
	}

    if(gBoolParams.find("hide.losses.newick")->second)
        if(!gBoolParams.find("write.newick")->second)
        {
            cout << "hide.losses.newick = 1 : forcing write.newick=1" << endl;
            gBoolParams["write.newick"] = true;
        }

    // create path prefix for output files
    if( gStringParams.find("output.prefix")->second != "none" ) 
        gPathPrefix = gStringParams.find("output.prefix")->second;
    outputDir = gStringParams.find("output.dir")->second;
    if( outputDir != "none" ) {
        if( !mkpath( outputDir ) ) {
            cout << "ERROR: " << outputDir << " is not a directory." 
                 << endl;
            exit(1);
        }
        gPathPrefix = outputDir + "/" + gPathPrefix;
    }

    if( gDoubleParams.find("DTLRecCoev.Temp")->second <= 0)
    {
        cerr << "Negative or null value of DTLRecCoev.Temp which must be a positive temperature : forcing back to default = 0.5" << endl;
        gDoubleParams["DTLRecCoev.Temp"] = 0.5;
    }



}


void printParameters()
{
    int lineLimit = 50;
	for (map<string,bool>::iterator it=gBoolParams.begin(); it!=gBoolParams.end(); ++it)
	{
    	cout << it->first ;
		for( unsigned j = 0; j < lineLimit - it->first.size(); j++)
			cout <<".";
    	cout << it->second << endl;
    }

	for (map<string,int>::iterator it=gIntParams.begin(); it!=gIntParams.end(); ++it)
	{
    	cout << it->first ;
		for( unsigned j = 0; j < lineLimit - it->first.size(); j++)
			cout <<".";
    	cout << it->second << endl;
    }

	for (map<string,double>::iterator it=gDoubleParams.begin(); it!=gDoubleParams.end(); ++it)
	{
    	cout << it->first ;
		for( unsigned j = 0; j < lineLimit - it->first.size() ; j++)
			cout <<".";
    	cout << it->second << endl;
    }

	for (map<string,string>::iterator it=gStringParams.begin(); it!=gStringParams.end(); ++it)
	{
    	cout << it->first ;
		for( unsigned j = 0; j < lineLimit - it->first.size() ; j++)
			cout <<".";
    	cout << it->second << endl;
    }


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
	
	try 
	{
        // fill global parameter variables
        map<string, string> params = 
            bpp::AttributesTools::parseOptions(args, argv);
        readParameters( params );

		bool randomtree = false;	
		bool dateAsBootstrap = false; // a bit harsh -> make option?
		bool verbose = gBoolParams.find("verbose")->second;
		bool superverbose = gBoolParams.find("superverbose")->second;

        if( verbose ) {
	        cout << "**********************************************"
                    "********************" << endl;
	        cout << "*                 DeCoSTAR, version " 
                 << version << "          *" << endl;
	        cout << "* Authors: W. Duchemin, C. Scornavacca, E. Jacox         "
                    "Created     02/03/16 *" << endl;
	        cout << "*                                           "
                    "Last Modif  14/03/16 *" << endl;
	        cout << "**********************************************"
                    "********************" << endl;
	        cout << endl;
        }

		if( verbose ) 
			printParameters();


        clock_t begin;
        clock_t end;
        double elapsed_secs;


		///////////////////////////////////////
		/// Reading data
		///////////////////////////////////////

		MySpeciesTree * speciesTree = getSpeciesTree(gStringParams.find("species.file")->second, dateAsBootstrap);
		int maxTS = processSpeciesTree( speciesTree ,  dateAsBootstrap, gBoolParams.find("dated.species.tree")->second , gBoolParams.find("with.transfer")->second, gBoolParams.find("with.transfer")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second, verbose );




		vector<string> geneFiles = readGeneDistributionsFile( gStringParams.find("gene.distribution.file")->second );
		
		vector <GeneFamily *> * GeneFamilyList = new vector <GeneFamily *>;

		readGeneDistributions( GeneFamilyList, speciesTree , geneFiles, gBoolParams.find("ale")->second ,  gBoolParams.find("already.reconciled")->second, gStringParams.find("char.sep")->second[0] ,verbose, superverbose);


        for(unsigned i = 0 ; i < GeneFamilyList->size();i++)
            GeneFamilyList->at(i)->setDefaultMaxLikelihood();//quick and dirty
            //GeneFamilyList->at(i)->setMaxLikelihood(); // in order to get a normalized score


        map < string, map <string , double> > adjacencyScores; // will stay empty, just here to ensure compat with deco's current version
		vector< pair <string,string > > adjacencies = readAdjacencies( gStringParams.find("adjacencies.file")->second );

		if(verbose)
			cout << "read " << adjacencies.size() << " adjacencies." << endl;


        string LoadPrefix = ""; // used only if load.save == true

        if(gBoolParams.find("load.save")->second) //we load a previous DeCo instance reconciled gene trees
        {
            if(verbose)
            {
                cout << "loading a save";
            }

            string LoadFolder = gStringParams.find("load.folder")->second;
            if(LoadFolder !="none") //name of folder where a previous instance can be found (used with load.save=1)..
            {
                LoadPrefix  += LoadFolder + "/";
                if(verbose)
                    cout << " from folder: "<< LoadFolder + "/";
            }

            string LoadFilePrefix = gStringParams.find("load.prefix")->second;
            if(LoadFilePrefix !="none") //prefix of file of a previous instance (used with load.save=1).3
            {
                LoadPrefix  += LoadFilePrefix;
                if(verbose)
                    cout << " with prefix: "<< LoadFilePrefix;
            }
            cout << endl;


            LoadGeneFamilyReconciledTrees( GeneFamilyList, LoadPrefix, speciesTree, superverbose ,true); // for the "one file" stuff

            for(unsigned i =0; i <  GeneFamilyList->size(); i++ )
                GeneFamilyList->at(i)->setRecScore(gDoubleParams.find("dupli.cost")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second ); // Counts events in the RecTree and uses the costs to compute the score





            if(verbose)
                cout << "reconciled tree loaded" << endl;




        }
        else // we compute an initial solution by considering each family / adjacency equivalence class family independantly
        {

            ///////////////////////////////////////
    		/// Computing MPR for each gene family
            ///////////////////////////////////////

            //declaring a bunch of default args

    		double gWeight = 0;	  ///< split weight
    		if(gBoolParams.find("try.all.amalgamation")->second)
    			gWeight = gDoubleParams.find("Topology.weight")->second;

		

    		for(int i = 0; i < GeneFamilyList->size(); i++)
    		{
    			GeneFamily * Gfamily = GeneFamilyList->at(i);
    
    			if(!gBoolParams.find("already.reconciled")->second)//only reconcile if it is not already done
    			{
    				if(!gBoolParams.find("try.all.amalgamation")->second) // if not all amalgamation are tried, we have to set the gene tree
    					Gfamily->makeUnrootedTree(randomtree);
    	
    				Gfamily->makeReconciliation(speciesTree,
    											gBoolParams.find("with.transfer")->second, gBoolParams.find("with.transfer")->second,  // transfer and transferLoss
    											gDoubleParams.find("dupli.cost")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second , // costs
    											maxTS, gWeight, gBoolParams.find("dated.species.tree")->second, //here, dated.species.tree means that species tree is subdivided
    											gBoolParams.find("try.all.amalgamation")->second);
    			}
    			else
                    Gfamily->setRecScore(gDoubleParams.find("dupli.cost")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second );
    


            }
        }

        map <string,int> LeafToSpMap; // used in scaffolding mode

        for(unsigned i =0; i <  GeneFamilyList->size(); i++ )
        {

            if(gBoolParams.find("bounded.TS")->second)
                GeneFamilyList->at(i)->setReconciledTreeTimeSlicesToBTS(speciesTree); // set the tree to a bounded time slice one (BTS)
            else if(gBoolParams.find("dated.species.tree")->second)
                GeneFamilyList->at(i)->setReconciledTreeTimeSlicesToTS(speciesTree); // set the tree to a complete time sliced one (TS)
            else
                GeneFamilyList->at(i)->setReconciledTreeTimeSlicesToNoTS();//no TS

            if(verbose)
            {
                cout << "Tree computed and reconciled for distribution " << i << ". Likelihood: " << GeneFamilyList->at(i)->getTreeLikelihood() << ". Reconciliation score: " << GeneFamilyList->at(i)->getRecScore() << " TS status " << GeneFamilyList->at(i)->getReconciledTreeTimeSliceStatus()<< endl;
                if(superverbose)
                    GeneFamilyList->at(i)->printRecTree();
            }

            //GeneFamilyList->at(i)->dumpStuff(); // saves memory by dumping CCP distrib and potential RecMat


            if(gBoolParams.find("scaffolding.mode")->second)
                fillLeafToSpMap( GeneFamilyList->at(i)->getRecTree(), LeafToSpMap);
        }




        // Create map of chromosome number expected by species
        map < string, int > speciesChrNb;
        map < string, int >::iterator itChr;
        if(gBoolParams.find("scaffolding.mode")->second)
        {
            if(gStringParams.find("chromosome.file")->second != "none"){
                speciesChrNb = storeSpeciesChrNumber(gStringParams.find("chromosome.file")->second);
                if(verbose){
                    cout<<endl<<endl<<"Chromosome number expected by species:"<<endl;
                    for(itChr=speciesChrNb.begin();itChr!=speciesChrNb.end();itChr++)
                    {
                        cout<<"\t"<<(*itChr).first<<"\t-> "<<(*itChr).second<<" chr"<<endl;
                    }
                }
            }
            else{
                cerr<<"ERROR: If scaffolding mode is ON, then you have to give as input a chromosome file (See param: chromosome.file)!"<<endl;
                exit(EXIT_FAILURE);
            }
        }



        ////////////////////////////////////////////////////
        /// Putting adjacencies into Equivalence Class
        ////////////////////////////////////////////////////

        
        vector < EquivalenceClassFamily > * ECFams = BasicCreateEquivalenceClassFamilies( adjacencies, GeneFamilyList, 
                                                                                            verbose, superverbose);

        //CreateEquivalenceClassFamilies( adjacencies, GeneFamilyList, 
        //                                                         ( gBoolParams.find("with.transfer")->second || gBoolParams.find("scaffolding.mode")->second ) ,
        //                                                         gBoolParams.find("all.pair.equivalence.class")->second , 
        //                                                         verbose, superverbose);


        // Dictionary that will contain c1 and c0 cost for EXT/EXT case for each species if using ARt-DeCo algorithm
        map<int,vector<float> > speciesC0C1;
        map<int,vector<float> >::iterator it1;
        map<int, map<string,int> > speGeneAdjNb;
        map<int,string> species_id_name;

        

        if(gBoolParams.find("scaffolding.mode")->second)
        {
            if(verbose)
                cout<<endl<<"Scaffolding mode enabled"<<endl;
            //ECFams = CreateEquivalenceClassFamilies( adjacencies, GeneFamilyList, speciesTree, speciesChrNb, speciesC0C1, species_id_name, speGeneAdjNb, gDoubleParams.find("ABreak.cost")->second, gBoolParams.find("with.transfer")->second , gBoolParams.find("all.pair.equivalence.class")->second, verbose, superverbose);

            //WARNING: Use species tree that can be changed !!!
            vector<MySpeciesNode*> extantSpecies = speciesTree->getLeaves();
            for(unsigned j =0; j < extantSpecies.size() ; j++)
            {
                string name=extantSpecies[j]->getName();
                int id=speciesTree->getRPO(extantSpecies[j]->getId());
                species_id_name[id]=name;
            }

            map < string, map <string , double> > * ASP = &adjacencyScores;

            bool problem = computeArtDeCoMaps(adjacencies, LeafToSpMap, speciesTree, speciesChrNb, speciesC0C1,
                                species_id_name, speGeneAdjNb, gDoubleParams.find("ABreak.cost")->second,verbose, superverbose,
                                ASP ,false);

            if(problem)
            {
                cerr << "scaffolding mode: less observed contigs than expected chromosome in some species."<<endl;
                //gBoolParams["scaffolding.mode"] = false;
            }
        }
        else
        {
            //cout<<endl<<"#################\n### DeCo mode ###\n#################"<<endl;
            //ECFams = CreateEquivalenceClassFamilies( adjacencies, GeneFamilyList, gBoolParams.find("with.transfer")->second , gBoolParams.find("all.pair.equivalence.class")->second, verbose, superverbose);
        }

        speciesChrNb.clear();
        adjacencies.clear();
        
        //LeafToGFMap.clear();
        LeafToSpMap.clear();



    	if(verbose)
    	{
            int totalAdj = 0;
    		for(unsigned i = 0 ; i < ECFams->size(); i++)
            {
    			ECFams->at(i).printMe(superverbose);
                totalAdj += ECFams->at(i).getNbAdj();
            }
            cout << "total number of adjs in ECFs : "<<totalAdj << endl;
    	}


        if(gBoolParams.find("load.save")->second) //we load a previous DeCo instance adjacency trees
        {
            string filename = LoadPrefix + ".adjacencyTrees.xml";
            bool ok = readECFamTreesOneFile(filename, ECFams, gBoolParams.find("always.AGain")->second, superverbose);

            if(!ok)
            {
                cerr << "Trying separate files reading" << endl;

            
                for(unsigned i = 0 ; i < ECFams->size(); i++)
                {
                    
                    //EquivalenceClassFamily * ECF;
                    //*ECF = ECFams->at(i);
    
                    LoadECFamTrees( &(ECFams->at(i)) , LoadPrefix, gBoolParams.find("always.AGain")->second, superverbose);
                }
            }

            if(verbose)
                cout << "adjacency trees loaded" << endl;
        }
        else
        {
    		////////////////////////////////////////////////////
    		/// Computing and backtracking Adjacency matrices
    		////////////////////////////////////////////////////



    		if(verbose)
    		{
    			begin = clock();
    			cout << "Computing Adj Matrices" << endl;
    		}

    		double gAdjWeight = gDoubleParams.find("Adjacency.weight")->second;
    		if(gAdjWeight == 0)
    			gAdjWeight = 0.000001; // in order to avoid division by 0
    

    		//computing matrices
    		ComputeEquivalenceClassFamilies(ECFams, GeneFamilyList , 
                                            adjacencyScores,speciesC0C1, speGeneAdjNb, // dummy arguments because of artdeco and ADseq
    													gDoubleParams.find("AGain.cost")->second, gDoubleParams.find("ABreak.cost")->second , 
    													false,//gBoolParams.find("use.boltzmann")->second, 
    													gBoolParams.find("substract.reco.to.adj")->second , 
    														gDoubleParams.find("dupli.cost")->second * gDoubleParams.find("Reconciliation.weight")->second / gAdjWeight, 
    														gDoubleParams.find("loss.cost")->second * gDoubleParams.find("Reconciliation.weight")->second / gAdjWeight, 
    														gDoubleParams.find("HGT.cost")->second * gDoubleParams.find("Reconciliation.weight")->second / gAdjWeight, 
    													verbose,superverbose, 
    													1,//gDoubleParams.find("boltzmann.temperature")->second,
                                                        gDoubleParams.find("absence.penalty")->second);


                                                


    		if(verbose)
    		{
    			end = clock();
    			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    			cout << "time elapsed for Adjacency Matrix computing:" << elapsed_secs << endl;
    
    			begin = clock();
    		}


    		//backtracking
    		vector <NECFsample * > * AllSamples;


   			backtrackInPlaceEquivalenceClassFamilies( ECFams, GeneFamilyList , 
   														false, //gBoolParams.find("use.boltzmann")->second, //should be false
   														verbose, superverbose,
   														gBoolParams.find("always.AGain")->second , gDoubleParams.find("C1.Advantage")->second,
                                                        false ); // <- last false to ignore number of backtrack in backtrack...


    		if(verbose)
    		{
    			end = clock();
    			elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    			cout << "time elapsed for Adjacency Matrix backtracking:" << elapsed_secs << endl;
    
                
    		}

        }


        //// Co-Event Building
        if(verbose)
        {
            cout <<  "reading Co-events" << endl;
            begin = clock();
        }

        vector <CoEvent> * CoEventSet = new vector <CoEvent>;
  
        for(unsigned i = 0 ; i < ECFams->size(); i++)
        {
            int gfam1 = ECFams->at(i).getGfamily1();
            int gfam2 = ECFams->at(i).getGfamily2();

            ReconciledTree * Rtree1 = GeneFamilyList->at(gfam1)->getRecTree();
            ReconciledTree * Rtree2 = GeneFamilyList->at(gfam2)->getRecTree();

            for(unsigned j = 0; j < ECFams->at(i).getNbEqClasses() ; j++)
            {
                //cout << "CoEventSet size " << CoEventSet->size()<< endl;
                PopulatesCoeventsFromAdjForest( ECFams->at(i).getAdjForest(j), CoEventSet, j ,  gfam1,  gfam2,  Rtree1,  Rtree2, gBoolParams.find("load.save")->second );
                //cout << "CoEventSet size " << CoEventSet->size()<< endl;
            }
        }

        /*0
        //July 2017  - testing some stuff ofn loading coevents
        map <int,int> coevCounts;
        map <int,int> coevCountsGenes;

        for(unsigned i = 0 ; i < CoEventSet->size(); i++)
        {
            int evt = CoEventSet->at(i).getEvent();
            if( coevCounts.find(evt) == coevCounts.end() )
            {
                coevCounts[evt] = 0;
                coevCountsGenes[evt]=0;
            }
            coevCounts[evt] += 1;
            coevCountsGenes[evt] += CoEventSet->at(i).getNumberOfGene() - CoEventSet->at(i).getNbConnexComponents();

            if( evt == 5 )
            {
                cout << "transferBack in sp: " << CoEventSet->at(i).getSpecies() << " nbG : " << CoEventSet->at(i).getNumberOfGene() << " nbC : " << CoEventSet->at(i).getNbConnexComponents() << endl;
            }
        }

        map<int,int>::iterator MyIt;

        for( MyIt = coevCounts.begin(); MyIt != coevCounts.end(); ++MyIt ) 
        {
            cout << "coevent of type" << MyIt->first << "nb : " << MyIt->second << " expected nb evt reimbursed : " << coevCountsGenes[MyIt->first] << endl;
        }
        */


        if(verbose)
        {
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

            cout << "time elapsed for Co-event building:" << elapsed_secs << endl;
            begin = clock();
        }

        //////////////////// COMPUTING THE SYSTEM SCORE
            double SystemScore = 0;
    
            double TopoScore = computeTopoScore( *GeneFamilyList, gDoubleParams.find("Topology.weight")->second);
    
            if(verbose)
                cout << "Topology part of the score: " <<  TopoScore << endl;

    
            double ReconScore = computeReconciliationScore( *GeneFamilyList, gDoubleParams.find("Reconciliation.weight")->second );
    
            if(verbose)
                cout << "Single reconciliation part of the score: " <<  ReconScore << endl;
    


    
            double AdjScore =  computeAdjacenciesScore( ECFams,gDoubleParams.find("AGain.cost")->second, gDoubleParams.find("ABreak.cost")->second , gDoubleParams.find("Adjacency.weight")->second );
    
    
            if(verbose)
                cout << "Adjacency part of the score: " <<  AdjScore << endl;
    
            double CoEventScore =  computeCoEventScore(*CoEventSet, gDoubleParams.find("dupli.cost")->second , 
            														gDoubleParams.find("loss.cost")->second ,
            														gDoubleParams.find("HGT.cost")->second ,
            														gDoubleParams.find("Reconciliation.weight")->second );
    


            if(verbose)
                cout << "CoEvent balancing of the score: " <<  CoEventScore << endl;
    
            SystemScore = TopoScore + ReconScore + AdjScore + CoEventScore ;
    
            if(verbose)
                cout << "Score of the system: " << TopoScore << " + " << ReconScore << " + " << AdjScore << " + " << CoEventScore << " = " << SystemScore  << endl;

        if(verbose)
        {
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

            cout << "time elapsed for score computing:" << elapsed_secs << endl;
            begin = clock();
        }


        //////////////////////////////
        /// Sample in the Gfam
        //////////////////////////////


        
        double temp = gDoubleParams.find("sampling.temp")->second ;

        int redrawAlgo = gIntParams["redraw.algo"];
        /*
        0 : only change reconciliation (not topologie)
        1 : use only the ccp
        2 : DLRecCoev
        */

        bool trySeveralECFSolutions = false;

        

        bool CanContinueSampling = gBoolParams.find("continue.sampling.while.accepted")->second;
        bool continueSampling = true;

        int tryLimit = gIntParams.find("sampling.nbMax")->second;
        if(tryLimit == 0)
            continueSampling = false;


        int totalAdj = 0;
        for(unsigned i = 0 ; i < ECFams->size(); i++)
        {
            //ECFams->at(i).printMe(superverbose);
            totalAdj += ECFams->at(i).getNbAdj();
        }
        cout << "total number of adjs in ECFs : "<<totalAdj << endl;


        double DupliCost = gDoubleParams.find("dupli.cost")->second;
        double HGTCost = gDoubleParams.find("HGT.cost")->second;
        double LossCost = gDoubleParams.find("loss.cost")->second;

        double Again  = gDoubleParams.find("AGain.cost")->second;
        double Abreak = gDoubleParams.find("ABreak.cost")->second;
        
        double TopoWeight = gDoubleParams.find("Topology.weight")->second;
        double RecWeight  = gDoubleParams.find("Reconciliation.weight")->second;
        double AdjWeight  = gDoubleParams.find("Adjacency.weight")->second;

        double probaFlat = gDoubleParams.find("proba.flat")->second;

        bool noCoev = gBoolParams.find("DTLRecCoev.NoCoev")->second;
        double DLRecCoevTemp = gDoubleParams.find("DTLRecCoev.Temp")->second;

        int optTurn = 0;
        while(continueSampling)
        {
            continueSampling = false;
            
            vector<int> randomIndexes;
            randomIndexes.reserve(GeneFamilyList->size());
            for(int i = 0; i < GeneFamilyList->size(); i++)
                randomIndexes.push_back(i);

            random_shuffle(randomIndexes.begin(),randomIndexes.end());

            //cout << "********************"<<endl;
            //for(int i = 0; i < randomIndexes.size(); i++)
            //    cout << " " << randomIndexes[i];
            //cout << endl;
            //cout << "********************"<<endl;
            for(int i = 0; i < randomIndexes.size(); i++)
            {
                int GfamIdToReset = randomIndexes[i];



                MyCladesAndTripartitions * BiasedCndT = NULL;
        
                if(gDoubleParams.find("proba.bias")->second > 0)
                {
                    if(verbose)
                        cout << "enter bias" << endl;
                    //// biasing the CandT of GF 0
                    int GFtoBias = GfamIdToReset;
        
                    vector < AdjTree * > * BiasingForest = new vector < AdjTree * >();
                    vector <bool> CompatList;
                    vector <int> GFindexes;
        
                    //1. detecting which ECFams bias which geneFamily
                    for(unsigned i = 0 ; i < ECFams->size();i++)
                    {
                        bool add = false;
                        int GfamIndex = 0;
                        if( ECFams->at(i).getGfamily1() == GFtoBias)
                        {
                            add = true;
                        }
                        else if( ECFams->at(i).getGfamily2() == GFtoBias )
                        {
                            GfamIndex = 1;
                            add = true;
                        }
        
                        if(add)
                        {
                            CompatList.push_back( !gBoolParams.find("load.save")->second ); 
                            GFindexes.push_back( GfamIndex );
        
                            for(unsigned EclassInd = 0 ; EclassInd <  ECFams->at(i).getNbEqClasses() ; EclassInd++)
                            {
                                vector <AdjTree * > * tmp =  ECFams->at(i).getAdjForest(EclassInd);
                                for(unsigned AtreeInd = 0 ; AtreeInd < tmp->size(); AtreeInd++)
                                {
                                    BiasingForest->push_back( tmp->at(AtreeInd) );
                                }
                            }
                        }
        
                    }
                    if(verbose)
                        cout << "gotten biasing forest: " << BiasingForest->size() << endl;
        
                    //2. biasing
                    BiasedCndT = GeneFamilyList->at(GFtoBias)->setBias( BiasingForest , GFindexes, CompatList);

                    if(verbose)
                        cout << "exit bias" << endl;
                }



                //bool accepted = updateGfam(GfamIdToReset, tryLimit ,redrawAlgo, trySeveralECFSolutions, temp,
                //                            GeneFamilyList, ECFams, &CoEventSet, speciesTree,
                //                            maxTS, gBoolParams.find("dated.species.tree")->second, gBoolParams.find("bounded.TS")->second, gBoolParams.find("with.transfer")->second, 
                //                            gDoubleParams.find("dupli.cost")->second , gDoubleParams.find("HGT.cost")->second, gDoubleParams.find("loss.cost")->second,
                //                            gDoubleParams.find("Topology.weight")->second, gDoubleParams.find("Reconciliation.weight")->second, gDoubleParams.find("Adjacency.weight")->second,
                //                            gDoubleParams.find("AGain.cost")->second, gDoubleParams.find("ABreak.cost")->second,
                //                            gDoubleParams.find("absence.penalty")->second, gBoolParams.find("all.pair.equivalence.class")->second, gBoolParams.find("substract.reco.to.adj")->second,
                //                            gBoolParams.find("always.AGain")->second , gDoubleParams.find("C1.Advantage")->second,
                //                            verbose, superverbose,
                //                            speciesC0C1,
                //                            speGeneAdjNb,
                //                            gDoubleParams.find("proba.bias")->second, BiasedCndT,
                //                            outputDir
                //                             );


                int nbtry = 0;
                bool accepted = false;
    
    
                

                //1. keeping the initial solution
                double normTreeLkhINIT = - log( GeneFamilyList->at(GfamIdToReset)->getNormalizedTreeLikelihood());
            
                GfamSave historic =  createHistoric(GfamIdToReset , 
                                                    GeneFamilyList, ECFams, &CoEventSet,
                                                    Again, Abreak,
                                                    DupliCost, HGTCost, LossCost,
                                                    RecWeight
                                                    );
                
                double historicScore = normTreeLkhINIT * TopoWeight + historic.RecScore * RecWeight + historic.AdjScore * AdjWeight + historic.CoEventScore ;
    
    
                

                while(nbtry < tryLimit)
                {
                    //accepted = tryGfamRedraw( GfamIdToReset, redrawAlgo, trySeveralECFSolutions, temp, historicScore, historic,
                    //                        GeneFamilyList, ECFams, &CoEventSet, speciesTree,
                    //                        maxTS, gBoolParams.find("dated.species.tree")->second, gBoolParams.find("bounded.TS")->second, gBoolParams.find("with.transfer")->second,
                    //                        DupliCost, HGTCost, LossCost,
                    //                        TopoWeight, RecWeight, AdjWeight,
                    //                        Again, Abreak,
                    //                        gDoubleParams.find("absence.penalty")->second, gBoolParams.find("all.pair.equivalence.class")->second, gBoolParams.find("substract.reco.to.adj")->second,
                    //                        gBoolParams.find("always.AGain")->second , gDoubleParams.find("C1.Advantage")->second,
                    //                        verbose, superverbose,
                    //                        speciesC0C1,
                    //                        speGeneAdjNb,
                    //                        gDoubleParams.find("proba.bias")->second, BiasedCndT,
                    //                        outputDir
                    //                         );
                    CoEventSet = new vector< CoEvent >;

                    ReDrawGeneFamily( GfamIdToReset, redrawAlgo, 
                                        GeneFamilyList, ECFams, CoEventSet, speciesTree,
                                        gBoolParams.find("dated.species.tree")->second, gBoolParams.find("bounded.TS")->second, gBoolParams.find("with.transfer")->second,
                                        DupliCost, HGTCost, LossCost, maxTS, TopoWeight,
                                        Again, Abreak, 
                                        gDoubleParams.find("absence.penalty")->second, gBoolParams.find("all.pair.equivalence.class")->second, gBoolParams.find("substract.reco.to.adj")->second,
                                        DupliCost * RecWeight / AdjWeight, LossCost * RecWeight / AdjWeight, HGTCost * RecWeight / AdjWeight,
                                        gBoolParams.find("always.AGain")->second , gDoubleParams.find("C1.Advantage")->second,
                                        verbose, superverbose ,
                                        trySeveralECFSolutions, temp, RecWeight, false,
                                        speciesC0C1,
                                        speGeneAdjNb, 
                                        probaFlat,
                                        gDoubleParams.find("proba.bias")->second, BiasedCndT,
                                        outputDir, nbtry, tryLimit,
                                        noCoev , DLRecCoevTemp,
                                        gStringParams.find("DTLRecCoevExecPath")->second
                                        );


                    //3. compute the new scores.
                    double topoScore = (- log( GeneFamilyList->at(GfamIdToReset)->getNormalizedTreeLikelihood())) * TopoWeight;
                    double recScore = RecWeight * GeneFamilyList->at(GfamIdToReset)->getRecScore();//only for the changed tree
                
                    double adjScore = 0;
                
                    map<int, EquivalenceClassFamily>::iterator it;
                    for(it = historic.ECFamsMap.begin() ; it != historic.ECFamsMap.end() ; ++it)
                    {
                        int ecId = it->first;
                        adjScore += ((double) ECFams->at(ecId).getNbAdjGain() * Again);
                        adjScore += ((double) ECFams->at(ecId).getNbAdjBreak() * Abreak);
                    }
                
                    adjScore *= AdjWeight;
                
                    double coEventScore = computeCoEventScore( *CoEventSet, DupliCost, LossCost, HGTCost, RecWeight);
                
                    if(verbose)
                        cout << "comparing scores for gfam " << GfamIdToReset << " : " << topoScore << " + " << recScore << " + " << adjScore << " + " <<  coEventScore << " = " << topoScore + recScore + adjScore + coEventScore << " <-> " << historicScore << endl; 
                
                    //4. accepting or no the new proposition
                    accepted = acceptProposition( historicScore ,
                                                    topoScore + recScore + adjScore + coEventScore,
                                                    temp);


                    if(accepted)
                        break;
                    else 
                    {
                        if(verbose)
                            cout << "not accepted" << endl;

                        CoEventSet->clear();
                        delete CoEventSet;
                    }
            
                    nbtry++;
                }
            
                if(verbose)
                {
                    if(accepted)
                        cout << "Accepted after try " << nbtry << endl;
                    else
                        cout << "Not accepted after try " << nbtry <<  endl;
                }
    
    
                
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
                    map<int, EquivalenceClassFamily>::iterator it;
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
            
                    //CoEventSet->clear();
                    //delete CoEventSet;
            
                    CoEventSet = historic.CoEventSet; // MEMWORK
            
            
                }
                else 
                {
                    
                    //cout << "rec nodes in historic" << historic.recTree->getNumberOfNodes() << endl;
                    historic.clearP();
                    delete historic.recTree;                  
            

                    historic.ECFamsMap.clear();
                    //cout << historic.CoEventSet->size() << endl;
                    //cout << "<<<<<" << endl;
            
                }



                ///////////////////////////////////////////////////////////

                int NEWtotalAdj = 0;
                for(unsigned j = 0 ; j < ECFams->size(); j++)
                {
                    //ECFams->at(i).printMe(superverbose);
                    NEWtotalAdj += ECFams->at(j).getNbAdj();
                }
                cout << "total number of adjs in ECFs : "<< NEWtotalAdj << endl;
                if( NEWtotalAdj != totalAdj )
                {
                    cout << "ERROR : wrong adj number!!" << endl;
                    exit(1);
                }




                if((accepted)&&(CanContinueSampling))
                    continueSampling = true;

                if(verbose)
                {
                    if(accepted)
                        cout << "accepted a solution for Gfam " << GfamIdToReset << endl;
                    else
                        cout << "did not accept a solution for Gfam " << GfamIdToReset << endl;

                }
            }

            if(continueSampling)
            {
                if(verbose)
                    cout << "reporting new solution of optimization turn " << optTurn << endl;
                // trees
                string tmpPrefix = "";
                if(outputDir !="")
                {
                    tmpPrefix += outputDir + "/";
                }
                tmpPrefix += "optTurn" + static_cast<ostringstream*>( &(ostringstream() << optTurn) )->str() ;
                if(!mkpath(tmpPrefix))
                {
                    cout << "ERROR: " << tmpPrefix << " is not a directory." << endl;
                    exit(1);
                }
                tmpPrefix += "/";

                if( gStringParams.find("output.prefix")->second != "none" )
                {
                    tmpPrefix += gStringParams.find("output.prefix")->second;
                }
                WriteGeneFamilyReconciledTrees(  GeneFamilyList, gBoolParams.find("write.newick")->second , gBoolParams.find("hide.losses.newick")->second,  tmpPrefix);
                for(unsigned i = 0 ; i < ECFams->size() ; i++)
                    WriteECFamTrees( &ECFams->at(i), gBoolParams.find("write.newick")->second, gBoolParams.find("hide.losses.newick")->second , tmpPrefix);

            }

            if(verbose)
            {
                if(continueSampling)
                {
                    cout << "A turn of optimization yielded change: continue optimization." << endl;
                }
                else if(CanContinueSampling)
                    cout << "A turn of optimization did not yield change: stop optimization." << endl;
            }
            optTurn++;
        }
        //5. score update
        SystemScore = 0;
    
        TopoScore = computeTopoScore( *GeneFamilyList, gDoubleParams.find("Topology.weight")->second);
        ReconScore = computeReconciliationScore( *GeneFamilyList, gDoubleParams.find("Reconciliation.weight")->second );
        AdjScore =  computeAdjacenciesScore( ECFams,gDoubleParams.find("AGain.cost")->second, gDoubleParams.find("ABreak.cost")->second , gDoubleParams.find("Adjacency.weight")->second );
        CoEventScore =  computeCoEventScore(*CoEventSet, gDoubleParams.find("dupli.cost")->second , 
                                                                    gDoubleParams.find("loss.cost")->second ,
                                                                    gDoubleParams.find("HGT.cost")->second ,
                                                                    gDoubleParams.find("Reconciliation.weight")->second );
    
        SystemScore = TopoScore + ReconScore + AdjScore + CoEventScore ;
    
        if(verbose)
            cout << TopoScore << " + " << ReconScore << " + " << AdjScore << " + " << CoEventScore << " = " << SystemScore  << endl;


        if(verbose)
        {
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

            cout << "time elapsed for sampling:" << elapsed_secs << "s."<< endl;// That's "<< elapsed_secs / nbtry << "s / try" << endl;
            begin = clock();
        }


        ////////////////////////////////////////////////////
        /// Reporting
        ////////////////////////////////////////////////////

        // trees
        WriteSpeciestree( speciesTree,  gBoolParams.find("write.newick")->second , gPathPrefix);



        //write all rec trees in one file
        string recFileName = gPathPrefix;
        string adjTreeFileName = gPathPrefix;
        
        if(( recFileName[ recFileName.size() - 1 ] != '/')&&(recFileName.size() > 0))
        {
            recFileName += ".";
            adjTreeFileName += ".";
        }

        recFileName += "reconciliations";
        adjTreeFileName += "adjacencyTrees";

        if( gBoolParams.find("write.newick")->second ) 
        {
            recFileName += ".newick" ;
            adjTreeFileName += ".newick" ;
        }
        else
        {
            recFileName += ".xml";
            adjTreeFileName += ".xml";
        }

        InitReconciledTreeFile( recFileName, gBoolParams.find("write.newick")->second );

        for(unsigned i = 0 ; i < GeneFamilyList->size() ; i++)
            AddOneGeneFamilyReconciledTreeToFile( GeneFamilyList->at(i), gBoolParams.find("write.newick")->second , gBoolParams.find("hide.losses.newick")->second, recFileName, i);

        FinishReconciledTreeFile(recFileName, gBoolParams.find("write.newick")->second);

        InitAdjacencyTreeFile( adjTreeFileName, gBoolParams.find("write.newick")->second );

        for(unsigned i = 0 ; i < ECFams->size() ; i++)
        {
            //WriteECFamTrees( &ECFams->at(i), gBoolParams.find("write.newick")->second, gBoolParams.find("hide.losses.newick")->second , gPathPrefix);
            AddECFToFile(adjTreeFileName, &ECFams->at(i),
                             gBoolParams.find("write.newick")->second, gBoolParams.find("hide.losses.newick")->second, 
                            gDoubleParams.find("AGain.cost")->second,
                            gDoubleParams.find("ABreak.cost")->second);

        }


        if(gBoolParams.find("write.adjacencies")->second)
            WriteAdjacencies( ECFams, gPathPrefix);

        if(gBoolParams.find("write.adjacencies")->second)// writng co-events.
        {
            DeCoOutputManager DOM;

            string filename = gPathPrefix;
            if( filename[ filename.size() - 1 ] != '/')
                filename += ".";
            filename += "coevents.txt";

            ofstream ofs;
            ofs.open(filename.c_str(),ofstream::out);

            DOM.WriteCoEventSet(ofs, CoEventSet);
            
            ofs.close();
        }

        if(verbose)
        {
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

            cout << "time elapsed for writing:" << elapsed_secs << endl;
            begin = clock();
        }


        ECFams->clear();
        delete ECFams;
        delete speciesTree;

        for(unsigned i = 0 ; i < GeneFamilyList->size();i++)
            delete GeneFamilyList->at(i);

        delete GeneFamilyList;


    }
	catch(exception & e)
	{
		cerr << e.what() << endl;
		exit(-1);
	}
	
	return 0;


}
