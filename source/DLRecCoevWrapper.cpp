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

This file contains a wrapper for the DLRecCoev algorithm

Created the: 16-01-2017
by: Wandrille Duchemin

Last modified the: 16-01-2017
by: Wandrille Duchemin

*/

#include "DLRecCoevWrapper.h"

// change to boost::filesystem?
static bool doMkdir( const char *path )
{
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
int mkPath( string path ) // mode_t mode )
{
    char *copypath = strdup(path.c_str());
    char *sp;
    int status = 0;
    char *pp = copypath;
    while( status == 0 && (sp = strchr(pp, '/')) != 0 ) {
        if (sp != pp) { 
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = doMkdir(copypath ); //, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    //if (status == 0)
        status = doMkdir(path.c_str()); //, mode);
    free(copypath);
    return status;
}


int DLRecCoevWrapper::MakePath()
{
    return mkPath( DirPath );
}

void DLRecCoevWrapper::writeCoevFile(GeneFamily * GuideF)
{

    map<int, vector<int> > guideEvents = GuideF->countEventTypePerSpecies();
    //map associating species RPO to the count of: Duplications; Losses; HGT
    string fileName = DirPath + "coev";

    //cout << fileName << endl;

    //for( map<int, vector<int> >::iterator it = guideEvents.begin(); it != guideEvents.end(); ++it )
    //{
    //    cout << it->first << " -> D:"<< it->second[0] << " L:" << it->second[1] << " T:" << it->second[2] << endl;
    //}

    ofstream ofs;
    ofs.open(fileName.c_str(),ofstream::out );


    map<int, vector<int> >::iterator it;

    vector <int>  Snodes = spTree->getNodesId();
    for(unsigned j = 0 ; j < Snodes.size(); j++)
    {
        int RPO = spTree->getRPO(Snodes[j]);
        it = guideEvents.find( RPO );
        if( it != guideEvents.end() )
        {
            vector< string > names = spTree->cloneSubtree(Snodes[j])->getLeavesNames(); // list of leaves under this node

            for(unsigned i = 0 ; i < 3 ; i ++)
            {
                for(unsigned k = 0 ; k < it->second[i] ; k++)
                {
                    ///print evt type + names
                    //if(i==0)
                    //    cout << "D";
                    //else if(i == 1)
                    //    cout << "L";
                    //else if(i == 2)
                    //    cout << "T";
                    //
                    //for(unsigned l = 0 ; l< names.size(); l++)
                    //    cout << " " << names[l];
                    //cout << endl;

                    if(i==0)
                        ofs << "D";
                    else if(i == 1)
                        ofs << "L";
                    else if(i == 2)
                        ofs << "T";
                    
                    for(unsigned l = 0 ; l< names.size(); l++)
                        ofs << " " << names[l];
                    ofs << endl;
                }

            }
        }
    }

    //GuideF->printRecTree();

}


void DLRecCoevWrapper::writeCoevFile()
{
    string fileName = DirPath + "coev";

    ofstream ofs;
    ofs.open(fileName.c_str(),ofstream::out );
    ofs.close();
}

void DLRecCoevWrapper::writeSpeciesFile()
{
    WriteSpeciestree( spTree, true , DirPath);
}

void DLRecCoevWrapper::setCommand()
{
    command += "time python ~/Documents/fooling_around/RecCoev/DTLCoevBoltzmann/";
    if(!Transfer)
        command += "DLCoevBoltzmann.py";
    else
        command += "DTLCoevBoltzmann.py";

    command += " -s " + DirPath + "speciesTree.newick";
    command += " -g " + GFamToRedraw->getSourceFile();

    if(!NoCoev)
        command += " -G " + DirPath + "coev";

    command += " -n " + static_cast<ostringstream*>( &(ostringstream() << nbSample) )->str();
    
    command += " -o " + outputFile;
    
    command += " --temperature " + static_cast<ostringstream*>( &(ostringstream() << temp) )->str();
    command += " -T " +  static_cast<ostringstream*>( &(ostringstream() << RelativeTopoWeight) )->str();
    command += " -d " +  static_cast<ostringstream*>( &(ostringstream() << DupCost) )->str();
    command += " -l " +  static_cast<ostringstream*>( &(ostringstream() << LossCost) )->str();
    command += " --de " +  static_cast<ostringstream*>( &(ostringstream() << ExtensionCost) )->str();
    command += " --le " +  static_cast<ostringstream*>( &(ostringstream() << ExtensionCost) )->str();

    //transfer specific commands
    if(Transfer)
    {
        if(!dated)
            command += " --undated";

        command += " -t " + static_cast<ostringstream*>( &(ostringstream() << HGTCost) )->str();
        command += " --te " +  static_cast<ostringstream*>( &(ostringstream() << ExtensionCost) )->str();
        command += " --inconsistent.scenarios.when.failure";
        command += " --ignore.coevents.consistency"; // as we don't read co-events from the reconciliation.
    }

    stringstream ss;
    string sep;
    ss << GFamToRedraw->getCharSep() ;
    ss >> sep;
    command += " --sep " + sep ;

    command += " --xml ";
    
    //cout << command <<endl;
}

void DLRecCoevWrapper::launchCommand()
{
    int returnValue = system( command.c_str() ) ;

    if( returnValue != 0)
    {
        cerr << command << endl;
        cerr << "returned " << returnValue << endl;
    }
}

void DLRecCoevWrapper::readTree(int index )
{
    cout << "before:" <<endl;
    GFamToRedraw->printRecTree();
    


    ifstream fileStream(outputFile.c_str());
    if( !fileStream.is_open() ) 
    {
        throw Exception("DLRecCoevWrapper::readTree : could not open reconciled tree file : " + outputFile);
        exit(1);
    }
    string RecGeneTreeTag = "recGeneTree";
    string startTag = "clade";
    int current = 0;
    while( goToNextTag(fileStream,RecGeneTreeTag) )  //goes to the next RecGeneTreeTag -> next reconciled gene tree
    {
        if(current == index)
        {
            if( goToNextTag(fileStream,startTag) ) // go to clade -> the root of the reconciled gene tree
            {
                ReconciledTree Rtree(fileStream, spTree, false);
                GFamToRedraw->setReconciliation(Rtree);
                break;
            }
        }
        current++;
    }


    if(current != index)
    {
        cerr << "!!ERROR!! did not find the recGeneTree with desired index (" << index << ")" << endl;
        exit(1);
    }

    GFamToRedraw->setRecScore(DupCost, HGTCost, LossCost);

    fileStream.close();

    cout << "after:" <<endl;
    GFamToRedraw->printRecTree();

}