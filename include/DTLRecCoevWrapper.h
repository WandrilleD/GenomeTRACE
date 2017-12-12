#ifndef DLRECCOEV_WRAPPER_H_
#define DLRECCOEV_WRAPPER_H_

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


#include <sys/stat.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>

#include "GeneFamily.h"
#include "DeCoUtils.h"

using namespace std;

static bool doMkdir( const char *path );
int mkPath( string path ); 

class DTLRecCoevWrapper
{
protected:

    string execPath;

    string DirPath;
    string command;
    string outputFile;

    MySpeciesTree * spTree;
    GeneFamily * GFamToRedraw;

    double temp;
    double RelativeTopoWeight;
    double DupCost;
    double LossCost;
    double HGTCost;
    double ExtensionCost;

    int nbSample;

    bool Transfer;
    bool dated;

    bool NoCoev;

public:
    DTLRecCoevWrapper()
    {DTLRecCoevWrapperAux();}

    DTLRecCoevWrapper(string path)
    {
        DirPath = path;
        outputFile = DirPath + "recOutput.xml";
        DTLRecCoevWrapperAux();
    }

    DTLRecCoevWrapper(string path, string Epath)
    {
        DTLRecCoevWrapperAux();
        DirPath = path;
        outputFile = DirPath + "recOutput.xml";
        execPath = Epath;

    }


    void DTLRecCoevWrapperAux()
    {
        execPath = "";
        temp = 1;
        RelativeTopoWeight = 1;
        DupCost = 2;
        HGTCost = 3;
        LossCost = 1;
        ExtensionCost = 0;
        nbSample = 1;
        Transfer = false;
        dated = false;
        NoCoev=false;
    }

    void setTemperature(double v){temp = v;};
    void setRelativeTopoWeight(double v){RelativeTopoWeight = v;}
    void setDupCost(double v){DupCost = v;}
    void setLossCost(double v){LossCost = v;}
    void setHGTCost(double v){HGTCost = v;}

    void setnbSample(int s){nbSample = s;}
    void setTransfer(bool t){Transfer = t;}
    void setDated(bool d){dated = d;}

    void setNoCoev(){NoCoev = !NoCoev;}

    void setSpTree(MySpeciesTree * t){spTree = t;}

    void setGFamToRedraw(GeneFamily * GFam){GFamToRedraw = GFam;}

    void writeCoevFile(GeneFamily * GuideF);
    void writeCoevFile(); // overload when no guide F

    void writeSpeciesFile();

    ~DTLRecCoevWrapper(){}

    int MakePath();

    void setCommand();

    void launchCommand();

    void readTree(int index = 0);
};


#endif