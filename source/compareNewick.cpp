
#include <iostream>
#include <boost/foreach.hpp>
#include "MyGeneTree.h"

// This can be replaced by MyTree::readTree
MyGeneTree* readTree( char *fileName ) {

    ifstream treePath( fileName, ios::in);  

    if(! treePath) 
        throw bpp::IOException ("Newick::read: failed to read from stream"); 
    
    if( treePath.eof() ) {
        cout << "empty file " << fileName << endl;
        exit(1);
    }

    string temp;
	getline(treePath, temp, '\n'); 


    size_t index = temp.find(";");
    if(index == string::npos) {
        cout << "no semi-colon in " << fileName << endl;
        exit(1);
    }

    string errString = "";
    MyGeneTree *geneTree = new MyGeneTree( temp.substr(0, index + 1), 
                                           errString );
    if( errString != "" ) {
        cout << "Error reading tree: " << errString << endl;
        exit(1);
    }


    return geneTree;
}

// binary trees only
string orderedNewick( MyGeneNode *node, string &least ) {
    string newick;
    if( node->isLeaf() ) {
        newick = least = node->getName();
    } else {
	    if( node->getNumberOfSons() != 2 ) 
            throw bpp::Exception("Non-binary tree");

        string sonLeast1, sonLeast2;
        string newick1 = orderedNewick(node->getSon(0), sonLeast1);
        string newick2 = orderedNewick(node->getSon(1), sonLeast2);
        if( sonLeast1.compare( sonLeast2 ) <= 0 ) {
            least = sonLeast1;
            newick = "(" + newick1 + "," + newick2 + ")";
        } else {
            least = sonLeast2;
            newick = "(" + newick2 + "," + newick1 + ")";
        }
    }
    return newick;

}

int main(int args, char ** argv) {
    if( args != 3 ) {
        cout << "Need two trees" << endl;
        exit(1);
    }
    MyGeneTree *tree1 = readTree( argv[1] );
    MyGeneTree *tree2 = readTree( argv[2] );

    cout << "node count: " << tree1->getNumberOfNodes() << endl;
    if( tree1->getNumberOfNodes() != tree2->getNumberOfNodes() ) {
        cout << "unequal node count: " << tree2->getNumberOfNodes() << endl;
        exit(1);
    } 

    string x,y;
    string newick1 = orderedNewick( tree1->getRootNode(), x );
    cout << newick1 << endl;
    string newick2 = orderedNewick( tree2->getRootNode(), y );

    cout << endl;
    cout << newick2 << endl;

    if( newick1.compare( newick2 ) == 0 ) {
        cout << "EQUAL" << endl;
    } else {
        cout << "NOT EQUAL" << endl;
    }

    delete tree1;
    delete tree2;
}
