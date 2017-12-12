#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        22-Feb-2017         #
##  Last modified:  22-Feb-2017         #
#########################################

import ete3
from ete3 import Tree, TreeNode

import sys

if __name__ == "__main__":

    help =  """
                Given a tree with one or more polytomies, this script returns the
                list of all trees (in newick format) resulting from the combination of
                all possible solutions of the multifurcated nodes.

                usage : python expandPolytomies.py -i fileIn [-o fileOut -n 5]
                            -i fileIn   : name of the input file with a multifurcated tree
                            -o fileOut  : (optional) name of the output file for solved trees (default is fileIn + ".expanded" )
                            -n 5        : (optional) maximum polytomy size (default is 5, maximum is 9)

    
                Please note that the number of of possible binary trees grows
                exponentially with the number and size of polytomies. Using this
                script with large multifurcations is not feasible:
             
                polytomy size: 3 number of binary trees: 3
                polytomy size: 4 number of binary trees: 15
                polytomy size: 5 number of binary trees: 105
                polytomy size: 6 number of binary trees: 945
                polytomy size: 7 number of binary trees: 10395
                polytomy size: 8 number of binary trees: 135135
                polytomy size: 9 number of binary trees: 2027025
               """


    OK = True

    nextKEY = None
    params = {
                            "-i" : None, #fileIn   : name of the input file with a multifurcated tree
                            "-o" : None, #fileOut  : (optional) name of the output file for solved trees (default is fileIn + ".expanded" )
                            "-n" : 5   ,  #   : (optional) maximum polytomy size (default is 5)
                            "--force" : False
    }

    flagArgs = ["--force"]

    for i in range(1,len(sys.argv)):

        if not nextKEY is None:
            params[nextKEY] = sys.argv[i]
            print "argument ",nextKEY,":", sys.argv[i]
            nextKEY = None
            continue

        if sys.argv[i] in params.keys():

            if sys.argv[i] in flagArgs:
                params[sys.argv[i]] = True
                print sys.argv[i],"flag activated"
            else:
                nextKEY = sys.argv[i]
            continue
        else:
            print "unknown argument", sys.argv[i]

    if params["-i"] is None:
        OK = False
        print "error: input file not given."

    if OK:
        
        ## treating positive float options
        for pname in []:
            try:
                params[pname] = float(params[pname])
                if params[pname] < 0:
                    print "error: ",pname ,"must be a positive number."
                    OK = False
            except:
                print "error:",pname,"must be a positive number."
                OK = False

        ## treating positive int options
        for pname in ["-n"]:
            try:
                params[pname] = int(params[pname])
                if params[pname] < 1:
                    print "error: ",pname ,"must be a positive integer."
                    OK = False
            except:
                print "error:",pname,"must be a positive number."
                OK = False


    if params["-n"] > 9:
        print "-n parameter set back to its maximum value : 9."
        params["-n"] = 9

    if not OK:
        print help
        exit(1)



    defaultOutputSuffix = ".expanded"
    if params["-o"] is None:
        params["-o"] = params["-i"] + defaultOutputSuffix


    print "reading input tree."


    dPsNbT = {3 : 3,
        4 : 15,
        5 : 105,
        6 : 945,
        7 : 10395,
        8 : 135135,
        9 : 2027025}

    tree = Tree(params["-i"])
    nbE = 1
    for n in tree.traverse():
        nbC = len(n.children)
        if nbC > 3:
            nbE *= dPsNbT[nbC]

    print "expanding will create", nbE, "topologies"
    if nbE > 1000000 : 
        print "More than 100000 topologies. Aborting to avoid potential memory hog."
        print "Use --force option to override."
        exit(1)

    print "expanding poytomies"

    try:
        TL = tree.expand_polytomies(polytomy_size_limit= params["-n"] )
    except ete3.coretype.tree.TreeError as e:
        se = str(e)
        se = se.replace("\\n","\n")
        print se
        #print e
        exit(1)

    print "generated",len(TL),"bifurcating trees"

    OUT = open(params["-o"],"w")
    for l in TL:
        OUT.write(l + "\n")
    OUT.close()