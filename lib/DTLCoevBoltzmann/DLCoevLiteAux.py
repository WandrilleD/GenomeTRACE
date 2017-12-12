#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        13-Jan-2017         #
##  Last modified:  22-Feb-2017         #
#########################################


from CCP_distribution_ETEbased import CCP_distribution
from ETE3_based_ReconciledTree import RecEvent, ReconciledTree


from ete3 import Tree
import numpy as np
from copy import deepcopy


def readSpTree(fileName , fullNames = False):
    """ 
    reads a file containing a species tree and returns the tree
    and structures for the quick association between the nodes of the tree and their post-order.
    If fullNames is True, all names present in the tree will be added to SpNameToPO (not only leaves).
    """
    format = 0
    if fullNames:
        format = 1
    tree = Tree(fileName, format = format)
    SpNameToPO = {}
    POToNode = []
    i = 0
    for node in tree.traverse("postorder"):
        if node.is_leaf() or fullNames:
            SpNameToPO[node.name] = i
        node.name = i
        POToNode.append(node)
        i += 1
    #print tree
    #print SpNameToPO
    return tree, SpNameToPO, POToNode


def setupPOandTSdict(spTree):
    """
    Takes:
        - spTree (ete3.Tree) : ultrametric species tree that has been subdivided 
                                and where nodes posses a timeSlice attribute

    Returns:
        (tuple):
                (ete3.Tree) : trees where nodes have an additionnal PO feature
                (dict) : SpNameToPO -> associates species name to postOrder
                (list) : POToNode   -> associates postOrder to Node object
                (list) : TStoPO   -> associates timeSlice to postOrder / or None if the tree is not dated
    """
    SpNameToPO = {}
    POToNode = []


    TStoPO = None
    maxTS = None

    if "timeSlice" in spTree.features:
        maxTS = spTree.timeSlice
        TStoPO = [ [] for i in range(maxTS + 1)]

    i = 0
    for node in spTree.traverse("postorder"):
        if node.name != "":
            SpNameToPO[node.name] = i

        node.add_feature("PO", i)
        POToNode.append(node)

        if not TStoPO is None:
            TStoPO[node.timeSlice].append(i)

        i += 1

    return spTree , SpNameToPO , POToNode , TStoPO


def trimSpTree(spTree , spList):
    """
    * changes spTree *
    Takes
        - spTree (ete3.TreeNode) : the species tree (without dead lineage)
        - spList (list) : list of species name

    Returns:
            (ete3.TreeNode) : trimmed species tree so that the root species is the LCA of the genes

    """
    try:
        t = spTree.get_common_ancestor(spList)
    except ValueError as e:
        print "!!ERROR!! while trimming the species tree."
        print "Some genes are not associated to a valid species.",
        se = str(e)
        b,j,l = se.partition("[")
        l = "[" + l
        l = eval(l)
        print len(l) , "invalid species (" ,l[0],
        i = 1
        while i < len(l) and i < 3:
            print ",",l[i],
            i += 1
        if len(l) > i:
            print  ",... )"
        else:
            print ")"

        exit(1)
    t.detach()
    return t

def NodeToPreTrimPO(spTree):
    dNodeToPreTrimPO = {}
    i = 0
    for node in spTree.traverse("postorder"):
        dNodeToPreTrimPO[node] = i
        i += 1
    return dNodeToPreTrimPO

def POLeafList(spTree):
    dPreTrimPOLeafList = {}
    dNodeToPO = {}
    i = 0
    for node in spTree.traverse("postorder"):
        dNodeToPO[node] = i
        dPreTrimPOLeafList[i] = set()
        if node.is_leaf():
            dPreTrimPOLeafList[i].add(node.name)
        else:
            for c in node.children:
                dPreTrimPOLeafList[i]  = dPreTrimPOLeafList[i].union( dPreTrimPOLeafList[ dNodeToPO[c] ] )

        i += 1
    return dPreTrimPOLeafList



class RecCaseCoEv:
    def __init__(self):
        """
        self.minScore (float): minimal score over all tag
        self.minScores (dict): keys are tag, values are the minimal score associated with this tag
        self.solutions (dict): keys are tag, values are lists of solution with this tag and the minimal score for this tag
        """
        self.minScore = np.inf
        self.minScores = {}
        
        self.solutions = {}

    def getTags(self):
        return self.minScores.keys()

    def getMinScore(self, tag=None):
        if tag is None:
            return self.minScore
        return self.minScores[tag]

    def addSolution(self, solution):
        """ presuming solution is a RecSolutionCoEv instance 
            Only keeps optimal solutions for each possible tag
        """
        key = solution.getHashableConstraint()

        if ( not self.minScores.has_key(key) ) or self.minScores[key] > solution.score:
            self.minScores[key] = solution.score
            self.solutions[key] = [ solution ]
            if self.minScore > solution.score:
                self.minScore = solution.score
        elif self.minScores[key] == solution.score:
            self.solutions[key].append( solution )
        

    def getAllOptSolutions(self):
        opt = []
        for k in self.solutions.keys():
            if self.minScores[k] == self.minScore:
                opt += self.solutions[k]
        return opt


    def chooseOneOptSolution(self, tag = None):
        """ 
            randomly choose a solution with the best score
            if tag is not None, only choose solutions with this tag
        """

        opt = []
        if tag is None:
            opt = self.getAllOptSolutions()
        else:
            opt = self.solutions[tag]

        #print opt, self.minScores, self.minScore
        r = np.random.randint(0, len(opt))
        #print r
        return opt[ r ]


    def getScore(self,tag = None):
        return self.getMinScore(tag = None)

class RecSolutionCoEv:
    def __init__(self, score = np.inf, component = [] , event = None, tag = 0, coev = False):
        """
        represent a solution for reconciliation, with a score, some components (cases it points to) and event
        but also the tag (which is linked to the number of duplication in this history) and wether or not a coevent is formed
        """
        self.score = score
        self.component = component
        self.event = event
        self.tag = tag
        self.coev = coev
        

    def getNbEvent(self):
        return self.tag

    def addEvent(self):
        self.tag += 1

    def __str__(self):
        return str(self.event) + "  " + str(self.score) + " " + str(self.component) + " " + str(self.tag)  + " " + str(self.coev)

    def getHashableConstraint(self):
        return self.tag

    def copy(self):
        return RecSolutionCoEv(self.score, self.component , self.event, self.tag , self.coev )

class RecMatrix:
    def __init__(self , caseClass = RecCaseCoEv):
        self.caseClass = caseClass
        self.mat = []

    def setupIndexes(self, sortedBip , sortedSp):

        self.dbipId = {}
        self.dSpId = {}

        for i,bip in enumerate(sortedBip):
            self.dbipId[bip] = i

            self.mat.append( [] )

            for j,sp in enumerate(sortedSp):
                self.dSpId[sp] = j

                self.mat[i].append( self.caseClass() )

    def __getitem__(self, bipId ):
        return self.mat[ self.dbipId[bipId] ]

    def has_case(self, bipId, spId ):
        return True
        #if self.mat.has_key(bipId):
        #    if self.mat[bipId].has_key(spId):
        #        return True
        #return False

    def set_case(self, bipId, spId, solution ):
        """ add a solution to a given case """
        #if not self.mat.has_key(bipId):
        #    self.mat[bipId] = {}
        #if not self.mat[bipId].has_key(spId):
        #    self.mat[bipId][spId] = self.caseClass()

        i = self.dbipId[bipId]
        j = self.dSpId[spId]

        self.mat[i][j].addSolution(solution)
        

    def add_case(self, bipId, spId, solution ):
        """ add the solution if the case is empty OR the score of the solution is not infinity """

        i = self.dbipId[bipId]
        j = self.dSpId[spId]

        if solution.score == 0.:
            return
        else:
            self.mat[i][j].addSolution(solution)

        return

    def get_case_tags(self, bipId, spId ):
        i = self.dbipId[bipId]
        j = self.dSpId[spId]

        return self.mat[i][j].getTags()


    def get_case_score(self, bipId, spId , tag = None):
        #print "get_case_score" , bipId, spId , tag, "->" ,self.has_case(bipId, spId)
        i = self.dbipId[bipId]
        j = self.dSpId[spId]

        return self.mat[i][j].getScore(tag)


class RecMatrixOld:
    def __init__(self , caseClass = RecCaseCoEv):
        self.caseClass = caseClass
        self.mat = {}

    def __getitem__(self, bipId ):
        return self.mat[bipId]

    def has_case(self, bipId, spId ):
        if self.mat.has_key(bipId):
            if self.mat[bipId].has_key(spId):
                return True
        return False

    def set_case(self, bipId, spId, solution ):
        """ add a solution to a given case """
        if not self.mat.has_key(bipId):
            self.mat[bipId] = {}
        if not self.mat[bipId].has_key(spId):
            self.mat[bipId][spId] = self.caseClass()

        self.mat[bipId][spId].addSolution(solution)
        

    def add_case(self, bipId, spId, solution ):
        """ add the solution if the case is empty OR the score of the solution is not infinity """
        if not self.has_case(bipId,spId):
            self.set_case(bipId, spId, solution)
            return
        else:
            if solution.score == np.inf:
                return
            else:
                self.mat[bipId][spId].addSolution(solution)
        return

    def get_case_tags(self, bipId, spId ):
        return self.mat[bipId][spId].getTags()


    def get_case_score(self, bipId, spId , tag = None):
        #print "get_case_score" , bipId, spId , tag, "->" ,self.has_case(bipId, spId)
        return self.mat[bipId][spId].getScore(tag)

