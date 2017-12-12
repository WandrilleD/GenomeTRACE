#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        02-Feb-2017         #
##  Last modified:  07-Feb-2017         #
#########################################


from ete3 import Tree, TreeNode
from copy import deepcopy

def getDistFromRootDic(spTree , checkUltrametric = True):
    """
    Takes:
        - spTree (ete3.Tree) : an ULTRAMETRIC species tree
        - checkUltrametric (bool) [default = True] : check if the tree is ultrametric and returns None

    Returns:
        (dict) : a dictionnary where keys are nodes and values are heights
        or
        (None) : if the tree isn't ultrametic and checkUltrametric is True
    """

    Dheight = {}
    

    for n in spTree.traverse():
        if n.is_root():
            Dheight[n] = n.dist
        else:
            Dheight[n] = n.dist + Dheight[ n.up ]

    if checkUltrametric:
        leaves = spTree.get_leaves()
        refD = Dheight[ leaves[0] ]
        #print [(l.name , Dheight[l] ) for l in leaves]
        for l in leaves[1:]:
            if Dheight[ l ] != refD:
                return None
    return Dheight




def subdivideSpTree(spTree):
    """
    Takes:
        - spTree (ete3.Tree) : an ULTRAMETRIC species tree

    Returns:
        (ete3.Tree) : subdivided species tree where all nodes have a timeSlice feature
        or
        None if the species tree is not ultrametric
    """
    newSpTree = deepcopy(spTree)

    featureName = "timeSlice"

    ##1/ getting distance from root.
    Dheight = getDistFromRootDic(newSpTree , checkUltrametric = True)

    if Dheight is None:
        print "!!ERROR!! : the species tree is not ultrametric"
        return None

    # we know that there is n-1 internal nodes (where n is the number of leaves)
    # hence the maximal timeSlice is n-1 (all leaves have timeSlice 0)

    ##2/assign timeSlice to nodes
    currentTS = len(newSpTree.get_leaves()) - 1


    for n,h in sorted(Dheight.iteritems(), key=lambda (k,v): (v,k)):
        n.add_feature(featureName, currentTS )

        if currentTS != 0:
            currentTS -= 1


    #print newSpTree.get_ascii(attributes=[featureName,"name"])

    ##3/subdivide according to timeSlice
    RealNodes = [i for i in  newSpTree.traverse()]

    for n in RealNodes:
        if n.is_root():
            continue

        nodeToAdd = n.up.timeSlice - n.timeSlice - 1

        while nodeToAdd > 0:
            parentNode = n.up
            
            n.detach()
            
            NullNode = TreeNode()
            NullNode.add_feature( featureName, parentNode.timeSlice - 1 )

            if "dead" in n.features:
                NullNode.add_feature("dead" , n.dead)

            parentNode.add_child(NullNode)
            NullNode.add_child(n)
            nodeToAdd -= 1 

    #print newSpTree.get_ascii(attributes=[featureName,"name"])
    return newSpTree

#def setupPOandTSdict(spTree):
#    """
#    Takes:
#        - spTree (ete3.Tree) : ultrametric species tree that has been subdivided 
#                                and where nodes posses a timeSlice attribute
#
#    Returns:
#        (tuple):
#                (ete3.Tree) : trees where nodes have an additionnal PO feature
#                (dict) : SpNameToPO -> associates species name to postOrder
#                (list) : POToNode   -> associates postOrder to Node object
#                (list) : TStoPO   -> associates timeSlice to postOrder / or None if the tree is not dated
#    """
#    SpNameToPO = {}
#    POToNode = []
#
#
#    TStoPO = None
#    maxTS = None
#
#    if "timeSlice" in spTree.features:
#        maxTS = spTree.timeSlice
#        TStoPO = [ [] for i in range(maxTS + 1)]
#
#    i = 0
#    for node in spTree.traverse("postorder"):
#        if node.name != "":
#            SpNameToPO[node.name] = i
#
#        node.add_feature("PO", i)
#        POToNode.append(node)
#
#        if not TStoPO is None:
#            TStoPO[node.timeSlice].append(i)
#
#        i += 1
#
#    return spTree , SpNameToPO , POToNode , TStoPO
#


def addDeadLineage(spTree):
    """ 
    Takes:
        - spTree (ete3.Tree) : species tree

    Returns:
        (ete3.Tree) : same tree with a dead lineage (name "-1") as outgroup
                      AND all nodes have a "dead" feature (bool that is True only for the dead lineage and the new root)
    """

    newSpTree = deepcopy(spTree)
    newSpTree.dist = 0.1

    for n in newSpTree.traverse():
        n.add_feature("dead",False)

    newRoot = TreeNode()
    newRoot.add_feature("dead",True)
    newRoot.dist = 0.0


    newRoot.add_child(newSpTree)
    rootHeight = newRoot.get_distance(newRoot.get_leaves()[0])

    deadLineage = TreeNode()
    deadLineage.add_feature("dead",True)    
    deadLineage.name = "-1"
    deadLineage.dist  = rootHeight

    newRoot.add_child(deadLineage)

    return newRoot


def setupRealNodePO(spTree):
    """
    setup a RealNodePO features on all nodes to has the PO of their closest descendant that is either a leaf or bifurcating
    """

    currentRPO = None
    for n in spTree.traverse("postorder"):
        if n.is_leaf() or len(n.children) > 1:
            currentRPO = n.PO

        n.add_feature("RealNodePO",currentRPO)

    return spTree

def getPOtoDead(POToNode):
    return [n.dead for n in POToNode]


def setupPossibleReceiverUndated(spTree, current = []):
    """
    Takes:
        - spTree (ete3.TreeNode) : undated, unsubdivided species tree

    Returns:
        (list) : list where the element at index i is the list of potential 
                    transfer receiver post-order of the node with post-order i
                    in an undated tree (ie. all other node excepted for i and its ancestors)

    """


    if current == []:
        current = [[] for i in range(spTree.PO + 1) ]

        current[ spTree.PO ] = range( spTree.PO )
    else:
        current[ spTree.PO ] = current[ spTree.up.PO ][:]
        current[ spTree.PO ].remove(spTree.PO)

    for c in spTree.get_children():
        current = setupPossibleReceiverUndated(c,current)

    return current


if __name__ == "__main__":


    ## testing height
    t = Tree("((a:0.1,b:0.1):0.1,c:0.2);")
    print "ultrametric tree : ",getDistFromRootDic(t)

    t = Tree("((a:0.1,b:0.3):0.1,c:0.2);")
    print "non ultrametric tree : ",getDistFromRootDic(t)
    print "non ultrametric tree without safeguard: ",getDistFromRootDic(t,False)

    ## testing subdivision
    t = Tree("((a:0.1,b:0.1):0.2,(c:0.2,d:0.2):0.1);")
    st = subdivideSpTree(t)
    print t.get_ascii(attributes=["name"])
    print st.get_ascii(attributes=["timeSlice","name"])

    st , SpNameToPO , POToNode , TStoPO  = setupPOandTSdict(st)
    print st.get_ascii(attributes=["PO","timeSlice","name"])
    print TStoPO

    ## testing dead lineage
    dt = addDeadLineage(t)
    #print dt.get_ascii(attributes=["dist","name","dead"])
    st = subdivideSpTree(dt)
    st , SpNameToPO , POToNode , TStoPO  = setupPOandTSdict(st)
    print st.get_ascii(attributes=["PO","timeSlice","name","dead"])
    print TStoPO

    t , SpNameToPO , POToNode , TStoPO  = setupPOandTSdict(t)
    print t.get_ascii(attributes=["PO","name"])
    print setupPossibleReceiverUndated(t)

    dt , SpNameToPO , POToNode , TStoPO  = setupPOandTSdict(dt)
    print dt.get_ascii(attributes=["PO","name","dead"])
    print setupPossibleReceiverUndated(dt)


    st = setupRealNodePO(st)
    print st.get_ascii(attributes=["PO","RealNodePO","name","dead"])