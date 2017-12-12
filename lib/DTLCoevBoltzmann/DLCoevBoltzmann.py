#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        13-Jan-2017         #
##  Last modified:  22-Feb-2017         #
#########################################


from CCP_distribution_ETEbased import CCP_distribution
from ETE3_based_ReconciledTree import RecEvent, ReconciledTree
from DLCoevLiteAux import RecCaseCoEv , RecSolutionCoEv, RecMatrix, readSpTree, trimSpTree, setupPOandTSdict, NodeToPreTrimPO

from ete3 import Tree
import numpy as np
from copy import deepcopy


## coevent per sp branch format: {evtcode: #evt}

DUPCODE = "D"
LOSSCODE = "L"
LEAFCODE = "C"
SPECIATIONCODE = "S"
SPECIATIONLOSSCODE = "SL"

#EVENTINDEXES = {DUPCODE : COEVDUPINDEX , LOSSCODE : COEVLOSSINDEX}


class RecCaseCoEvBoltzmann:
    def __init__(self):
        """
        self.totalScore (float): total sum of scores over all tags
        self.sumScores (dict): keys are tag, values are the sum of scores associated with this tag
        self.solutions (dict): keys are tag, values are lists of solution with this tag and score for this tag
        """
        self.totalScore = 0.

        self.sumScores = {}
    
        self.solutions = {}


    def getScore(self, tag = None):
        if tag is None:
            return self.totalScore
        return self.sumScores[tag]

    def getTags(self):
        return self.sumScores.keys()


    def addSolution(self, solution):
        """ 
            presuming solution is a RecSolutionCoEv instance 
            Only keeps all non impossible 
            (ie. solution with 0 boltzmann proba ) 
            solutions for each possible tag
        """
        key = solution.getHashableConstraint()

        if solution.score == 0:
            return #ignore solution with score == 0 (impossible solution)

        if not self.sumScores.has_key(key) :
            self.sumScores[key] = 0
            self.solutions[key] = []
            
        self.solutions[key].append( solution )
        self.sumScores[key] += solution.score
        self.totalScore += solution.score
        

    def getAllSolutions(self):
        s = []
        for k,v in self.solutions.items():
            s += v
        return s

    def chooseOneSolution(self, tag = None):
        """ 
            randomly choose a solution according to its score
            if tag is not None, only choose solutions with this tag
        """

        opt = []
        total = 0.
        if tag is None:
            opt = self.getAllSolutions()
            total = self.totalScore
        else:
            opt = self.solutions[tag]
            total = self.sumScores[tag]

        #print opt, self.minScores, self.minScore
        r = np.random.random() * total

        #print r , total , len(opt)

        i = -1
        while r > 0:
            i +=1
            r -= opt[i].score

        return opt[ i ]



class BoltzmannDLCoevReconciliator:
    def __init__(self, tree_distFile , spTreeFile, geneNameToSpeciesNameAssociation , costs , sep = "_" , aleFile = False , trim = False):
        """ reconciliator able to accept some events in a guide trees

            Takes:
             - tree_distFile: file containing trees in newock format
             - spTreeFile:
             - geneNameToSpeciesNameAssociation : dictionnary associating leaf names to species name
             - costs : dictionnary associating a cost to different events:

                        COSTS = {"SPE" : exp(- SPECOST / Temperature ),
                                "DUP" : exp(- DUPCOST / Temperature ),
                                "LOSS" : exp(- LOSSCOST / Temperature ),
                                "TOPO" :  TOPOWEIGHT / Temperature,
                                "LOSS_EXTENSION" : 1.,
                                "DUP_EXTENSION"  : 1.,
                                }

             - sep = "_" : separator between sp name and gene name in the gene trees (in case geneNameToSpeciesNameAssociation was an empty dictionnary)
             - aleFile = False : wether the tree_distFile contains trees or is a .ale file
        """
        #self.spTree, self.spNameToPO, self.spPOToNode = readSpTree(spTreeFile , True)



        tree = Tree(spTreeFile, True)

        ## presuming that a list if trees is given
        self.ccp_dist = CCP_distribution()
        IN = open(tree_distFile,"r")

        f1 = None
        f2 = None

        if not aleFile:
            f1 = self.ccp_dist.read_from_treelist_handle
            f2 = self.ccp_dist.read_from_ale_handle
        else:
            f1 = self.ccp_dist.read_from_ale_handle
            f2 = self.ccp_dist.read_from_treelist_handle


        ok = f1(IN)
        if ok != 0:
            #PB. trying to read as ale?
            print "wrong tree distribution format detected. trying alternative reader."
            IN.close()
            IN = open(tree_distFile,"r")
            ok = f2(IN)
            if ok != 0:
                print "failed to read tree distribution file",tree_distFile,"as a .ale file or a list of trees."
                exit(1)

        IN.close()
        self.ccp_dist.set_root_split()


        if trim:

            dNodeToPreTrimPO = NodeToPreTrimPO(tree)

            gList = self.ccp_dist.dleaf_id.values()
            spList = [ g.partition(sep)[0] for g in gList]
            self.spTree = trimSpTree(tree, spList)

            print "trimmed sp tree to",len(self.spTree),"leaves"
        else:
            self.spTree = tree

        #3.setting PO and indexing for easy access later on
        self.spTree , self.spNameToPO , self.spPOToNode , self.TStoPO = setupPOandTSdict(self.spTree)

        ## gestion of ids before and after POs...
        if trim:
            self.POtoPreTrimPO = { po :dNodeToPreTrimPO[n] for po,n in enumerate(self.spPOToNode)}
        else:
            self.POtoPreTrimPO = None

        #print self.spTree.get_ascii(attributes=["PO","RealNodePO","timeSlice","name","dead"])


        #association of leaf and species
        if len(geneNameToSpeciesNameAssociation) == 0:
            geneNameToSpeciesNameAssociation = associateGtoSwithSep(self.ccp_dist.dleaf_id.values() , sep)
        self.AssociateGeneToSpecies(geneNameToSpeciesNameAssociation)

        ## other setup:
        self.costs = costs ## dictionnary

        self.guideCoEvents = {}
        ## dict:
        ##     keys : spId
        ##     values: {evtcode: #evt}

        

    def addGuideCoEvents(self, spId, eventCode , eventNb):
        """ adds guide events. eventCode is supposed to be one of "D" or "L" """
        if not self.guideCoEvents.has_key(spId):
            self.guideCoEvents[spId] = {}

        if not self.guideCoEvents[spId].has_key(eventCode):
            self.guideCoEvents[spId][eventCode] = 0
        self.guideCoEvents[spId][eventCode] += eventNb

    def getCoEvents(self,spId):
        """ returns a dictionnary {evtcode: #evt} """
        return self.guideCoEvents.get(spId , {})

    def getNbAuthorizedCoEvent(self,spId, eventCode):
        """ returns the number of guide coevent of a specific type in a specific species"""
        tmp = self.getCoEvents(spId)
        if not tmp.has_key(eventCode):
            return 0
        return tmp[eventCode]

    def AssociateGeneToSpecies(self, geneNameToSpeciesNameAssociation):
        """ internally associating genes to species """
        self.LeafIdToSpeciesPOAssociation = {}
        for leafId in self.ccp_dist.dleaf_id.keys():
            sp = geneNameToSpeciesNameAssociation.get( self.ccp_dist.dleaf_id[leafId] , None )
            if sp is None:
                print "ERROR : leaf",self.ccp_dist.dleaf_id[leafId] , "is not associated to any species."
                exit(1)
            spId = self.spNameToPO.get(sp,None)
            if spId is None:
                print "ERROR : leaf",self.ccp_dist.dleaf_id[leafId] , "is associated to the species",sp,"which is absent from the given species tree."
                exit(1)
            leafBip = self.ccp_dist.get_bip_from_leafset(set([leafId]))
            self.LeafIdToSpeciesPOAssociation[leafBip] = spId

    def setupSortedIndexes(self):
        """ internally setup indexes of gene tree bipatitions and species post-ordered nodes """
        ## iterating over all bip ids and sorting them by length
        bipLength = {}
        for k in self.ccp_dist.dbip_set.keys():
            l = len(self.ccp_dist.dbip_set[k])
            if not bipLength.has_key(l):
                bipLength[l] = []
            bipLength[l].append(k)

        self.sortedBip = []
        lengths = bipLength.keys()
        lengths.sort()
        for currentSize in lengths:
            self.sortedBip += bipLength[currentSize]
    
        self.sortedSpId = range(len(self.spPOToNode))


    def setupRecMat(self):
        """ internally setup reconciliation matrix and the indexes to parse it """
        self.RMAT = RecMatrix(caseClass = RecCaseCoEvBoltzmann)
        self.setupSortedIndexes()
        self.RMAT.setupIndexes( self.sortedBip , self.sortedSpId)

    def fillRecMat(self):
        """ fill the reconciliation matrix. self.setupRecMat() MUST have been called """
        for bipId in self.sortedBip:
            #print "treating bip" , bipId ,"->", self.ccp_dist.dbip_set[bipId]
            for spId in self.sortedSpId:
                self.treatCase(bipId,spId)

    def treatCase(self, bipId, spId):
        """ filling a case of the matrix """

        leafGene = self.LeafIdToSpeciesPOAssociation.has_key(bipId)

        spChildren = [c.PO for c in self.spPOToNode[spId].get_children()]
        spNbChildren = len( spChildren )
        LeafSp = ( spNbChildren == 0 )

        ## first case : both are leaves --> trying leaf association
        if LeafSp and leafGene:
            self.CurrentCase( bipId, spId )

        if not leafGene:
            ## iterate over possible partition of the current bip
            for part,part_count in self.ccp_dist.ddip_count[bipId].items():
                topologyScore =  ( part_count / float(self.ccp_dist.dbip_count[bipId]) )**self.costs["TOPO"]
                #print "topologyScore",topologyScore , "<-" , "(", part_count, "/", float(self.ccp_dist.dbip_count[bipId]) , ")**", self.costs["TOPO"]
                ## speciation
                if spNbChildren == 2:
                    self.SpeciationCase(bipId, spId, part, topologyScore)

                ## duplication
                self.DuplicationCase(bipId, spId, part, topologyScore)

        if spNbChildren == 2: ## speciationLoss
            self.SpeciationLossCase(bipId , spId)

    def CurrentCase(self, bipId, spId):
        """ clade where the event is 'leaf' (if bipId and spId are both leaves) """
        score = 0.
        if self.LeafIdToSpeciesPOAssociation[bipId] == spId:
            score = 1.
        self.RMAT.set_case(bipId,spId , RecSolutionCoEv(score , [] , LEAFCODE) )

    def SpeciationLossCaseAux(self, bipId, spId, childBearingSp, lossBearingSp):
        """ auxiliary function for the speciationLoss case """

        ## score without any co-event
        scoreBase = self.costs["LOSS"] * self.costs["SPE"] 
        ## score in case of co-event
        scoreExtension = self.costs["LOSS_EXTENSION"] * self.costs["SPE_EXTENSION"]
        scoreChild =  self.RMAT.get_case_score( bipId , childBearingSp )

        components = [( bipId , childBearingSp )]

        ## no co-event
        score = scoreBase * scoreChild
        self.RMAT.add_case(bipId,spId , RecSolutionCoEv(score , components , SPECIATIONLOSSCODE) )

        ## shunt in case of infinite solution.
        if scoreChild == 0.:
            return

        ## co-event
        eventCode = LOSSCODE
        Authorized = self.getNbAuthorizedCoEvent(lossBearingSp, eventCode)
        if Authorized > 0:

            score = scoreExtension * scoreChild

            coevent = True
            recSolution = RecSolutionCoEv(score , components , SPECIATIONLOSSCODE, 0, coevent)

            self.RMAT.add_case(bipId,spId , recSolution )                    
            #print "adding co speciation loss (" , bipId , "," , spId , ")" 

        return


    def SpeciationLossCase(self, bipId, spId):
        """ speciation loss case """
        spChildren = [c.PO for c in self.spPOToNode[spId].get_children()]

        self.SpeciationLossCaseAux( bipId, spId, spChildren[0], spChildren[1])
        self.SpeciationLossCaseAux( bipId, spId, spChildren[1], spChildren[0])
        return

    def DuplicationCase(self, bipId, spId, part, topologyScore):
        """ 
        duplication case, the presence of a split in the gene tree 
        means that a score for the topology is needed
        """

        ## simple duplication whithout co-event
        scoreChild0 = self.RMAT.get_case_score( part[0] , spId )
        scoreChild1 = self.RMAT.get_case_score( part[1] , spId )

        components = [  (part[0],spId ),(part[1],spId ) ]

        score = self.costs["DUP"] * topologyScore * scoreChild0 * scoreChild1

        if score == 0: ## solution impossible anyway, just put dummy 
            self.RMAT.add_case(bipId,spId , RecSolutionCoEv(score , components, DUPCODE ) )
            return

        TagsChild0 = self.RMAT.get_case_tags( part[0], spId )
        TagsChild1 = self.RMAT.get_case_tags( part[1], spId )

        ## number of guide co-events
        eventCode = DUPCODE
        Authorized = self.getNbAuthorizedCoEvent(spId, eventCode)

        ## scores without and with co-event
        scoreBase = self.costs["DUP"] * topologyScore
        scoreExtension = self.costs["DUP_EXTENSION"] * topologyScore

        for tag0 in TagsChild0:
            for tag1 in TagsChild1:

                components = [  (part[0],spId , tag0),(part[1],spId , tag1)  ]

                scoreChild0 = self.RMAT.get_case_score( part[0], spId , tag0)
                scoreChild1 = self.RMAT.get_case_score( part[1], spId , tag1)

                nbCoDup = 0
                ## here, tags are simply the number of duplication
                
                nbCoDup = max(tag0, tag1)

                isPossible = ( nbCoDup <= Authorized )
                if not isPossible:
                    continue

                #print nbCoDup, "<>" , Authorized

                ## it is possible to associate the two
                coevent = False # no coevent at this node
                score = scoreBase * scoreChild0 * scoreChild1
                self.RMAT.add_case(bipId,spId , RecSolutionCoEv(score , components , eventCode, nbCoDup , coevent) )
                #print "added dup", bipId,spId, nbCoDup ,"->",score
                ###NB: this works only because the duplication are the only non terminal event treated here

                if nbCoDup < Authorized:
                    ## there is room for at least one more co-duplication
                    coevent = True
                    score = scoreExtension * scoreChild0 * scoreChild1
                    self.RMAT.add_case(bipId,spId , RecSolutionCoEv(score , components , eventCode, nbCoDup + 1, coevent) )
                    #print "added co-dup", bipId,spId , nbCoDup + 1,"->",score

    def SpeciationCase(self, bipId, spId, part, topologyScore):
        """
        speciation case, the presence of a split in the gene tree 
        means that a score for the topology is needed
        """
        spChildren = [c.PO for c in self.spPOToNode[spId].get_children()]

        CaseComponents = [ ( (part[0],spChildren[0]),(part[1],spChildren[1]) ) , # c0 in s0 , c1 in s1
                      ( (part[0],spChildren[1]),(part[1],spChildren[0]) ) ] # converse
        for components in CaseComponents: ## for each case
            score =  self.costs["SPE"] * topologyScore 
            for c in components: ## for each component of the case
                score *= self.RMAT.get_case_score( c[0] , c[1] )
            self.RMAT.add_case(bipId,spId ,  RecSolutionCoEv(score , components , "S") )
        return


    def backtrackRecMat(self, bipId, spId , solutionTag = None ):
        """ 
        !recursive function!

        Takes:
            - bipId (int) : bipartition id
            - spId  (int) : species id
            - solutionIndex [ default= None ] : None or Tag (should be an int)

        Returns:
            (tuple): reconciled subtree
        """
        ## 1. chosing a solution
        
        chosenSolution = self.RMAT[bipId][spId].chooseOneSolution(solutionTag)

        ## 2. backtracking that solution
        childSubTrees = []
        childSubTreesCoevents = []
        for component in chosenSolution.component:
            ## component is a tuple of size either 2 or 1
            childBipId = component[0]
            childspId = component[1]
            childTag = None
            if len(component)>2:
                childTag = component[2]
            sub  = self.backtrackRecMat(childBipId , childspId , childTag)
            childSubTrees.append( sub )


        ##3. creating the reconciliation event associated with the reconciliation
        
        timeSlice = None #no timeslice in this version
        species = spId
        additionnal = {}
        
        #if self.spPOToNode[spId].is_leaf():
        #    species = self.spPOToNode[spId].name
        #
        #el
        if not self.POtoPreTrimPO is None:
            species = self.POtoPreTrimPO[spId]

        if chosenSolution.coev:
            additionnal["coevent"] = True

        event = RecEvent(chosenSolution.event , species , timeSlice , additionnal)

        ##4. creating a new node if need be
        currentNode = None
        if len(childSubTrees) != 1 : ## more or less than 1 child subtree -> n-furcation in gene tree -> create new node
            currentNode = ReconciledTree()

            name = bipId
            if self.ccp_dist.is_leaf(bipId):
                name = self.ccp_dist.get_leaf_name(bipId)
            currentNode.setName(name)
            for child in childSubTrees:
                currentNode.add_child(child)

        else: ##only one subtree -> we add the event at the beginning of it 
            currentNode = childSubTrees[0]

        currentNode.addEvent(event , append = False)
        return currentNode


def readGuideEvents(Gfile, sptree, spNameToPO):
    """
    Each line is an event. It begins with "D", "T" or "L" (to indicate wether the event is a duplication, a transfer reception or a loss)
                                                            and is followed by the species it is in as the list (separated by spaces) of the leaves it is an ancestor to.
                                    example: D 1 2 3
                                            signifies a duplication in the species that is the LCA of the leaves a, b and c.
    """

    IN = open(Gfile,"r")
    l = IN.readline()
    guideEvents = []
    while l != "":
        sl = l.strip().split()
        evtType = sl[0]
        if not evtType in ["D","T","L"]:
            print "unknown co-event type",evtType
        else:
            leaves = []
            species = None
            for l in sl[1:]:
                nodes = sptree.get_leaves_by_name(l)
                if len(nodes) == 0:
                    print "no species called",l
                leaves.append(nodes[0])
            #print leaves, spNameToPO
            if len(leaves) == 0:
                print "unassigned coevent : no valid leaves."
            elif len(leaves) == 1:
                species = leaves[0].PO
            else:
                species = leaves[0].get_common_ancestor(leaves[:]).PO
            if not species is None:
                guideEvents.append( tuple([ species , evtType ]) )
        l = IN.readline()
    IN.close()
    #print guideEvents
    return guideEvents


#def readGuideEvents(Gfile, sptree, spNameToPO):
#    """
#    Each line is an event. It begins with "D" or "L" (to indicate wether the event is a duplication or a loss)
#                                                            and is followed by the species it is in as the list (separated by spaces) of the leaves it is an ancestor to.
#                                    example: D 1 2 3
#                                            signifies a duplication in the species that is the LCA of the leaves a, b and c.
#    """
#
#    IN = open(Gfile,"r")
#    l = IN.readline()
#    guideEvents = []
#    while l != "":
#        sl = l.strip().split()
#        evtType = sl[0]
#        if not evtType in ["D","L"]:
#            print "unknown co-event type",evtType
#        else:
#            leaves = []
#            species = None
#            for l in sl[1:]:
#                PO = spNameToPO.get(l,-1)
#                print spNameToPO
#                print sptree
#                nodes = sptree.get_leaves_by_name(PO)
#                print l,PO , nodes
#                if len(nodes) == 0:
#                    print "no species called",l
#                leaves.append(nodes[0])
#            #print leaves, spNameToPO
#            if len(leaves) == 0:
#                print "unassigned coevent : no valid leaves."
#            elif len(leaves) == 1:
#                species = leaves[0].PO
#            else:
#                species = leaves[0].get_common_ancestor(leaves[:]).PO
#            if not species is None:
#                guideEvents.append( tuple([ species , evtType ]) )
#        l = IN.readline()
#    IN.close()
#    #print guideEvents
#    return guideEvents

def associateGtoSwithSep(LeafList , sep ="_"):
    assoc = {}
    for l in LeafList:
        sp , j , jbis = l.partition(sep)
        assoc[l] = sp
    return assoc



if __name__ == "__main__":


    import sys

    help = """ 
    this script reconciles a gene distribution with a species tree
    by computing a score including topology, Duplication, Loss, and co-events with a provided guide tree
    in a pseudo likelihood framework such that a more parsimonious tree will have more chance to be sampled.
    usage:
        python DLCoevBoltzmann.py  -s species_tree -g gene_tree_distribution [options]
    options:
        -s : name of the species tree file (newick)
        -g : name of the gene tree distribution file (newick). leaf names MUST be in the format : speciesName, followed by geneName, separated by a character (option --sep, '_' by default)
        -G : (optional) name of a file containing potential co-events in a guide tree 
                                    Each line is an event. It begins with "D" or "L" (to indicate wether the event is a duplication or a loss)
                                                            and is followed by the species it is in as the list (separated by spaces) of the leaves it is an ancestor to.
                                    example: D 1 2 3
                                            signifies a duplication in the species that is the LCA of the leaves a, b and c.

        -a : (optional) use to specify that the gene tree distribution file is, in fact, a .ale file (will try to detect automatically anyway)
        -n : (optional) number of trees to sample (default: 1)

        -d : (optional) cost of a single duplication (default: 2)
        -l : (optional) cost of a single loss (default: 1)
        -T : (optional) weight of the topology (default: 1)
        --temperature : (optional) a higher temperature will mean that less parsimonious trees will have a higher chance to be sampled (default: 0.1)
        --de : (optional)  cost of the extension of a co-duplication (default: 1)
        --le : (optional)  cost of the extension of a co-loss (default: 0.5)

        -o : (optional) name of the file to print the generated tree in (default is stdout)
        --sep : (optional) separator between species name and gene name in the gene tree leaves (default '_')
        --seed : (optional) seed to use in the randomisation (by default, time() is used)
        --topoOnly : (optional) if present, only topology (without reconciliation annotation) will be outputted
        --xml : (optional) if present, the tree will be written in an XML format

        --trim.sp.tree : (optional) trim the species tree above the LCA of the species which have a leaf associated to them
    """
        #--guide : (optional) name of the guide tree file (recPhyloXML)

    OK = True

    nextKEY = None
    params = {
        "-s" : None,#name of the species tree file (newick)
        "-g" : None,#name of the gene tree distribution file (newick)
        "-G" : None,#(optional) name of a file containing potential co-events in a guide tree 

        "-a" : False,#(optional) use to specify that the gene tree distribution file is, in fact, a .ale file
        "-n" : 1, #(optional) number of trees to sample (default: 1)

        "-d" : 2. ,#(optional) cost of a single duplication (default: 2)
        "-l" : 1.,#(optional) cost of a single loss (default: 1)
        "-T" : 1.,#(optional) weight of the topology (default: 1)
        "--temperature" : 0.1,#(optional) a higher temperature will mean that less parsimonious trees will have a higher chance to be sampled (default: 1)
        "--de" : 1., #(optional)  cost of the extension of a co-duplication (default: 1)
        "--le" : 0.5, #(optional)  cost of the extension of a co-loss (default: 0.5)

        "-o" : None, #(optional) name of the file to print the generated tree in (default is stdout)
        "--sep" : '_',#(optional) separator between species name and gene name in the gene tree leaves (default '_')
        "--seed": None, #(optional) seed to use in the randomisation (by default, time() is used)
        "--topoOnly" : False, #(optional) if True, only topology (without reconciliation annotation) will be outputed
        "--xml" : False, #(optional), if True, the tree will be written in an XML format.
        "--trim.sp.tree" : False #(optional) trim the species tree above the LCA of the species which have a leaf associated to them
    }

    flagArgs = ["-a", "--topoOnly", "--xml", "--trim.sp.tree"]

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

    if params["-s"] is None:
        OK = False
        print "error: species tree not given."
    elif params["-g"] is None :
        OK = False
        print "error: gene tree not given."


    if OK:
        
        ## treating positive float options
        for pname in ["-d","-l","-T","--de","--le","--temperature"]:
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


    if not OK:
        print help
        exit(1)



    if params["--seed"] is None:
        np.random.seed()
    else :
        np.random.seed(int(params["--seed"])) ## put a fixed seed to ensure reproductibility

    spFilename = params["-s"]

    distrFilename = params["-g"] 

    SPECOST = 0.

    Temperature = params["--temperature"]
    

    COSTS = {"SPE" : np.exp(- SPECOST / Temperature ) ,
            "DUP" : np.exp(- params["-d"] / Temperature ),
            "LOSS" : np.exp(- params["-l"] / Temperature ),
            "TOPO" : params["-T"] / Temperature,
            "LOSS_EXTENSION" : np.exp(- params["--le"] / Temperature ),
            "DUP_EXTENSION"  : np.exp(- params["--de"] / Temperature ),
            "SPE_EXTENSION"  : np.exp(- SPECOST / Temperature ),
            }

    ### association between genes and species
    
    geneNameToSpeciesNameAssociation = {}

    Rec = BoltzmannDLCoevReconciliator(distrFilename , spFilename , geneNameToSpeciesNameAssociation , COSTS, params['--sep'], params['-a'], params["--trim.sp.tree"])

    print "trees read and bipartitions computed"
    ##adding some guide events
    if not params["-G"] is None:
        coevents = readGuideEvents(params["-G"], Rec.spTree, Rec.spNameToPO)
        for e in coevents:
            Rec.addGuideCoEvents( spId = e[0] , eventCode= e[1] , eventNb=1)

    print "potential co-events added"

    Rec.setupRecMat()
    Rec.fillRecMat()

    print "reconciliation matrix computed"



    ##loop to check for inconsistencies
    exponentialLimit = 309
    exponentialLimit -= int( - np.log( min(COSTS.values()) ) / np.log(10.) )

    #print "exponentialLimit" , exponentialLimit
    potentialOverflow = False
    for bipId in Rec.sortedBip:
            #print "treating bip" , bipId ,"->", Rec.ccp_dist.dbip_set[bipId]

        for spId in Rec.sortedSpId:
            s = Rec.RMAT[bipId][spId].totalScore 
            if s !=0 and int( - np.log(s) / np.log(10.) ) > exponentialLimit:
                print "!!WARNING!! some scores are very close to overflow"
                print "Try with a higher temperature to avoid problems."
                potentialOverflow = True
                break
        if potentialOverflow:
            break
        

    rootBip = Rec.ccp_dist.get_root_bip_id() 



    totalAtRoot = 0.
    RootP = []
    RootList = []
    for j,dsol in enumerate(Rec.RMAT[rootBip]):
        sp = Rec.sortedSpId[j]
        RootList.append(sp)
        RootP.append(dsol.totalScore)
        totalAtRoot += dsol.totalScore

    
    OUT = sys.stdout
    if not params["-o"] is None:
        OUT = open(params["-o"],"w")

    speciesNames = {}#{v:k for k,v in Rec.spNameToPO.items()}
    #print speciesNames

    if totalAtRoot == 0:
        print "No non impossible solution found"
        print "This can be caused by an overflow problem. Try with a higher temperature."
    else:
        print "sampling" ,params["-n"], "reconciliations"

        if not params["--xml"]:
            print "format newick nodes: name|evt0Species.evt0Code|evt1Species.evt1Code ... "
        #print "internal node name are their bipartition name in the CCP distribution."
        for n in range(params["-n"]):
            r = np.random.random() * totalAtRoot
            i = -1
            while r > 0:
                i += 1
                r -= RootP[i]
            sp = RootList[i]
            RT = Rec.backtrackRecMat(rootBip, sp)

            ### uncommenting this will add a series of counts for each event at the beginning of the line
            #DE = RT.countEvents(True)
            #OUT.write( ",".join([str(DE.get(e,0)) for e in ["D","SL","coD","coSL"]]) + " " )

            if not params["--xml"]:
                OUT.write( RT.getTreeNewick(topoOnly=params['--topoOnly']) + "\n" )
            else:
                OUT.write( RT.getTreeRecPhyloXML( speciesNames = speciesNames, topoOnly=params['--topoOnly']) + "\n" ) 

    OUT.close()