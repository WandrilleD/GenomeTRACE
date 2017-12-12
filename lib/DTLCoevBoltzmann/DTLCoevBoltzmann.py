#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        02-Feb-2017         #
##  Last modified:  08-Feb-2017         #
#########################################


from CCP_distribution_ETEbased import CCP_distribution
from ETE3_based_ReconciledTree import RecEvent, ReconciledTree
from DLCoevLiteAux import RecCaseCoEv , RecSolutionCoEv, RecMatrix, setupPOandTSdict, NodeToPreTrimPO, trimSpTree, POLeafList, RecMatrixOld
from DTLCoevAux import getDistFromRootDic, subdivideSpTree, addDeadLineage, getPOtoDead, setupPossibleReceiverUndated, setupRealNodePO

import ete3
from ete3 import Tree
import numpy as np
from copy import deepcopy

#from memory_profiler import profile
## coevent per sp branch format: {evtcode: #evt}

DUPCODE = "D"
LOSSCODE = "L"
LEAFCODE = "C"
SPECIATIONCODE = "S"
SPECIATIONLOSSCODE = "SL"

NULLCODE = "N"

TRANSFERCODE = "T"
TRANSFERLOSSCODE = "TL"

SPECIATIONOUTCODE = "So"
SPECIATIONOUTLOSSCODE = "SoL"

BIFOUTCODE = "Bo"

#EVENTINDEXES = {DUPCODE : COEVDUPINDEX , LOSSCODE : COEVLOSSINDEX}


class RecCaseCoEvBoltzmann:
    def __init__(self):
        """
        self.totalScore (float): total sum of scores over all tags
        
        #self.sumScores (dict): keys are tag, values are the sum of scores associated with this tag
        #self.solutions (dict): keys are tag, values are lists of solution with this tag and score for this tag
        
        self.sumScores (list): keys are tag, values are the sum of scores associated with this tag
        
        self.solutions (list): keys are tag, values are lists of solution with this tag and score for this tag
        """
        self.totalScore = 0.

        self.sumScores = []
        self.solutions = []

        #self.sumScores = {}
        #self.solutions = {}



    def getScore(self, tag = None):
        if tag is None:
            return self.totalScore
        return self.sumScores[tag]

    def getTags(self):
        #return self.sumScores.keys()
        return range(len(self.sumScores))


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

        while len(self.sumScores) <= key:
            self.sumScores.append(0.)
            self.solutions.append([])
        #if not self.sumScores.has_key(key) :
        #    self.sumScores[key] = 0
        #    self.solutions[key] = []

        

        self.solutions[key].append( solution )
        self.sumScores[key] += solution.score
        self.totalScore += solution.score
        
        #print "addSolution" , solution.event , solution.tag , self.solutions.keys()

    def getAllSolutions(self):
        s = []
        #for k,v in self.solutions.items():
        for k,v in enumerate(self.solutions):
            s += v
        return s

    def chooseOneSolution(self, tag = None, excludeFromTL = None):
        """ 
            randomly choose a solution according to its score
            if tag is not None, only choose solutions with this tag

            if excludeFromTL is not None, will exclude event marked as TL and SoL toward the specified species from the choice
        """

        opt = []
        total = 0.
        if tag is None:
            opt = self.getAllSolutions()
            total = self.totalScore
        else:
            opt = self.solutions[tag]
            total = self.sumScores[tag]

        if not excludeFromTL is None:

            newopt = []
            total = 0.
            for o in opt:
                if o.event in [TRANSFERLOSSCODE, SPECIATIONOUTLOSSCODE]:
                    if o.component[0][1] == excludeFromTL: ## component of a TL of SoL has one member whose second object is the id of the receiver
                        continue
                newopt.append(o)
                total += o.score

            opt = newopt

        #print opt, self.minScores, self.minScore
        r = np.random.random() * total

        #print r , total , len(opt)

        i = -1
        while r > 0:
            i +=1
            r -= opt[i].score

        return opt[ i ]



class BoltzmannDTLCoevReconciliator:
    def __init__(self, tree_distFile , spTreeFile, geneNameToSpeciesNameAssociation , costs , sep = "_" , aleFile = False, dated = True, trim = False ):
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
             - aleFile = False : whether the tree_distFile contains trees or is a .ale file
             - dated = True : whether the species tree is dated or not (contdition subdivision)
             - trim = False : if True, trim the species tree at the LCA of the species that are associated to a gene
        """

        #self.spTree, self.spNameToPO, self.spPOToNode = readSpTree(spTreeFile , True)
        #self.TStoPO  = None

        ##for now expecting a fully ultrametric tree

        #1.reading
        try:
            tree = Tree(spTreeFile, True)
        except ete3.parser.newick.NewickError as e:
            print "!!ERROR!! while reading species tree file :" , e
            exit(1)



        ## presuming that a list if trees is given
        self.ccp_dist = CCP_distribution()

        try:
            IN = open(tree_distFile,"r")
        except IOError as e:
            print "!!ERROR!! while reading gene tree file :" , e
            exit(1)

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

            #dNodeToPreTrimPO = NodeToPreTrimPO(tree)
            dPreTrimPOLeafList = POLeafList(tree)
            #print dPreTrimPOLeafList

            gList = self.ccp_dist.dleaf_id.values()
            spList = [ g.partition(sep)[0] for g in gList]
            self.spTree = trimSpTree(tree, spList)


        else:
            self.spTree = tree



        #1.5 adding the dead
        self.has_dead_lineage = True ####### SHOULD ALWAYS BE TRUE WHEN UNDATED BECAUSE OF MEMORY OPTIM !!! ######
        if self.has_dead_lineage:
            self.spTree = addDeadLineage(self.spTree)

        if dated:
            #2.subdividing
            self.spTree = subdivideSpTree(self.spTree)
            if self.spTree is None:
                print "use the --undated option for non ultrametric trees"
                exit(1) ## happens when the species tree isn't ultrametric
        
        #3.setting PO and indexing for easy access later on
        self.spTree , self.spNameToPO , self.spPOToNode , self.TStoPO = setupPOandTSdict(self.spTree)

        if not dated: ## in the undated version, the possible transfer receiver are compued differently
            self.undatedPossibleReceiver = setupPossibleReceiverUndated(self.spTree)

        if self.has_dead_lineage:
            self.spPOtoDead = getPOtoDead(self.spPOToNode)

        self.spTree = setupRealNodePO(self.spTree)

        #print self.spTree.get_ascii(attributes=["PO","RealNodePO","timeSlice","name","dead"])
        
        ## gestion of ids before and after POs...
        if trim:

            dPostTrimPOLeafList = POLeafList(self.spTree)
            #print dPostTrimPOLeafList

            self.POtoPreTrimPO = {}
            for postPO in dPostTrimPOLeafList:

                for prePO in dPreTrimPOLeafList:

                    if dPostTrimPOLeafList[postPO] == dPreTrimPOLeafList[prePO]:
                        self.POtoPreTrimPO[postPO] = prePO

            #print self.POtoPreTrimPO
            #self.POtoPreTrimPO = { po :dNodeToPreTrimPO[n] for po,n in enumerate(self.spPOToNode)}
            #print self.POtoPreTrimPO
        else:
            self.POtoPreTrimPO = None









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
        """ adds guide events. eventCode is supposed to be one of "D", "T" or "L" """
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
        realSpId  = self.spPOToNode[spId].RealNodePO  ##we have to deal with the articial nodes of the subdivision...

        tmp = self.getCoEvents(realSpId)
        return tmp.get(eventCode,0)

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
    
        self.sortedSpId = range(len(self.spPOToNode) - self.has_dead_lineage) ## when there is a dead lineage, the root is artificial and we don't want to compute it


    def setupRecMat(self):
        """ internally setup reconciliation matrix and the indexes to parse it """
        #self.RMAT = RecMatrixOld(caseClass = RecCaseCoEvBoltzmann)
        self.RMAT = RecMatrix(caseClass = RecCaseCoEvBoltzmann)
        self.setupSortedIndexes()
        self.RMAT.setupIndexes( self.sortedBip , self.sortedSpId)



    #@profile
    def fillRecMat(self):
        """ fill the reconciliation matrix. self.setupRecMat() MUST have been called """
        nbBip = len(self.sortedBip)
        for i,bipId in enumerate(self.sortedBip):
            #print "treating bip" , bipId ,"->", self.ccp_dist.dbip_set[bipId]
            if self.TStoPO is None: ## species tree is not subdivided
                for spId in self.sortedSpId:
                    self.treatCase(bipId,spId)
                DTbtoSelf = self.TransferBackCase( bipId, None)
                self.TransferLossCase( bipId, None, DTbtoSelf )
            else: ##species tree is subdivided
                nbTS = len(self.TStoPO)
                if self.has_dead_lineage:
                    nbTS -= 1 ##when there is a dead lineage, the root is artificial. 

                for TS in range(nbTS):
                    for spId in self.TStoPO[TS]:
                        self.treatCase(bipId,spId,TS)
                        #print "treat",bipId,spId
                    DTbtoSelf = self.TransferBackCase( bipId, TS)
                    self.TransferLossCase( bipId, TS, DTbtoSelf )
            print "computed clade",i,"/",nbBip,"\r",
        print "computed clade",nbBip,"/",nbBip



    def TransferBackCase(self, bipId, timeSlice):
        """
        generate transfer back solutions towards all species at the given timeSlice

        Takes:
            - bipId (int) : current gene clade
            - timeSlice (int) : current timeSlice / None if undated

        Returns:
            (dict): keys are species id , values are the sums of the costs of doing transferBack to the given species

        """
        if not self.has_dead_lineage:
            print "!!ERROR!! : this version of the algorithm necessarily transfers trough dead lineages. please allow dead lineages."
            exit(1)
        
        Species = None
        if not timeSlice is None:
            Species = self.TStoPO[timeSlice]
        else:
            Species = range( len(self.undatedPossibleReceiver)-self.has_dead_lineage ) ##in case of dead linage, the root is artificial and we ignore it

        Solutions = []
        CoSolutions = {} # : spId , v : Solution
        
        deadGiver = None

        ###get the solution scores BEFORE we add Tb event
        for receiverSp in Species: ##different iteration when undated
            score = self.RMAT.get_case_score( bipId , receiverSp )

            evtCode = TRANSFERLOSSCODE

            if self.spPOtoDead[receiverSp]: 
                deadGiver = receiverSp
                continue ## no transfer back to the dead
            else:
                score *= self.costs[TRANSFERCODE]
            


            Solutions.append( RecSolutionCoEv(score, [(bipId , receiverSp)] , evtCode, 0 , 0 ) )
            ##coevent goes there
            
            ## checking potential co-events
            nbCoev = self.getNbAuthorizedCoEvent(receiverSp,"T")
            if nbCoev > 0: ## co-event possible -> we place it
                score = self.RMAT.get_case_score( bipId , receiverSp ) * self.costs["TRANSFER_EXTENSION"]
                CoSolutions[receiverSp] =  RecSolutionCoEv(score , [(bipId , receiverSp)] , evtCode, 0, 1)


        DTbtoSelf = {}
        
        receiversIndex = range(len(Solutions))

        for j in receiversIndex:

            solution = Solutions[j]
            Rsol = solution.copy()

            self.RMAT.add_case( bipId, deadGiver , Rsol )

            receiverSp = Rsol.component[0][1]
            DTbtoSelf[receiverSp] = Rsol.score

            if CoSolutions.has_key( receiverSp ):
                ## there is a co-event to form in the receiving species!
                Rsol = CoSolutions[receiverSp].copy()
                
                self.RMAT.add_case( bipId, deadGiver , Rsol )

                DTbtoSelf[receiverSp] += Rsol.score

        return DTbtoSelf



    def TransferLossCase(self, bipId, timeSlice , DTbtoSelf ):
        """
        generate transfer loss solutions accross all species at the given timeSlice
        TOWARD the dead lineage, so all this is only speciationOutLoss

        NB : forcing all transfers to go through the dead lineage does not change the model (if the necessary adjustment to prevent transfer looping are made)
                and makes the algorithm to go from |G|*|S|**2 to |G|*|S|*2 in time and memory

        Takes:
            - bipId (int) : current gene clade
            - timeSlice (int) : current timeSlice / None if undated
            - DTbtoSelf (dict) : keys are species id , values are the sums of the costs of doing transferBack to the given species

        """
        if not self.has_dead_lineage:
            print "!!ERROR!! : this version of the algorithm necessarily transfers trough dead lineages. please allow dead lineages."
            exit(1)

        
        Species = None
        if not timeSlice is None:
            Species = self.TStoPO[timeSlice]
        else:
            Species = range( len(self.undatedPossibleReceiver)-self.has_dead_lineage ) ##in case of dead linage, the root is artificial and we ignore it

        Solution = None
        
        ###get the solution for the only possible receiver -> the dead species
        for receiverSp in Species: ##different iteration when undated
            if self.spPOtoDead[receiverSp]: ##transferLoss to the dead, transfer free
                score = self.RMAT.get_case_score( bipId , receiverSp )

                evtCode = SPECIATIONOUTLOSSCODE

                Solution =  RecSolutionCoEv(score, [(bipId , receiverSp)] , evtCode, 0 , 0 ) 
                break
            else:
                continue
    

        for giverSp in Species:

            scoreModif = 1.

            if self.spPOtoDead[giverSp]: ##transferLoss from the dead -> loss is free
                continue


            scoreModif *= self.costs["LOSS"]

            TbToSelfCost = DTbtoSelf[giverSp]


            Rsol = Solution.copy()            
            Rsol.score -= TbToSelfCost ## modifying the cost to ignore the cost of doing a transfer back toward the same species
            Rsol.score *= scoreModif ## adding the LOSS to the cost

            self.RMAT.add_case( bipId, giverSp , Rsol )

            ## checking potential co-events -> co-losses in the iver species
            nbCoev = self.getNbAuthorizedCoEvent(giverSp,"L")

            if nbCoev > 0 : # co-losses
                Rsol = Solution.copy()
                
                Rsol.score -= TbToSelfCost ## modifying the cost to ignore the cost of doing a transfer back toward the same species
                Rsol.score *= self.costs["LOSS_EXTENSION"]

                Rsol.coev += 2

                self.RMAT.add_case( bipId, giverSp , Rsol )

        return

    def treatCase(self, bipId, spId, timeSlice = None):
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

                ## transfer
                self.TransferCase( bipId, spId, part, topologyScore, timeSlice )

        if spNbChildren == 2: ## speciationLoss
            self.SpeciationLossCase(bipId , spId)
        elif spNbChildren == 1:
            self.NullCase(bipId , spId)




    def iter_receivers(self, bipId, spId , timeSlice):
        """
        iterate across all case with clade bipId and species different from spId with the same timeSlice.

        Takes:
            - bipId (int) : current gene clade
            - spId (int) : current species id
            - timeSlice (int) : time slice of spId (None in undated case)

        Yields:
            (tuple)
                (float) : score
                (tuple) : bipId , spId
        """
        
        receiver = None
        if not timeSlice is None:
            receiver = self.TStoPO[timeSlice]
        else:
            receiver = self.undatedPossibleReceiver[spId]

        for sId in receiver:
            if sId != spId:
                yield self.RMAT.get_case_score( bipId , sId) , (bipId,sId)



    def generateTransferSolution(self, bipId , spId , stayingId, leavingId , topologyScore, timeSlice):
        """
        generates the different solutions when a transfer split bipId in stayingId and leavingId
        and stayingId stays in spId while leavingId is transfered to another species (with the same timeSlice) 

        Takes:
            - bipId (int) : current gene clade
            - spId (int) : current species id
            - stayingId (int) : id of the clade child that stays in spId
            - leavingId (id) : id of the clade child that leaves spId
            - topologyScore (float) : part of the score corresponding to this split
            - timeSlice (int) : time slice of spId (None in undated case)

        """

        #some non impossible solutions exists
        #as child0 stays, we have to take the tag into account
        TagsChild = self.RMAT.get_case_tags( stayingId , spId )
        for tag in TagsChild:
            stayComponent = ( stayingId , spId , tag )
            baseScore = topologyScore * self.RMAT.get_case_score( *stayComponent )

            #iterate accross possible receivers to only keep the dead one:
            receiver = None
            if not timeSlice is None:
                receiver = self.TStoPO[timeSlice]
            else:
                receiver = self.undatedPossibleReceiver[spId]
    

            deadReceiver = None
            for sId in receiver:
                if self.has_dead_lineage:
                    if self.spPOtoDead[sId]:
                        deadReceiver = sId
                        break
                else:
                    print "!!ERROR!! : this version of the algorithm necessarily transfers trough dead lineages. please allow dead lineages."
                    exit(1)

            scoreReceiver = self.RMAT.get_case_score( leavingId , deadReceiver)

            ## this is a dead lineage -> a transfer to it costs nothing and is called a SpeciationOut
            evtCode = SPECIATIONOUTCODE

            score = baseScore * scoreReceiver * self.costs[evtCode]

            recSolution = RecSolutionCoEv(score , [stayComponent , ( leavingId, deadReceiver ) ] , evtCode, tag, False)
            self.RMAT.add_case(bipId,spId , recSolution )

            ## checking potential co-events -> no: only co-events in non-dead lineages
            #nbCoev = self.getNbAuthorizedCoEvent(ReceiverComponent[1],"T")
            #if nbCoev > 0: ## co-event possible -> we place it
            #    score = baseScore * scoreReceiver * self.costs["TRANSFER_EXTENSION"]
            #    recSolution = RecSolutionCoEv(score , [stayComponent , ReceiverComponent ] , evtCode, tag, True)
            #    self.RMAT.add_case(bipId,spId , recSolution )


    def TransferCase(self, bipId, spId, part, topologyScore, timeSlice ):
        """
        case of a transfer where 1 child stays in the tree 

        Takes:
            - bipId (int) : current gene clade
            - spId (int) : current species id
            - part (tuple) : ids of the two clade child in the gene tree
            - topologyScore (float) : part of the score corresponding to this split
            - timeSlice (int) : time slice of spId (None in undated case)

        """

        if self.has_dead_lineage:
            if self.spPOtoDead[spId]:
                return ## no simple transfer from the dead. -> we prefer a bifurcationOut event

        scoreChild0 = self.RMAT.get_case_score( part[0] , spId )
        scoreChild1 = self.RMAT.get_case_score( part[1] , spId )

        if scoreChild0 != 0: ## child0 stays in spId
            self.generateTransferSolution(bipId , spId , part[0], part[1] , topologyScore, timeSlice)
        if scoreChild1 != 0:## child1 stays in spId
            self.generateTransferSolution(bipId , spId , part[1], part[0] , topologyScore, timeSlice)
        
        ## if no possible solutions -> just put a dummy solution?

    def NullCase(self, bipId , spId):
        """ case where nothing happened that the clade just goes to the unique children of the species """
        spChildrenId = self.spPOToNode[spId].children[0].PO

        ##shunt for impossible solution
        if self.RMAT.get_case_score( bipId , spChildrenId) == 0.0:
            self.RMAT.add_case(bipId,spId , RecSolutionCoEv(0.0 ,  [ ( bipId , spChildrenId ) ] , NULLCODE ) )
            return

        ## as a null event stays in the same species branch , we have to keep the tag of the children and iterate through the tags of our children
        TagsChild = self.RMAT.get_case_tags( bipId , spChildrenId )

        #print "NULL",bipId , spId, TagsChild

        for tag in TagsChild:
            component = ( bipId , spChildrenId , tag )
            score = self.RMAT.get_case_score( *component )
            recSolution = RecSolutionCoEv(score , [component] , NULLCODE, tag, False)
            self.RMAT.add_case(bipId,spId , recSolution )


    def CurrentCase(self, bipId, spId):
        """ clade where the event is 'leaf' (if bipId and spId are both leaves) """
        score = 0.
        if self.LeafIdToSpeciesPOAssociation[bipId] == spId:
            score = 1.
        self.RMAT.set_case(bipId,spId , RecSolutionCoEv(score , [] , LEAFCODE, 0) )

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

        #print "SL ->",Authorized,"coev in ",lossBearingSp
        if Authorized > 0:

            score = scoreExtension * scoreChild

            coevent = True
            recSolution = RecSolutionCoEv(score , components , SPECIATIONLOSSCODE, 0, coevent)

            self.RMAT.add_case(bipId,spId , recSolution )                    
            #print "adding co speciation loss (" , bipId , "," , spId , ")" , score

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


        if self.has_dead_lineage:
            if self.spPOtoDead[spId]:
                ## if we are in the dead, the duplication is a bifurcationOut and is free
                score = topologyScore * scoreChild0 * scoreChild1
                self.RMAT.add_case(bipId,spId , RecSolutionCoEv(score , components, BIFOUTCODE ) )
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


    def backtrackRecMat(self, bipId, spId , solutionTag = None , TLbefore = None, timeConsistencyBehavior = 0):
        """ 
        !recursive function!

        Takes:
            - bipId (int) : bipartition id
            - spId  (int) : species id
            - solutionIndex [ default= None ] : None or Tag (should be an int)
            - TLbefore (int or None)[ default= None ] : indicates if a TL-like event was done before and gives the id of the giver species
                                                    because, according to the rules, TL events are computed on the solutions that don't have TL events to themselves
                                                    and thus we have to exclude these from the backtrack if we just did a TL-like event from this species
            - timeConsistencyBehavior (int) [default= 0] : defines the behavior towards time timeConsistencyBehavior
                                                            0 : ignores time consistency issues (if the species tree is ultrametric or there is no transfer)
                                                            1 : detect time inconsistent scenarios and abort backtrack if they are detected
        Returns:
            (ete3.TreeNode): reconciled subtree
            or
            (None): in case a time inconsistent tree was found
        """

        ## 0. pretreatment for timeConsistencyBehavior
        if timeConsistencyBehavior >0 and self.prohibitedReceiver is None:
            self.prohibitedReceiver = set()


        ## 1. chosing a solution

        chosenSolution = self.RMAT[bipId][spId].chooseOneSolution(solutionTag , TLbefore)


        isTLlike = None
        if chosenSolution.event in [TRANSFERLOSSCODE, SPECIATIONOUTLOSSCODE]:
            isTLlike = spId

        ## 1.5 updating the list of prohibited species
        if timeConsistencyBehavior>0:
            for anc in self.spPOToNode[spId].get_ancestors(): ## prohibited species are the ancestor of the current one
                ancId = anc.PO
                if ancId in self.prohibitedReceiver:
                    break##if it is already there, it means that its ancestors are too
                else:
                    self.prohibitedReceiver.add(ancId)
            #print bipId , spId, 'prohibitedReceiver',self.prohibitedReceiver

        ## 2. backtracking that solution
        childSubTrees = []
        childSubTreesCoevents = []
        for component in chosenSolution.component:
            ## component is a tuple of size either 2 or 1
            childBipId = component[0]
            childspId = component[1]
        
            ## 2.5 checking time consistensy
            if timeConsistencyBehavior >0:
                if childspId in self.prohibitedReceiver:
                    #print "found prohibited sp:",childspId
                    return None

            childTag = None
            if len(component)>2:
                childTag = component[2]
            sub  = self.backtrackRecMat(childBipId , childspId , childTag,isTLlike , timeConsistencyBehavior )

            ## bypass for time inconsistency
            if sub is None:
                return None

            childSubTrees.append( sub )


        ##3. creating the reconciliation event associated with the reconciliation
        
        timeSlice = None #no timeslice in this version
        species = spId
        additionnal = {}
        #if self.spPOToNode[spId].name != '' :
        #    species = self.spPOToNode[spId].name

        if chosenSolution.coev:
            additionnal["coevent"] = chosenSolution.coev

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

def associateGtoSwithSep(LeafList , sep ="_"):
    assoc = {}
    for l in LeafList:
        sp , j , jbis = l.partition(sep)
        assoc[l] = sp
    return assoc



def checkReconciliationCoEventCompatibility(reconciliator , recTree , alreadySeen  ={}):
    """
    *RECURSIVE*

    Takes:
        - reconciliator (BoltzmannDTLCoevReconciliator)
        - recTree ( ETE3_based_ReconciledTree.ReconciledTree)
        - alreadySeen  (dict) [default={}] : keys are species id, 
                                                values are list (nbCoDup, nbCoTransfer, nbCoLoss)
                                                        which contains the number of co-events in the ancestors of recTree

    Returns:
        (bool) : True if recTree's co-events don't violate the ones in the reconciliator
    """

    events = recTree.getEvents()
    for i,e in enumerate(events):
        if e.isCoev():

            #determining the type of coevent
            sp = e.species
            evt = e.eventCode
            coev = e.additionnalInfo["coevent"]

            hasD = ( evt == "D" )
            hasT = ( evt[0] == "T" )
            hasL = ( evt[-1] == "L" )
            if hasT and hasL:
                if coev == 1:
                    hasL = False
                elif coev == 2:
                    hasT = False

            #print e , "->" , coev , evt , "(" , hasD , hasT, hasL,")"

            # D or L -> fairly simple
            if hasD or hasL:
                if not alreadySeen.has_key(sp):
                    alreadySeen[sp] = [0,0,0]

                nEv = None
                nAuth = None
                if hasD:
                    alreadySeen[sp][0] += 1
                    nAuth = reconciliator.getNbAuthorizedCoEvent(sp , "D")
                    nEv = alreadySeen[sp][0]

                elif hasL:
                    ## losses can be complex to assess, and will never cause problem here anyway...
                    nEv = 1
                    nAuth = 1## dummy because we don't really care...

                    #lossBearingSp = sp
                    #if evt == "SL":
                    #    ##we have to find the children that doesn't have 
                    #alreadySeen[sp][2] += 1
                    #nAuth = reconciliator.getNbAuthorizedCoEvent(sp , "L")
                    #nEv = alreadySeen[sp][2]

                if nEv > nAuth : 
                    print "more co",evt,"than authorized(",nEv,"<>",nAuth,")" 
                    return False 

            if hasT:
                ## for transfer event, we have to redefine the species as they are defined by species of arrival.
                if evt == TRANSFERLOSSCODE:## -> we go to the next event
                    sp = events[i+1].species
                else: ## simple transfer -> go look in children
                    for children in recTree.get_children():
                        sChild = children.getEvents()[0].species 
                        if sChild != sp:
                            sp = sChild
                            break

                #print "co transfer from" , e.species, "to",sp

                if not alreadySeen.has_key(sp):
                    alreadySeen[sp] = [0,0,0]

                alreadySeen[sp][1] += 1
                nAuth = reconciliator.getNbAuthorizedCoEvent(sp , "T")
                nEv = alreadySeen[sp][1]

                if nEv > nAuth :
                    print "more co",evt,"than authorized (",nEv,"<>",nAuth,")" 
                    return False 



            #print "alreadySeen" , alreadySeen

    for children in recTree.get_children():
        if not checkReconciliationCoEventCompatibility(reconciliator , children , alreadySeen.copy() ):
            return False

    return True




def formatRecTreeForXML(recTree , reconciliator):
    """
    *RECURSIVE*

    recursively performs a fe modification to prepare XML export of the recTree
    """
    i = 0
    while i < len(recTree.getEvents()):
        e = recTree.eventRecs[i]

        if e.eventCode == TRANSFERCODE:
            # all transfer are seen as So  + TL from the dead
            sp = e.species
            for children in recTree.get_children():
                sChild = children.getEvents()[0].species 
                if sChild != sp:
                    sChild
                    ##new transfer back event in the children
                    children.eventRecs.insert(0 , RecEvent( "Tb" , sChild , ts = e.timeSlice , additionnalInfo = e.additionnalInfo) )
                    break

            recTree.eventRecs[i].eventCode = SPECIATIONOUTCODE

        elif e.eventCode == TRANSFERLOSSCODE:
            sp = e.species
            spCh = recTree.eventRecs[i+1].species

            recTree.eventRecs[i].eventCode = SPECIATIONOUTLOSSCODE

            inTheDead = False

            
            if reconciliator.has_dead_lineage:
                if reconciliator.spPOtoDead[sp]: ## we are in the dead-> just replace TL by Tb
                    inTheDead = True
                    recTree.eventRecs[i].eventCode = "Tb"
                    recTree.eventRecs[i].species = spCh

            if not inTheDead:
                ##adding new Tb Event
                recTree.eventRecs.insert( i+1 , RecEvent( "Tb" , spCh , ts = e.timeSlice , additionnalInfo = e.additionnalInfo) )


        sp = e.species
        N = reconciliator.spPOToNode[sp]

        sp = N.RealNodePO        

        realN = reconciliator.spPOToNode[sp]
        if realN.name != '':
            sp = realN.name
        else:
            if not reconciliator.POtoPreTrimPO is None:
                sp = reconciliator.POtoPreTrimPO[sp]

        recTree.eventRecs[i].species = sp

        i += 1

    for children in recTree.get_children():
        formatRecTreeForXML(children , reconciliator)
    return 


if __name__ == "__main__":


    import sys

    help = """ 
    this script reconciles a gene distribution with a species tree
    by computing a score including topology, Duplication, Transfer, Loss, and co-events with a provided guide tree
    in a pseudo likelihood framework such that a more parsimonious tree will have more chance to be sampled.
    usage:
        python DTLCoevBoltzmann.py  -s species_tree -g gene_tree_distribution [options]
    options:
        -s : name of the species tree file (newick)
        -g : name of the gene tree distribution file (newick). leaf names MUST be in the format : speciesName, followed by geneName, 
                                                                            separated by a character (option --sep, '_' by default)
        -G : (optional) name of a file containing potential co-events in a guide tree 
                                    Each line is an event. It begins with "D", "T" or "L" (to indicate wether the event is a duplication,
                                    a transfer reception or a loss) and is followed by the species it is in as the 
                                    list (separated by spaces) of the leaves it is an ancestor to.
                                    example: D 1 2 3
                                            signifies a duplication in the species that is the LCA of the leaves a, b and c.

        --undated : (optional) use to specify that the species tree is not ultrametric
        -a : (optional) use to specify that the gene tree distribution file is, in fact, a .ale file (will try to detect automatically)
        -n : (optional) number of trees to sample (default: 1)

        -d : (optional) cost of a single duplication (default: 2)
        -t : (optional) cost of a single transfer (default: 3)
        -l : (optional) cost of a single loss (default: 1)
        -T : (optional) weight of the topology (default: 1)

        --temperature : (optional) a higher temperature will mean that less parsimonious trees will have a higher chance to be sampled 
                                                                                                                        (default: 0.1)
        --de : (optional)  cost of the extension of a co-duplication (default: 1)
        --le : (optional)  cost of the extension of a co-loss (default: 0.5)
        --te : (optional)  cost of the extension of a co-transfer (default: 1.5)

        -o : (optional) name of the file to print the generated tree in (default is stdout)
        --sep : (optional) separator between species name and gene name in the gene tree leaves (default '_')
        --seed : (optional) seed to use in the randomisation (by default, time() is used)
        --topoOnly : (optional) if present, only topology (without reconciliation annotation) will be outputted
        --xml : (optional) if present, the tree will be written in an XML format

        --ignore.time.consistency : (optional) when '--undated' is used, some generated scenarios will be time inconsistent 
                                                    (ie. transfer to an ancestor). By default, these time inconsistent scenarios will
                                                    be detected and omitted from the results. Use this option to keep them.
        --max.inconsistency.sample.overhead : (optional) means we won't sample more than ten time than -n in order to find some 
                                                        time consistent scenarios. (default: 10)


        --ignore.coevents.consistency : (optional) in the presence of transfer, some solutions might inlcude cases where a guide event
                                                 is linked to a node AND one of its ancestor, which is not a desired result. 
                                                 By default, such scenarios are detected when sampling, and eliminated. 
                                                 Activate this option if you want to include them.

        --trim.sp.tree : (optional) trim the species tree above the LCA of the species which have a leaf associated to them


    """
        #--guide : (optional) name of the guide tree file (recPhyloXML)

    OK = True

    nextKEY = None
    params = {
        "-s" : None,#name of the species tree file (newick)
        "-g" : None,#name of the gene tree distribution file (newick)
        "-G" : None,#(optional) name of a file containing potential co-events in a guide tree 

        "--undated" : False, #(optional) use to specify that the species tree is not ultrametric
        "-a" : False,#(optional) use to specify that the gene tree distribution file is, in fact, a .ale file
        "-n" : 1, #(optional) number of trees to sample (default: 1)

        "-d" : 2. ,#(optional) cost of a single duplication (default: 2)
        "-t" : 3. ,#(optional) cost of a single transfer (default: 3)
        "-l" : 1.,#(optional) cost of a single loss (default: 1)
        "-T" : 1.,#(optional) weight of the topology (default: 1)

        "--temperature" : 0.1,#(optional) a higher temperature will mean that less parsimonious trees will have a higher chance to be sampled (default: 1)
        "--de" : 1., #(optional)  cost of the extension of a co-duplication (default: 1)
        "--le" : 0.5, #(optional)  cost of the extension of a co-loss (default: 0.5)
        "--te" : 1.5, #(optional)  cost of the extension of a co-transfer reception (default: 1.5)

        "-o" : None, #(optional) name of the file to print the generated tree in (default is stdout)
        "--sep" : '_',#(optional) separator between species name and gene name in the gene tree leaves (default '_')
        "--seed": None, #(optional) seed to use in the randomisation (by default, time() is used)
        "--topoOnly" : False, #(optional) if True, only topology (without reconciliation annotation) will be outputed
        "--xml" : False, #(optional), if True, the tree will be written in an XML format.
        
        "--ignore.time.consistency" : False, #(optional) when '--undated' is used, some generated scenarios will be time inconsistent (ie. transfer to an ancestor).
                                            #    By default, these time inconsistent scenarios will be detected and omitted from the results. Use this option to keep them.
        "--max.inconsistency.sample.overhead" : 10., #(optional) means we won't sample more than ten time than -n in order to find some time consistent scenarii
                                                    # only useful when --undated and not --ignore.time.consistency

        "--ignore.coevents.consistency" : False, #(optional) in the presence of transfer, some solutions might inlcude cases where a guide event is linked to a node AND one of its ancestor
                                                # which is not a desired result. As such, such scenarios are detected when sampling, and eliminated. activate this option if you want to include them
        
        "--trim.sp.tree" : False ,#(optional) trim the species tree above the LCA of the species which have a leaf associated to them

        "--UNTIL" : False, #(hidden) : if True, the sampling will continue until the desired number of correct reconciliation is achieved. This has the potential for a (rare) infinite loop
        "--inconsistent.scenarios.when.failure" : False
    }

    flagArgs = ["-a", "--topoOnly", "--xml", "--undated", "--ignore.time.consistency", "--ignore.coevents.consistency", "--trim.sp.tree", "--UNTIL","--inconsistent.scenarios.when.failure"]

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
        for pname in ["-d","-l","-t","-T","--de","--le","--te","--temperature", "--max.inconsistency.sample.overhead"]:
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
            "SPE_EXTENSION"  : 1.,

            TRANSFERCODE : np.exp(- params["-t"] / Temperature),

            "TRANSFER_EXTENSION" : np.exp(- params["--te"] / Temperature),

            SPECIATIONOUTCODE : 1.,## speciation toward the dead
            }

    ### association between genes and species
    
    geneNameToSpeciesNameAssociation = {}

    Rec = BoltzmannDTLCoevReconciliator(distrFilename , spFilename , geneNameToSpeciesNameAssociation ,
                                                 COSTS, params['--sep'], params['-a'] , not params["--undated"] , params["--trim.sp.tree"] )

    print "trees read and bipartitions computed"
    ##adding some guide events
    if not params["-G"] is None:
        coevents = readGuideEvents(params["-G"], Rec.spTree, Rec.spNameToPO)
        #print "coev", coevents
        for e in coevents:
            Rec.addGuideCoEvents( spId = e[0] , eventCode= e[1] , eventNb=1)

    print "potential co-events added"

    Rec.setupRecMat()
    Rec.fillRecMat()

    print "reconciliation matrix computed"
    

    #for bipId in Rec.sortedBip:
    #        #print "treating bip" , bipId ,"->", self.ccp_dist.dbip_set[bipId]
    #        if Rec.TStoPO is None: ## species tree is not subdivided
    #            for spId in Rec.sortedSpId:
    #                print "c(",bipId,",",spId,") =", Rec.RMAT[bipId][spId].totalScore
    #        else: ##species tree is subdivided
    #            nbTS = len(Rec.TStoPO)
    #            if Rec.has_dead_lineage:
    #                nbTS -= 1 ##when there is a dead lineage, the root is artificial. 
    #
    #            for TS in range(nbTS):
    #                for spId in Rec.TStoPO[TS]:
    #                    
    #                    print "c(",bipId,",",spId,") =", Rec.RMAT[bipId][spId].totalScore
    #                    for sol in Rec.RMAT[bipId][spId].getAllSolutions():
    #                        print "\t",sol



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




    ## gestion of the time consistency
    timeConsistencyBehavior = 0
    if params["--undated"] and not params["--ignore.time.consistency"]:
        timeConsistencyBehavior = 1



    rootBip = Rec.ccp_dist.get_root_bip_id() 

    totalAtRoot = 0.
    RootP = []
    RootList = []
    #for sp,dsol in Rec.RMAT[rootBip].items():
    for j,dsol in enumerate(Rec.RMAT[rootBip]):
        sp = Rec.sortedSpId[j]
        RootList.append(sp)
        RootP.append(dsol.totalScore)
        totalAtRoot += dsol.totalScore
    

    #print RootList
    #print RootP

    OUT = sys.stdout
    if not params["-o"] is None:
        OUT = open(params["-o"],"w")

    speciesNames = {v:k for k,v in Rec.spNameToPO.items()}

    if totalAtRoot == 0:
        print "no non impossible solution found"
    else:
        print "sampling" ,params["-n"], "reconciliations"
        if not params["--xml"]:
            print "format newick nodes: name|evt0Species.evt0Code|evt1Species.evt1Code ... "
        #print "internal node name are their bipartition name in the CCP distribution."
        nbValid = 0
        nbTentative = 0
        maxNbTentative = int( params["--max.inconsistency.sample.overhead"] * params["-n"] )
        while nbValid < params["-n"]:

            if not params["--UNTIL"]:
                if nbTentative == maxNbTentative:
                    print "sampled more than", params["--max.inconsistency.sample.overhead"], "times the sample size"
                    print "to find time consistent or coevent consistent scenarios. Aborting to avoid a potential infinite loop."
                    print "Valid results are still reported."
                    if params["--inconsistent.scenarios.when.failure"]:
                        print "to ensure the correct number of reconciliation in the output, time inconsistent solutions will be added."
                        print "the ouput will contain",nbValid,"time consistent reconciliations"
                        print "and",params["-n"] - nbValid,"reconciliations whose time-consistency have not been tested."
                        params["--ignore.coevents.consistency"] = True
                        timeConsistencyBehavior = 0
                        nbTentative = 0
                    else:
                        break


            r = np.random.random() * totalAtRoot
            i = -1
            while r > 0:
                i += 1
                r -= RootP[i]
            sp = RootList[i]

            if timeConsistencyBehavior>0:
                Rec.prohibitedReceiver = None

            RT = Rec.backtrackRecMat(rootBip, sp, timeConsistencyBehavior = timeConsistencyBehavior)

            nbTentative+=1
            if RT is None:
                ##time inconsistency detected
                continue
                        

            if not params["--ignore.coevents.consistency"]:## we have to check co-event coherence
                if not checkReconciliationCoEventCompatibility(Rec , RT , {}):
                    continue ## 
            

            nbValid+=1

            ### uncommenting this will add a series of counts for each event at the beginning of the line
            #DE = RT.countEvents(True)
            #OUT.write( ",".join([str(DE.get(e,0)) for e in ["D","SL","T","TL"]]) + " " )
            
            #print RT.getTreeNewick(topoOnly=params['--topoOnly']), "->",
            formatRecTreeForXML(RT , Rec)
            #print RT.getTreeNewick(topoOnly=params['--topoOnly'])


            if not params["--xml"]:
                OUT.write( RT.getTreeNewick(topoOnly=params['--topoOnly']) + "\n" )
            else:
                OUT.write( RT.getTreeRecPhyloXML( speciesNames = speciesNames, topoOnly=params['--topoOnly']) + "\n" ) 

        if nbTentative != nbValid:
            print "backtracked",nbTentative,"times to get",nbValid,"trees"
    OUT.close()