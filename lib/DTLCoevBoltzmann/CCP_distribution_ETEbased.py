#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################################
##  Author:         Wandrille Duchemin  #
##  Created:        23-Jul-2015         #
##  Last modified:  13-Jan-2017         #
#########################################

from ete3 import Tree,TreeNode
import random



class CCP_distribution:
    """
    Distribution of CCP
    This class is linked to the ete3.TreeNode class for tree representation purposes

    Attributes:
        - dleaf_id (dict): keys are leaves ids (int), values are leaves names (str)
        - nb_observation (int): number of trees composing the distribution
        - dbip_count (dict): keys are bipartition ids, values are biparition count (number of times the bipartition was seen in the tree distribution)
        - dbip_bls (dict): keys are bipartition ids, values are the sum of the length of the branch leading to the bipatition 
        - dbip_set (dict): keys are bipartition ids, values are sets containing the  id of the leaves composing the bipartition
        - constructor_string (str): string describing the modelled tree in newick format
        - ddip_count (dict): 
                keys are bipartition ids, 
                values are (dict):
                    with keys (tuple): two clades
                    values (int): count of the observation where the bipartition resolve in the two clades
    """
    def __init__(self):
        
        self.dleaf_id = {}
        self.nb_observation = 0  #number of trees composing the distribution
        self.dbip_count = {} #keys are bipartition ids, values are biparition count (number of times the bipartition was seen in the tree distribution)
        self.dbip_bls = {} #keys are bipartition ids, values are the sum of the length of the branch leading to the bipatition 
        self.dbip_set = {} #keys are bipartition ids, values are sets containing the  id of the leaves composing the bipartition
        self.constructor_string = ""
        self.ddip_count = {}
        return

    def is_leaf(self, bipId):
        """ return True is the clade bipId is of size 1 """
        return len(self.dbip_set[bipId]) == 1

    def get_leaf_name(self, bipId):
        """ returns the name of the leaf is the clade bipId is of size 1 , returns None otherwise"""
        if self.is_leaf(bipId):
            return self.dleaf_id[ tuple( self.dbip_set[bipId] )[0] ]
        return None

    def read_from_ale_handle(self,filehandle):
        """
        Takes:
            - filehandle (file): reading handle of a file containing CCPs information in the .ale format
            
        Returns:
            0 -> no problem
            1 -> problem in the CCP file format
        """
        
        header_lines = [
        "#constructor_string",
        "#observations",
        "#Bip_counts",
        "#Bip_bls",
        "#Dip_counts",
        "#last_leafset_id",
        "#leaf-id",
        "#set-id",
        "#END"
        ]
        current_header_index = 0
        
        def check_header(l,header_lines,current_header_index):
            if l.strip() != header_lines[current_header_index]:
                return False
            return True
        
        l = filehandle.readline()
        
        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1
            
        
        self.constructor_string = filehandle.readline().strip()
        
        ####
        
        l = filehandle.readline()

        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1


        self.nb_observation = int(filehandle.readline().strip())
        
        
        ####
        
        l = filehandle.readline()##  "#Bip_counts\n"

        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1


        l = filehandle.readline()
        
        while not l.startswith("#"): ## format: "bip_id bip_count\n"
            sl = [i.strip() for i in l.split()]
            
            self.dbip_count[int(sl[0])] = int(sl[1])
        
            l = filehandle.readline()
            
        ## l is supposed to be: "#Bip_bls"
        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1



        l = filehandle.readline()
        
        while not l.startswith("#"): ## format: "bip_id sum_blen\n"
            sl = [i.strip() for i in l.split()]
            
            self.dbip_bls[int(sl[0])] = float(sl[1])
        
            l = filehandle.readline()
        
        
        ## l is supposed to be: "#Dip_counts"
        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1

        l = filehandle.readline()
        while not l.startswith("#"): ## format: "bip_id bip_id1 bip_id2 count\n"
            sl = [i.strip() for i in l.split()]
            
            if not self.ddip_count.has_key(int(sl[0])):
                self.ddip_count[int(sl[0])] = {}
        
            self.ddip_count[int(sl[0])][(int(sl[1]),int(sl[2]))] = int(sl[3])
        
            l = filehandle.readline()
        
        
        ## l is supposed to be: "#last_leafset_id"
        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index],l
            return 1
        current_header_index+=1

        l = filehandle.readline()
        
        l = filehandle.readline()
        ## l is supposed to be: "#leaf-id"
        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1
        
        l = filehandle.readline()
        while not l.startswith("#"): ## format: "leaf_name leaf_id\n"
            sl = [i.strip() for i in l.split()]
            
            self.dleaf_id[int(sl[1])] = sl[0]

            l = filehandle.readline()

        
        ## l is supposed to be: "#set-id"
        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1

        l = filehandle.readline()
        while not l.startswith("#"): ## format: "bip_id : leaf_id1 leaf_id2 ...\n"
            sl = [i.strip() for i in l.split()]
            
            self.dbip_set[int(sl[0])] = set([int(i) for i in sl[2:]])
            
            l = filehandle.readline()

        ## l is supposed to be: "#END"
        if not check_header(l,header_lines,current_header_index):
            print "error in the CCP file format",header_lines[current_header_index]
            return 1
        current_header_index+=1
        
        ## checking if leaves were absent of the .ale ...
        for k in self.dbip_set.keys():
            if not self.dbip_count.has_key(k):
                self.dbip_count[k] = self.nb_observation

        return 0


    def set_root_split(self):
        """ setup a root split that contains all leaves """

        root_split = self.get_total_split()

        root_splitId = -1
        self.dbip_set[root_splitId] = set(self.dleaf_id.keys())
        self.ddip_count[root_splitId] = root_split
        self.dbip_count[root_splitId] = sum( root_split.values() )

    def get_root_bip_id(self):
        if not self.dbip_count.has_key( -1 ):
            print "ERROR : asked for root bip while it is not set : returning None"
            return None
        return -1


    def number_of_leaves(self):
        """
        Returns:
            (int): number of leaves of the tree
        """
        return len(self.dleaf_id)

    def get_bip_from_leafset(self,leafset,default = None):
        """
        Takes:
            - leafset (set): set of leaves id
            - default (None): default value to return if there is no such bipartition in the tree
        
        Returns:
            (int): Bip id that has leafset for leaves
        """
        for k in self.dbip_set.keys():
            if leafset == self.dbip_set[k]:
                return k
        return default
    
    def add_bip(self,leafset):
        """
        Takes:
            leafset (set): set of leaf associated with the bip
            
        Returns:
            (int): bip id
        """
        bip_id = self.get_bip_from_leafset(leafset)
        if bip_id is None:##the clade is unknown
        
            bip_id = len(self.dbip_count) + 1 ##new id
            
            self.dbip_count[bip_id] = 0
            self.dbip_bls[bip_id] = 0.
            self.dbip_set[bip_id] = leafset
        

        return bip_id

    def get_leaf_id(self,leaf_name):
        """
        Adds the leaf if not present
        
        Takes:
            leaf_name (str): name of the leaf
            
        Returns:
            (int): id of the leaf
        """
        
        if not leaf_name in self.dleaf_id.values():
            new_id = len(self.dleaf_id) + 1
            self.dleaf_id[new_id] = leaf_name
            return new_id
        
        rev = {v:k for k,v in self.dleaf_id.items()}
        
        return rev[leaf_name]

    def add_tree_branch_to_distribution(self,node):
        """
        add the clades and bipartitions of a branch to the distribution
        
        Takes:
            - node (ete3.TreeNode): Node in a tree

        Returns:
            (None)
        """
        
        BIPS = []
        SUBBIPS = []
        LENGTH = 0
        
        isSubRoot = False
        i = 0
        for n in node.iter_ancestors():
            if i > 0 :
                break
            if n.is_root():
                isSubRoot = True
                break
            i += 1

        if isSubRoot:
            return ##if this is the root of the tree, do nothing as CCPS are computed on unrooted trees
        elif node.is_root():
            
            for subroot in node.get_children():
                
                SUBBIPS.append([])
                
                LENGTH += subroot.dist##adding the length of both subroots
                
                if subroot.is_leaf():
                    BIPS.append(self.add_bip(set([self.get_leaf_id( subroot.name )])))##adding the bip
                    continue
                
                for c in subroot.get_children():
                    SUBBIPS[-1].append(set([self.get_leaf_id(i) for i in c.get_leaf_names()]))
                
                
                BIPS.append(self.add_bip(set.union(*SUBBIPS[-1]))) ##adding the bipartition
            
            
        else:#node is neither root nor subroot

            ##getting first bip: what is under the node

            LENGTH = node.dist##getting the length

            SUBBIPS.append([])

            if node.is_leaf():
                BIPS.append(self.add_bip(set([self.get_leaf_id(node.name)])))##adding the bip

            else:
                for c in node.get_children():
                    SUBBIPS[-1].append(set([self.get_leaf_id(i) for i in c.get_leaf_names()]))
            
            
                BIPS.append(self.add_bip(set.union(*SUBBIPS[-1])))##adding the bip
            ##getting second Bip: what is above the node
 
            #starting with the brother

            #bro = node.go_brother()
            bro = node.get_sisters()[0]
            
            SUBBIPS.append([])
            SUBBIPS[-1].append(set([self.get_leaf_id(i) for i in bro.get_leaf_names()]))
            
            rest_of_leaves = set(self.dleaf_id.keys())
            
            rest_of_leaves.difference_update(self.dbip_set[BIPS[0]])
            rest_of_leaves.difference_update(SUBBIPS[-1][-1])
            
            SUBBIPS[-1].append(rest_of_leaves)
            
            BIPS.append(self.add_bip(set.union(*SUBBIPS[-1])))
        
        for i,bip in enumerate(BIPS):
            
            self.dbip_count[bip] += 1 ##adding the bip count
            self.dbip_bls[bip] += LENGTH ##adding the bip length

            if len(SUBBIPS[i]) !=0: ##the bip is not a leaf
                if not self.ddip_count.has_key(bip):
                    self.ddip_count[bip] = {}
        
                dip = [self.add_bip(s) for s in SUBBIPS[i]]
                dip.sort()
                
                dip_id = tuple(dip)
            
                if not self.ddip_count[bip].has_key(dip_id):
                    self.ddip_count[bip][dip_id] = 0
            
                if bip in dip_id:
                
                    print "PROBLEM"
                    print branch,"-> adding dip", bip,"->",dip_id            
            
                self.ddip_count[bip][dip_id] += 1
        
    def add_tree_to_distribution(self,tree):
        """
        Add the bipartition of a tree to the CCP distribution
        
        Takes:
            - tree (ete3.Tree): phylogenetic tree
            
        """
        
        if len(tree.children) == 3: 
            ## special unrroted case where the tree begin by a trifurcation ...
            ## we artificially remove the trifurcation to avoid future problems
            a = TreeNode()
            b = tree.children[1]
            c = tree.children[2]
            b.detach()
            c.detach()
            tree.add_child(a)
            a.add_child(b)
            a.add_child(c)
            #print " special rerooting "

        for i in tree.traverse():
            if len(i.children) > 2:
                print "multifurcation detected! Please provide bifurcating trees."
                print "exiting now"
                exit(1)



        if self.nb_observation == 0:##no tree has been observed yet: add all the leaves
            for l in tree.get_leaf_names():
                self.get_leaf_id(l)##adds the leaves to the CCP
        
        
        for node in tree.traverse("postorder"): ##for each branch of the tree
            self.add_tree_branch_to_distribution(node)
        
        self.nb_observation += 1

        return
    
    def read_from_treelist_handle(self,filehandle):
        """
        Takes:
            - filehandle (file): reading handle of a file containing trees (1 per line, newick)
            
        Returns:
            0 -> no problem
            1 -> error in the format
        
        """

        for l in filehandle.readlines():
            try:
                tree = Tree(l.strip())#Node(nwk=l.strip())
            
            except:
                return 1
            
            self.add_tree_to_distribution(tree)
        
            if self.constructor_string == "":
                self.constructor_string = l.strip()

        return 0

    def __str__(self):
        s = "leaves id" + "\n"
        s += str(self.dleaf_id) + "\n"
        
        s += "clades" + "\n"
        s += str(self.dbip_set) + "\n"
        
        s += "clades split\n"
        s += str(self.ddip_count) + "\n"
        
        return s

    def write_ale_file(self,filehandle):
        """
        Writes the CCPs in the same format as ALEobserve
        
        Takes:
            - filehandle (file): writing handle
        """
        
        header_lines = [
        "#constructor_string",
        "#observations",
        "#Bip_counts",
        "#Bip_bls",
        "#Dip_counts",
        "#last_leafset_id",
        "#leaf-id",
        "#set-id",
        "#END"
        ]
        cur_header = 0
        
        
        ##constructor string
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1
        
        filehandle.write(self.constructor_string + "\n")
        
        ##nb observations
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1
        
        filehandle.write(str(self.nb_observation) + "\n")
        
        ## count of bipartition
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1
        
        for k in self.dbip_count.keys():
            #print self.dbip_set[k], "->" , len(self.dbip_set[k])
            if len(self.dbip_set[k])> 1:
                #print "plop", str(k) + "\t" + str(self.dbip_count[k]) 
                filehandle.write(str(k) + "\t" + str(self.dbip_count[k]) + "\n")
        
        
        ## bipartition br length
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1

        for k in self.dbip_bls.keys():
            filehandle.write(str(k) + "\t" + str(self.dbip_bls[k]) + "\n")

        ##dip count
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1

        for k in self.ddip_count.keys():
            for k2 in self.ddip_count[k].keys():
                filehandle.write( str(k) + "\t" + str(k2[0]) + "\t" + str(k2[1]) + "\t" + str(self.ddip_count[k][k2]) + "\n")
        
        
        ##max leafset id
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1
    
        filehandle.write(str(max(self.dbip_count.keys())) + "\n")
        
        ##leaf-id
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1
    
        for k in self.dleaf_id.keys():
            filehandle.write(str(self.dleaf_id[k]) + "\t" + str(k) + "\n")
            
        ##set-id
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1
    
        for k in self.dbip_set.keys():
            filehandle.write(str(k) + "\t:\t" + "\t".join([str(i) for i in self.dbip_set[k]]) + "\n")
        
        ##END
        filehandle.write(header_lines[cur_header] + "\n")
        cur_header += 1
    
        return
    
    def get_total_split(self):
        """
        return the bipartitions that concern the whole leafset

        Returns:
            (dict): keys are tuples of bip id, values are associated counts
        """
        
        root_split = {}
        
        bips = self.dbip_count.keys()
        
        ##to each bip correspond a complementary bip that form the all leafset
        
        nb_leaves = self.number_of_leaves()
        
        while len(bips) >0:
        
            if len(bips) == 1:
                print self.dbip_set[ bips[0] ]
                print "there is a problem in your tree(s) (n-furcation at the root maybe?)."
                exit(1)
            for i in range(len(bips)):
                for j in range(i+1,len(bips)):
                    if len(self.dbip_set[bips[i]]) + len(self.dbip_set[bips[j]]) == nb_leaves:
                        if self.dbip_set[bips[i]].isdisjoint(self.dbip_set[bips[j]]):##found the complementary
                            root_split[(min(bips[i],bips[j]),max(bips[i],bips[j]))] = self.dbip_count[bips[j]]
                            
                            bips.pop(j)
                            bips.pop(i)
                            break
                break
        
        return root_split
    
    def number_of_amalgamable_tree_in_bip(self,bip= None):
        """
        Takes:
            - bip (int): id of the bipartition to get the number from
            
        Returns:
            (int): number of amalgamable subtree given the bipartition bip
        """
        
        if bip == None:##"root", actually an arbitrary leaf
            for b in self.dbip_set.keys():
                if len(self.dbip_set[b]) == self.number_of_leaves() -1 : ##all leaves but one type of clade
                    bip = b


        if len(self.dbip_set[bip])  == 1:##leaf
            return 1
        
        nb_tree = 0
        
        for child_bips in self.ddip_count[bip].keys():

            nb_st1 = self.number_of_amalgamable_tree_in_bip(child_bips[0])
            nb_st2 = self.number_of_amalgamable_tree_in_bip(child_bips[1])
            
            nb_tree += (nb_st1 * nb_st2)
    
        #print bip,len(self.ddip_count[bip]),nb_tree
    
        return nb_tree

    
    def get_ccp(self,bip,child_bips,default_proba = 0.):
        """
        Takes:
            - bip (int): set id
            - child_bips (tuple): tuple of two bip ids
            - default_proba (float) (default: 0.): default probability to give to a split absent from the ccp distribution
        
        Returns:
            (float): conditionnal clade probability associated with that split
        """
        
        if not self.dbip_count.has_key(bip):
            return default_proba
        
        total = 1. * self.dbip_count[bip]
        
        if not self.ddip_count.get(bip,{}).has_key(child_bips):
            return default_proba
        
        return 1. * self.ddip_count[bip][child_bips] / total

    def get_all_ccps(self,bip):
        """
        Takes:
            - bip (int): set id
        
        Returns:
            (dict): keys are tuples of bip ids, values are corresponding conditionnal clade probability
        """
        
        ccps = {}
        
        for child_bips in self.ddip_count.get(bip,{}).keys():
            ccps[child_bips] = self.get_ccp(bip,child_bips)
        
        return ccps
    
    def get_max_ccp(self,bip):
        """
        Takes:
            - bip (int): set id
        
        Returns:
            (tuple): tuple of ids corresponding to the more frequent observed split of bip
            (None) if bip isn't split (i.e. bip isn't observed or is only a leaf)
            
            Split with the same probability will be chosen randomly
        """
        
        
        ccps = self.get_all_ccps(bip)
        
        
        #print "bip:",bip,"ccps:",ccps
        
        
        if len(ccps) == 0:
            return None
        
        max_ccp = 0.
        lmax_id = []
        
        for k,ccp in ccps.items():
            if ccp > max_ccp:
                max_ccp = ccp
                lmax_id = []
                lmax_id.append(k)
            elif ccp == max_ccp:
                lmax_id.append(k)
            
        return random.choice(lmax_id)

    def sample_ccp(self,bip):
        """
        Takes:
            - bip (int): set id
        
        Returns:
            (tuple): tuple of ids corresponding an observed split of bip
            (None) if bip isn't split (i.e. bip isn't observed or is only a leaf)
            
            Splits are chosen with a probability equals to their ccp
        """
        
        
        
        ccps = self.get_all_ccps(bip)
        
#        print bip,"->",ccps
        
        
        if len(ccps) == 0:
            return None
        
        res = random.random()
        
        for k in ccps.keys():
            if ccps[k] > res:
                return k
            
            res -= ccps[k]
        
        return None ##should never go there


    def get_tree_from_CCP(self,method,bip=None, node=None):
        """
        RECURSIVE
        build a tree from the CCP distribution
        
        Takes:
            - method (function) : function that takes a bipartition id and returns a tuple of children bipartition ids
            - bip (int): bip id
            - node (ete3.TreeNode): current tree
            
        Returns:
            (ete3.TreeNode): phylogenetic tree drawn from the CCP distribution
        """
        
        DIP = []
        BLEN = []
        
        if node == None:##nothing is created yet
            root_bip = None ##bip that contains all leaves but one
            leaf_bip = None ##bip that contains only 1 leaf
            
            for bip in self.dbip_set.keys():
                if len(self.dbip_set[bip]) == self.number_of_leaves() -1 : ##all leaves but one type of clade
                    root_bip = bip
                    
                    one_leaf_set = set(self.dleaf_id.keys())  - self.dbip_set[bip] ##complementary leaf set. Only one leaf
                    
                    leaf_bip = self.get_bip_from_leafset(one_leaf_set)
                    
                    break
            
            leaf_bip_count = self.dbip_count[leaf_bip]
            leaf_bip_blen = self.dbip_bls[leaf_bip] /leaf_bip_count
            
            ##as this is the root, this length will be divided between both of the root children.
            DIP = [leaf_bip,root_bip]
            BLEN = [leaf_bip_blen/2.,leaf_bip_blen/2.]
            
            ##... creating the root node
            node = Tree()
        
        
        elif len(self.dbip_set[bip]) == 1: ##the current bip is a leaf
        
            leaf_id = [i for i in self.dbip_set[bip]][0]
            leaf_name = self.dleaf_id[leaf_id]
            
            node.name = leaf_name
            
            return node
        
        else: ##bipartition node that is not the root: draw bipatition using the method function
            
            DIP = method(bip) ##choosing a split of the clade
            
            for d in DIP:
                BLEN.append(self.dbip_bls[d] *1. / self.dbip_count[d])
            
            
        ##for each new clade in dip, we create a child node
        
        for i,d in enumerate(DIP):
            new = TreeNode(dist=BLEN[i])
            #new = node.newnode()##creating new node

            node.add_child(new)
            #node.link_child(new,newlen=BLEN[i])##linking it as child and giving it its length
            
            self.get_tree_from_CCP(method,d, new)#RECURSION
        
        return node
    
    def get_tree_likelihood(self,tree,default_proba = 10.**-6):
        """
        RECURSIVE
        
        Takes:
            - tree (tree2.Node): a phylogenetic tree
            - default_proba (float) (default: 10.**-6): default probability to give to a split absent from the ccp distribution
        
        Returns:
            (float): likelihood of the tree according to the ccp distribution
        """

        prob = 1.

        ##end of the recursion
        if tree.is_leaf():
            
            return prob

        SUBBIPS = []
        SUBBIPS_id = []

        for c in tree.get_children():
            SUBBIPS.append(set([self.get_leaf_id(i) for i in c.get_leaf_names()]))
            SUBBIPS_id.append(self.add_bip(SUBBIPS[-1]))

            prob *= self.get_tree_likelihood(c,default_proba) ##recursion on the children

        if not tree.is_root():
            BIP = self.add_bip(set.union(*SUBBIPS))

            SUBBIPS_id.sort()
            
            prob *= self.get_ccp(BIP,tuple(SUBBIPS_id),default_proba)



        else:
            ##actual root of the whole tree: we have to give the probability of that first split
            c = self.dbip_count.get(SUBBIPS_id[0], self.nb_observation * default_proba)
            prob *= c *1./self.nb_observation
        
        
        
        return prob

if __name__ == "__main__":
    
    import sys



    help_str = """
    
    a simple CCP manipulation software.
    programs and options:
        observe             -> Creates a CCP distribution from a set of trees. Requires a file containing a list of trees (newick format) and the name of an output file. Can take an additionnal integer as burnout.
        ml                  -> Creates the Maximum Likelihood tree of a CCP distribution. Requires a CCP distribution file (created with the observe option for example) and the name of an output file.
        sample              -> Creates a sample of trees by drawing them from the CCP distribution. Requires a CCP distribution file (created with the observe option for example) and the name of an output file as well as an integer specifying how many trees to draw.
        count_amalgamable   -> Print the number of amalgamable trees in a CCP distribution. Requires a CCP distribution file (created with the observe option for example).
        proba               -> Compute the likelihoods of a set of trees. Requires a CCP distribution file (created with the observe option for example), a file containing a list of trees (newick format) and the name of an output file.
    
    """


    if len(sys.argv) < 2:
        print help_str
        exit(1)

    prg_type = sys.argv[1]

    if prg_type == "observe":
        
        if len(sys.argv) < 4:
            print help_str
            exit(1)


        treelist_filename = sys.argv[2]
        ale_filename = sys.argv[3]
        burnin = 0
        if len(sys.argv)>4:
            burnin = int(sys.argv[4])
        
        ccp_dist = CCP_distribution()
        
        IN = open(treelist_filename,"r")
        
        for i in range(burnin):
            IN.readline() ##discarding the firsts lines of the treelist
        
        
        ccp_dist.read_from_treelist_handle(IN)
        
        IN.close()
        
        OUT = open(ale_filename,"w")
        
        ccp_dist.write_ale_file(OUT)
        
        OUT.close()
        
#        print ccp_dist


    elif prg_type in ["ml","sample"]:
        
        if len(sys.argv) < (4 + int(prg_type == "sample")):
            print help_str
            exit(1)

        ale_filename = sys.argv[2]
        out_file = sys.argv[3]
        
        ccp_dist = CCP_distribution()
        
        IN = open(ale_filename,"r")
        
        if ccp_dist.read_from_ale_handle(IN) == 1:
            print "error while reading ale file"
            exit(1)
        
        IN.close()

        #print ccp_dist
        #print "plop"
        #print ccp_dist.dleaf_id
        #print ccp_dist.dbip_set
        #print ccp_dist.dbip_count
        
        nb_tree = 1
        method = ccp_dist.get_max_ccp
        
        if prg_type == "sample":
            nb_tree = int(sys.argv[4])
            method = ccp_dist.sample_ccp
        
        OUT = open(out_file,"w")
        
        
#        print ccp_dist
        
        for i in range(nb_tree):
            t = ccp_dist.get_tree_from_CCP(method)

            
            
            OUT.write(t.write() + "\n")
            
        OUT.close()
    
    elif prg_type == "count_amalgamable":
        
        if len(sys.argv) < 3:
            print help_str
            exit(1)

        
        ale_filename = sys.argv[2]
        
        ccp_dist = CCP_distribution()
        
        IN = open(ale_filename,"r")
        
        if ccp_dist.read_from_ale_handle(IN) == 1:
            print "error while reading ale file"
            exit(1)
        
        IN.close()

        print "number of amalgamable trees:", ccp_dist.number_of_amalgamable_tree_in_bip()
    
    elif prg_type == "proba":
        help_proba_str= """
        Takes:
            - ale_filename: name of an ale file
            - tree_filename: name of a file containing newick trees (1 per line)
            - out_filename: name of a file to write in
            
        For each tree in tree_filename, the program computes its likelihood according to the CCPs in ale_filename and outputs the results in a line in out_filename
        """
        
        if len(sys.argv)<5:
            print help_proba_str
            exit(0)
        
        ale_filename = sys.argv[2]
        tree_filename = sys.argv[3]
        out_filename = sys.argv[4]

        IN = open(ale_filename,"r")
        
        ccp_dist = CCP_distribution()
        
        if ccp_dist.read_from_ale_handle(IN) == 1:
            print "error while reading ale file"
            exit(1)
        
        IN.close()

        IN = open(tree_filename,"r")
        OUT = open(out_filename,"w")
        
        for l in IN.readlines():
            
            tree = Tree(l.strip())
            
            
            if len(tree.children) >2: ##the root is not a bifurcation
                ##creating a bifurcation; this does not matter for CCP computation as rooting doesn't change anything
                print "please give a bifurcating, rooted tree (the rooting won't change the likelihood)."
#                new = tree.newnode()
#                
#                
#                children = [tree.go_child(i) for i in range(2)]
#                for ch in children:
#                    tree.rm_child(ch)
#                    new.link_child(ch,newlen=ch.lg())
#                
#                tree.link_child(new,newlen=0.)
            
#            print tree.arborescence_ASCII()

            proba = ccp_dist.get_tree_likelihood(tree)

            OUT.write(str(proba) + "\n")
        
        IN.close()
        OUT.close()
        
