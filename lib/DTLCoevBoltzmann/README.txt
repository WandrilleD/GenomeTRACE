D(T)LRecCoevBoltzmann is a python software made to reconciles a gene distribution with a species tree
 by computing a score including topology, Duplication, Transfer, Loss, and co-events with a provided guide tree
 in a pseudo likelihood framework such that a more parsimonious tree will have more chance to be sampled.

For more information on the reconciliation algorithm between a gene distribution and a species tree, <TERA algo>
For more information on the pseudo likelihood framework, see an application on a different dynamic programming example : <DeClone>
The XML output mentioned in the options is described here : http://phylariane.univ-lyon1.fr/recphyloxml/

In order for the executable to work properly, you need to have installed the libraries:
 - numpy (numpy.org)
 - ete3 (http://etetoolkit.org/)


This folder contains several python script that are used by two main executables:
DLCoevBoltzmann.py and DTLCoevBoltzmann.py

both have very similar functionnality and behavior:
the first one reconciles using only duplication and loss events,
the second one also incorporates transfers.

To access the help of each software, just execute them without any argument:

./DLCoevBoltzmann.py

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





./DTLCoevBoltzmann.py

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




#### SMALL TUTORIAL ####

## DL reconciliation ##


the folder testdata/ contains some data to test the software.

For instance, use the command:
./DLCoevBoltzmann.py -s testdata/sptree.txt -g testdata/toy_chain.treelist10 -o testdata/DLoutput

to sample 1 (the default sample size) reconciliation of the tree distribution in testdata/toy_chain.treelist10 with the species tree testdata/sptree.txt. The results are in testdata/DLoutput


Alternatively, you can use:
./DLCoevBoltzmann.py -s testdata/sptree.txt -g testdata/toy_chain.treelist10 -o testdata/DLoutput -n 1000 --temperature 1

To sample 1000 reconciliation, with a temperature of 1 (the default is 0.1). A higher temperature means than reconciliations deviating from parsimony will be less penalized and thus will have a higher chance to be sampled.


You can than analyse the results of the sample. For instance, if you want to know how many sampled reconciliation contains a duplication event in species 2 (which is the ancestor of leaves a and b).

The command:
grep "|2.D" testdata/DLoutput
testdata/DLoutput:12

So on 1000 sampled reconciliation, onlty 12 have a duplication in this species (NB: results may differ a bit as sampling is relying on radnomness).


## using coevents ##

Aside from classical reconciliation, D(T)LRecCoev can also use guide events : specific events with whom the tree to reconcile can form co-events (with the idea that forming a co-event is less costly than forming the event independently). Guide events can be specified in a file that is given to the software using the option '-G', as shown here:

./DLCoevBoltzmann.py -s testdata/sptree.txt -g testdata/toy_chain.treelist10 -o testdata/DLoutputCoev -n 1000 --temperature 1 -G testdata/coev


The file testdata/coev contains the line:
D a b

Which means that there is a guide Duplication (D) in the common ancestor the species a, b.
Looking at reconciled trees, the common ancestor the species a, b and c is named 2.

using the same command as before:

grep "|2.D" testdata/DLoutputCoev
testdata/DLoutputCoev:14

We can see that the number of reconciliations containing a duplication event in species 2 did not change much.
However, this is only the "independent" dupications: duplications that do not form co-events. To count these, use:

grep "|2.coD" testdata/DLoutputCoev
testdata/DLoutputCoev:32

One can see that there are 32 additionnal sampled duplication in this species, thanks to the formation of co-duplication with the guide event (NB : here this effectively triples the number of time we see a duplication here).


Tip:
    using the command:

    grep -c -P "\|2\.(co)*D" testdata/DLoutput*

    you can see the sum of both duplication and co-duplication for both sample (without and with coevents)

    testdata/DLoutput:12
    testdata/DLoutputCoev:46



## DTL reconciliation ##

The introduction of lateral gene transfer compexifies the reconciliation algorithm and, apart from the new costs,
an important parameter is whether or not the species tree is ultrametric or not.

An ultrametric species tree can be harder to obtains, but will always yield "time consistent" scenarios.
A "time consistent" scenario is a reconciliation where a lineage is never (directly or indirectly) transfered to one of its ancestor (which is not feasible without time travelling genes).
By default, DTLCoevBoltzmann considers the species tree to be ultrametric:

./DTLCoevBoltzmann.py -s testdata/sptree.txt -g testdata/toy_chain.treelist10 -o testdata/DTLoutput -n 1000 --temperature 1


will consider testdata/sptree.txt to e ultrametric and reconcile the distribution testdata/toy_chain.treelist10 with it to produce a sample of size 1000.


If, however, the species tree is not ultrametric, you can use the --undated option and DTLCoevBoltzmann will accomodate for this.
By default, the scenarios that are time inconsistent are discarded and DTLCoevBoltzmann will continue sampling until the wanted number of time consistent scenarios is reached.


./DTLCoevBoltzmann.py -s testdata/sptree.txt -g testdata/toy_chain.treelist10 -o testdata/DTLoutputUNDATED -n 1000 --temperature 1 --undated

An additionnal line of output should appear:

backtracked 3556 times to get 1000 trees

Meaning that 2556 scenarios where not time consistent and were discarded.
NB:
    - here the high number of time inconsistent scenarios were mainly due to a high temperature. Try with default temperature to compare

    - you can keep time inconsistent scenarios by using option --ignore.time.consistency

    - By default, DTLCoevBoltzmann will not sample more than 10 times the desired sample size in order to have the correct number of time consistent scenarios. This factor corresponds to the --max.inconsistency.sample.overhead option.



Coevents works much in the same way as in the DL version, with the exception that a few cases where where a guide event is linked to a node AND one of its ancestor, which is not a desired result. 
This only happen when descendant and ancestor are in the same species AND separated by a transfer (meaning that a gene in species a was transfered to another species b, then transfered again to species a ). 
By default, such scenarios are detected when sampling, and eliminated (there is, however, an option to keep them). 
