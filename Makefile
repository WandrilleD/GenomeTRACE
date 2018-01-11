CXX=g++
ODIR=obj
BDIR=bin
SDIR=source
IDIR=include



CPPFLAGS=-g -O3 -I$(IDIR) -I$(HOME)/local/bpp/dev/include
#-std=c++11
LDFLAGS=-L$(HOME)/local/bpp/dev/lib
LIBS=-lm -lbpp-core -lbpp-seq -lbpp-phyl



W_FILES = GenomeTRACE.cpp WmodifUtils.cpp DTLRecCoevWrapper.cpp DeCoUtils.cpp EquivalenceClassFamily.cpp CoEvent.cpp DeCoOutputManager.cpp MultiRootEquivalenceClass.cpp EquivalenceClass.cpp AdjMatrix.cpp AdjMatrix_CostFunctions.cpp MyMatrixAdjSolution.cpp AdjTree.cpp GeneFamily.cpp ReconciledTree.cpp CladeReconciliation.cpp ReconciliationEvent.cpp DTLMatrix.cpp MyCladesAndTripartitions.cpp MyCladesAndTripartitionsHelper.cpp CladesAndTripartitions.cpp MyMatrix.cpp MyGeneTree.cpp MySpeciesTree.cpp DTLGraph.cpp MyMatrixV.cpp MyMatrixT.cpp XMLUtils.cpp
W_OBJS = $(subst .cpp,.o,$(W_FILES))
W_OBJS := $(addprefix $(ODIR)/, $(W_OBJS))
W_SRCS = $(addprefix $(SDIR)/, $(W_FILES))


BTS_FILES = BiasedTopoSampler.cpp MyCladesAndTripartitions.cpp CladesAndTripartitions.cpp MyGeneTree.cpp
BTS_OBJS = $(subst .cpp,.o,$(BTS_FILES))
BTS_OBJS := $(addprefix $(ODIR)/, $(BTS_OBJS))
BTS_SRCS = $(addprefix $(SDIR)/, $(BTS_FILES))


_EXES= GenomeTRACE

EXES=$(patsubst %,$(BDIR)/%,$(_EXES))


.PHONY: clean all

all: $(EXES)

debug: CPPFLAGS = -g
debug: clean all

## SPECIFIC RULES
$(BDIR)/GenomeTRACE: $(W_OBJS) $(BPP_LIB)

$(BDIR)/BiasedTopoSampler: $(BTS_OBJS) $(BPP_LIB)


#general rules

$(ODIR)/%.o: $(SDIR)/%.cpp
	mkdir -p $(ODIR)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -c -MD -o $@ $<
	@cp $(ODIR)/$*.d $(ODIR)/$*.P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	    -e '/^$$/ d' -e 's/$$/ :/' < $(ODIR)/$*.d >> $(ODIR)/$*.P; \
	rm -f $(ODIR)/$*.d

$(BDIR)/%: $(ODIR)/%.o
	mkdir -p bin 
	$(CXX) -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LIBS)



all: $(EXES)


clean:
		rm -f $(ODIR)/*.o *~ $(SDIR)/*~ core  $(BDIR)/*
		rm -rf $(ODIR) $(BDIR)
