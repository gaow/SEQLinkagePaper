CXX = g++ 
CFLAGS=-O3 -I ./ -I./libsrc -I./merlin -I./pdf -I./clusters -D_FILE_OFFSET_BITS=64 -D__ZLIB_AVAILABLE__ -Wall -std=c++11
BINDIR = demo 
MERLIN = demo/demo.exe
EXECUTABLES = $(MERLIN)
# MERLIN File Set
MERLINBASE = merlin/AssociationAnalysis merlin/FastAssociation \
 merlin/AnalysisTask merlin/Conquer \
 merlin/ConquerHaplotyping merlin/DiseaseModel \
 merlin/ParametricLikelihood merlin/GenotypeInference \
 merlin/Houdini \
 merlin/KongAndCox merlin/Manners merlin/MerlinBitSet \
 merlin/MerlinCluster merlin/MerlinCore \
 merlin/MerlinError merlin/MerlinFamily merlin/MerlinIBD \
 merlin/InformationContent merlin/MerlinCache merlin/MerlinModel \
 merlin/MerlinKinship merlin/MerlinKinship15 \
 merlin/MerlinHaplotype merlin/MerlinMatrix merlin/MerlinParameters \
 merlin/MerlinPDF merlin/MerlinSimulator merlin/MerlinSimwalk2 \
 merlin/MerlinSort merlin/NPL-ASP merlin/NPL-QTL \
 merlin/Magic merlin/Mantra merlin/Parametric merlin/QtlModel \
 merlin/Tree \
 merlin/TreeBasics merlin/TreeIndex merlin/TreeManager \
 merlin/TreeInfo merlin/TreeFlips merlin/VarianceComponents
MERLINHDR = $(MERLINBASE:=.h) merlin/TreeNode.h Core.hpp
MERLINSRC = $(MERLINBASE:=.cpp) demo/LibmerlinDemo.cpp Core.cpp
MERLINOBJ = $(MERLINSRC:.cpp=.o)

# Files for dealing with clustered markers
CLUSTERS = clusters/HaploFamily clusters/HaploGraph \
 clusters/HaploSet clusters/HaploTree clusters/Likelihood \
 clusters/SparseLikelihood clusters/Unknown
CLUSTERCPP = $(CLUSTERS:=.cpp)
CLUSTERHDR = $(CLUSTERS:=.h)
CLUSTEROBJ = $(CLUSTERS:=.o)
# Utility Library File Set
LIBFILE = libsrc/lib-goncalo.a
LIBMAIN = libsrc/BasicHash libsrc/Error libsrc/FortranFormat \
 libsrc/GenotypeLists libsrc/InputFile libsrc/IntArray libsrc/Hash \
 libsrc/LongArray libsrc/Kinship libsrc/KinshipX  libsrc/MapFunction \
 libsrc/MathCholesky libsrc/MathDeriv libsrc/MathFloatVector \
 libsrc/MathGenMin libsrc/MathGold libsrc/MathMatrix libsrc/MathStats \
 libsrc/MathNormal libsrc/MathSVD libsrc/MathVector \
 libsrc/MemoryInfo libsrc/MiniDeflate \
 libsrc/Parameters libsrc/Pedigree libsrc/PedigreeAlleleFreq \
 libsrc/PedigreeDescription libsrc/PedigreeFamily libsrc/PedigreeGlobals \
 libsrc/PedigreePerson libsrc/QuickIndex libsrc/Random libsrc/Sort \
 libsrc/StringArray libsrc/StringBasics libsrc/StringMap \
 libsrc/StringHash libsrc/TraitTransformations
LIBPED = libsrc/PedigreeLoader libsrc/PedigreeTwin libsrc/PedigreeTrim
LIBSRC = $(LIBMAIN:=.cpp) $(LIBPED:=.cpp)
LIBHDR = $(LIBMAIN:=.h) libsrc/Constant.h \
 libsrc/MathConstant.h libsrc/PedigreeAlleles.h libsrc/LongInt.h
LIBOBJ = $(LIBSRC:.cpp=.o)
# PDF Library File Sets
PDFLIB = pdf/libpdf.a
PDFFILES = pdf/PDF pdf/PDFfont pdf/PDFinfo pdf/PDFpage \
 pdf/PDFchartbasics pdf/PDFchartbar pdf/PDFlinechart \
 pdf/PDFhistogram \
 pdf/PDFchartaxis pdf/PDFchartlegend pdf/PDFchartmarker \
 pdf/PDFchartline pdf/PDFchartobject
PDFSRC = $(PDFFILES:=.cpp)
PDFHDR = $(PDFFILES:=.h)
PDFOBJ = $(PDFFILES:=.o)

# make everything
all : $(EXECUTABLES)

$(EXECUTABLES) : $(BINDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

# dependencies for executables
$(MERLIN) : $(LIBFILE) $(PDFLIB) $(MERLINOBJ) $(CLUSTEROBJ)
	$(CXX) $(CFLAGS) -o $@ $(MERLINOBJ) $(CLUSTEROBJ) $(PDFLIB) $(LIBFILE) -lm -lz

$(LIBFILE) : $(LIBOBJ) $(LIBHDR)
	ar -cr $@ $(LIBOBJ)
	ranlib $@

$(PDFLIB) : $(PDFOBJ)
	ar -cr $@ $(PDFOBJ)
	ranlib $@

$(MERLINOBJ) : $(MERLINHDR) $(CLUSTERHDR) $(LIBHDR)

$(CLUSTEROBJ) : $(CLUSTERHDR) $(MERLINHDR) $(LIBHDR)

$(LIBOBJ) : $(LIBHDR)

$(PDFOBJ) : $(PDFHDR)

clean :
	rm -f */*.a */*.o $(EXECUTABLES)

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c

#.cpp.X.o : 
#	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -D__CHROMOSOME_X__

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp 

.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)
