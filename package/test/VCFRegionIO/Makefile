CXX = g++ 
CC = gcc
UMICHLIB=../../source/umich
CORELIB=../../source
CFLAGS=-O3 -I./ -I$(CORELIB) -I $(UMICHLIB)/general -I$(UMICHLIB)/klib -I$(UMICHLIB)/vcf -Wall
CXXFLAGS= $(CFLAGS) -std=c++11
BINDIR = ./ 
VCF = demo.exe
EXECUTABLES = $(VCF)
# VCF File Set
VCFBASE = $(UMICHLIB)/vcf/VcfFile $(UMICHLIB)/vcf/VcfFileReader \
	$(UMICHLIB)/vcf/VcfGenotypeField $(UMICHLIB)/vcf/VcfGenotypeFormat \
	$(UMICHLIB)/vcf/VcfGenotypeSample $(UMICHLIB)/vcf/VcfRecord \
	$(UMICHLIB)/vcf/VcfHeader $(UMICHLIB)/vcf/VcfRecordField \
	$(UMICHLIB)/vcf/VcfRecordGenotype $(UMICHLIB)/vcf/VcfRecord \
	$(UMICHLIB)/vcf/VcfRecordInfo $(UMICHLIB)/vcf/VcfSubsetSamples
VCFHDR = $(VCFBASE:=.h) $(CORELIB)/VcfRegionIO.hpp
VCFSRC = $(VCFBASE:=.cpp) demo.cpp 
VCFOBJ = $(VCFSRC:.cpp=.o)

# Utility Library File Set
LIBFILE = $(UMICHLIB)/general/lib-vcf.a
LIBMAIN = $(UMICHLIB)/general/InputFile $(UMICHLIB)/general/Tabix $(UMICHLIB)/general/FortranFormat \
 $(UMICHLIB)/general/GenotypeLists $(UMICHLIB)/general/InputFile $(UMICHLIB)/general/IntArray $(UMICHLIB)/general/Hash \
 $(UMICHLIB)/general/MemoryInfo $(UMICHLIB)/general/MiniDeflate \
 $(UMICHLIB)/general/Parameters $(UMICHLIB)/general/QuickIndex $(UMICHLIB)/general/Random $(UMICHLIB)/general/Sort \
 $(UMICHLIB)/general/StringArray $(UMICHLIB)/general/StringBasics $(UMICHLIB)/general/StringMap \
 $(UMICHLIB)/general/StringHash $(UMICHLIB)/libsrc/TraitTransformations \
 $(UMICHLIB)/klib/bgzf $(UMICHLIB)/klib/knetfile $(UMICHLIB)/general/PhoneHome \
 $(UMICHLIB)/general/GzipHeader $(UMICHLIB)/general/BgzfFileType $(UMICHLIB)/general/BgzfFileTypeRecovery \
 $(UMICHLIB)/general/GzipFileType $(UMICHLIB)/general/UncompressedFileType $(UMICHLIB)/general/FileType

LIBPED = $(UMICHLIB)/general/PedigreeLoader $(UMICHLIB)/general/PedigreeTwin $(UMICHLIB)/general/PedigreeTrim
LIBSRC = $(LIBMAIN:=.cpp) $(LIBPED:=.cpp)
LIBHDR = $(LIBMAIN:=.h) $(UMICHLIB)/general/Constant.h \
 $(UMICHLIB)/general/MathConstant.h $(UMICHLIB)/general/PedigreeAlleles.h $(UMICHLIB)/general/LongInt.h
LIBOBJ = $(LIBSRC:.cpp=.o)

# make everything
all : $(EXECUTABLES)

$(EXECUTABLES) : $(BINDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

# dependencies for executables
$(VCF) : $(LIBFILE) $(VCFOBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(VCFOBJ) $(LIBFILE) -lm -lz

$(LIBFILE) : $(LIBOBJ) $(LIBHDR)
	ar -cr $@ $(LIBOBJ)
	ranlib $@

$(VCFOBJ) : $(VCFHDR) $(LIBHDR)

$(LIBOBJ) : $(LIBHDR)

clean :
	rm -f ./*.o $(UMICHLIB)/*/*.a $(UMICHLIB)/*/*.o $(EXECUTABLES) $(CORELIB)/*.o

.c.o :
	$(CC) $(CFLAGS) -o $@ -c $*.c

.cpp.o : 
	$(CXX) $(CXXFLAGS) -o $@ -c $*.cpp 

.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)
