
LINKEDLIBS := $(shell $(ROOTSYS)/bin/root-config --noauxlibs --glibs) -lTreePlayer -lMinuit -lRootAuth -lSpectrum -lSpectrumPainter -lfftw3 -lmxml -lm 

all: egg2ROOT 

CXXFLAGS := $(CXXFLAGS) -I/usr/include/ -I$(ROOTSYS)/include/ -I../monarch/ -I../cicada

egg2ROOT: monarch.o cicada.o egg2ROOT.cpp cROOTFile.o
	g++ -o egg2ROOT $(CXXFLAGS) egg2ROOT.cpp monarch.o cicada.o cROOTFile.o $(LINKEDLIBS) 


cROOTFile.o: cROOTFile.cpp cROOTFile.h
	g++ -c cROOTFile.cpp $(CXXFLAGS)

monarch.o: ../monarch/monarch.c ../monarch/monarch.h
	gcc -c ../monarch/monarch.c 

cicada.o: ../cicada/cFFTransform.c ../cicada/cFFTransform.h
	gcc -c -I../monarch/ ../cicada/cFFTransform.c -o cicada.o

clean:
	rm *.o
	rm egg2ROOT
