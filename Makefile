CXXFLAGS=-Wall # -pedantic -ansi
ROOT_LIB:=`root-config --libs --glibs`
ROOT_FLAGS:=`root-config --cflags --ldflags`
ROOT_INCLUDE:=`root-config --incdir`

DEPS= interface/setOutputTree.h interface/METzCalculator.h interface/analysisUtils.h interface/setInputTree.h interface/METzCalculator_Run2.h interface/pseudodataNtuple.h interface/setInputPseudodata.h interface/PUWeight.h interface/readJSONFile.h 
DEPS_OBJ= lib/setOutputTree.o lib/METzCalculator.o lib/analysisUtils.o lib/setInputTree.o lib/METzCalculator_Run2.o lib/pseudodataNtuple.o lib/setInputPseudodata.o lib/PUWeight.o lib/readJSONFile.o ${CMSSW_BASE}/lib/${SCRAM_ARCH}/libBaconAnaDataFormats.so

CC = g++
CFLAGS = -Wall

lib/%.o: src/%.cc $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(ROOT_LIB) $(ROOT_FLAGS)

all: produceWWNtuples.exe produceWWpseudodata.exe PUWeightCalculator.exe hadronicWstudies.exe

PUWeightCalculator.exe: bin/PUWeightCalculator.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

produceWWNtuples.exe: bin/produceWWNtuples.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

produceWWpseudodata.exe: bin/produceWWpseudodata.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

hadronicWstudies.exe: bin/hadronicWstudies.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

clean:
	rm -f lib/*.o
	rm -f lib/*.d
	rm -f *.exe
