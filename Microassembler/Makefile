BAMTOOLS_DIR := ../bamtools/
ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
CXX := g++
CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
#CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -pg -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/src/ -I${ABS_BAMTOOLS_DIR}/build/src/ -L$(ABS_BAMTOOLS_DIR)/build/src/

all: Microassembler

Microassembler: Microassembler.cc Microassembler.hh align.cc util.cc Mer.hh Ref.hh ReadInfo.hh ReadStart.hh Transcript.hh Edge.cc Edge.hh ContigLink.hh Node.cc Node.hh Path.cc Path.hh ContigLink.cc Graph.cc Graph.hh
	$(CXX) $(CXXFLAGS) $(INCLUDES) Microassembler.cc Edge.cc Node.cc Graph.cc Path.cc ContigLink.cc align.cc util.cc -o Microassembler -lbamtools -lz
