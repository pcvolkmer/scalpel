BAMTOOLS_DIR := ./bamtools-2.3.0/
ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
SAMTOOLS_DIR := ./samtools-1.1/
ABS_SAMTOOLS_DIR := $(realpath $(SAMTOOLS_DIR))
BCFTOOLS_DIR := ./bcftools-1.1/
ABS_BCFTOOLS_DIR := $(realpath $(BCFTOOLS_DIR))
CXX := g++
CXXFLAGS := -Wno-deprecated -Wall -O3 -fexceptions -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -L$(ABS_BAMTOOLS_DIR)/lib/

all: bamtools samtools bcftools Microassembler

.PHONY : Microassembler
Microassembler:
	cd Microassembler; make; cd ../

.PHONY : samtools
samtools:
	cd $(ABS_SAMTOOLS_DIR); make; cd ../

.PHONY : bcftools
bcftools:
	cd $(ABS_BCFTOOLS_DIR); make; cd ../

.PHONY : bamtools
bamtools:
	mkdir $(ABS_BAMTOOLS_DIR)/build; cd $(ABS_BAMTOOLS_DIR)/build; cmake ..; make; cd ../../

.PHONY : cleanbamtools
cleanbamtools:
	cd $(ABS_BAMTOOLS_DIR)/build; make clean; cd ../../

#.PHONY : clean
clean:
	rm -rf $(ABS_BAMTOOLS_DIR)/build; rm -f Microassembler/Microassembler; cd $(ABS_BCFTOOLS_DIR); make clean-all; cd ..; cd $(ABS_SAMTOOLS_DIR); make clean-all; cd ..;
