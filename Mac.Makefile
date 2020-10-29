.EXPORT_ALL_VARIABLES:

.PHONY: clean all docs

BIN_DIR = $(HOME)/bin
LIB_DIR = $(HOME)/lib
TARTSYS=/home/gamma20/packages/anaroot_v4.5.38
#TARTSYS=/home/gamma20/exp/anaroot
#TARTSYS=/home/wimmer/mercurius/anaroot

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)
ROOTGLIBS   := $(shell root-config --glibs)
ROOTINC     := -I$(shell root-config --incdir)

CPP         = g++
CFLAGS	    = -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC -D_FILE_OFFSET_BITS=64 -MMD

INCLUDES    = -I./inc -I$(TARTSYS)/include
BASELIBS    = -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR) -L$(TARTSYS)/lib  -lXMLParser
LIBS  	    =  $(BASELIBS) -lCommandLineInterface -lHiCARI -lanaroot -lananadeko -lanacore -lanabrips -lanaloop -lBigRIPS

LFLAGS        = -g -fPIC -shared
CFLAGS         += -Wno-overloaded-virtual -Wno-unused-variable -Wno-write-strings
LFLAGS      += -dynamiclib -single_module -undefined dynamic_lookup

CLICFLAGS   = -g2 -O2 -fPIC
CLILFLAGS   = -g -fPIC -dynamiclib -single_module -undefined dynamic_lookup


LIB_O_FILES = build/Gretina.o build/GretinaDictionary.o build/Settings.o build/SettingsDictionary.o build/RunInfo.o build/RunInfoDictionary.o build/Trace.o build/TraceDictionary.o build/HiCARI.o build/HiCARIDictionary.o
BRLIB_O_FILES = build/Settings.o build/SettingsDictionary.o build/RunInfo.o build/RunInfoDictionary.o build/PPAC.o build/PPACDictionary.o build/FocalPlane.o build/FocalPlaneDictionary.o build/Beam.o build/BeamDictionary.o 
O_FILES = build/RawHistograms.o build/CalHistograms.o build/Calibration.o build/UnpackedEvent.o
MO_FILES = build/BuildEvents.o build/MergeHistograms.o
HO_FILES = build/RawHistograms.o build/CalHistograms.o build/MergeHistograms.o


USING_ROOT_6 = $(shell expr $(shell root-config --version | cut -f1 -d.) \>= 6)
ifeq ($(USING_ROOT_6),1)
	EXTRAS =  GretinaDictionary_rdict.pcm HiCARIDictionary_rdict.pcm SettingsDictionary_rdict.pcm RunInfoDictionary_rdict.pcm TraceDictionary_rdict.pcm PPACDictionary_rdict.pcm FocalPlaneDictionary_rdict.pcm BeamDictionary_rdict.pcm
endif

all: $(LIB_DIR)/libCommandLineInterface.so $(LIB_DIR)/libHiCARI.so $(LIB_DIR)/libBigRIPS.so  $(EXTRAS) HFC Unpack Calibrate MakeMode2 Raw_histos Cal_histos BigRIPSTree Merge Merge_histos Gated_histos 

SimCalculate: SimCalculate.cc $(LIB_DIR)/libHiCARI.so $(O_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@ 

Sim_histos: Sim_histos.cc $(LIB_DIR)/libHiCARI.so $(HO_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(HO_FILES) -o $(BIN_DIR)/$@ 

Unpack: Unpack.cc $(O_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@ 

Raw_histos: Raw_histos.cc $(LIB_DIR)/libHiCARI.so $(HO_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(HO_FILES) -o $(BIN_DIR)/$@ 

MakeMode2: MakeMode2.cc $(O_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@ 

Calibrate: Calibrate.cc $(O_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@ 

Cal_histos: Cal_histos.cc $(LIB_DIR)/libHiCARI.so $(HO_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(HO_FILES) -o $(BIN_DIR)/$@ 

BigRIPSTree: BigRIPSTree.cc $(LIB_DIR)/libBigRIPS.so
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) -o $(BIN_DIR)/$@ 

Merge: Merge.cc $(LIB_DIR)/libBigRIPS.so $(LIB_DIR)/libHiCARI.so $(MO_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(MO_FILES) -o $(BIN_DIR)/$@ 

Merge_histos: Merge_histos.cc $(LIB_DIR)/libHiCARI.so $(LIB_DIR)/libBigRIPS.so $(HO_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(HO_FILES) -o $(BIN_DIR)/$@ 

Gated_histos: Gated_histos.cc $(LIB_DIR)/libHiCARI.so $(LIB_DIR)/libBigRIPS.so $(HO_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(HO_FILES) -o $(BIN_DIR)/$@ 

HFC:
	@cd hfc; $(MAKE)
	@echo "GEB_HFC compiled"

$(LIB_DIR)/libHiCARI.so: $(LIB_O_FILES)
	@echo "Making $@"
	@$(CPP) $(LFLAGS) -o $@ $^ -lc

$(LIB_DIR)/libBigRIPS.so: $(BRLIB_O_FILES)
	@echo "Making $@"
	@$(CPP) $(LFLAGS) -o $@ $^ -lc

$(LIB_DIR)/libCommandLineInterface.so: build/CommandLineInterface.o  
	@echo "Making $@"
	@$(CPP) $(CLILFLAGS) -o $@ $^ -lc

build/UnpackedEvent.o: src/UnpackedEvent.cc inc/UnpackedEvent.hh $(LIB_O_FILES)
	@echo "Compiling $@"
	@mkdir -p $(dir $@)
	@$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

build/Calibration.o: src/Calibration.cc inc/Calibration.hh $(LIB_O_FILES)
	@echo "Compiling $@"
	@mkdir -p $(dir $@)
	@$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

build/RawHistograms.o: src/RawHistograms.cc inc/RawHistograms.hh $(LIB_O_FILES)
	@echo "Compiling $@"
	@mkdir -p $(dir $@)
	@$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

build/CalHistograms.o: src/CalHistograms.cc inc/CalHistograms.hh $(LIB_O_FILES)
	@echo "Compiling $@"
	@mkdir -p $(dir $@)
	@$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

build/%.o: src/%.cc inc/%.hh
	@echo "Compiling $@"
	@mkdir -p $(dir $@)
	@$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

build/%Dictionary.o: build/%Dictionary.cc
	@echo "Compiling $@"
	@mkdir -p $(dir $@)
	@$(CPP) $(CFLAGS) $(INCLUDES) -fPIC -c $< -o $@

build/%Dictionary.cc: inc/%.hh inc/%LinkDef.h
	@echo "Building $@"
	@mkdir -p $(dir $@)
	@rootcint -f $@ -c $(INCLUDES) $(ROOTCFLAGS) $(SWITCH) $(notdir $^)

build/%Dictionary_rdict.pcm: build/%Dictionary.cc
	@echo "Confirming $@"
	@touch $@

%Dictionary_rdict.pcm: build/%Dictionary_rdict.pcm 
	@echo "Placing $@"
	@cp build/$@ $(LIB_DIR)

build/CommandLineInterface.o: src/CommandLineInterface.cc inc/CommandLineInterface.hh
	@echo "Building $@"
	@mkdir -p $(dir $@)
	@$(CPP) $(CLICFLAGS) $(INCLUDES) -c $< -o $@


doc:	doxy-config
	doxygen doxy-config


clean:
	@echo "Cleaning up"
	@rm -rf build doc
	@rm -f inc/*~ src/*~ *~
	@rm -f scripts/*~  scripts/*_C.*
	@cd hfc; make clean
