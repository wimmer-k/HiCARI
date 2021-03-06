.EXPORT_ALL_VARIABLES:

.PHONY: clean all docs

BIN_DIR = $(HOME)/bin
LIB_DIR = $(HOME)/lib

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)
ROOTGLIBS   := $(shell root-config --glibs)
ROOTINC     := -I$(shell root-config --incdir)

CPP         = g++
CFLAGS	    = -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC -D_FILE_OFFSET_BITS=64 -MMD -std=c++11

INCLUDES    = -I./inc 
BASELIBS    = -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR) 
LIBS  	    =  $(BASELIBS) -lCommandLineInterface -lHRArray

LFLAGS	    = -g -fPIC -shared
CFLAGS 	    += -Wno-overloaded-virtual -Wno-unused-variable -Wno-write-strings
LFLAGS      += -dynamiclib -single_module -undefined dynamic_lookup

CLICFLAGS   = -g2 -O2 -fPIC
CLILFLAGS   = -g -fPIC -dynamiclib -single_module -undefined dynamic_lookup

LIB_O_FILES = build/Gretina.o build/GretinaDictionary.o build/GammaSim.o build/GammaSimDictionary.o build/Miniball.o build/MiniballDictionary.o build/ZeroDeg.o build/ZeroDegDictionary.o build/MINOS.o build/MINOSDictionary.o build/Settings.o build/SettingsDictionary.o build/TrackSettings.o build/TrackSettingsDictionary.o

USING_ROOT_6 = $(shell expr $(shell root-config --version | cut -f1 -d.) \>= 6)
ifeq ($(USING_ROOT_6),1)
        EXTRAS=GretinaDictionary_rdict.pcm GammaSimDictionary_rdict.pcm MiniballDictionary_rdict.pcm ZeroDegDictionary_rdict.pcm MINOSDictionary_rdict.pcm SettingsDictionary_rdict.pcm TrackSettingsDictionary_rdict.pcm 
endif


O_FILES = build/SimHistograms.o build/RawHistograms.o build/CalHistograms.o build/Calibration.o build/Tracking.o build/UnpackedEvent.o 
HO_FILES = build/SimHistograms.o build/RawHistograms.o build/CalHistograms.o 

all: $(LIB_DIR)/libCommandLineInterface.so $(LIB_DIR)/libHRArray.so $(EXTRAS) SimCalculate Sim_histos

SimCalculate: SimCalculate.cc $(LIB_DIR)/libHRArray.so $(O_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@ 

Sim_histos: Sim_histos.cc $(LIB_DIR)/libHRArray.so $(HO_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(HO_FILES) -o $(BIN_DIR)/$@ 


$(LIB_DIR)/libHRArray.so: $(LIB_O_FILES) 
	@echo "Making $@"
	@$(CPP) $(LFLAGS) -o $@ $^ -lc

$(LIB_DIR)/libCommandLineInterface.so: build/CommandLineInterface.o  
	@echo "Making $@"
	@$(CPP) $(CLILFLAGS) -o $@ $^ -lc

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
	@rootcint -f $@ -c $(INCLUDES) $(ROOTCFLAGS) $(notdir $^)

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


