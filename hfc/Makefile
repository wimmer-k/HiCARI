CC = g++

BIN_DIR = $(HOME)/bin
LIB_DIR = $(HOME)/lib

FLAG =  -O3 -g
MFILE_LIB =
RUNFILE = GEB_HFC

INCLUDE =
LIBS =

OBJFILES = GEB_HFC.o HFC.o

$(BIN_DIR)/$(RUNFILE): GEB_HFC.cpp $(OBJFILES)
	$(CC) $(FLAG) $(OBJFILES) -o $(BIN_DIR)/$(RUNFILE) $(LIBS)

GEB_HFC.o: GEB_HFC.cpp HFC.h
	$(CC) $(FLAG) -c $<

HFC.o: HFC.cpp HFC.h
	$(CC) $(FLAG) -c $<


clean:
	rm -f $(BIN_DIR)/$(RUNFILE) $(OBJFILES)
