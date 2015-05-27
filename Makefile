CC = g++
INCLUDE = -I../ -I./ -I./HepMC/include -I./fastjet/include
CFLAGS = -Wall -g $(shell root-config --cflags) $(INCLUDE)
LINKER = g++

LINKERFLAGS = -undefined dynamic_lookup $(shell root-config --libs) -lEG -lGenVector -lMinuit -headerpad_max_install_names

SRCDIR = src
OBJDIR = obj
BINDIR = bin
LNKDIR = links

SOURCES = $(wildcard $(SRCDIR)/*.cxx) $(patsubst $(LNKDIR)/%_linkdef.h,$(SRCDIR)/%_dict.cxx,$(wildcard $(LNKDIR)/*.h))
BINARIES = $(wildcard $(BINDIR)/*.cc)
OBJECTS = $(patsubst %,$(OBJDIR)/%.o,$(basename $(notdir $(SOURCES)))) 
GENHEADERS = generator/GenParticle_p5.h generator/GenEvent_p5.h generator/GenVertex_p5.h generator/McEventCollection_p5.h
INIHHEADERS = inih/INIReader.h
LIBS = -L./HepMC/lib/ -lHepMC -lHepMCfio -L./fastjet/lib/ -lfastjet -lfastjettools -lfastjetplugins -lsiscone -lsiscone_spherical
EXE = analysis.exe

.PRECIOUS: $(OBJDIR)/%.o $(SRCDIR)/%_dict.cxx

.PHONY: all

all: $(EXE)

# General rule for making object files
$(SRCDIR)/%_dict.cxx: $(LNKDIR)/%_linkdef.h
	@rootcling -f $@ -c $(GENHEADERS) $^
	@cp $(SRCDIR)/McEventCollection_p5_dict_rdict.pcm .

$(OBJDIR)/%.o: $(BINDIR)/%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	$(CC) $(CFLAGS) $< -c -o $@

%.exe: $(OBJECTS) $(OBJDIR)/%.o
	$(CC) -Wall -o $@ $(LINKERFLAGS) $(LIBS) $^
	@install_name_tool -change libHepMC.4.dylib @executable_path/HepMC/lib/libHepMC.4.dylib $@
	@install_name_tool -change libHepMCfio.4.dylib @executable_path/HepMC/lib/libHepMCfio.4.dylib $@

clean:  
	rm -v -f \
        $(OBJDIR)/*.o *.exe \
	$(SRCDIR)/*_dict.cxx \
	*.pcm
	@echo "Done"
