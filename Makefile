CC = g++
INCLUDE = -I../ -I./ -I./HepMC/include -I./fastjet/include
CFLAGS = -Wall -g $(shell root-config --cflags) $(INCLUDE)
LINKER = g++

LINKERFLAGS = -undefined dynamic_lookup $(shell root-config --libs) -lEG -lGenVector -lMinuit

SRCDIR = src
OBJDIR = obj
BINDIR = bin
LNKDIR = links

SOURCES = $(wildcard $(SRCDIR)/*.cxx) $(patsubst $(LNKDIR)/%_linkdef.h,$(SRCDIR)/%_dict.cxx,$(wildcard $(LNKDIR)/*.h))
BINARIES = $(wildcard $(BINDIR)/*.cc)
OBJECTS = $(patsubst %,$(OBJDIR)/%.o,$(basename $(notdir $(SOURCES)))) $(patsubst %,$(OBJDIR)/%.o,$(basename $(notdir $(BINARIES)))) 
GENHEADERS = generator/GenParticle_p5.h generator/GenEvent_p5.h generator/GenVertex_p5.h generator/McEventCollection_p5.h
LIBS = -L./HepMC/lib/ -lHepMC -lHepMCfio -L./fastjet/lib/ -lfastjet

.PHONY: all

all: $(OBJECTS)
	@echo $(SOURCES)
	@echo $(OBJECTS)
	@echo "Making the program"
	$(CC) -Wall -o reader.exe $(LINKERFLAGS) $(LIBS) $^
	@cp $(SRCDIR)/McEventCollection_p5_dict_rdict.pcm .
	@install_name_tool -change libHepMC.4.dylib @executable_path/HepMC/lib/libHepMC.4.dylib reader.exe
	@install_name_tool -change libHepMCfio.4.dylib @executable_path/HepMC/lib/libHepMCfio.4.dylib reader.exe

# General rule for making object files
$(SRCDIR)/%_dict.cxx: $(LNKDIR)/%_linkdef.h
	@rootcling -f $@ -c $(GENHEADERS) $^

$(OBJDIR)/%.o: $(BINDIR)/%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	$(CC) $(CFLAGS) $< -c -o $@


clean:  
	rm -v -f \
        $(OBJDIR)/*.o *.exe \
	$(SRCDIR)/*_dict.cxx \
	*.pcm
	@echo "Done"
