_processor = $(shell uname -p)
_os = $(shell uname)
MOL_VERSION = 0.0.6
MOL_INCLUDE = "\"../mol.$(MOL_VERSION).h\""
CPPFLAGS := -D _MOL_VERSION_="\"$(MOL_VERSION)\"" $(CPPFLAGS)
CPPFLAGS := -D _MOL_INCLUDE_=$(MOL_INCLUDE) $(CPPFLAGS)
CPPFLAGS := -D _GIT_VERSION_="\"$(shell git describe --always)\"" $(CPPFLAGS)

# library archive file
LIB_FILE = libmol.$(MOL_VERSION).a

# library installation dir
LIB_INSTALL_DIR = $(HOME)/lib
# header files installation dir
HEADER_INSTALL_DIR = $(HOME)/include

ifeq ($(JSON), none)
	CPPFLAGS := -D _NO_JANSSON_ $(CPPFLAGS)
endif

ifeq ($(_os), Darwin)

CPPFLAGS := -D _DARWIN_ -arch x86_64 $(CPPFLAGS)
_os_ver = $(shell sw_vers -productVersion | cut -f 1-2 -d .)
ifeq ($(_os_ver), 10.6) #lion
CPPFLAGS := -D _DARWIN_SNOW_LEOPARD_ $(CPPFLAGS)
endif

endif

ifeq ($(MAKECMDGOALS), openmp)
	CFLAGS := -fopenmp $(CFLAGS)
endif

ifeq ($(MAKECMDGOALS), mol.debug)
	CPPFLAGS := -D _DEBUG_ $(CPPFLAGS)
	CFLAGS := -g $(CFLAGS)
endif
ifeq ($(MAKECMDGOALS), mol.profile)
	CFLAGS := -pg $(CFLAGS)
endif

# file utils
RM = /bin/rm -f
CP = /bin/cp

# archiver and indexer
AR = /usr/bin/ar
RANLIB = /usr/bin/ranlib

# library object files
LIB_OBJS = mol.$(MOL_VERSION)/mem.o \
		   mol.$(MOL_VERSION)/myhelpers.o \
		   mol.$(MOL_VERSION)/prms.o \
		   mol.$(MOL_VERSION)/io.o \
		   mol.$(MOL_VERSION)/icharmm.o \
		   mol.$(MOL_VERSION)/bond.o \
		   mol.$(MOL_VERSION)/atom.o \
		   mol.$(MOL_VERSION)/atom_group.o \
		   mol.$(MOL_VERSION)/_atom_group_copy_from_deprecated.o \
		   mol.$(MOL_VERSION)/init.o \
		   mol.$(MOL_VERSION)/protein.o \
		   mol.$(MOL_VERSION)/pdb.o \
		   mol.$(MOL_VERSION)/ms.o \
		   mol.$(MOL_VERSION)/octree.o \
		   mol.$(MOL_VERSION)/matrix.o \
		   mol.$(MOL_VERSION)/sasa.o \
		   mol.$(MOL_VERSION)/potential.o \
		   mol.$(MOL_VERSION)/energy.o \
                   mol.$(MOL_VERSION)/benergy.o \
                   mol.$(MOL_VERSION)/hbond.o \
                   mol.$(MOL_VERSION)/hbond_probev2.o \
                   mol.$(MOL_VERSION)/nbenergy.o \
		   mol.$(MOL_VERSION)/minimize.o   \
		   mol.$(MOL_VERSION)/compare.o \
		   mol.$(MOL_VERSION)/subag.o \
		   mol.$(MOL_VERSION)/lbfgs.o \
		   mol.$(MOL_VERSION)/version.o \
                   mol.$(MOL_VERSION)/mol2.o \
		   mol.$(MOL_VERSION)/rotamer.o \
		   mol.$(MOL_VERSION)/rmsd.o \
		   mol.$(MOL_VERSION)/gbsa.o \
		   mol.$(MOL_VERSION)/sdf.o \
		   mol.$(MOL_VERSION)/json.o \
		   mol.$(MOL_VERSION)/rigid_body.o \

# library header files
LIB_MAIN_HEADER = mol.$(MOL_VERSION).h
LIB_HEADERS = mol.$(MOL_VERSION)/mem.h \
			  mol.$(MOL_VERSION)/myhelpers.h \
			  mol.$(MOL_VERSION)/prms.h \
			  mol.$(MOL_VERSION)/io.h \
			  mol.$(MOL_VERSION)/icharmm.h \
			  mol.$(MOL_VERSION)/bond.h \
			  mol.$(MOL_VERSION)/atom.h \
			  mol.$(MOL_VERSION)/atom_group.h \
			  mol.$(MOL_VERSION)/_atom_group_copy_from_deprecated.h \
			  mol.$(MOL_VERSION)/init.h \
			  mol.$(MOL_VERSION)/protein.h \
			  mol.$(MOL_VERSION)/pdb.h \
			  mol.$(MOL_VERSION)/ms.h \
			  mol.$(MOL_VERSION)/octree.h \
			  mol.$(MOL_VERSION)/matrix.h \
			  mol.$(MOL_VERSION)/sasa.h \
			  mol.$(MOL_VERSION)/potential.h \
			  mol.$(MOL_VERSION)/energy.h \
                          mol.$(MOL_VERSION)/benergy.h \
                          mol.$(MOL_VERSION)/hbond.h \
                          mol.$(MOL_VERSION)/nbenergy.h \
			  mol.$(MOL_VERSION)/minimize.h \
			  mol.$(MOL_VERSION)/compare.h \
			  mol.$(MOL_VERSION)/subag.h \
		    mol.$(MOL_VERSION)/lbfgs.h \
			  mol.$(MOL_VERSION)/version.h \
			  mol.$(MOL_VERSION)/mol2.h \
			  mol.$(MOL_VERSION)/vector.h \
			  mol.$(MOL_VERSION)/rotamer.h \
			  mol.$(MOL_VERSION)/rmsd.h \
			  mol.$(MOL_VERSION)/gbsa.h \
			  mol.$(MOL_VERSION)/sdf.h \
			  mol.$(MOL_VERSION)/json.h \
			  mol.$(MOL_VERSION)/enums.h \
			  mol.$(MOL_VERSION)/yeti.h \
			  mol.$(MOL_VERSION)/rigid_body.h \

# compiler flags
#CFLAGS := $(CFLAGS) -ffast-math -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -std=c11 -Winline -pedantic
CFLAGS := $(CFLAGS) -ffast-math -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -std=c99 -Winline -pedantic
ifneq ($(MAKECMDGOALS), mol.debug)
#	CFLAGS := $(CFLAGS) -O3 -DNDEBUG -flto
	CFLAGS := $(CFLAGS) -O3 -DNDEBUG
#	CFLAGS := $(CFLAGS) -O0 -flto -DNDEBUG
endif

# create library
$(LIB_FILE): $(LIB_OBJS)
	$(AR) rc $(LIB_FILE) $(LIB_OBJS)
	$(RANLIB) $(LIB_FILE)

mol.$(MOL_VERSION)/hbond_probev2.c:
	gperf --output-file=mol.$(MOL_VERSION)/hbond_probev2.c -m 20 mol.$(MOL_VERSION)/hbond_probev2.gperf

# set version.c to be a phony target
.PHONY: mol.$(MOL_VERSION)/version.c
# compile C files into object files.
%.o: mol.$(MOL_VERSION)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<
#$(CC) $(CPPFLAGS) $(CFLAGS) -c -dD -E $<

$(LIB_OBJS): $(LIB_HEADERS) $(LIB_MAIN_HEADER)

# additional rules
all: $(LIB_FILE)
mol.debug: $(LIB_FILE)
mol.profile: $(LIB_FILE)
openmp: all
clean:
	$(RM) $(LIB_OBJS) $(LIB_FILE)
install:
	mkdir -p $(HEADER_INSTALL_DIR)
	mkdir -p $(LIB_INSTALL_DIR)
	-cp $(LIB_FILE) $(LIB_INSTALL_DIR)
	cp $(LIB_MAIN_HEADER) $(HEADER_INSTALL_DIR)
	rm -rf $(HEADER_INSTALL_DIR)/mol.$(MOL_VERSION)/
	mkdir $(HEADER_INSTALL_DIR)/mol.$(MOL_VERSION)/
	$(CP) $(LIB_HEADERS) $(HEADER_INSTALL_DIR)/mol.$(MOL_VERSION)/
