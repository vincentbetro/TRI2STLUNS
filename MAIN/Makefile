
include ../Makefile.options

LIBRARY_PATH += -L../GEOMETRY -L../IO -L../UTIL

LIBS += -lgeometry -lio -lutility -loctree -lSpaceField -lcgns -lhdf5
#LIBS += -lgeometry -lio -lutility -loctree -lcgns /usr/local/hdf5-1.6.4-icc/lib/libhdf5.a -ltet

LDFLAGS = $(LIBRARY_PATH) $(LIBS)

OPTC_SRCS = TRI2STL_WRAP.cpp \
			TRI2STL.cpp \

OPTSRCS = $(OPTC_SRCS)
OPTOBJECTS = $(OPTC_SRCS:.cpp=.o)

all: TRI2STLUNS

TRI2STLUNS: geometry io utility $(OPTOBJECTS)
	$(LINK) -o $@ $(OPTOBJECTS) $(LDFLAGS)
	cp $@ ../$@.$(EXE_SUFFIX)

geometry:
	cd ../GEOMETRY; $(MAKE)

io:
	cd ../IO; $(MAKE)

utility:
	cd ../UTIL; $(MAKE)

clean:
	/bin/rm -f *.o
	/bin/rm -f TRI2STLUNS

include ../suffix.rules
