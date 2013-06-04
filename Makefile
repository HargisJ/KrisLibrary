SRCS= $(wildcard *.cpp)
LIBNAME= include
INCDIR= .
DEFINES= 
include Makefile.config
include Makefile.template

MAKE = make -j 2
SUB_LIBS = . planning optimization statistics geometry meshing spline image robotics GLdraw math math3d camera  utils

default: KrisLibrary docs

PQP:
	cd geometry/PQP; $(MAKE)

all: lib PQP
	 cd utils; $(MAKE)
	 cd math; $(MAKE)
	 cd math3d; $(MAKE)
	 cd camera; $(MAKE)
	 cd GLdraw; $(MAKE)
	 cd optimization; $(MAKE)
	 cd geometry; $(MAKE)
	 cd meshing; $(MAKE)
	 cd image; $(MAKE)
	 cd statistics; $(MAKE)
	 cd robotics; $(MAKE)
	 cd planning; $(MAKE)
	 cd spline; $(MAKE)

KrisLibrary: all 
	 $(AR) rcs $(LIBDIROUT)/libKrisLibrary.a $(addsuffix /$(OBJDIR)/*.o,$(SUB_LIBS)) $(GEOMETRY_EXTRAOBJS) $(OPTIMIZATION_EXTRAOBJS)
	 $(RANLIB) $(LIBDIROUT)/libKrisLibrary.a

docs:
	 doxygen doxygen.conf

