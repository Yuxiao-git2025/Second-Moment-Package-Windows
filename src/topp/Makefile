CMD	= topp
EXE	= $(CMD).exe
CC	= gcc
#FC	= g77
FC	= gfortran
#FC	= f77
SRCS	= $(CMD).f \
	  direct1.f refract.f \
	  tiddid.f  \
	  ttime.f vmodel.f 
CSRCS	= hypot_.c 
OBJS	= $(SRCS:%.f=%.o) $(CSRCS:%.c=%.o)
INCLDIR	= ./include
CFLAGS	= -O -I$(INCLDIR)

# Flags for GNU g77 compiler
#FFLAGS	= -O -I$(INCLDIR) -fno-silent -Wall -implicit

# Flags for SUN f77 compiler
FFLAGS	= -I$(INCLDIR)

# Embed the GNU Fortran/C runtimes so MATLAB can launch topp.exe even when
# it was started before MinGW was added to the Windows PATH.
LDFLAGS	= -static-libgfortran -static-libgcc -static-libquadmath

# Following line needed on HP-UX (hasn't been tested, though).
#LDFLAGS	= +U77

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@

# Extensive lint-like diagnostic listing (SUN f77 only)
hypoDD.lst: $(SRCS)
	f77 -e -Xlist -c $(SRCS)

clean:
	-del /Q $(EXE) $(CMD) *.o core a.out *.fln junk 2>NUL || exit 0

# Include-file dependencies

dist.o		: $(INCLDIR)/geocoord.inc
redist.o	: $(INCLDIR)/geocoord.inc
sdc2.o		: $(INCLDIR)/geocoord.inc
setorg.o	: $(INCLDIR)/geocoord.inc

cluster1.o	: $(INCLDIR)/hypoDD.inc
dtres.o		: $(INCLDIR)/hypoDD.inc
getdata.o	: $(INCLDIR)/hypoDD.inc
hypoDD.o	: $(INCLDIR)/hypoDD.inc
lsfit_lsqr.o	: $(INCLDIR)/hypoDD.inc
lsfit_svd.o	: $(INCLDIR)/hypoDD.inc
partials.o	: $(INCLDIR)/hypoDD.inc
refract.o	: $(INCLDIR)/hypoDD.inc
resstat.o	: $(INCLDIR)/hypoDD.inc
skip.o		: $(INCLDIR)/hypoDD.inc
tiddid.o	: $(INCLDIR)/hypoDD.inc
trialsrc.o	: $(INCLDIR)/hypoDD.inc
ttime.o		: $(INCLDIR)/hypoDD.inc
vmodel.o	: $(INCLDIR)/hypoDD.inc
weighting.o	: $(INCLDIR)/hypoDD.inc

atoangle_.o	: $(INCLDIR)/compat.h
atoangle_.o	: $(INCLDIR)/f77types.h
datetime_.o	: $(INCLDIR)/f77types.h
rpad_.o		: $(INCLDIR)/f77types.h
sscanf3_.o	: $(INCLDIR)/compat.h
sscanf3_.o	: $(INCLDIR)/f77types.h

