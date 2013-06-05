# Makefile for DEC geo

SHELL = /bin/sh
.IGNORE:

# System-specific settings

F77 =		f77
F77FLAGS =	-O -I/home/u/sjplimp/bin/GEO
CC =		cc
CCFLAGS =	-O -I/home/u/sjplimp/bin/GEO -DFFT_DEC
LINK =		f77
LINKFLAGS =	-O -L/home/u/sjplimp/bin/GEO
USRLIB =	-lmpi -ldxml
SYSLIB =
SIZE =		size

# Link rule

$(EXE):	$(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -o $(EXE)
	$(SIZE) $(EXE)

# Compilation rules

.f.o:
	@cp $< $*.F
	$(F77) $(F77FLAGS) -c $*.F
	@rm $*.F

.c.o:
	$(CC) $(CCFLAGS) -c $<

# Individual dependencies

$(OBJ):	$(INC)
