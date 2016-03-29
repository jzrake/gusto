# if there is no Makefile.in then use the template
# --------------------------------------------------
ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = $(PWD)/Makefile.in
else
MAKEFILE_IN = $(PWD)/Makefile.in.template
endif
include $(MAKEFILE_IN)



ifeq ($(HAVE_HDF5), 1)
HDF5_I = -I$(HDF5_HOME)/include
HDF5_L = -L$(HDF5_HOME)/lib -lhdf5
endif


SRC = gusto.c ser.c quartic.c srmhd_c2p.c
OBJ = $(SRC:.c=.o)
APP = gusto
HEADERS = gusto.h


default : $(APP)


ser.c : ser.h

%.c : %.jin.c build.py
	$(JINJA2) > $@

%.h : %.jin.h build.py
	$(JINJA2) > $@

%.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) $< $(HDF5_I) -c

$(APP) : $(OBJ) $(COW_LIB)
	$(CC) $(CFLAGS) $^ $(HDF5_L) $(FFTW_L) $(RNPL_L) $(CLIBS) -o $@

clean :
	$(RM) $(APP) $(OBJ)

.FORCE :
