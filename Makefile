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

SRC = gusto.c mesh.c ser.c quartic.c srmhd_c2p.c chkpt.c utils.c physics.c
OBJ = $(SRC:.c=.o)
APP = gusto
HEADERS = gusto.h

default : $(APP)

gusto.h : ser.h

%.c : %.jin.c build.py
	$(JINJA2) > $@

%.h : %.jin.h build.py
	$(JINJA2) > $@

%.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) $< $(HDF5_I) -c

$(APP) : $(OBJ)
	$(CC) $(CFLAGS) $^ $(HDF5_L) $(CLIBS) -o $@

clean :
	$(RM) $(APP) $(OBJ)
