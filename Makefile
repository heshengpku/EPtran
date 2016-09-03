#-----------------------------------------------------
# Makefile to compile the EPtran_driver system.
#-----------------------------------------------------

# Define compilers and flags.
# Compilers and flags

include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=EPtran_driver

LLIB=EPtran_lib

OBJECTS = EPtran_use_parameter.o\
          EPtran_use_transport.o\
          EPtran_to_tglf.o\
          EPtran_read_input.o \
          EPtran_comp_eq_plasma.o \
          EPtran_comp_alpha_slowing.o\
          EPtran_write_output.o\
          EPtran_transport.o\
          EPtran_mainsub.o 

.SUFFIXES : .o .f90 .f

all: $(LLIB).a $(EXEC)

$(EXEC): $(LLIB).a $(EXEC).o
	$(FC) $(FFLAGS) $(FOMP) -o $(EXEC) $(EXEC).o $(LLIB).a

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) $(FOMP) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXEC) $(LLIB).a
