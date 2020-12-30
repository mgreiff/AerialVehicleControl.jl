IDIR =../include
CC=gcc
CC_quiet=@gcc

LINK_LAPACK = true

override CFLAGS += -I$(IDIR) -ansi -pedantic -Wall -Werror -O3 -g

ODIR = obj
LDIR = ../lib
DDIR = ../data

LIBS = -lm
ifeq ($(LINK_LAPACK), true)
$(info Compiling: Using LAPACK...)
LIBS            += -llapack -lblas
override CFLAGS += -DLINK_LAPACK
else
$(info Compiling: Using no external  libreries...)
endif

_DEPS  = cont_types.h cont_math.h cont_main.h cont_matrix_math.h
_DEPS += cont_attitude_utils.h
_DEPS += cont_attitude_FSF_SO3_continuous.h
_DEPS += cont_attitude_FSF_SO3_robust.h
_DEPS += cont_attitude_FSF_SU2_continuous.h
_DEPS += cont_attitude_FSF_SU2_discontinuous.h
_DEPS += cont_attitude_FSF_SU2_robust.h
_DEPS += cont_attitude_FOF_SO3_continuous.h
_DEPS += cont_power_distribution.h
_DEPS += cont_attitude_reference_generator.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ  = cont_types.o cont_math.o cont_main.o cont_matrix_math.o
_OBJ += cont_attitude_utils.o
_OBJ += cont_attitude_FSF_SO3_continuous.o
_OBJ += cont_attitude_FSF_SO3_robust.o
_OBJ += cont_attitude_FSF_SU2_continuous.o
_OBJ += cont_attitude_FSF_SU2_discontinuous.o
_OBJ += cont_attitude_FSF_SU2_robust.o
_OBJ += cont_attitude_FOF_SO3_continuous.o
_OBJ += cont_power_distribution.o
_OBJ += cont_attitude_reference_generator.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

run: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

tests: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

#.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ tests run $(DDIR)/*.csv trace.txt

debug:
	valgrind --tool=memcheck --leak-check=yes ./run
