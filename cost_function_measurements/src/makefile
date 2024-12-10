ARCHNAME=$(shell uname -m)
HOSTNAME = $(shell cat /etc/hostname)

# ---------- Variables ----------

DO_OPTIMIZE = 1
DO_DEBUGINFO = 0
DO_ASSERTS = 0
DO_PROFILE_GPROF = 0
DO_SANITIZERS = 0

#CC       := clang
#CXX      := clang++
CC       := gcc
CXX      := g++

CXXSTD = -std=c++20

#CFLAGS := -Wall -Wextra
#CXXFLAGS := $(CXXSTD) -Wall -Wextra
CFLAGS := -Wall #-Wextra #-Winline
CXXFLAGS := $(CXXSTD) #-Wall -Wextra #-Winline

ifeq (1, ${DO_OPTIMIZE})
CFLAGS := ${CFLAGS} -fno-stack-protector -O3
CXXFLAGS := $(CXXFLAGS) -fno-stack-protector -O3
else
$(warning WARNING: Optimization disabled!)
CFLAGS := ${CFLAGS} -O0
CXXFLAGS := $(CXXFLAGS) -O0
endif

ifeq (1, ${DO_DEBUGINFO})
$(warning WARNING: Compiling with debug info: -g3!)
CFLAGS := ${CFLAGS} -ggdb3
CXXFLAGS := ${CXXFLAGS} -ggdb3
else
CFLAGS := ${CFLAGS} -g0
CXXFLAGS := ${CXXFLAGS} -g0
endif

ifeq (1, ${DO_ASSERTS})
CFLAGS := ${CFLAGS} -DDEBUG
CXXFLAGS := ${CXXFLAGS} -DDEBUG
else
$(warning WARNING: Asserts are disabled!)
CFLAGS := ${CFLAGS} -DNDEBUG
CXXFLAGS := ${CXXFLAGS} -DNDEBUG
endif

ifeq (1, ${DO_PROFILE_GPROF})
CFLAGS := ${CFLAGS} -pg
CXXFLAGS := ${CXXFLAGS} -pg
endif

DEPFLAGS = -MMD -MP

ifeq (1, ${DO_SANITIZERS})
SANITIZERS=-fsanitize=address,undefined
endif

# build directory
BUILDDIR = ./build/0
ifdef HDI
	BUILDDIR = ./build/$(HOSTNAME)
else ifdef ODI
	BUILDDIR = ./build/$(ARCHNAME)
endif


TARGETS  = $(BUILDDIR)/main_experiment.out
MAINOBJS = $(TARGETS:.out=.o)
GMOBJS   = $(BUILDDIR)/cb_glob.o \
	   $(BUILDDIR)/GenRandIntVec.o \
	   $(BUILDDIR)/line.o \
	   $(BUILDDIR)/meas_eval.o \
	   $(BUILDDIR)/gminfra/prime_table_1_1.o
OBJS     = $(GMOBJS) \
	   $(BUILDDIR)/hashtable_statistics.o \
	   $(BUILDDIR)/dfinfra/cmeasure.o \
	   $(BUILDDIR)/dfinfra/CrystalTimer.o \
	   $(BUILDDIR)/dfinfra/csv_writer.o \
	   $(BUILDDIR)/dfinfra/cbind_to_hw_thread.o


# ---------- Build Config ----------
$(info ******* BUILD CONFIG ********)
$(info host: ${HOSTNAME})
$(info arch: $(ARCHNAME))
$(info build dir: $(BUILDDIR))
$(info OPT  DEBUG  ASSERT  PROF  SAN)
$(info ${DO_OPTIMIZE}    ${DO_DEBUGINFO}      ${DO_ASSERTS}       ${DO_PROFILE_GPROF}     ${DO_SANITIZERS}) 
$(info *****************************)


# ---------- Targets & Rules ----------
.PHONY: all clean printenv

all: $(TARGETS)

# Linking
$(TARGETS): %.out: %.o $(OBJS)
	@echo "Linking $@"
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(SANITIZERS) $^ -o $@ 

# Compilation
$(BUILDDIR)/%.o: %.cc
	@echo "Compiling $@ (C++ compiler)"
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(SANITIZERS) $(DEPFLAGS) $(IPATHS) -c -o $@ $<

$(BUILDDIR)/%.o: %.c
	@echo "Compiling $@ (C compiler)"
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(DEPFLAGS) $(IPATHS) -c -o $@ $<

# include all dependency files (*.d) in the build directory recursively
DEPFILES = $(shell find $(BUILDDIR) -type f -name '*.d')
-include $(DEPFILES)
# below is not recursive!
#-include $(BUILDDIR)/*.d

clean:
	-rm -rf $(BUILDDIR)

printenv:
	@echo "HOSTNAME         = ${HOSTNAME}"
	@echo "DO_OPTIMIZE      = ${DO_OPTIMIZE}"
	@echo "DO_DEBUGINFO     = ${DO_DEBUGINFO}"
	@echo "DO_ASSERTS       = ${DO_ASSERTS}"
	@echo "DO_PROFILE_GPROF = ${DO_PROFILE_GPROF}"
	@echo "DO_SANITIZERS    = ${DO_SANITIZERS}"
	@echo "CC               = $(CC)"
	@echo "CFLAGS           = $(CFLAGS)"
	@echo "CXX              = $(CXX)"
	@echo "CXXFLAGS         = $(CXXFLAGS)"
	@echo "DEPFLAGS         = $(DEPFLAGS)"
	@echo "TARGETS          = $(TARGETS)"
	@echo "MAINOBJS         = $(MAINOBJS)"
	@echo "OBJS             = $(OBJS)"
	@echo "TAGFLAGS         = $(TAGFLAGS)"
	@echo "TAGFILES         = $(TAGFILES)"
	

# # Some makefile documentation
#
# ## Automatic variables
# - $@
#   The file name of the target of the rule.
#
# - $<
#   The name of the first prerequisite.
#
# - $?
#   The names of all the prerequisites that are newer than the target, with spaces between them. 
#
# - $^
#   The names of all the prerequisites, with spaces between them.
#
# - $(@D)
#   The directory part of the file name of the target, with the trailing slash removed.
#   If the value of ‘$@’ is dir/foo.o then ‘$(@D)’ is dir. This value is . if ‘$@’ does not contain a slash. 
