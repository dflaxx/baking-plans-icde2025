CC        = gcc
#CFLAGS    = -Wall -Wextra -Winline -O3
CFLAGS    = -Wall -Wextra -Winline -O0 -g3

CXX       = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -Weffc++ -Winline -O3
#CXXFLAGS = -std=c++20 -Wall -Wextra -Weffc++ -Winline -O0 -g3 -fsanitize=address
DEPFLAGS = -MMD -MP


CTAGS      = ctags
TAGFLAGS   = --c-kinds=+p --fields=+iaS --extras=+q
TAGFILES   = $(shell find . -type f -regex '.*\(hh\|cc\|h\|c\)' -not -path './.git/*')


TARGETS  = \

MAINOBJS = $(TARGETS:.out=.o)
OBJS     = cbind_to_hw_thread.o \
	   cmeasure.o \
	   CrystalTimer.o \
	   csv_writer.o

.PHONY: all clean printenv

all: $(TARGETS) $(OBJS)

# Linking
$(TARGETS): %.out %.o $(OBJS)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) $^ -o $@ 

# Compilation
%.o: %.cc
	@echo "Compiling $@"
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

%.o: %.c
	@echo "Compiling $@"
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

-include *.d

clean:
	-rm -f $(TARGETS) $(MAINOBJS) $(OBJS) *.d

printenv:
	@echo "CXX      = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "DEPFLAGS = $(DEPFLAGS)"
	@echo "TAGFLAGS = $(TAGFLAGS)"
	@echo "TAGFILES = $(TAGFILES)"
	@echo "TARGETS  = $(TARGETS)"
	@echo "MAINOBJS = $(MAINOBJS)"
	@echo "OBJS     = $(OBJS)"
	

