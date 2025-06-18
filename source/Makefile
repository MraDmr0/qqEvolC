ifeq ($(COMP),)
COMP=gnu
endif

COMMONDIR=./common

ifeq ($(COMP),gnu)
CC=gcc
CXX=g++                  # <--- importante!
CCFLAGS=-g -O2 -march=native -Wall -I$(COMMONDIR) -DNDEBUG
LDFLAGS=
LIBS=
endif

EXE=main
OBJS=main.o algorithms.o envelopes.o potentials.o 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(CCFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LIBS)

%.o: $(COMMONDIR)/%.c
	$(CC) $(CCFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CCFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

.PHONY: clean
clean:
	del /Q $(EXE) $(OBJS) *.o *.png *~ 2>nul || exit 0
