# common definition for Makefile
# for GNU c++
CCACHE=$(shell eval ls /usr/local/bin/ccache 2>/dev/null)

CXX = $(CCACHE) g++ -x c++
CC = $(CCACHE) gcc
LD = g++
CP = cp -f
AR = ar r
MKDIR=mkdir -p
RM=rm -f
CFLAGS = -O3 -fomit-frame-pointer -DNDEBUG -msse2 -mfpmath=sse -march=core2
CFLAGS_WARN=-Wall -Wextra -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wfloat-equal -Wpointer-arith #-Wswitch-enum -Wstrict-aliasing=2
CFLAGS_ALWAYS = -D_FILE_OFFSET_BITS=64 -fno-operator-names
LDFLAGS = -s -lm -lzm $(LIB_DIR)
AS = nasm
AFLAGS = -f elf -D__unix__

BIT=-m32
ifeq ($(shell uname -m),x86_64)
BIT=-m64
endif
ifeq ($(shell uname -s),Darwin)
BIT=-m64
endif

ifeq ($(DBG),on)
CFLAGS += -O0 -g3 -UNDEBUG
LDFLAGS += -g3
endif

.SUFFIXES: .cpp
.cpp.o:
	$(CXX) -c $< -o $@ $(CFLAGS) $(CFLAGS_WARN) $(CFLAGS_ALWAYS) $(INC_DIR) $(BIT)

.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) $(CFLAGS_WARN) $(CFLAGS_ALWAYS) $(INC_DIR) $(BIT)

INC_DIR= -I../src -I../../xbyak -I../include
LIB_DIR= -L../lib
