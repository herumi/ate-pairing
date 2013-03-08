# common definition for Makefile
# for GNU c++
#CCACHE=$(shell eval ls /usr/local/bin/ccache 2>/dev/null)

#CXX = g++
#CC = gcc
#LD = g++
CP = cp -f
AR = ar r
MKDIR=mkdir -p
RM=rm -f
CFLAGS = -O3 -fomit-frame-pointer -DNDEBUG -msse2 -mfpmath=sse -march=native
CFLAGS_WARN=-Wall -Wextra -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wfloat-equal -Wpointer-arith #-Wswitch-enum -Wstrict-aliasing=2
CFLAGS_ALWAYS = -D_FILE_OFFSET_BITS=64 -fno-operator-names
LDFLAGS = -s -lm -lzm $(LIB_DIR) -lgmp -lgmpxx
AS = nasm
AFLAGS = -f elf -D__unix__

# for only 64-bit
BIT=-m64
#BIT=-m32
#ifeq ($(shell uname -m),x86_64)
#BIT=-m64
#endif
#ifeq ($(shell uname -s),Darwin)
#BIT=-m64
#endif
ifeq ($(shell uname -o),Cygwin)
# install mingw64-x86_64-gcc-g++
CXX=x86_64-w64-mingw32-g++
LD=x86_64-w64-mingw32-g++
AR=x86_64-w64-mingw32-ar r
#LDFLAGS+=-L/usr/x86_64-w64-mingw32/sys-root/mingw/lib
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
