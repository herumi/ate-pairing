include ../common.mk

TARGET=../lib/libzm.a
all:$(TARGET)

OBJ=zm2.o zm.o
clean:
	$(RM) *.o $(TARGET)

../lib/libzm.a: $(OBJ)
	$(MKDIR) ../lib
	-$(RM) $@
	$(AR) $@ $(OBJ)

zm.o: zm.cpp ../include/zm.h
zm2.o: zm2.cpp ../include/zm2.h ../include/bn.h


