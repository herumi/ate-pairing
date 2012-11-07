include ../common.mk

HEADER= ../include/zm.h ../include/bn.h
TARGET=test_zm bn ../lib/libzm.a

ifeq ($(BIT),-m64)
TARGET += bench_test sample
endif

all:$(TARGET)

.SUFFIXES: .cpp

test_zm: test_zm.o ../lib/libzm.a
	$(CXX) -o $@ $< $(LDFLAGS) $(BIT)

bn: bn.o ../lib/libzm.a
	$(CXX) -o $@ $< $(LDFLAGS) $(BIT)

sample: sample.o ../lib/libzm.a
	$(CXX) -o $@ $< $(LDFLAGS) $(BIT)

bench_test: bench.o ../lib/libzm.a
	$(CXX) -o $@ $< $(LDFLAGS) $(BIT)

../lib/libzm.a:
	$(MAKE) -C ../src ../lib/libzm.a

clean:
	$(RM) *.o $(TARGET)
	$(MAKE) -C ../src clean

test: $(TARGET)
	./test_zm && ./bn

test_zm.o: test_zm.cpp $(HEADER)
bn.o: bn.cpp $(HEADER)
bench.o: bench.cpp $(HEADER)

