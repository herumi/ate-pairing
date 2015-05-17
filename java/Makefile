include ../common.mk
LDFLAGS=-lgmp -lgmpxx
ifeq ($(UNAME_S),Darwin)
  JAVA_INC=-I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers/
  LIB_SUF=dylib
else
  JAVA_INC=-I/usr/lib/jvm/default-java/include
#JAVA_INC=-I/usr/lib/jvm/java-7-openjdk-amd64/include
  LIB_SUF=so
  CFLAGS+=-z noexecstack
  LIB+=-lrt
endif
CFLAGS+= -shared -fPIC $(JAVA_INC)

PACKAGE_NAME=mcl.bn254
PACKAGE_DIR=$(subst .,/,$(PACKAGE_NAME))

TARGET=../bin/libbn254_if_wrap.$(LIB_SUF)
JAVA_EXE=cd ../bin && LD_LIBRARY_PATH=./:$(LD_LIBRARY_PATH) java -classpath ../java
all: $(TARGET)

bn254_if_wrap.cxx: bn254_if.i bn254_if.hpp
	$(MKDIR) $(PACKAGE_DIR)
	swig -java -package $(PACKAGE_NAME) -outdir $(PACKAGE_DIR) -c++ -Wall bn254_if.i

$(TARGET): bn254_if_wrap.cxx
	$(MKDIR) ../bin
	$(PRE)$(CXX) $< -o $@ $(CFLAGS) $(LDFLAGS) -I../include ../src/zm2.cpp ../src/zm.cpp -I../../xbyak

%.class: %.java
	javac $<

BN254Test.class: BN254Test.java $(TARGET)

jar:
	jar cvf bn254.jar mcl

test: BN254Test.class $(TARGET)
	$(JAVA_EXE) BN254Test

clean:
	rm -rf *.class $(TARGET) $(PACKAGE_DIR)/*.class bn254_if_wrap.cxx

