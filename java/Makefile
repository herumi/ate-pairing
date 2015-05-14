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

PACKAGE_NAME=mcl.bn
PACKAGE_DIR=$(subst .,/,$(PACKAGE_NAME))

TARGET=../bin/libbn_if_wrap.$(LIB_SUF)
JAVA_EXE=cd ../bin && LD_LIBRARY_PATH=./:$(LD_LIBRARY_PATH) java -classpath ../java
all: $(TARGET)

bn_if_wrap.cxx: bn_if.i bn_if.hpp
	$(MKDIR) $(PACKAGE_DIR)
	swig -java -package $(PACKAGE_NAME) -outdir $(PACKAGE_DIR) -c++ -Wall bn_if.i

$(TARGET): bn_if_wrap.cxx
	$(MKDIR) ../bin
	$(PRE)$(CXX) $< -o $@ $(CFLAGS) $(LDFLAGS) -I../include ../src/zm2.cpp ../src/zm.cpp

%.class: %.java
	javac $<

BnTest.class: BnTest.java $(TARGET)

jar:
	jar cvf bn.jar mcl

test: BnTest.class $(TARGET)
	$(JAVA_EXE) BnTest

clean:
	rm -rf *.class $(TARGET) $(PACKAGE_DIR)/*.class bn_if_wrap.cxx

