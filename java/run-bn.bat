@echo off
echo [[compile BnTest.java]]
%JAVA_DIR%\bin\javac BnTest.java

echo [[run BnTest]]
pushd ..\bin
%JAVA_DIR%\bin\java -classpath ..\java BnTest %1 %2 %3 %4 %5 %6
popd
