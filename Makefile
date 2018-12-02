CC=g++
FLAGS=-std=c++11 -lblas
testVector:
	$(CC) -I . $(FLAGS) test/maths/testVector.cpp -o bin/test/testVector

testMatrix:
	$(CC) -I . $(FLAGS) test/maths/testMatrix.cpp -o bin/test/testMatrix
