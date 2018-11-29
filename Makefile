CC=gcc
FLAGS=-std=c++11

testVector:
	$(CC) -I . $(FLAGS) test/maths/testVector.cpp -o bin/testVector
