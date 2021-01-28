lsh:main.o Functions.o Classes.o
	g++ main.o Functions.o Classes.o -o lsh

main.o:main.cpp Classes.h Functions.h
	g++ -c main.cpp

Functions.o:Functions.cpp Classes.h Functions.h
	g++ -c Functions.cpp

Classes.o:Classes.cpp Classes.h Functions.h
	g++ -c Classes.cpp

run:
	./lsh


clean:
	rm *.o lsh
