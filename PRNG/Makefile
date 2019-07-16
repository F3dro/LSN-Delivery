CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o BASES/random.o
	$(CC) BASES/random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : BASES/random.cpp BASES/random.h
	$(CC) -c BASES/random.cpp -o BASES/random.o $(CFLAGS)
clean :
	rm *.o main.exe BASES/seed.out
