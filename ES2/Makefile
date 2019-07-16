CC = g++
CFLAGS = -Wall -O3 --std=c++11

EXECUTABLES := main.exe main2.exe
all : $(EXECUTABLES)

%.exe : %.o ../PRNG/BASES/random.o
	$(CC) ../PRNG/BASES/random.o $< -o $@
%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)
random.o : ../PRNG/BASES/random.cpp ../PRNG/BASES/random.h
	$(CC) -c ../PRNG/BASES/random.cpp -o ../PRNG/BASES/random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
