CC = g++
CFLAGS = -Wall -O03 --std=c++11

EXECUTABLES := SA.exe
PAREXEC := ParSA.pexe

all : $(EXECUTABLES)
par : $(PAREXEC)

%.exe : %.o /home/marco/LSN/PRNG/BASES/random.o
	$(CC) /home/marco/LSN/PRNG/BASES/random.o $< -o $@

%.pexe : %.o /home/marco/LSN/PRNG/BASES/random.o
	mpicxx /home/marco/LSN/PRNG/BASES/random.o $< -o $@

%.o : %.cpp
	mpicxx -c $< -o $@ $(CFLAGS)
random.o : /home/marco/LSN/PRNG/BASES/random.cpp /home/marco/LSN/PRNG/BASES/random.h
	mpicxx -c /home/marco/LSN/PRNG/BASES/random.cpp -o /home/marco/LSN/PRNG/BASES/random.o $(CFLAGS)
clean :
	rm *.o Genetic_Alg.exe seed.out
