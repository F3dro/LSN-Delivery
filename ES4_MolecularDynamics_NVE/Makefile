CC = g++
CFLAGS = -Wall -O03 --std=c++11

EXECUTABLES := MolDyn_NVE.exe
all : $(EXECUTABLES)

%.exe : %.o /home/marco/LSN/PRNG/BASES/random.o
	$(CC) /home/marco/LSN/PRNG/BASES/random.o $< -o $@
%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)
random.o : /home/marco/LSN/PRNG/BASES/random.cpp /home/marco/LSN/PRNG/BASES/random.h
	$(CC) -c /home/marco/LSN/PRNG/BASES/random.cpp -o /home/marco/LSN/PRNG/BASES/random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
