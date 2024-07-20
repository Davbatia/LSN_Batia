CC = g++
CFLAGS = -Wall -O3 --std=c++11
AFLAGS = -larmadillo
MPATH = ../

main.exe : main.o random.o Metropolis.o
	$(CC) main.o random.o Metropolis.o -o main.exe $(AFLAGS)
main.o : main.cpp $(MPATH)metropolis.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : $(MPATH)random.cpp $(MPATH)random.h
	$(CC) -c $(MPATH)random.cpp -o random.o $(CFLAGS)
Metropolis.o : $(MPATH)Metropolis.cpp $(MPATH)Metropolis.h
	$(CC) -c $(MPATH)Metropolis.cpp -o Metropolis.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
