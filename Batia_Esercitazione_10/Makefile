C = mpic++ -w
CC = g++
CFLAGS = -Wall -O3 --std=c++11
AFLAGS = -larmadillo

main.exe : main.o random.o
	$(C) main.o random.o -o main.exe $(AFLAGS)
main.o : main.cpp Genetic_Algorithm.h
	$(C) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(C) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
