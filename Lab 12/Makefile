all : run clean

CC = g++
CFLAGS = -c
# DEBUG_FLAG = -g
# DEBUG = gdb

run: compile
	$(DEBUG) ./run

compile: main
	$(CC) $(DEBUG_FLAG) -o run main.o gauss_seidel.o display.o

main:	display gauss_seidel
	$(CC) $(CFLAGS) main.cpp

display :
	$(CC) $(CFLAGS) display.cpp

gauss_seidel:
	$(CC) $(CFLAGS) gauss_seidel.cpp

clean:
	rm *.o run
