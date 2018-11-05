all : run

CC = g++
CFLAGS = -c
# DEBUG_FLAG = -g
# DEBUG = gdb

run: compile
	$(DEBUG) ./run

compile: main
	$(CC) $(DEBUG_FLAG) -o run main.o parabolic_forward_diff.o gauss_seidel.o parabolic_backward_diff.o parabolic_CN_scheme.o

main:	forward backward gauss_seidel CN_scheme
	$(CC) $(CFLAGS) main.cpp

forward :
	$(CC) $(CFLAGS) parabolic_forward_diff.cpp

backward :
	$(CC) $(CFLAGS) parabolic_backward_diff.cpp

gauss_seidel:
	$(CC) $(CFLAGS) gauss_seidel.cpp

CN_scheme:
	$(CC) $(CFLAGS) parabolic_CN_scheme.cpp

clean:
	rm *.o run
