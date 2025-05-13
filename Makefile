CFLAGS := -g -O3 -Wall -Iinclude
CC := gcc
main: src/main.c
	$(CC) $(CFLAGS) src/*.c -lm -o main
