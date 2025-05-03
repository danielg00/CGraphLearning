main: main.c
	gcc -g -Wno-unused-result -O0 -Wall main.c linalg.c -o main
	./main
