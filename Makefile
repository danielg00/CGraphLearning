main: main.c
	gcc -g -Wno-unused-result -O0 -Wall main.c linalg.c io.c graph.c -o main
	./main
