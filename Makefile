main: main.c
	gcc -g -Wno-unused-result -O0 -Wall main.c linalg.c io.c graph.c score_functions.c utils.c -lm -o main
	./main

.PRECIOUS: main
