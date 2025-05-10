main: main.c
	clang -g -Wno-unused-result -O0 -Wall main.c linalg.c io.c graph.c score_functions.c -lm -o main
	./main
