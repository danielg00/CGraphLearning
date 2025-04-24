main: main.c
	gcc  -Wno-unused-result -O3 -Wall main.c linalg.c -o main
	./main
