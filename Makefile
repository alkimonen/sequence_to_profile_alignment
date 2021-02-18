all: build
build: seq2profileAlignment.c
	gcc -g -o seq2profileAlignment seq2profileAlignment.c.c
clean:
	rm -fr seq2profileAlignment.c seq2profileAlignment.c.o *~