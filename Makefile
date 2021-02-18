all: build
build: seq2profileAlignment.c
	gcc -g -o alignSeqToProfile seq2profileAlignment.c
clean:
	rm -fr alignSeqToProfile seq2profileAlignment.o *~
