HOME := $(shell pwd)

.PHONY: all compile clean install uninstall

all: compile

install:
	cd lcptools && make STATS=1 PREFIX=. install

uninstall:
	cd lcptools && make uninstall PREFIX=.

compile: vg_lcp

vg_lcp: vg_lcp.o vcf_parser.o fasta_reader.o
	$(CXX) -o vg_lcp vg_lcp.o vcf_parser.o fasta_reader.o -L$(HOME)/lcptools/lib -llcptoolsS -Wl,-rpath,$(HOME)/lcptools/lib

vg_lcp.o: vg_lcp.cpp
	$(CXX) -std=c++11 -Wall -Wextra -g -I$(HOME)/lcptools/include -c vg_lcp.cpp -o vg_lcp.o

vcf_parser.o: vcf_parser.cpp vcf_parser.h
	$(CXX) -std=c++11 -Wall -Wextra -g -I$(HOME)/lcptools/include -c vcf_parser.cpp -o vcf_parser.o

fasta_reader.o: fasta_reader.cpp fasta_reader.h
	$(CXX) -std=c++11 -Wall -Wextra -g -I$(HOME)/lcptools/include -c fasta_reader.cpp -o fasta_reader.o

clean:
	rm -f vg_lcp vg_lcp.o vcf_parser.o fasta_reader.o
