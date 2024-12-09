HOME := $(shell pwd)

.PHONY: all compile clean install uninstall

all: compile

install:
	cd lcptools && make STATS=1 PREFIX=. install

uninstall:
	cd lcptools && make uninstall PREFIX=.

compile:
	$(CXX) -std=c++11 -Wall -Wextra -g -I$(HOME)/lcptools/include -c vg_lcp.cpp -o vg_lcp.o
	$(CXX) -o vg_lcp vg_lcp.o -L$(HOME)/lcptools/lib -llcptools -Wl,-rpath,$(HOME)/lcptools/lib

clean:
	rm -f vg_lcp vg_lcp.o