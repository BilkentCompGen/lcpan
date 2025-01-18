# programs
TARGET := lcp_vg
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)

GXX := gcc
CXXFLAGS = -O3 -g -Wall -Wextra -Wpedantic

# object files that need lcptools
LCPTOOLS_CXXFLAGS := -I$(CURRENT_DIR)/lcptools/include
LCPTOOLS_LDFLAGS := -L$(CURRENT_DIR)/lcptools/lib -llcptools -Wl,-rpath,$(CURRENT_DIR)/lcptools/lib -lz

$(TARGET): $(OBJS)
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -o $@ $^ $(LCPTOOLS_LDFLAGS) -lm
	rm $(OBJS)

%.o: %.c
	$(GXX) $(CXXFLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

install:
	cd lcptools && make PREFIX=. install

clean:
	rm -f $(TARGET) $(OBJS)
