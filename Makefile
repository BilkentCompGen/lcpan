# programs
TARGET := lcpan
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)

CC := cc
CFLAGS ?= -O3 -Wall -Wextra -Wpedantic
THREAD_FLAGS ?= -pthread

# object files that need lcptools
LCPTOOLS_CXXFLAGS := -I./lcptools/include
LCPTOOLS_LDFLAGS := -L./lcptools/lib -llcptools -Wl,-rpath,./lcptools/lib -lz

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(THREAD_FLAGS) $(LCPTOOLS_CXXFLAGS) -o $@ $^ $(LCPTOOLS_LDFLAGS) -lm
	@mkdir -p bin
	mv $(TARGET) bin
	rm *.o

%.o: %.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(THREAD_FLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@

install:
	@echo "Installing lcptools"
	cd lcptools && make PREFIX=. install

clean:
	rm -f $(TARGET) $(OBJS)

profile: PROF_FLAGS = -g -fno-omit-frame-pointer -fno-optimize-sibling-calls
profile: clean $(TARGET)
	rm *.o

.PHONY: profile install clean $(TARGET)