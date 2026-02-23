# programs
TARGET := lcpan
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)

CC := gcc
CFLAGS ?= -O3 -Wall -Wextra -Wpedantic
THREAD_FLAGS ?= -pthread

# object files that need lcptools
LCPTOOLS_CXXFLAGS := -I$(CURRENT_DIR)/lcptools/include
LCPTOOLS_LDFLAGS := -L$(CURRENT_DIR)/lcptools/lib -llcptools -Wl,-rpath,$(CURRENT_DIR)/lcptools/lib -lz

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(LCPTOOLS_CXXFLAGS) -o $@ $^ $(LCPTOOLS_LDFLAGS) -lm $(THREAD_FLAGS)
	@mkdir -p bin
	mv $(TARGET) bin
	rm *.o

%.o: %.c
	$(CC) $(CFLAGS) $(PROF_FLAGS) $(LCPTOOLS_CXXFLAGS) -c $< -o $@ $(THREAD_FLAGS)

install: install-lcptools
	@chmod +x lcpan-merge.sh

install-lcptools:
	@echo "Installing lcptools"
	cd lcptools && \
	make install PREFIX=.

clean:
	rm -f $(TARGET) $(OBJS)

profile: PROF_FLAGS = -g -fno-omit-frame-pointer -fno-optimize-sibling-calls
profile: clean $(TARGET)
	rm *.o

.PHONY: profile install install-lcptools clean $(TARGET)
