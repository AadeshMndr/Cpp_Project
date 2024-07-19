# Compiler
CC := g++

# Compiler flags
CFLAGS := -Wall -Wextra -std=c99

# Libraries
LDFLAGS := -lraylib -lGL -lm -lpthread -ldl -lrt -lX11

# Source files
SRCS := main.cpp

# Object files
OBJS := $(SRCS:.cpp=.o)

# Executable name
TARGET := output

# Default target
all: $(TARGET)

# Compile source files into object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link object files into executable
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# Clean up object files and executable
clean:
	rm -f $(TARGET) $(OBJS)
	rm -rf *.dSYM
	rm -rf *.gch
	rm -rf *.txt
	rm -rf *.o
	rm -rf *.out
