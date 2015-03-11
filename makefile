# Compiler.
CC = gcc

# Comiler flags.
#  -std=c99 Compiles using the C99 standard.
#  -Wall    Turns on most compiler warnings.
CFLAGS += -std=c99 -Wall

# Libraries to link.
#  -lm Math library.
LDLIBS += -lm

# Source files for the project.
SRC_FILES := $(wildcard *.c)

# The build target executable.
TARGET = mtrx

all: $(TARGET)

$(TARGET): $(SRC_FILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_FILES) $(LDLIBS)

clean:
	@- $(RM) $(TARGET)

distclean: clean
