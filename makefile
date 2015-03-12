# Compiler.
CC := gcc

# Comiler flags.
#  -std=c99 Compiles using the C99 standard.
#  -I       Indicates locations to look for headers.
CFLAGS += -std=c99 -I src

# Libraries to link.
#  -lm Math library.
LDLIBS += -lm

# Source files for the project.
SRC_FILES := $(wildcard src/*.c)

# Example files.
EXAMPLE_FILES := $(wildcard example/*.c)

# The build target executable.
TARGET := mtrx

# Items for building and running the unit tests.
TEST_TARGET := testit
TEST_FILES := $(wildcard tests/*.c)
TEST_FLAGS := -D_GNU_SOURCE

.PHONY: all clean distclean test

all: clean $(TARGET)

$(TARGET): $(SRC_FILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_FILES) $(EXAMPLE_FILES) $(LDLIBS)

clean:
	@- $(RM) $(TARGET)

distclean: clean

test:
	$(CC) $(CFLAGS) -o $(TEST_TARGET) $(SRC_FILES) $(TEST_FILES) $(LDLIBS) $(TEST_FLAGS)
	./$(TEST_TARGET)
	@- $(RM) $(TEST_TARGET)
