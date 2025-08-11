CC = gcc
CFLAGS = -w -O2 -fPIC -I./src
LDFLAGS = -lm

SRC_DIR = src/image_compression/qr_decomposition
TEST_DIR = tests/image_compression/qr_decomposition

MAIN_SRC = src/image_compression/main.c
SVD_SRC = src/image_compression/svd.c

TARGET = prueba
LIB = libimage_compression.a

SRC_ALL = $(wildcard $(SRC_DIR)/*.c) $(SVD_SRC) $(MAIN_SRC)

SRC_LIB = $(filter-out $(MAIN_SRC),$(SRC_ALL))

OBJ_LIB = $(SRC_LIB:.c=.o)
OBJ_MAIN = $(MAIN_SRC:.c=.o)

# Tests
TEST_SRCS = $(wildcard $(TEST_DIR)/test_*.c)
TEST_BINS = $(patsubst $(TEST_DIR)/%.c,$(TEST_DIR)/%,$(TEST_SRCS))

all: $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB): $(OBJ_LIB)
	ar rcs $@ $^

$(TARGET): $(OBJ_MAIN) $(LIB)
	$(CC) $(CFLAGS) $(OBJ_MAIN) -o $@ $(LIB) $(LDFLAGS)

$(TEST_DIR)/test_%: $(TEST_DIR)/test_%.c $(LIB)
	$(CC) $(CFLAGS) $< -o $@ $(LIB) $(LDFLAGS)

run-tests: $(TEST_BINS)
	@for t in $(TEST_BINS); do \
		echo "Running $$t..."; \
		$$t || exit 1; \
	done

clean:
	rm -f $(TARGET) $(TEST_BINS) $(OBJ_LIB) $(OBJ_MAIN) $(LIB)

.PHONY: all run-tests clean