CC = gcc
CFLAGS = -w -O2 -fPIC
LDFLAGS = -lm

SRCS = $(wildcard $(SRC_DIR)/*.c) svd.c
SRC_DIR = src/image_compression/qr_decomposition
MAIN_SRC = src/image_compression/main.c
ALL_SRCS = $(SRCS) $(MAIN_SRC)
TARGET = prueba

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean