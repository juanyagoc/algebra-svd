CC = gcc
CFLAGS = -w -O2 -fPIC
LDFLAGS = -lm

SRCS = householder.c tools.c qr_decomp.c main.c
TARGET = qr_decomposition

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean