CC = gcc
CFLAGS = -Wall -Wextra -O2 -fPIC
LDFLAGS = -lm

SRCS = householder.c tools.c qr_decomp.c main.c
OBJS = $(SRCS:.c=.o)
TARGET = qr_decomposition

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

householder.o: householder.c householder.h tools.h
tools.o: tools.c tools.h
qr_decomp.o: qr_decomp.c qr_decomp.h householder.h
main.o: main.c qr_decomp.h tools.h

.PHONY: all clean