CC = g++
CFLAGS  = -O3 -Wall -Winline -Wshadow -std=c++17
TARGET = mgsolve

all: main.cpp 
	$(CC) $(CFLAGS) main.cpp Grid.cpp -o $(TARGET)