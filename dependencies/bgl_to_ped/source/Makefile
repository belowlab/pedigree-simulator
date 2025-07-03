CC=g++
CFLAGS=-O3

all: bgl_to_ped ped_to_bgl search

bgl_to_ped:
	$(CC) -o $@ $@.cpp $(CFLAGS) $(LIBS)

ped_to_bgl:
	$(CC) -o $@ $@.cpp $(CFLAGS) $(LIBS)

search:
	$(CC) -o $@ $@.cpp $(CFLAGS) $(LIBS)

clean:
	rm bgl_to_ped ped_to_bgl search
