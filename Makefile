CC=mpicc
SRC=parallel_levinson.c
EXEC=par_levinson
CFLAGS=-O3

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(SRC)

clean:
	rm -f $(EXEC)
