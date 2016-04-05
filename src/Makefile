CC = gcc
CCFLAGS = -g -I/usr/local/include
LINK = gcc
LINKFLAGS = -L/usr/local/lib

all:	cartesianRTE

%.o:	%.c
	$(CC) $(CCFLAGS) -c $< -o $@

cartesianRTE:	main.o visit_writer.o	
	$(LINK) $(LINKFLAGS) main.o visit_writer.o -o cartesianRTE -lgsl

clean:
	rm *.o cartesianRTE