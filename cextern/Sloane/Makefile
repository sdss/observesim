all: creconstruct

creconstruct: creconstruct.o cicosagen.o cdot.o
	cc creconstruct.o cicosagen.o cdot.o -lm -o creconstruct
creconstruct.o: creconstruct.c
	cc -c creconstruct.c
cicosagen.o: cicosagen.c
	cc -c cicosagen.c
cdot.o: cdot.c
	cc -c cdot.c

clean:
	rm -f *.o $(OBSERVESIM_DIR)/bin/creconstruct
