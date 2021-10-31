# makefile

program: main.o
	mpicc -std=c99 main.o -o matrixMain
	mpirun -n 3 ./matrixMain

main.o: main.c matrixFunctions.h
	mpicc -std=c99 -c main.c

time: timing.o
	mpicc -std=c99 timing.o -o time		
	time mpirun -n 5 ./time

timing.o: timing.c matrixFunctions.h
	mpicc -std=c99 -c timing.c

test:
	time  mpirun -n 10 ./time 100 100 
	time  mpirun -n 10 ./time 1000 1000
	time  mpirun -n 10 ./time 10000 10000
	time   mpirun -n 10 ./time 20000 20000
	time   mpirun -n 10 ./time 30000 30000
	time   mpirun -n 10 ./time 40000 40000

timeGJ: timingGJ.o
	mpicc -std=c99 timingGJ.o -o timeGJ		
	mpirun -n 5 ./timeGJ 100 100
	

timingGJ.o: timingGJ.c matrixFunctions.h
	mpicc -std=c99 -c timingGJ.c

testGJ:
	time  mpirun -n 10 ./timeGJ 100 100 
	time  mpirun -n 10 ./timeGJ 1000 1000
	time  mpirun -n 10 ./timeGJ 10000 10000
	time   mpirun -n 10 ./timeGJ 20000 20000
	time   mpirun -n 10 ./timeGJ 30000 30000
	time   mpirun -n 10 ./timeGJ 40000 40000

file: fileread.o
	mpicc -std=c99 fileread.o -o fileMain
	mpirun -n 3 ./fileMain

fileread.o: fileread.c matrixFunctions.h
	mpicc -std=c99 -c fileread.c


clean:
	rm *.o
	rm time