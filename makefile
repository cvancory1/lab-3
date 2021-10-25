# makefile

program: main.o
	mpicc -std=c99 main.o -o matrixMain
	mpirun -n 5 ./matrixMain

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


# commitTest:
# 	git add timing.txt
# 	git commit -m"send timing"
# 	git push

clean:
	rm *.o
	rm time
	rm matrixMain