
/*

Chloe VanCory
Lab 3

purpose: used for matrix timing outputs .. Currently lab 3 version 



Kalyn - comment out the 
TO RUN : 
  make time 
  make test 

* look at the makefile if you want to change the parameters. Look under the test
  program makes two matricies of size argv[1] , argv[2]




*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "matrixFunctions.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  world = MPI_COMM_WORLD;

  // passing the container for the result in second param
  MPI_Comm_size(world, &worldSize);
  MPI_Comm_rank(world, &rank);
  MPI_Get_processor_name(name, &nameLen);

  if(argc != 3){
    puts("3 arguments need to be passed in to run ... or else will error");
    return 1;

  }

  /* make array of size */ 
  int rows = atoi(argv[1]);
  int cols = atoi(argv[2]);

  Matrix A, B; 

  if (rank == ROOT) {
    initMatrix(&A, rows, cols);
    initMatrix(&B, rows, cols);

    int count = 0;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        // ACCESSA(A,i,j) = rand() % 100;
        // ACCESSB(B,i,j) = rand() % 100;

        ACCESS(A, i, j) = count;
        ACCESS(B, i, j) = count;
        count++;
      }
    }

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        // printf("%4d ", ACCESS(A, i, j));
        // printf("%4d ", ACCESS(B, i, j));
      }
      // printf("\n");
    }
  } else if (rank != ROOT) {
    // for all other nodes define what A is for them, prevents seg faults
    A.data = NULL;
    A.rows = rows;
    A.cols = cols;

    B.data = NULL;
    B.rows = rows;
    B.cols = cols;
  }

  
  Matrix result= matrixAdd(A,B);

  // if(rank == ROOT){
  //   printMatrx(C);
  // }

  result = matrixSub(A,B);

  
  //*** MULTIPLICATION PART */ 
  Matrix F, G;
  int F_rows, F_cols, G_rows, G_cols;
  F_rows = atoi(argv[1]);
  F_cols = atoi(argv[2]);

  G_rows = atoi(argv[1]);
  G_cols = atoi(argv[2]);


  if (rank == ROOT) {
    initMatrix(&F, F_rows, F_cols);
    initMatrix(&G, G_rows, G_cols);
    int count = 1;
    for (int i = 0; i < F.cols * F.rows; i++) {
      F.data[i] = count++;
    }
    count = 1;
    for (int i = 0; i < G.cols * G.rows; i++) {
      G.data[i] = count++;
    }

  } else {
    F.data = NULL;
    F.rows = F_rows;
    F.cols = F_cols;

    G.data = NULL;
    G.rows = G_rows;
    G.cols = G_cols;
  }

  result = multiply(F ,G );


  MPI_Finalize();
  return 0;
}
