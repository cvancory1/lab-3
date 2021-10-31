
/*

Chloe VanCory & Kalyn Howes
Lab 3
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

//   if(argc != 3){
//     puts("3 arguments need to be passed in to run ... or else will error");
//     return 1;

//   }

  /* make array of size */ 
  int rows = atoi(argv[1]);
  int cols = atoi(argv[2]);

  MatrixD A, B; 
  A.rows = rows;
  A.cols = cols;

  B.rows = rows;
  B.cols = 1;

  if(rank ==ROOT){
      A.data = malloc(A.rows * A.cols * sizeof(double));
      B.data = malloc(B.rows * B.cols * sizeof(double));
    srand(time(0));
    for(int i =0; i <A.rows*A.cols;i++){
        A.data[i] = rand() % 100;
        // printf("A.data[i=%f", A.data[i]);
    }

     for(int i =0; i <B.rows*B.cols;i++){
        B.data[i] = rand() % 100;
        // printf("B.data[i=%f", B.data[i]);

    }
  }else{
      A.data = NULL;
      B.data =NULL;

  }

//   gauss_j(A,B);

    double t1 = MPI_Wtime();
    gauss_j(A,B);
    double t2 = MPI_Wtime();
    if (rank == 0) {
     printf("Gauss-Jordan Elimination for a %d by %d matrix took %f seconds\n", A.rows, A.cols, t2 - t1);
    }
  


  MPI_Finalize();
  return 0;
}
