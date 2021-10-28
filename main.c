
/*

Chloe VanCory
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

  

  // int i;
  // initial set up of creating and filling Matrix A & B
  int rows = 5;
  int cols = 5;
  MatrixD A, B; 

  if (rank == ROOT) {
    A.rows = rows;
    A.cols = cols;
    A.data = malloc(A.rows * A.cols * sizeof(double));

    B.rows = rows;
    B.cols = cols;
    B.data = malloc(B.rows * B.cols * sizeof(double));

    int count = 1;
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
        // printf("%4f ", ACCESS(A, i, j));
        // printf("%4f ", ACCESS(B, i, j));

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

  // printf("rank = %d \n", rank);
  Vector b;
  b.length = rows;

  if(rank == ROOT){
    b.data = malloc(rows * sizeof(double));
    int count=1;
    for(int i =0; i< rows ; i++){
      b.data[i] = count++;
    }
  }




  /* when calling make sure there are no 0's off the main diag*/
  // gauss_j(A,b);




  /*Testing Kalyns Matrix */

  MatrixD C;
  C.rows = 3;
  C.cols = 3;
  if (rank == ROOT) {
   
    C.data = malloc(C.rows * C.cols * sizeof(double));

    
    int hardCoded [] = {1,4,2 ,1,2,3,2,1,3};
    for (int i = 0; i < C.cols * C.rows; i++) {
      C.data[i] = hardCoded[i];
      // printf("i = %d data = %f\n", i , C.data[i]);
    }

  } else if (rank != ROOT) {
    // for all other nodes define what A is for them, prevents seg faults
    C.data = NULL;
  }


  puts("");

  MatrixD vectorB;
  vectorB.rows = 3;
  vectorB.cols = 1;

  if(rank ==0 ){
   vectorB.data = malloc( vectorB.rows * vectorB.cols * sizeof(double));
    int vectorHardCoded [] = {11,11,13};
    for (int i = 0; i < vectorB.rows ; i++) {
      vectorB.data[i] = vectorHardCoded[i];
        // printf("i = %d data = %f\n", i , vectorB.data[i]);
    }

  }else{
    vectorB.data =NULL;
  }
  
 

  gauss_j(C,vectorB);





  MPI_Finalize();
  return 0;
}
