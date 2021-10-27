
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
  Matrix A, B; 

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

  printf("rank = %d \n", rank);
  Vector b;
  b.length = rows;

  // if(rank == ROOT){
    b.data = malloc(rows * sizeof(double));
  
    int count=1;
    for(int i =0; i< rows ; i++){
      b.data[i] = count++;
      
    }

  // }

    gauss_j(A,b);

  

  // //  do not try to print unless you are the root since the other nodes have
  // //  empty Matricies
  // Matrix D = matrixAdd(A,B);
  // if(rank ==ROOT){
  //   printf("testing adding A +B  matricies\n");
  //   printf("A:\n");
  //   printMatrx(A);

  //   printf("B:\n");
  //   printMatrx(B);
  //   printf("result :\n");

  //   printMatrx(D);
  // }

  // Matrix E = matrixSub(A,B);
  // if(rank ==ROOT){
  //   printf("testing adding A -B  matricies ....will equal the 0 matrix \n");
  //   printMatrx(E);
  // }

  //  testing if unmatched matricies for addition and subtraction
  // if(rank ==ROOT){
  //   printf("\ntesting with two matricies that cannot be added or subtracted\n");
  //   A.cols = B.cols +1;
  //   matrixAdd(A,B);
  // }
  

  // Matrix aT = transpose(A);
  // if (rank == ROOT) {
  //   printf("transpose before\n");
  //   printMatrx(A);

  //   printf("transpose after\n");
  //   printMatrx(aT);
  // }

  // // if(argc!= 2){
  // //   puts("error- multplication requres two arguments");
  // //   return 1; 
  // // }

  // Matrix F, G;
  // int F_rows, F_cols, G_rows, G_cols;
  // F_rows = atoi(argv[1]);
  // F_cols = atoi(argv[2]);
  // G_rows = atoi(argv[1]);
  // G_cols = atoi(argv[2]);


  // if (rank == ROOT) {
  //   initMatrix(&F, F_rows, F_cols);
  //   initMatrix(&G, G_rows, G_cols);
  //   int count = 1;
  //   for (int i = 0; i < F.cols * F.rows; i++) {
  //     F.data[i] = count++;
  //   }
  //   count = 1;
  //   for (int i = 0; i < G.cols * G.rows; i++) {
  //     G.data[i] = count++;
  //   }

  // } else {
  //   F.data = NULL;
  //   F.rows = F_rows;
  //   F.cols = F_cols;

  //   G.data = NULL;
  //   G.rows = G_rows;
  //   G.cols = G_cols;
  // }

  // if (rank == ROOT) {
  //   printf("Testing matrix multiplication F 5x6 , G 6x5 \n");
  //   printf("-------\n");
  //   printf("printing F \n");
  //   printMatrx(F);
  //   printf("\nprinting G \n");
  //   printMatrx(G);
  //   printf("-------\n");
  // }

  // Matrix C = multiply(F, G);
  // if (rank == ROOT){
  //   printf("\nprinting multplication  \n");
  //   printMatrx(C);
  //   printf("-------\n");
  // }



  // // free(A.data);
  // // free(aT.data);
  // // free(B.data);
  // // free(C.data);
  // // free(D.data);
  // // free(E.data);
  // // free(F.data);
  // // free(G.data);




  MPI_Finalize();
  return 0;
}
