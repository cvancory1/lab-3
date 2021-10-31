
/*

  Chloe VanCory & Kalyn Howes 
  COSC 420 
  10.31.21
  Lab 3 - Algorithm 1 Gauss-Jordan Elimination

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
  MPI_File fh;


  char * fname = "output.txt";
  MPI_File_open(
    world,  // comm
    fname,  // filename
    MPI_MODE_CREATE | MPI_MODE_RDWR, // mode
    MPI_INFO_NULL, // info structure
    &fh );

 
  // MatrixD A, B; 
//   A.rows = 5;
//   A.cols = 5;
//   B.rows = 5;
//   B.cols = 1;

//   if (rank == ROOT) {
//     A.data = malloc(A.rows * A.cols * sizeof(double));
//     B.data = malloc(B.rows * B.cols * sizeof(double));

//     int count = 1;
//     for (int i = 0; i < A.rows* A.cols; i++) {
//       A.data[i] = count++;
//     }

//     count = 1;
//     for(int i =0; i< B.cols* B.rows ;i++){
//       B.data[i] = count++;

//     }

//   } else  {
//     A.data = NULL;
//     B.data = NULL;
//   }

// if(rank ==ROOT){
//     printf("Testing GJ\n");
//     printf("Printing Matricies\n");
//     printMatrxD(A);
//     puts("");
//     printMatrxD(B);
//     puts("");

//   }
//   gauss_j(A,B);


//   if(rank ==ROOT){
//     printf("RESULT MATRIX\n");
//     printMatrxD(result);
//   }





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


  // puts("");

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

  if(rank ==ROOT){
    printf("Testing GJ\n");
    printf("Printing Matricies\n");
    printMatrxD(C);
    puts("");
    printMatrxD(vectorB);
    puts("");

  }

    //  result =
   gauss_j(C,vectorB);


  // if(rank ==ROOT){
  //   printf("RESULT MATRIX\n");
  //   printMatrxD(result);
  // }


  // MatrixD E, vectorC ;
  // E.rows = 5;
  // E.cols = 5;
  // vectorC.rows= 5;
  // vectorC.cols= 1;

  // if(rank == ROOT){
  //   E.data = malloc(E.cols * E.rows *sizeof(double));
  //   int hardCoded[] = {1,2,3,6,5,
  //                       0,3,2,2,3,
  //                       -1,3,4,3,5,
  //                       0,3,2,4,3,
  //                       1,6,3,8,-3};
  //   for(int i =0; i < E.cols * E.rows; i++ ){
  //     E.data[i] = hardCoded[i];
  //   }

  //   int hardCoded2[] = {192,46,23,384,960 };
  //   for(int i =0; i < vectorC.rows*vectorC.cols ;i++){
  //     vectorC.data[i] = hardCoded2[i];

  //   }

  // }else{
  //   E.data = NULL;
  // }

  // if(rank ==ROOT){
  //   // printf("Testing GJ\n");
  //   // printf("Printing Matricies\n");
  //   // printMatrxD(E);
  //   // puts("");
  //   // printMatrxD(vectorC);
  //   // puts("");

  // }

  //  gauss_j(E,vectorC);
  



  MPI_Finalize();
  return 0;
}
