
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define ROOT 0

#define INDEX(i, j, n, m) i *m + j
#define ACCESS(A, i, j) A.data[INDEX(i, j, A.rows, A.cols)]

// declare the MPI set ups in the global scope
MPI_Comm world;
int worldSize, rank;
char name[MPI_MAX_PROCESSOR_NAME];
int nameLen;


typedef struct Matrix {
  int rows;
  int cols;
  double *data;
} Matrix;

typedef struct Vector {
  double *data;
  int length;
} Vector;


void initMatrix(Matrix *A, int rows, int cols) {
  A->rows = rows;
  A->cols = cols;
  A->data = malloc(A->rows * A->cols * sizeof(int));
}

  // // to error check 
  // char *arrbuf = bufArr(recvbufA, sendcts[rank]);
  // printf("Rank %d received %s\n", rank, arrbuf);
  // arrbuf = bufArr(recvbufB, sendcts[rank]);
  // printf("\nRank %d received %s\n", rank, arrbuf);
char *bufArr(int *arr, int n) {
  char *buf = malloc(n * (4 + 1) + 1);
  buf[0] = '\0';

  const char *fmt = "%4d ";
  char tmp[4 + 1 + 1];  // same width as fmt plus a null term

  for (int i = 0; i < n; i++) {
    sprintf(tmp, fmt, arr[i]);
    strcat(buf, tmp);
  }

  return buf;
}

void printArr(Matrix A) {
  for (int i = 0; i < A.rows * A.cols; i++) {
    printf("%4d ", A.data[i]);
  }
  printf("\n");
}

// void printARR(int* A) {
//   for (int i = 0; i < A.rows * A.cols; i++) {
//     printf("%4d ", A[i]);
//   }
//   printf("\n");
// }

void printMatrx(Matrix A) {
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.cols; j++) {
      printf("%2f ", ACCESS(A, i, j));
    }
    printf("\n");
  }
}

/* To automate this garbage and store the results in a nice struct */
typedef struct SGData {
  int *cnts;
  int *displs;
} SGData;

/*
Take the dimensions of a matrix and do the calculations to either
scatter or gather the rows
 */
SGData getSGCounts(int rows, int cols, int worldSize) {
  SGData temp;
  temp.cnts = malloc(worldSize * sizeof(int));
  temp.displs = malloc(worldSize * sizeof(int));

  int minSend = rows / worldSize;
  for (int i = 0; i < worldSize; i++) {
    temp.cnts[i] = minSend;
  }

  for (int i = 0; i < rows % worldSize; i++) {
    temp.cnts[i]++;
  }

  for (int i = 0; i < worldSize; i++) {
    temp.cnts[i] *= cols;
    // printf("sendcount=%d\n",temp.cnts[i]);
  }

  temp.displs[0] = 0;
  for (int i = 1; i < worldSize; i++) {
    temp.displs[i] = temp.displs[i - 1] + temp.cnts[i - 1];
  }

  return temp;
}

/*

  2 dummy functions call this. Addition and subtraction of matricies use the same logic 
*/
Matrix matrix_add_sub(Matrix A, Matrix B, int operationType) {
  
  SGData AB_counts = getSGCounts(A.rows , A.cols , worldSize);

  /*  make containers for the nodes to receive their data after the scatter */
  int *recvbufA = malloc(AB_counts.cnts[rank] * sizeof(int));
  int *recvbufB = malloc(AB_counts.cnts[rank] * sizeof(int));

  /*  scatter two matricies A & B  to add*/
  MPI_Scatterv(A.data,         // sendbuf
               AB_counts.cnts,        // sendcnts
               AB_counts.displs,         // displacements
               MPI_INT,        // datatype
               recvbufA,       // recvbuf
               AB_counts.cnts[rank],  // recvcnt
               MPI_INT, ROOT, world);

  MPI_Scatterv(B.data,         // sendbuf
               AB_counts.cnts,        // sendcnts
               AB_counts.displs,         // displacements
               MPI_INT,        // datatype
               recvbufB,       // recvbuf
               AB_counts.cnts[rank],  // recvcnt
               MPI_INT, ROOT, world);

  /* each node computes a row of  Matrix C */ 
  int *sumArr = malloc(AB_counts.cnts[rank] * sizeof(int));
  if (operationType == 0) {
    for (int i = 0; i < AB_counts.cnts[rank]; i++) {
      sumArr[i] = recvbufA[i] + recvbufB[i];
      // printf("rank = %d sumArr[%d] = %d\n", rank, i, sumArr[i]);
    }
  } else if (operationType == 1) {
    for (int i = 0; i < AB_counts.cnts[rank]; i++) {
      sumArr[i] = recvbufA[i] - recvbufB[i];
      // printf("rank = %d sumArr[%d] = %d\n", rank, i, sumArr[i]);
    }
  }

  /* Matrix C gets made for the gather */ 
  Matrix C;
  if (rank == ROOT) {
    initMatrix(&C, A.rows, A.cols);
  } else {
    C.data = NULL;
    C.rows = A.rows;
    C.cols = A.cols;
  }

  MPI_Gatherv(
      sumArr,         // sendbuf - address of send buffer,
      AB_counts.cnts[rank],  // sendcut-number of elements in send buffer( array )
      MPI_INT,        // stype - data type of send buff elements
      C.data,         // rbuf - address of receive containter
      AB_counts.cnts,        // recvcount- arr of size amount being received from each
                      // process
      AB_counts.displs,         // displs
      MPI_INT,        // data type of recv buffer
      ROOT, world);


  return C;
}

/*
  dummy function for main to access and add 2 matricies
  call matrix_add_sub with a 1 to denote addition 
*/ 
Matrix matrixAdd(Matrix A, Matrix B) {
  Matrix C;
  C.data = NULL;
  C.cols = 0;
  C.rows = 0;
  if (A.rows != B.rows && rank == ROOT) {
    printf(
        "ERROR - rows of A and B do not match.. returning an empty matrix \n");
    return C;
  } else if (A.cols != B.cols && rank == ROOT) {
    printf(
        "ERROR - cols of A and B do not match.. returning an empty matrix \n");
    return C;
  } else {
    return matrix_add_sub(A, B, 0);
  }
}

/*
  dummy function for main to access and subtract 2 matricies
  call matrix_add_sub with a 1 to denote subtraction 
*/ 

Matrix matrixSub(Matrix A, Matrix B) {
  Matrix C;
  C.data = NULL;
  C.cols = 0;
  C.rows = 0;

  if (A.rows != B.rows && rank == ROOT) {
    printf(
        "ERROR - rows of A and B do not match.. returning an empty matrix \n");
    return C;
  } else if (A.cols != B.cols && rank == ROOT) {
    printf(
        "ERROR - cols of A and B do not match.. returning an empty matrix \n");
    return C;
  } else {
    return matrix_add_sub(A, B, 1);
  }
}

Matrix transpose(Matrix A) {
  int rows = A.rows;
  int cols = A.cols;
  // printf("in transpose function \n%d - rows %d cols %d \n",rank,rows,cols);

  // transpose of a 1x1
  if ((rows == 1 && cols == 1)) {
    return A;
  } else if (A.data == NULL) {
    Matrix C;
    C.data = NULL;
    C.rows = A.cols;
    C.cols = A.rows;
    return C;
  }

  Matrix B;
  initMatrix(&B, cols, rows);

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      ACCESS(B, i, j) = ACCESS(A, j, i);
      // printf("1= %d 2= %d\n", ACCESS(B, i, j), ACCESS(A, j, i));
    }
  }
  return B;
}

/*
  scatters rows of matrix A , and bcasts matrix B to all other nodes 
  calcs the multplicatation for all nodes afterwards
*/
Matrix multiply(Matrix A, Matrix B) {
  Matrix C;

  int A_rows = A.rows;
  int A_cols = A.cols;
  int B_rows = B.rows;
  int B_cols = B.cols;

  if (A.cols != B.rows && rank == ROOT) {
    printf(
        "ERROR - cols of A and B do not match.. returning an empty matrix \n");
    C.data = NULL;
    C.rows = 0;
    C.cols = 0;
    return C;
  }

  /*calculate B transpose before scattering*/
  Matrix bT = transpose(B);

  if (rank == ROOT) {
    initMatrix(&C, A.rows, B.cols);
  } else {
    C.data = NULL;
    C.rows = A.rows;
    C.cols = B.cols;
  }
  SGData A_counts = getSGCounts(A.rows , A.cols , worldSize);

  /*  container for the nodes to receive their a row from Matrix A*/ 
  int *recvbufA = malloc(A_counts.cnts[rank] * sizeof(int));

    MPI_Scatterv(A.data,          // sendbuf
                 A_counts.cnts,        // sendcnts
                 A_counts.displs,         // displacements
                 MPI_INT,         // datatype
                 recvbufA,        // recvbuf
                 A_counts.cnts[rank],  // recvcnt
                 MPI_INT, ROOT, world);

    /* make container for bT which gets broadcasted by all other nodes*/
    if (rank != ROOT) {
      bT.data = malloc(bT.rows * bT.cols * sizeof(int));
    }
    MPI_Bcast(bT.data, bT.rows * bT.cols, MPI_INT, ROOT, world);

    int sum = 0;
    int *sumArr = malloc(C.cols * sizeof(int));
    int rowOffset = 0;
    for (int i = 0; i < C.rows; i++) {
      for (int j = 0; j < C.cols; j++) {
        sum += recvbufA[j] * bT.data[ j + rowOffset];
      }
      sumArr[i] = sum;
      rowOffset += bT.rows;
      sum =0;
      
    }

    SGData C_counts = getSGCounts(C.rows , C.cols , worldSize);
    MPI_Gatherv(sumArr,         // sendbuf - address of send buffer,
              C_counts.cnts[rank],  // sendcut-number of elements in send buffer, 
              MPI_INT,        // stype - data type of send buff elements
              C.data,         // rbuf - address of receivecontainter
              C_counts.cnts,        // rount- arr of amount being received  from each
              C_counts.displs,         // displs MPI_INT,  //   data type of recv
              MPI_INT, ROOT, world);


  // free (C_counts.displs);
  // free (C_counts.cnts);
  return C;
}


//TODO : add error checking if a non-root node calls this function 
void gauss_j(Matrix A ,Vector b ){
  Vector l; 
  l.length = A.rows;
  l.data = malloc(A.rows * sizeof(double));
  printf("l - length = %d \n", l.length);

  if(A.cols < A.rows){
    puts(" A's cols must be greater or equal to the rows");
    // return 1;
  }
  
  // todo: is there an edge case where the top left column is 0 ? 
  for(int k =0 ; k < A.rows ; k++){ // for each pivot row
      double denom;
      if(rank == ROOT ){
        denom = ACCESS(A,k,k);
        // printf("k =%f denom = %f ", k, denom );
      }
    for(int i =0 ; i < A.rows ; i++){ // for 1, ... , n 
      if(rank == ROOT ){
        l.data[i] = ACCESS(A,i,k) / denom; 
        printf("l.data[i] = %f ACCESS(A,i,k)= %f denom=%f \n", l.data[i], ACCESS(A,i,k), denom);

      }
    }
    // puts(" here\n");

    MPI_Bcast(l.data , l.length , MPI_INT , ROOT, world);
    printf("rank  = %d \n", rank);


    // for(int i = 0 ; i < l.length ;i++){
    //   printf("rank  = %d data=%d ", rank, l.data[i]);
    // }

  }

}