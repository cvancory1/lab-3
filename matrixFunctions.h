
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
  int *data;
} Matrix;

// data = double
typedef struct MatrixD {
  int rows;
  int cols;
  double *data;
} MatrixD;

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
char *bufArr(double *arr, int n) {
  char *buf = malloc(n * (10 + 1) + 1);
  buf[0] = '\0';

  const char *fmt = "%0.4f ";
  char tmp[10 + 1 + 1];  // same width as fmt plus a null term

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
      printf("%2d ", ACCESS(A, i, j));
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
SGData getGaussCounts(int rows, int cols, int worldSize) {
  SGData temp;
  temp.cnts = malloc(worldSize * sizeof(double));
  temp.displs = malloc(worldSize * sizeof(double));

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

  2 dummy functions call this. Addition and subtraction of matricies use the
  same logic
*/
Matrix matrix_add_sub(Matrix A, Matrix B, int operationType) {
  SGData AB_counts = getSGCounts(A.rows, A.cols, worldSize);

  /*  make containers for the nodes to receive their data after the scatter */
  int *recvbufA = malloc(AB_counts.cnts[rank] * sizeof(int));
  int *recvbufB = malloc(AB_counts.cnts[rank] * sizeof(int));

  /*  scatter two matricies A & B  to add*/
  MPI_Scatterv(A.data,                // sendbuf
               AB_counts.cnts,        // sendcnts
               AB_counts.displs,      // displacements
               MPI_INT,               // datatype
               recvbufA,              // recvbuf
               AB_counts.cnts[rank],  // recvcnt
               MPI_INT, ROOT, world);

  MPI_Scatterv(B.data,                // sendbuf
               AB_counts.cnts,        // sendcnts
               AB_counts.displs,      // displacements
               MPI_INT,               // datatype
               recvbufB,              // recvbuf
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
      sumArr,  // sendbuf - address of send buffer,
      AB_counts
          .cnts[rank],  // sendcut-number of elements in send buffer( array )
      MPI_INT,          // stype - data type of send buff elements
      C.data,           // rbuf - address of receive containter
      AB_counts.cnts,  // recvcount- arr of size amount being received from each
                       // process
      AB_counts.displs,  // displs
      MPI_INT,           // data type of recv buffer
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
  SGData A_counts = getSGCounts(A.rows, A.cols, worldSize);

  /*  container for the nodes to receive their a row from Matrix A*/
  int *recvbufA = malloc(A_counts.cnts[rank] * sizeof(int));

  MPI_Scatterv(A.data,               // sendbuf
               A_counts.cnts,        // sendcnts
               A_counts.displs,      // displacements
               MPI_INT,              // datatype
               recvbufA,             // recvbuf
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
      sum += recvbufA[j] * bT.data[j + rowOffset];
    }
    sumArr[i] = sum;
    rowOffset += bT.rows;
    sum = 0;
  }

  SGData C_counts = getSGCounts(C.rows, C.cols, worldSize);
  MPI_Gatherv(
      sumArr,               // sendbuf - address of send buffer,
      C_counts.cnts[rank],  // sendcut-number of elements in send buffer,
      MPI_INT,              // stype - data type of send buff elements
      C.data,               // rbuf - address of receivecontainter
      C_counts.cnts,        // rount- arr of amount being received  from each
      C_counts.displs,      // displs MPI_INT,  //   data type of recv
      MPI_INT, ROOT, world);

  // free (C_counts.displs);
  // free (C_counts.cnts);
  return C;
}

// QUESTION - how does this work if we used the lab2 version where all nodes are
// already working
// how could we use that scattered version where each takes a piece
// returns the vector scaled by value
void scalarByVector(double *arr, int size, int scalar) {
  for (int i = 0; i < size; i++) {
    arr[i] = arr[i] * scalar;
  }
}

// TODO : add error checking if a non-root node calls this function
void gauss_j(MatrixD A, MatrixD B) {
  /*  Error checking  */
  if (A.cols < A.rows) {
    puts(" A's cols must be greater or equal to the rows");
    // return 1;
  }

  /*
      Calc - sendcounts to scatter matrix A
      Allocate for local container A
      scatter A
  */

  // printf("A.rows= %d A.cols= %d \n",A.rows, A.cols );
  SGData A_counts = getGaussCounts(A.rows, A.cols, worldSize);
  SGData B_counts = getGaussCounts(B.rows, B.cols, worldSize);

  // printf("rank = %d sendcounts = %d\n", rank , A_counts.cnts[rank] );
  MatrixD localA;
  localA.rows = A_counts.cnts[rank] / A.rows;
  localA.cols = A.cols;
  localA.data = malloc(A_counts.cnts[rank] * sizeof(double));

  // printf("LOCAL rank = %d  A.rows= %d A.cols= %d  A_counts.cnts[rank]
  // =%d\n",rank, localA.rows,  localA.cols, A_counts.cnts[rank] );

  MatrixD localB;
  localB.cols = B.cols;
  localB.rows = B_counts.cnts[rank];
  localB.data = malloc(B_counts.cnts[rank] * sizeof(double));
  // printf("LOCAL rank = %d  B.rows= %d B.cols= %d  B_counts.cnts[rank]
  // =%d\n",rank, localB.rows,  localB.cols, B_counts.cnts[rank] );

  // if (rank == 0) {
  //   printArr(A);
  // }

  /*  scatters A by rows so each processor can calc and update A in every
   * iteration of outermost for loop  */

  MPI_Scatterv(A.data,               // sendbuf
               A_counts.cnts,        // sendcnts
               A_counts.displs,      // displacements
               MPI_DOUBLE,           // datatype
               localA.data,          // recvbuf
               A_counts.cnts[rank],  // recvcnt
               MPI_DOUBLE, ROOT, world);

  // error check what everyone is recv
  // for(int i =0 ;i< localA.rows * localA.cols ;i++){
  // printf("rank = %d  localA[%d]= %f\n", rank , i, localA.data[i] );
  // }
  // printf("rank = %d  A_counts.cnts[rank]= %d\n", rank , A_counts.cnts[rank]);

  // arr that holds which row of l to send out
  // TODO: bcast this to all nodes OR just have the root check every time

  MPI_Scatterv(B.data,               // sendbuf
               B_counts.cnts,        // sendcnts
               B_counts.displs,      // displacements
               MPI_DOUBLE,           // datatype
               localB.data,          // recvbuf
               B_counts.cnts[rank],  // recvcnt
               MPI_DOUBLE, ROOT, world);

  // printf("rank = %d  B_counts.cnts[rank]= %d\n", rank , B_counts.cnts[rank]);
  // printf("rank = %d  B_counts.cnts[rank]= %d\n", rank , B_counts.displs[rank]);


  // char* buf = bufArr(localB.data,  B_counts.cnts[rank] );
  // printf("Rank %d has pivot row: %s\n", rank, buf);

  /*  using the sendcounts from above create an array which stores which
     proccessor is in charge of a
      corresponding row k . ex:  pivotTracker[3] = 2 ... so row 3 is being calc
     by proccessor 2  */

  // created so all other nodes "know" which rows each processor is in charge of
  int *pivotTracker = malloc(A.rows * sizeof(double));
  if (rank == ROOT) {
    int index = 0;
    int total = A_counts.cnts[index];
    pivotTracker[index] = index;

    for (int i = 0; i < A.rows; i++) {
      total -= A.rows;
      pivotTracker[i] = index;
      if (total == 0) {
        index++;
        total = A_counts.cnts[index];
      }
    }
    if (rank == ROOT) {
      for (int i = 0; i < A.rows; i++) {
        // printf("pivotTracker=%d\n", pivotTracker[i]);
      }
    }
  }
  MPI_Bcast(pivotTracker, A.rows, MPI_DOUBLE, ROOT, world);

  // // /*    begin Alg    */
  int localCount = 0;  // counts how many times a node will be the pivot row
  for (int k = 0; k < A.rows; k++) {  // for each pivot row
    // rank is going to be sending l
    double *pivotRow = malloc(localA.cols * sizeof(double));
    double *pivotRow_B = malloc(localB.cols * sizeof(double));

    MPI_Barrier(world);
    if (rank == pivotTracker[k]) {
      // printf("Enter rank = %d pivot tracker = %d , k = %d localCount = %d\n",
      // rank ,pivotTracker[k] , k, localCount);
      for (int i = 0; i < localA.cols; i++) {
        pivotRow[i] = ACCESS(localA, localCount, i);
        // printf("Enter rank = %d  pivotRow[i] = %f , k = %d \n", rank
        // ,ACCESS(localA, localCount, i), k );
      }
      // char* buf = bufArr(pivotRow, localA.cols);
      // printf("It %d Rank %d has pivot row: %s\n", k, rank, buf); free(buf);


      for(int i =0 ;i< localB.cols;i++){
        pivotRow_B[i] = ACCESS(localB, localCount, i);
      }

      // char* buf = bufArr(pivotRow_B, localB.cols);
      // printf("It %d Rank %d has pivot row: %s\n", k, rank, buf); 


      localCount++;
      if(localCount == localA.rows) localCount =0;
    }

    MPI_Bcast(pivotRow, localA.cols, MPI_DOUBLE, pivotTracker[k], world);
    MPI_Bcast(pivotRow_B, localB.cols, MPI_DOUBLE, pivotTracker[k], world);
    //   // char *arrbuf = bufArr(pivotRow, localA.cols);
    //   // printf("it %d Rank %d received %s\n", k, rank, arrbuf);

    MPI_Barrier(world);
    double denom = pivotRow[k];

    // printf("it %d Enter rank = %d denom = %f \n",k, rank, denom);

    double *l =  malloc(localA.rows *  sizeof(double));  // everyone allocate their local version of l
    for (int i = 0; i < localA.rows; i++) {  // compute l
      l[i] = ACCESS(localA, i, k) / denom;
      // printf("Enter rank = %d l[i] = %f num = %f denom = %f \n",rank, l[i],
      // ACCESS(localA, i , k ),  denom);
    }

    // error checking- L vector to ensure it is correct after every iteration
    // if(rank ==ROOT){
    //   for(int i =0; i < localA.rows ;i++){ // compute l
    //     printf("Enter rank = %d l[i] = %f num = %f denom = %f \n",rank, l[i],
    //     ACCESS(localA, i , k ),  denom);
    //   }

    // }


    // if (k == 0 && rank == 1) {
    MPI_Barrier(world);

    if (rank != pivotTracker[k] ) {
      for (int row = 0; row < localA.rows; row++) {
        for (int cols = 0; cols < localA.cols; cols++) {
          ACCESS(localA, row, cols) = ACCESS(localA, row, cols) - (l[row] * pivotRow[cols]);
          // if(rank == 1 ){
          // printf("sum= %f  \n", ACCESS(localA, row , cols ) );
          // printf( "rank = %d k= %d ACCESS(localA, row , cols )= %f  - l[i]=%f pivotRow[i]=%f \n",rank, k, ACCESS(localA, row, cols), l[row], pivotRow[cols]);
          // printf( "rank = %d Ans= %f\n",rank, ACCESS(localA, row, cols) );
          // }
        }
      }

      for (int col = 0; col < localB.cols; col++) {
        localB.data[col] = localB.data[col] - (l[col] * pivotRow_B[col]);
        // printf("k = %d Rank = %d localB[col]=%f l[col]=%f localB.data[col]=%f \n", k, rank, localB.data[col], l[col], pivotRow_B[col]);
        // printf("Rank = %d Ans= %f\n", rank, localB.data[col]);
      }
    }
    MPI_Barrier(world);
    // sleep(1);
    puts("");
  }

  // char *arrbuf = bufArr(localB.data,localB.rows);
  // printf(" Rank %d received %s\n",  rank, arrbuf);




  /* scale each row of A and b - divide the pivots by themselves to equal 1 */
  for(int i =0 ; i< localA.rows ;i++){
    double pivot;
    pivot = localA.data [ rank + i];
    // printf(" Rank = %d pivot= %f index = %d \n" ,rank, pivot, rank + i);
    ACCESS(localA, i, rank + i) /= pivot;


    for(int j =0 ; j< localB.cols ; j++){
      // printf("arr= %f pivot=%f\n",ACCESS(localB, i , j) , pivot );
      ACCESS(localB, i , j) /=pivot;

    }

  }


  /* gather the b vector back to the root */

  MatrixD result; 
  result.rows = B.rows;
  result.cols = B.cols;
  if(rank == ROOT){
    result.data = malloc(result.rows * result.cols * sizeof(double));
  }else{
    result.data= NULL;
  }

  MPI_Gatherv(
      localB.data,  // sendbuf - address of send buffer,
      B_counts.cnts[rank],  // sendcut-number of elements in send buffer( array )
      MPI_DOUBLE,          // stype - data type of send buff elements
      result.data,           // rbuf - address of receive containter
      B_counts.cnts,  // recvcount- arr of size amount being received from each process
      B_counts.displs,  // displs
      MPI_DOUBLE,           // data type of recv buffer
      ROOT, world);



  // printf("row=%d cols=%d ",result.rows , result.cols);

  if(rank == ROOT){
    for(int i=0; i < result.rows * result.cols ;i++ ){
      for(int j=0; j < result.cols ;j++ ){
        printf("arr=%f \n", ACCESS(result, i , j ) );
      }
    }
    
  }


  


    

}