///#include "stdafx.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<time.h>
// CPU Calculation to check the output
#define MAX_RANGE 10001
#define MASTER 0


void out_results(int *tst, int Nw, int Nh)
{
	int i, j;
	for (i = 0; i < Nh; i++)
	{
		for (j = 0; j < Nw; j++)
			printf("%d ", tst[i*Nw + j]);
		printf("\n");
	}
}

void initialize_data(int **A, int **Ap,int N)
{
	int i, j, p,rank;
	srand(time(NULL));
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	/*Master node will have full memory allocation*/
	//if (rank == 0)
	{
		*A = (int *)malloc(N * N * sizeof(int));
		*Ap = (int *)calloc(N * N ,sizeof(int));
	}
	
	for (i = 0; i < N ; i++)
		for (j = 0; j < N ; j++)
			(*A)[i + N*j] =  rand()% MAX_RANGE;
/*	if (rank == 0)
	{   
		printf("input matrix \n");
		out_results(*A,N,N);
	}
	*/
}
void scatter_data(int *A,int N)
{
	int i, p, rank,compute_size;
	int *offsets, *chunk_sizes;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	offsets = (int *)malloc(p*sizeof(int)); /*offsets from which data is sent of input A from master*/
	chunk_sizes = (int *)malloc(p*sizeof(int)); /*sizes of the data chunks sent of input A from master*/
	
	compute_size = ((N - 2) / p);
	for (i = 0; i<p; ++i) {
		 /*data chunks that would be computed at each node*/
		chunk_sizes[i] = (compute_size + 2)*N + (i==0)*((N-2)%p)*N;
		offsets[i] =  i*compute_size*N+(i!=0)*((N-2)%p)*N;
	}
	if (((N - 2) / p) > 0)
	/*Scatter data from the master to the worker*/
	MPI_Scatterv(A, chunk_sizes, offsets, MPI_INT, A, (compute_size+2)*N, MPI_INT,
		MASTER, MPI_COMM_WORLD);
		
	free(chunk_sizes);
	free(offsets);
}
void mask_operation(int *A, int N, int * Ap)
{
	int i, j, p, rank, compute_size;
	int *res ;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	/*Values converted to float before computation to avoid accumulation from overflowing as MAX_RANGE 10000 */
	/*Offset added to avoid first row overwrite for master*/
	if (rank == 0)
	{
		compute_size = ((N - 2) / p)+((N-2)%p);/*data chunks + leftovers that would be computed at each node*/
		res = Ap +N;
	}
	else
	{
		compute_size = ((N - 2) / p);/*data chunks that would be computed at each node*/
		res = Ap;
	}
	for (i = 1; i < compute_size+1; i++)
	{
		for (j = 1; j < N - 1; j++)
		{
			res[(i-1)*N + j] =(int)( ((float) A[(i - 1) * N + j - 1] + (float)A[(i - 1)*N + j] + (float)A[(i - 1)*N + j + 1] + (float)A[i *N + j - 1] +
				2 * ((float)A[i*N + j]) + (float)A[i*N + j + 1] + (float)A[(i + 1)*N + j - 1] + (float)A[(i + 1)*N + j] + (float)A[(i + 1)*N + j + 1]) / 10);
		}
	}
	
}
void gather_results(int *Ap, int N) {
	int i, p, rank, compute_size;
	int *offsets, *chunk_sizes;
	int *res;
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	offsets = (int *)malloc(p*sizeof(int));/*offsets from which data is sent of input A from master*/
	chunk_sizes = (int *)malloc(p*sizeof(int));
	compute_size = ((N - 2) / p);
	for (i = 0; i<p; ++i) {
		 /*data chunks that would be computed at each node*/
		chunk_sizes[i] = (compute_size + 2)*N + (i==0)*((N-2)%p)*N;
		offsets[i] =  i*compute_size*N+(i!=0)*((N-2)%p)*N+N;
	}
	/*Offset added to avoid first row overwrite for master*/
	if (rank == 0)
		res = Ap +N;
	else
		res = Ap;
	// Gather the geenerated outputs
	if (((N - 2) / p) > 0) 	
	MPI_Gatherv(res, compute_size*N, MPI_INT, Ap, chunk_sizes, offsets, MPI_INT,
		MASTER, MPI_COMM_WORLD);
/*	if (rank == 0)
	{
		printf("Output Matrix\n");
		out_results(Ap, N,N);
	}*/
	free(chunk_sizes);
	free(offsets);
}
int main(int argc, char** argv) {

	int * A, *Ap,*tst;
	int p, rank;
	clock_t begin,end;
	double time_spent;
	MPI_Init(&argc, &argv);
	int N = atoi(argv[1]);
	begin = clock();
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	initialize_data(&A, &Ap,N);
	scatter_data(A, N);
	mask_operation(A, N, Ap);
	gather_results(Ap, N);
	if (rank == MASTER) {

		end = clock();
		time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("time spent %f sec\n",  time_spent);
	}
	MPI_Finalize();
	
	free(A);
	free(Ap);
	return 0;
}
