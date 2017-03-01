
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


#define pi 3.14159265358979323846  /* pi */

/*MPI_Type Enum*/
enum {
	MASTER = 0,
	FROM_MASTER = 1,
	FROM_WORKER = 2,
};

double estimate_g(double lower_bound, double upper_bound, long long int N)
{
	int		num_tasks;	/* number of tasks     */
	int		task_id;	/* number of processes  */
	int		num_workers;	/* number of worker tasks */
	int		source;		/* rank of sender       */
	int		dest;		/* rank of receiver     */
	int		mtype;		/* message type */
	MPI_Status	status;		/* return status for receive */
	int		rows;
	double		x, y;		/* First boundary condition, Second boundary condition */
	double		sum1 = 0, sum = 0;
	int		i;
	clock_t begin, end;
	double time_spent;

	/* Find out process rank  */
	MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

	/* Find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

	

	num_workers = num_tasks - 1;
	int average_row = N / num_workers;
	int left_overs = N % num_workers;

	if (task_id == MASTER) {
		//begin = clock();
		//Send matrix data to worker tasks
		mtype = FROM_MASTER;
		for (dest = 1; dest < num_tasks; dest++) {
			rows = (dest <= left_overs) ? average_row + 1 : average_row;
			printf("Sending %d rows to task %d\n", rows, dest);
			MPI_Send(&lower_bound, 1, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&upper_bound, 1, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&N, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
		}
		printf("Sent all the data!\n");
		
	}
	else {
		//Worker task
		mtype = FROM_MASTER;
		
		MPI_Recv(&lower_bound, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&upper_bound, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&N, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
		double result, final_val = 0, sum = 0;
		int j;
		srand(time(NULL)*task_id);
		for (j = 0; j < rows; j++)
		{
			result = exp(-1 * pow(2 * ((rand() / (RAND_MAX / (upper_bound - lower_bound))) + lower_bound), 2));
			final_val += result;
		}
		sum1 = (8 * sqrt(2 * pi))*((upper_bound - lower_bound) / N)*final_val;
	}
	return sum1;
}

void collect_results(double *result)
{
	int		num_tasks;	/* number of tasks     */	
	int		task_id;	/* number of processes  */

	double		sum1 = 0, sum = 0;
	int		i;

	/* Find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

	/* Find out process rank  */
	MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

	sum1 = *result;

	//Wait for results from workers
	MPI_Reduce(&sum1, &sum, num_tasks, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (task_id == MASTER)
	{
		printf("\n\n HERE Value of integration is:%lf\n", sum);
	}
}

int main(int argc, char *argv[]) {
	int		num_tasks;	/* number of tasks     */
	int		task_id;	/* number of processes  */
	clock_t begin, end;
	double time_spent;



	//MPI
	/* Start up MPI */
	MPI_Init(&argc, &argv);

	/* Find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

	/* Find out process rank  */
	MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

	float lower_bound = atof(argv[1]);
	float upper_bound = atof(argv[2]);
	long long int N = atof(argv[3]);
	
	begin = clock();

	double result = estimate_g(lower_bound, upper_bound, N);
	collect_results(&result);

	MPI_Finalize();

	if (task_id == MASTER) {

		end = clock();
		time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
		printf("time spent %f sec\n",  time_spent);
	}
	return 0;
}