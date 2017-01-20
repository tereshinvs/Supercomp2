#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include <omp.h>
double global_end_data[1000][1000];
double* arr;
#include "myfunctions1.h"

int main(int argc, char* argv[]) {
	int procs_rank;
	int procs_size;
	double start;
	double end;	
	
	int data_length = atoi(argv[1]);
	int data_width = atoi(argv[2]);
	float iteration_accuracy = atof(argv[3]);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procs_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &procs_rank);
	//printf("%d of %d total\n",procs_rank,procs_size);


	Initialization(data_length, data_width, procs_rank, procs_size);
	
	start = MPI_Wtime();
	int count = 0;
	int coun = 0;
	double current_accuracy = iteration_accuracy + 2;
	while (current_accuracy >= iteration_accuracy){
		arr = Iteration(count);
		current_accuracy = Iter_accuracy();
		count++;
	}
	
	for(int i = 0; i < field_width - 1; ++i){
	for(int j = 0; j < field_length -1; ++j){
	global_end_data[x_bias + j][y_bias + i] = arr[i * field_length + j];
	printf("%d %d %f\n", x_bias + j, y_bias + i, global_end_data[x_bias + j][y_bias + i]+1);
	}
	}
	double total_accuracy = Total_accuracy();
	MPI_Barrier(MPI_COMM_WORLD);

	end = MPI_Wtime();
	double final_time = end - start;
	
	//ResultCollection (global_end_data,arr,data_length,field_length);
	if(procs_rank == 0){
	//printf("%f",global_end_data[55][50]);
	for (int i =0; i < data_length - 1; ++i) {
		for (int j = 0; j < data_width - 1; ++j) {
		//printf("%f ",global_end_data[i][j]);
		}
		//printf("\n");
		}
	}
	
		//ResultCollection (global_end_data,array_data,data_length,field_length);
		//printf("It takes %f seconds to perform program\n",final_time);
		//printf("Final mistake: %f \n",total_accuracy);

	MPI_Finalize();
	//printf("%f",global_end_data[55][50]);
	return 0;
}