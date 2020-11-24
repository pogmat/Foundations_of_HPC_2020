#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv)
{
	long int buffer;
	int total_nodes, my_id;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &total_nodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

	if (my_id == 0) {
		buffer = 1;
		MPI_Send(&buffer, 1, MPI_LONG, 1, 10, MPI_COMM_WORLD);
		MPI_Recv(&buffer, 1, MPI_LONG, total_nodes - 1, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("The factorial of %d is %ld\n", total_nodes - 1, buffer);
	} else {
		MPI_Recv(&buffer, 1, MPI_LONG, my_id - 1, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		buffer *= my_id;
		MPI_Send(&buffer, 1, MPI_LONG, (my_id == total_nodes - 1) ? 0 : my_id + 1, 10, MPI_COMM_WORLD);
	}


	MPI_Finalize();

	return 0;
}
