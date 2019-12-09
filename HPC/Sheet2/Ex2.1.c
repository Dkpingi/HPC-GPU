#include <mpi.h>
#include <stdio.h>
int main ( int argc, char **argv )
{
 int i,j, size, rank;
 MPI_Status status;
 double starttime, endtime, signal;

 MPI_Init ( &argc, &argv ); //Initialisation of MPI
 MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
 MPI_Comm_size ( MPI_COMM_WORLD, &size );
 int m = 100000;

 signal = 666.0;
 starttime = MPI_Wtime();
 i = rank;
/**
for (j = 0; j < m; j++ ) {
	if (i % 2 == 0)
        {
	MPI_Send ( &signal, 1, MPI_DOUBLE, (rank+1)%size, 0, MPI_COMM_WORLD );
	//printf ( "%d: Sent to %d\n", rank, (rank+1)%size );
	MPI_Recv ( &signal, 1, MPI_DOUBLE, (rank-1 + size)%size, 0, MPI_COMM_WORLD, &status );
	//printf ( "%d: Received from %d\n", rank, (rank-1 + size)%size);
	i++;
	}
        else
	{
	MPI_Recv ( &signal, 1, MPI_DOUBLE, (rank-1 + size)%size, 0, MPI_COMM_WORLD, &status );
	//printf ( "%d: Received from %d\n", rank, (rank-1 + size)%size);
	MPI_Send ( &signal, 1, MPI_DOUBLE, (rank+1)%size, 0, MPI_COMM_WORLD );
	//printf ( "%d: Sent to %d\n", rank, (rank+1)%size );
	i--;
	}
}
**/

for (j = 0; j < 2*m; j++ ) {
	if (i % 2 == 0)
        {
	MPI_Send ( &signal, 1, MPI_DOUBLE, (rank+1)%size, 0, MPI_COMM_WORLD );
	//printf ( "%d: Sent to %d\n", rank, (rank+1)%size );
	}
        else
	{
	MPI_Recv ( &signal, 1, MPI_DOUBLE, (rank-1 + size)%size, 0, MPI_COMM_WORLD, &status );
	//printf ( "%d: Received from %d\n", rank, (rank-1 + size)%size);
	}
	if (j == m-1)
	{
	i++;
	}
}
	
	
 endtime = MPI_Wtime();
 double runtime = endtime-starttime;
 printf ( "%d: Time passed: %.17g\n", rank, runtime );
 printf ( "%d: AVGTimeperMessage passed: %.17g\n", rank, runtime/m );
if (rank != 0)
{
	signal = runtime;
	MPI_Send ( &signal, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
} 
else 
{
	for (j = 1; j < size; j++)
	{
		MPI_Recv ( &signal, 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, &status );
		runtime += signal;
	}
 	printf ( "Final Time passed: %.17g\n", runtime );
 	printf ( "Final AVGTimeperMessage passed: %.17g\n", runtime/(m*size) );
	
}

 MPI_Finalize(); //Deinitialisation of MPI
 return 0;
}
