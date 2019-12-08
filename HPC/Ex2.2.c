#include <mpi.h>
#include <stdio.h>
int main ( int argc, char **argv )
{
 int size, rank;
 MPI_Status status;
 double starttime, endtime, signal;

 MPI_Init ( &argc, &argv ); //Initialisation of MPI
 MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
 MPI_Comm_size ( MPI_COMM_WORLD, &size );
 int m = 1000;
 int j;
 signal = 666.0;
 starttime = MPI_Wtime();
void barrier(MPI_Comm comm)
{
	int rank,size;
	MPI_Comm_rank ( comm, &rank );
	MPI_Comm_size ( comm, &size );
	int i;
	if (rank == 0)
	{
		for (i = 1; i < size; i++)
		{
		MPI_Recv ( &signal, 1, MPI_DOUBLE, i, 0, comm, &status );
		}
		for (i = 1; i < size; i++)
		{
		MPI_Send ( &signal, 1, MPI_DOUBLE, i, 0, comm );
		}
	}
	else
	{
		MPI_Send ( &signal, 1, MPI_DOUBLE, 0, 0, comm );
		MPI_Recv ( &signal, 1, MPI_DOUBLE, 0, 0, comm, &status );
	}
}
for (j = 0; j < m; j++ ) {
//MPI_Barrier(MPI_COMM_WORLD);
//printf ( "%d: At Barrier %d \n",rank, j);
barrier(MPI_COMM_WORLD);
//printf ( "%d: Past Barrier %d \n", rank, j);
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
 	printf ( "Final AVGTimeperMessage passed: %.17g\n", runtime/m );
	FILE *fp;
	fp = fopen("Ex2.2res.txt","a");
	fprintf (fp, "%d %.17g\n", size, runtime );
	fclose(fp);
}

 MPI_Finalize(); //Deinitialisation of MPI
 return 0;
}
