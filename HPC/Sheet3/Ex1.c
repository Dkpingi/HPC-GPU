#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
int main ( int argc, char **argv )
{
 int size, rank;
 MPI_Status status;
 double starttime, endtime;

 MPI_Init ( &argc, &argv ); //Initialisation of MPI
 MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
 MPI_Comm_size ( MPI_COMM_WORLD, &size );
 int m = 10000;
 int j;
 int N = atoi(argv[1]);
 char *signal = (char*) malloc(N*sizeof(char));
 starttime = MPI_Wtime();
for (j = 0; j < m; j++ ) {
if(rank == 0)
{
MPI_Send ( signal, N, MPI_CHAR, 1, 0, MPI_COMM_WORLD );
//printf("%d: Sent", rank);
MPI_Recv ( signal, N, MPI_CHAR, 1, 0, MPI_COMM_WORLD, &status );
//printf("%d: Received", rank);
}
else
{
MPI_Recv ( signal, N, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status );
//printf("%d: Sent", rank);
MPI_Send ( signal, N, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
//printf("%d: Received", rank);
}

}

	
 endtime = MPI_Wtime();
 double runtime = endtime-starttime;
 free(signal);
if (rank == 0) 
{
 	printf ( "Final Time passed: %.17g\n", runtime );
 	printf ( "Final AVGTimeperMessage passed: %.17g\n", runtime/m );
	FILE *fp;
	fp = fopen("Ex1res.txt","a");
	fprintf (fp, "%li %.17g\n", N*sizeof(char), runtime/m );
	fclose(fp);
}

 MPI_Finalize(); //Deinitialisation of MPI
 return 0;
}
