#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
int main ( int argc, char **argv )
{
 int size, rank;
 MPI_Status status;
 MPI_Request req;
 double starttime, endtime;

 MPI_Init ( &argc, &argv ); //Initialisation of MPI
 MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
 MPI_Comm_size ( MPI_COMM_WORLD, &size );
 int m = 10000;
 int j;
 int N = atoi(argv[1]);
 long int messagesize = N*sizeof(char);
 char *signal = (char*) malloc(messagesize);
 starttime = MPI_Wtime();
for (j = 0; j < m; j++ ) {
if(rank == 0)
{
MPI_Send ( signal, N, MPI_CHAR, 1, 0, MPI_COMM_WORLD );
}
else
{
MPI_Recv ( signal, N, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status );
}

}
 endtime = MPI_Wtime();
 double runtimebl = endtime-starttime;

 starttime = MPI_Wtime();
for (j = 0; j < m; j++ ) {
if(rank == 0)
{
MPI_Isend ( signal, N, MPI_CHAR, 1, 0, MPI_COMM_WORLD ,&req);
MPI_Request_free(&req);
}
else
{
MPI_Recv ( signal, N, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status );
}

}
 endtime = MPI_Wtime();
 double runtimenonbl = endtime-starttime;

 free(signal);
if (rank == 0) 
{
 	printf ( "Final Time passed: %.17g\n", runtimebl );
 	printf ( "Bandwidth: %.17g\n", m*messagesize/runtimebl );
	FILE *fp;
	fp = fopen("Ex2resbl.txt","a");
	fprintf (fp, "%li %.17g\n", messagesize, m*messagesize/runtimebl );
	fclose(fp);

 	printf ( "Final Time passed: %.17g\n", runtimenonbl );
 	printf ( "Bandwidth: %.17g\n", m*messagesize/runtimenonbl );
	fp = fopen("Ex2resnonbl.txt","a");
	fprintf (fp, "%li %.17g\n", messagesize, m*messagesize/runtimenonbl );
	fclose(fp);
}

 MPI_Finalize(); //Deinitialisation of MPI
 return 0;
}
