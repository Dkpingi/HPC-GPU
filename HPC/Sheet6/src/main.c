#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <chTimer.hpp>

//stole this from gpu course
void printOutMatrix (double *matrix, int width) {
	int i;
	for (i = 0; i < width*width; i++) {
		printf ("%4.2g\t", matrix[i%width + (i/width) * width]);
		if ((i+1) % width == 0) printf ("\n");
		}
	printf ("\n");
}
//this also (converted to double)
bool MatrixCompare ( double* P, double* Q, long  matWidth)
{
	long i;

	for ( i = 0; i < matWidth * matWidth; i++ ) {
		//if ( P[i] != Q[i] )
		// Holger 09.04.2014 floating point calculations might have small errors depending on the operation order
		if ( fabs ( ( P[i]-Q[i] ) / ( P[i]+Q[i] ) ) > 1E-05 )
			return false;
	}
	return true;
}

void sequentialRelaxation(double **grid,int gridSize,int nIterations)
{
int i,j,t;
double* old_grid = *grid; 
double* new_grid = (double*) malloc(gridSize*gridSize*sizeof(double));
memcpy(new_grid,old_grid,gridSize*gridSize*sizeof(double));

for(t = 0;t<nIterations;t++)
{
	for(i = 1;i<gridSize-1;i++)
	{
		for(j=1;j<gridSize-1;j++)
		{
			new_grid[i*gridSize + j] = (24.0/100.0)*(old_grid[i*gridSize + j-1] 
	        	 		     + (old_grid[i*gridSize + j]/6.0) + old_grid[i*gridSize + j+1]
					     + old_grid[(i+1)*gridSize + j] + old_grid[(i-1)*gridSize + j]);
		}
	}
	double* temp = old_grid;
	old_grid = new_grid;
	new_grid = temp;
}
*grid = old_grid;
}
int main ( int argc, char **argv )
{

 //double sendTime = 0.0;
 //double recTime = 0.0;

double seqStart = 0.0,seqEnd = 0.0;
double parStart = 0.0,parEnd = 0.0;
 MPI_Init ( &argc, &argv ); //Initialisation of MPI
 //printf("blubb???\n");
 int i,j,t,size,rank;
 //double starttime, endtime;
 MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
 MPI_Comm_size ( MPI_COMM_WORLD, &size );
 
 int gridSize = atoi(argv[1]);
 int nIterations = atoi(argv[2]);
 ChTimer seqTimer,parTimer;
/// Sequential part on rank 0 for comparison ///
double* seqgrid = (double*) malloc(gridSize*gridSize*sizeof(double));
if(rank == 0)
{
 for(i=0;i<gridSize;i++)
 {
	for(j=0;j<gridSize;j++)
	{	
		if(i == 0 && j > floor(gridSize/4.0) && j < ceil(gridSize*3.0/4.0))
		{		
			seqgrid[i*gridSize + j] = 127.0;
		}
		else
		{
			seqgrid[i*gridSize + j] = 0.0;
		}
	}
 }
seqStart = MPI_Wtime();
sequentialRelaxation(&seqgrid,gridSize,nIterations);
seqEnd = MPI_Wtime();
}
MPI_Barrier(MPI_COMM_WORLD);
 //starttime = MPI_Wtime();
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

 ///
 /// Initialize Grid
 ///
 int nRows = floor((float)gridSize/(float)size) + (rank < gridSize%size);
 int nRowsMax = ceil((float)gridSize/(float)size);
 double* grid_old = (double*) malloc(nRowsMax*gridSize*sizeof(double));
 double* grid_new = (double*) malloc(nRowsMax*gridSize*sizeof(double));
 //printf("blubb\n");
 for(i=0;i<nRowsMax;i++)
 {
	int arrayrow = rank+i*size;
	for(j=0;j<gridSize;j++)
	{	
		if(arrayrow == 0 && j > floor(gridSize/4.0) && j < ceil(gridSize*3.0/4.0))
		{		
			grid_old[i*gridSize + j] = 127.0;
			grid_new[i*gridSize + j] = 127.0;
		}
		else
		{
			grid_old[i*gridSize + j] = 0.0;
			grid_new[i*gridSize + j] = 0.0;
		}
	}
 }
MPI_Barrier(MPI_COMM_WORLD);
 if(rank == 0)
 {
	parStart = MPI_Wtime();
 }	
double* upsignal   =   (double*) malloc(gridSize*sizeof(double));
double* downsignal =   (double*) malloc(gridSize*sizeof(double));
double* up         =   (double*) malloc(gridSize*sizeof(double));
double* down       =   (double*) malloc(gridSize*sizeof(double));
 for(t=0;t<nIterations;t++)
 {	
	MPI_Request reqs,reqs2,reqr,reqr2;
	if(rank <= 1)
	{
		memcpy(upsignal, grid_old + gridSize,gridSize*sizeof(double));
	}
	else
	{
		memcpy(upsignal, grid_old,gridSize*sizeof(double));
	}
		
	memcpy(downsignal, grid_old,gridSize*sizeof(double));

	MPI_Isend ( upsignal, gridSize, MPI_DOUBLE, (rank+size-1)%size, 0, MPI_COMM_WORLD,&reqs);		
	MPI_Isend ( downsignal, gridSize, MPI_DOUBLE, (rank + 1)%size, 0, MPI_COMM_WORLD,&reqs2);	
		
	MPI_Irecv ( down, gridSize, MPI_DOUBLE, (rank+1)%size, 0, MPI_COMM_WORLD, &reqr );
	MPI_Irecv ( up, gridSize, MPI_DOUBLE, (rank+size-1)%size, 0, MPI_COMM_WORLD,&reqr2 );
	
	for(i=0;i<nRowsMax;i++)
	{
		int arrayrow = rank+i*size;
		//sendTime -= MPI_Wtime();
		MPI_Wait(&reqs,MPI_STATUS_IGNORE);
		MPI_Wait(&reqs2,MPI_STATUS_IGNORE);
                //sendTime += MPI_Wtime();
		memcpy(downsignal, grid_old + (i+1)*gridSize,gridSize*sizeof(double));
		if(rank <= 1)
		{
			if(i+2 < nRowsMax)
			{
				memcpy(upsignal, grid_old + (i+2)*gridSize,gridSize*sizeof(double));
			}

		}
		else
		{
			memcpy(upsignal, grid_old + (i+1)*gridSize,gridSize*sizeof(double));
		}

		MPI_Isend ( upsignal, gridSize, MPI_DOUBLE, (rank+size-1)%size, i+1, MPI_COMM_WORLD,&reqs);
		MPI_Isend ( downsignal, gridSize, MPI_DOUBLE, (rank + 1)%size, i+1, MPI_COMM_WORLD,&reqs2);

		//recTime -= MPI_Wtime();
		MPI_Wait(&reqr,MPI_STATUS_IGNORE);
		MPI_Wait(&reqr2,MPI_STATUS_IGNORE);
		//recTime += MPI_Wtime();	

		

		for(j=1;j<gridSize-1;j++)
 		{
			if(rank == 0)
			{
				if (i+1 < nRows && (i+1)*size != gridSize-1)
				{
					grid_new[(i+1)*gridSize + j] = (24.0/100.0)*(grid_old[(i+1)*gridSize + j-1] 
								     + (grid_old[(i+1)*gridSize + j]/6.0) 
								     + grid_old[(i+1)*gridSize + j+1] + up[j] + down[j]);
				}
			}
			else if(arrayrow > 0 && arrayrow < gridSize-1)
			{

					grid_new[i*gridSize + j] = (24.0/100.0)*(grid_old[i*gridSize + j-1] 
								 + (grid_old[i*gridSize + j]/6.0) 
								 + grid_old[i*gridSize + j+1] + up[j] + down[j]);
			}

		}

		MPI_Irecv ( down, gridSize, MPI_DOUBLE, (rank+1)%size, i+1, MPI_COMM_WORLD, &reqr );
		MPI_Irecv ( up, gridSize, MPI_DOUBLE, (rank+size-1)%size, i+1, MPI_COMM_WORLD, &reqr2 );

	}
	
	double* temp = grid_old;
	grid_old = grid_new;
	grid_new = temp;
}
MPI_Barrier(MPI_COMM_WORLD);

if(rank == 0)
{
	parEnd = MPI_Wtime();
}

double* buffer = (double*) malloc(nRowsMax*gridSize*sizeof(double));
double* grid = (double*) malloc(gridSize*gridSize*sizeof(double));
 if(rank == 0)
 {

	for(i = 0;i<nRows;i++)
	{	
		if(rank + i*size < gridSize)
		{
		//printf("%d\n",(rank + i*size)*gridSize);
		memcpy(grid + (rank + i*size)*gridSize,grid_old+i*gridSize,gridSize*sizeof(double));
		}
	}	
	for(i = 1; i<size;i++)
	{
		int nRows_i = floor(gridSize/size) + (i < gridSize%size);
		
		//printf("rec:%d %d\n",i,nRows_i);
		MPI_Recv ( buffer, nRows_i*gridSize, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		//printf("%d\n",i);
		for(j = 0;j<nRows_i;j++)
		{
			if(i+j*size < gridSize)
			{
				//printf("%d memcpy\n",i);
				//printf("%d\n",(i + j*size)*gridSize);
				memcpy(grid + (i + j*size)*gridSize,buffer+j*gridSize,gridSize*sizeof(double));
			}
		}
		//printf("%d\n",i);

	}
 } 
 else
 {
	//printf("send:%d %d\n", rank, nRows);
	MPI_Send ( grid_old, nRows*gridSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	//printf("sent:%d %d\n", rank, nRows);
 }

fflush(NULL);
MPI_Barrier(MPI_COMM_WORLD);
double seqTime = (seqEnd - seqStart)*1e6/nIterations;
double parTime = (parEnd - parStart)*1e6/nIterations;
double Speedup = seqTime/parTime;
double Efficiency = Speedup/size;





if(rank == 0)
{

	printf("matrixsize: %d,processes: %d\n",gridSize,size);
	printf("nIterations: %d\n",nIterations);
	printf("seqTimer(ms): %f,parTimer(ms): %f\n",seqTime,parTime);
	printf("Speedup: %f,Efficiency: %f\n\n",Speedup,Efficiency);


bool resultOk = MatrixCompare(grid, seqgrid, 
                static_cast<long>(gridSize));


if (!resultOk) {
	std::cout << "\033[31m***" << std::endl
                      << "*** Error - The two matrices are different!!!" << std::endl
                      << "***\033[0m" << std::endl;
	if(gridSize < 16)
	{
	printOutMatrix(grid,gridSize);
	printOutMatrix(seqgrid,gridSize);
	}
        exit(-1);
        }
}
//sendTime /= nRowsMax*nIterations;
//recTime/=nRowsMax*nIterations;
//printf("rank %d avg sendWaitTime(ms): %g\n",rank,sendTime*1e6);
//printf("rank %d avg recWaitTime(ms): %g\n",rank,recTime*1e6);
//printf("yay!");
 //endtime = MPI_Wtime();
 //double runtime = endtime-starttime;

 MPI_Finalize(); //Deinitialisation of MPI
 
 return 0;
}
