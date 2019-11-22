#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>

#define bool int
#define false 0
#define true -1

long calcTime ( struct timeval* tp1, struct timeval* tp2 )
{
	//calculate time passed measured with gettimeofday
	//returns total usecs
	long usecs;

	usecs = (tp2->tv_sec - tp1->tv_sec) * 1E6 + tp2->tv_usec - tp1->tv_usec;

	return usecs;
}


void printTime ( char *label, struct timeval *tp1, struct timeval *tp2, long bytes, long mflops )
{
	//calculate and print out passed time measured with gettimeofday
	long usecs;

	usecs = calcTime (tp1, tp2);
	if ( bytes != 0 )
		printf ( "%s: %4ld usecs passed (%.2lf MB/s)\n", label, usecs, (double) bytes/usecs );
	else if ( mflops != 0 )
		printf ( "%s: %4ld usecs passed (%ld MFLOP, %.2lf MFLOP/s)\n", label, usecs, mflops, (double) mflops/usecs*1E06 );
	else
		printf ( "%s: %4ld usecs passed\n", label, usecs );
	fflush (stdout);
}

bool MatrixCompare ( float* P, float* Q, long  matWidth)
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

void MatrixMulOnHost(float* M, float* N, float* P, long matWidth)
{  
	long i, j, k;
	
	for ( i = 0; i < matWidth; ++i) {
		for ( j = 0; j < matWidth; ++j) {
			float sum = 0;
			for ( k = 0; k < matWidth; ++k) {
				float a = M[i * matWidth + k];
				float b = N[k * matWidth + j];
				//printf ("P[%ld][%ld] += M[%ld][%ld] * N[%ld][%ld]\n", j, i, k, i, j, k );
				sum += a * b;
			}
			P[i * matWidth + j] = sum;
		}
	}
}

long findMin ( long a, long b )
{
	if ( a <= b )
		return a;
	else
		return b;
}

void MatrixMulOnHostBlocked(float* M, float* N, float* P, long matWidth, long blockSize)
{
	long ii, jj, kk, i, j, k;
	float temp;

	//printf ("matWidth = %ld, blockSize = %ld\n", matWidth, blockSize );
	for (ii = 0; ii < matWidth; ii += blockSize) {
		for (jj = 0; jj < matWidth; jj += blockSize) {
			for (kk = 0; kk < matWidth; kk += blockSize) {
				for (i = ii; i < findMin(ii+blockSize, matWidth); i++) {
					for (j = jj; j < findMin(jj+blockSize, matWidth); j++) {
						temp = 0;
						for (k = kk; k < findMin(kk+blockSize, matWidth); k++) {
							//if ( j == 1 && i == 2 ) {
							//printf ("P[%ld][%ld] += M[%ld][%ld] * N[%ld][%ld]\n", j, i, k, i, j, k );
							//printf ("P[%ld][%ld] += %.2f * %.2f\n", j, i,  M[i * matWidth + k], N[k * matWidth + j] );
							//}
							temp += M[i * matWidth + k] * N[k * matWidth + j];
						}
						P[ i * matWidth + j] += temp;
						//printf ("P[%ld][%ld]=%.2f\n", j, i, P[ i * matWidth + j] );
					}
				}
			}
		}
	}  
}



void printOutMatrix (float *matrix, int width) {
	int i;
	for (i = 0; i < width*width; i++) {
		printf ("%4.2f\t", matrix[i%width + (i/width) * width]);
		if ((i+1) % width == 0) printf ("\n");
		}
	printf ("\n");
}

/*int main (int argc, char **argv)
{	
	int optDebug = false; //option for debug output
	float *M, *N, *P;
	float *Q;
	long matWidth; 			//size of the matrix in elements per dimension (matrix = square of Width * Width)	
	long matrixSize;		//size of the matrix in bytes
	struct timeval tp1, tp2;

	//command line parsing
	int factor = 1; 
	char *pos = NULL;	
	if (argc < 3) {
	  printf ("Usage: %s <problem size{k,M,G}>\n", argv[0]);
		exit (0);
	}
	
	if (argc == 4)
		if (strcmp (argv[3], "-debug") == 0) optDebug = true;
	
	pos = strrchr (argv[1], 'k');
	if (pos != NULL) {
		factor = 1024;
		*pos = '\0'; //terminate input string here
	}
	pos = strrchr (argv[1], 'M');
	if (pos != NULL) {
		factor = 1024*1024;
		*pos = '\0'; //terminate input string here
	}
	pos = strrchr ( argv[1], 'G' );
	if ( pos != NULL ) {
		factor = 1024*1024*1024;
		*pos = '\0'; //terminate input string here
	}
	matWidth = atol ( argv[1] );
	matWidth *= factor;	

	matrixSize = matWidth * matWidth * sizeof ( float );
	printf ( "Matrix size = %d x %d elements (total %ld MB or %ld)\n", matWidth, matWidth, 3 * matrixSize/1024/1024, matrixSize );
	
	//allocate memory
	M = ( float* ) malloc ( matrixSize );
	N = ( float* ) malloc ( matrixSize );
	P = ( float* ) malloc ( matrixSize );
	Q = ( float* ) malloc ( matrixSize );
	if ( M == NULL || N == NULL || P == NULL || Q == NULL ) {
		printf ("malloc failed\n");
		return 1;
	}
		
	// initialize matrices
	long i;
	for (i = 0; i < matWidth*matWidth; i++) {
		M[i] = i % matWidth;
		N[i] = (i % matWidth == i / matWidth) ? 1 : 0 ;		
		P[i] = 0;
	}

	long mflops = ( long ) matWidth * matWidth * matWidth * 2 / 1E06;
	
	/////////////////////////////////////
	// calculate on host
	/////////////////////////////////////  	
	gettimeofday (&tp1, NULL);
	MatrixMulOnHost (M, N, Q, matWidth);
	gettimeofday (&tp2, NULL);
	printTime ( "host calculation (naive)", &tp1, &tp2, 0, mflops );

	if (optDebug) {
		printf ("Matrix Q ===============\n");
	        printOutMatrix (Q, matWidth);
	}

	long blockSize;
	char output[50];
	for (blockSize = 1; blockSize <= 32; blockSize++) {
		for (i = 0; i < matWidth*matWidth; i++) // initialize P
			P[i] = 0;

	 	gettimeofday (&tp1, NULL);
 		MatrixMulOnHostBlocked (M, N, P, matWidth, blockSize );
		gettimeofday (&tp2, NULL);
		sprintf ( output, "host calculation (blocked by %dx%d)", blockSize, blockSize );
		printTime ( output, &tp1, &tp2, 0, mflops );

		if (optDebug) {
			printf ("Matrix P ===============\n");
	                printOutMatrix (P, matWidth);
		}
		//todo Compare P with reference result Q
		if ( !MatrixCompare ( P, Q, matWidth) )
			printf ("Matrix comparison failed!\n");
	}

	/////////////////////////////////////
	// debug output
	/////////////////////////////////////  	
	if ( optDebug ) {
		printf ("Matrix M ===============\n");
		printOutMatrix (M, matWidth);
		printf ("Matrix N ===============\n");
		printOutMatrix (N, matWidth);
		printf ("Matrix P ===============\n");
		printOutMatrix (P, matWidth);
		printf ("Matrix Q ===============\n");
		printOutMatrix (Q, matWidth);
	}
	
	free(M);
	free(N);
	free(P);
	free(Q);
	return 0;
}*/
