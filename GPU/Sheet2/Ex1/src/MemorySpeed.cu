/*
 *
 * nullKernelAsync.cu
 *
 * Microbenchmark for throughput of asynchronous kernel launch.
 *
 * Build with: nvcc -I ../chLib <options> nullKernelAsync.cu
 * Requires: No minimum SM requirement.
 *
 * Copyright (c) 2011-2012, Archaea Software, LLC.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions 
 * are met: 
 *
 * 1. Redistributions of source code must retain the above copyright 
 *    notice, this list of conditions and the following disclaimer. 
 * 2. Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in 
 *    the documentation and/or other materials provided with the 
 *    distribution. 
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <stdio.h>

#include "chTimer.h"

__global__
void
NullKernel()
{
}

int
main(int argc,char** argv)
{
    const int cIterations = 1000;
    printf( "Measuring memory transfer... " ); //fflush( stdout );
 
    chTimerTimestamp start, stop;
    size_t freemem, totalmem;

    int N = atoi(argv[1]);
    printf("\n Size: %li \n",N*sizeof(float));
    void *dmem;
    cudaMalloc (&dmem,N*sizeof ( float ) ); // Allocate GPU memory 
    void *hmem;
    hmem = (void *) malloc ( N*sizeof ( float ) );     // Allocate CPU memory 
    cudaMemGetInfo(&freemem,&totalmem);

    printf("%li KB free of total %li KB\n",freemem/1024,totalmem/1024); 
    chTimerGetTime( &start );
    for ( int i = 0; i < cIterations; i++ ) {
    // Transfer data from host to device 
    cudaMemcpy ( dmem, hmem, N*sizeof ( float ), cudaMemcpyHostToDevice ); 
    //cudaMemcpy ( hmem, dmem, N*sizeof ( float ), cudaMemcpyDeviceToHost ); 

    }   
 
    //free ( hmem );      // Free host buffer
    cudaThreadSynchronize();
    chTimerGetTime( &stop );
    {
        double microseconds = 1e6*chTimerElapsedTime( &start, &stop );
        double usPerLaunch = microseconds / (float) cIterations;

        printf( "%.2f us\n", usPerLaunch );
	FILE *fp;
	fp = fopen("MemoryMallocH2D.txt","a");
	fprintf (fp, "%li %.5g\n", N*sizeof(float), usPerLaunch );
	fclose(fp);
    }
    chTimerGetTime( &start );
    for ( int i = 0; i < cIterations; i++ ) {
    // Transfer data from host to device 
    //cudaMemcpy ( dmem, hmem, N*sizeof ( float ), cudaMemcpyHostToDevice ); 
    cudaMemcpy ( hmem, dmem, N*sizeof ( float ), cudaMemcpyDeviceToHost ); 

    }   
 
    //free ( hmem );      // Free host buffer
    cudaThreadSynchronize();
    chTimerGetTime( &stop );
    {
        double microseconds = 1e6*chTimerElapsedTime( &start, &stop );
        double usPerLaunch = microseconds / (float) cIterations;

        printf( "%.2f us\n", usPerLaunch );
	FILE *fp;
	fp = fopen("MemoryMallocD2H.txt","a");
	fprintf (fp, "%li %.5g\n", N*sizeof(float), usPerLaunch );
	fclose(fp);
    }
    free(hmem);
    cudaMallocHost(&hmem,N*sizeof ( float ) );

    chTimerGetTime( &start );
    for ( int i = 0; i < cIterations; i++ ) {
    // Transfer data from host to device 
    cudaMemcpy ( dmem, hmem, N*sizeof ( float ), cudaMemcpyHostToDevice ); 
    //cudaMemcpy ( hmem, dmem, N*sizeof ( float ), cudaMemcpyDeviceToHost ); 

    }   
 
    //free ( hmem );      // Free host buffer
    cudaThreadSynchronize();
    chTimerGetTime( &stop );
    {
        double microseconds = 1e6*chTimerElapsedTime( &start, &stop );
        double usPerLaunch = microseconds / (float) cIterations;

        printf( "%.2f us\n", usPerLaunch );
	FILE *fp;
	fp = fopen("MemoryMallocHostH2D.txt","a");
	fprintf (fp, "%li %.5g\n", N*sizeof(float), usPerLaunch );
	fclose(fp);
    }
    chTimerGetTime( &start );
    for ( int i = 0; i < cIterations; i++ ) {
    // Transfer data from host to device 
    //cudaMemcpy ( dmem, hmem, N*sizeof ( float ), cudaMemcpyHostToDevice ); 
    cudaMemcpy ( hmem, dmem, N*sizeof ( float ), cudaMemcpyDeviceToHost ); 

    }   
 
    //free ( hmem );      // Free host buffer
    cudaThreadSynchronize();
    chTimerGetTime( &stop );

    {
        double microseconds = 1e6*chTimerElapsedTime( &start, &stop );
        double usPerLaunch = microseconds / (float) cIterations;

        printf( "%.2f us\n", usPerLaunch );
	FILE *fp;
	fp = fopen("MemoryMallocHostD2H.txt","a");
	fprintf (fp, "%li %.5g\n", N*sizeof(float), usPerLaunch );
	fclose(fp);
    }
    cudaFree ( dmem );  // Free device buffer
    cudaFree(hmem);



    return 0;
}
