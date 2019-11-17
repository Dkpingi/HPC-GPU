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
Init()
{
}

int
main(int argc,char** argv)
{
    (void) argc; //so compiler stops being annoying
    const int cIterations = 100000000;
    printf( "Measuring memory transfer... " ); //fflush( stdout );
 
    chTimerTimestamp start, stop;

    long int size = atoi(argv[1])*sizeof(char);
    char *dmem = NULL; cudaMalloc((void**) dmem,size);
    char *dmem2 = NULL; cudaMalloc((void**) dmem2,size);
    char *hmem = (char*) malloc(size);
    for (int i = 0; i<size; i++)
    {
    hmem[i] = 's';
    }
    cudaMemcpy(hmem,dmem,size,cudaMemcpyHostToDevice);
    chTimerGetTime( &start );
    for ( int i = 0; i < cIterations; i++ ) {
    cudaMemcpy ( dmem, dmem2, size, cudaMemcpyDeviceToDevice ); 

    }   
 
    cudaThreadSynchronize();
    chTimerGetTime( &stop );
    {
        double bandwidth = (size*cIterations)/chTimerElapsedTime( &start, &stop );
	printf( "Message Size: %li Byte\n", size );
        printf( "Bandwidth: %.5g Byte/s\n", bandwidth );
	FILE *fp;
	fp = fopen("MemorySpeedD2D.txt","a");
	fprintf (fp, "%li %.5g\n", size, bandwidth );
	fclose(fp);
    }

    cudaFree(dmem);  // Free device buffer
    cudaFree(dmem2);



    return 0;
}
