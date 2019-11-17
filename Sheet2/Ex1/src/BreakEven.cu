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
__device__ long int dummyFlag;
__global__
void
WaitKernel(long int waitTime)
{


clock_t clockstart = clock();
long int timeDiff = clock() - clockstart;

//printf("timeDiff: %li \n waitTime: %li \n",timeDiff,waitTime);
while (timeDiff < waitTime)
{
//printf("(insideloop)timeDiff: %li \n waitTime: %li \n",timeDiff,waitTime);

timeDiff = clock() - clockstart;
//printf("(after new calc)timeDiff: %li \n waitTime: %li \n",timeDiff,waitTime);

if (threadIdx.x == 1000000)
{
	dummyFlag = timeDiff;
}
}
//printf("(endofkernel)timeDiff: %li \n waitTime: %li \n",timeDiff,waitTime);
}

int
main(int argc,char **argv)
{
    int device;
    cudaGetDevice(&device);

    struct cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);
    int clockrate = props.clockRate;
    printf("clockrate: %d\n",clockrate); 
    long int waitTimeCl = atoi(argv[1]);
    const int cIterations = 1e5;
    printf( "Measuring Asynchronous launch time... " ); fflush( stdout );
    //printf( "dummyFlag = %.5g\n",dummyFlag ); fflush( stdout );
    chTimerTimestamp start, stop;

    chTimerGetTime( &start );
    for ( int i = 0; i < cIterations; i++ ) {
        WaitKernel<<<1,1>>>(waitTimeCl);
    }
    cudaThreadSynchronize();
    chTimerGetTime( &stop );

    {
        double microseconds = 1e6*chTimerElapsedTime( &start, &stop );
        double usPerLaunch = microseconds / (float) cIterations;
	double AsyncTime = usPerLaunch;

        printf( "%.2f us\n", usPerLaunch );
	printf( "%li clocks per sec\n", CLOCKS_PER_SEC );
	FILE *fp;
	fp = fopen("BreakEven.txt","a");
	fprintf (fp, "%li %.5g\n", waitTimeCl, AsyncTime );
	fclose(fp);
    }

    return 0;
}
