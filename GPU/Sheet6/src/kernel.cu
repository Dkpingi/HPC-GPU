/**************************************************************************************************
 *
 *       Computer Engineering Group, Heidelberg University - GPU Computing Exercise 06
 *
 *                 Gruppe : TODO
 *
 *                   File : kernel.cu
 *
 *                Purpose : Reduction
 *
 **************************************************************************************************/
#include <cstdio>
//
// Reduction_Kernel
//
__global__ void
reduction_Kernel(int numElements, float* dataIn, float* dataOut)
{
	int elementId = blockIdx.x * blockDim.x + threadIdx.x;
	extern __shared__ float blockSum[];
	//printf("Thread %d start", threadIdx.x);
	if (elementId < numElements)
	{
		//printf("ThreadId: %d\n", threadIdx.x);
		//printf("ElementID: %d\n", elementId);
		//printf("NumElements: %d\n", numElements);
	        blockSum[threadIdx.x] = dataIn[elementId];
		//printf("thread just wrote: %d\n", threadIdx.x);
	} 
	else
	{
		//printf("%d\n", threadIdx.x);
		//printf("%d\n", elementId);
		//printf("%d\n", numElements);
		blockSum[threadIdx.x] = 0.0;
	}
	__syncthreads();
	
	for(int s = 1;s< blockDim.x; s*= 2)
	{
		//printf("4\n");
		int index = 2*s*threadIdx.x;
		if(index < blockDim.x)
		{
			//printf("5\n");
			blockSum[index] += blockSum[index + s];
		}
		__syncthreads();
	}
		
	if(threadIdx.x == 0)
		*(dataOut + blockIdx.x) = blockSum[0];
	
	//printf("6\n");
}

void reduction_Kernel_Wrapper(dim3 gridSize, dim3 blockSize, int numElements, float* dataIn, float* dataOut) {
	float* buffer = NULL;
	cudaMalloc(&buffer,gridSize.x*sizeof(float));
	//printf("7\n");
	reduction_Kernel<<< gridSize, blockSize, blockSize.x*sizeof(float) >>>(numElements, dataIn, buffer);
	unsigned int nthreads = 1;
	while(nthreads < gridSize.x)
	{
		nthreads*=2;
	}
	cudaDeviceSynchronize();
	//printf("8\n");
	reduction_Kernel<<< 1, nthreads, nthreads*sizeof(float) >>>(gridSize.x, buffer, dataOut);
	//printf("9\n");
	cudaFree(buffer);

}
