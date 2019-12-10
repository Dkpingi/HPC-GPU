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
	int tid = threadIdx.x;
	if (elementId < numElements)
	{
	        blockSum[tid] = dataIn[elementId];
	} 
	else
	{
		blockSum[tid] = 0.0;
	}
	__syncthreads();
	
	for(int s = 1;s< blockDim.x; s*= 2)
	{
		if(tid%(2*s) == 0)
		{
			blockSum[tid] += blockSum[tid + s];
		}
		__syncthreads();
	}
		
	if(tid == 0)
		*(dataOut + blockIdx.x) = blockSum[0];

}

void reduction_Kernel_Wrapper(dim3 gridSize, dim3 blockSize, int numElements, float* dataIn, float* dataOut) {
	float* buffer = NULL;
	cudaMalloc(&buffer,gridSize.x*sizeof(float));
	reduction_Kernel<<< gridSize, blockSize, blockSize.x*sizeof(float) >>>(numElements, dataIn, buffer);
	cudaDeviceSynchronize();

	unsigned int nthreads = 1;
	while(nthreads < gridSize.x)
	{
		nthreads*=2;
	}
	reduction_Kernel<<< 1, nthreads, nthreads*sizeof(float) >>>(gridSize.x, buffer, dataOut);


}


__global__ void
opt_reduction_Kernel(int numElements, float* dataIn, float* dataOut)
{
	int elementId = blockIdx.x * 2 * blockDim.x + threadIdx.x;
	extern __shared__ float blockSum[];
	int tid = threadIdx.x;
	if (elementId < numElements)
	{
		blockSum[tid] = dataIn[elementId] + dataIn[elementId + blockDim.x];
	}
	else
	{
		blockSum[tid] = 0.0;
	}

	__syncthreads();	
	for(unsigned int s = blockDim.x/2;s>32; s >>= 1)
	{
		if(tid < s)
		{
			blockSum[tid] += blockSum[tid + s];
		}
		__syncthreads();
	}
 	if( tid < 32 && blockDim.x >= 64) blockSum[tid] += blockSum[tid + 32];  
 	if ( tid < 16 && blockDim.x >= 32) blockSum[tid] += blockSum[tid + 16];  
 	if ( tid <  8 && blockDim.x >= 16) blockSum[tid] += blockSum[tid + 8];   
 	if ( tid <  4 && blockDim.x >=  8) blockSum[tid] += blockSum[tid + 4];   
 	if ( tid <  2 && blockDim.x >=  4) blockSum[tid] += blockSum[tid + 2];   
 	if ( tid <  1 && blockDim.x >=  2) blockSum[tid] += blockSum[tid + 1];
	
			
	if(tid == 0)
		*(dataOut + blockIdx.x) = blockSum[0];
	
}

void opt_reduction_Kernel_Wrapper(dim3 gridSize, dim3 blockSize, int numElements, float* dataIn, float* dataOut) {
	float* buffer = NULL;
	cudaMalloc(&buffer,gridSize.x*sizeof(float));
	blockSize.x /= 2;
	opt_reduction_Kernel<<< gridSize, blockSize, blockSize.x*sizeof(float) >>>(numElements, dataIn, buffer);
	cudaDeviceSynchronize();
	unsigned int nthreads = 1;
	while(nthreads < gridSize.x/2)
	{
		nthreads*=2;
	}
	opt_reduction_Kernel<<< 1, nthreads, nthreads*sizeof(float) >>>(gridSize.x, buffer, dataOut);

}
