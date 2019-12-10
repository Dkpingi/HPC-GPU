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
// Simple reduction kernel in global memory
// Unoptimized version
// Works, same as CPU
__global__ void
reduction_Kernel(int numElements, float *d_out, float *d_in)
{
	int myId = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;
	if (myId < numElements)
	{
		// do reduction in global mem
		// Go down from 1024 ... 512 ... 256 ... 128 ....... 1
		for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
		{
			if (tid < s)
			{
				d_in[myId] += d_in[myId + s];
			}
			__syncthreads();
		}

		// At the end, only thread 0 is out
		if (tid == 0)
		{
			d_out[blockIdx.x] = d_in[myId];
		}
	}
}

// Simple reduction wrapper
void reduction_Kernel_Wrapper(dim3 gridSize, dim3 blockSize, int numElements, float *dataIn, float *dataOut)
{
	float *rs; // tmp
	cudaMalloc(&rs, gridSize.x * sizeof(float));

	reduction_Kernel<<<gridSize, blockSize>>>(numElements, rs, dataIn);
	cudaDeviceSynchronize();
	// Only start 1 block
	//reduction_Kernel<<<gridSize, blockSize>>>(numElements, dataOut, rs);
	reduction_Kernel<<<1, gridSize>>>(numElements, dataOut, rs);

	cudaFree(rs);
}

// Optimized reduction kernel
// We used technique:
// SEQUENTIAL ADDRESSING NONDIVERGENT
__global__ void
opt_reduction_Kernel(int numElements, float *dataIn, float *dataOut)
{
	// sdata is allocated in the kernel call: 3rd arg to <<<b, t, shmem>>>
	extern __shared__ float sdata[];

	int myId = threadIdx.x + blockDim.x * blockIdx.x;
	int tid = threadIdx.x;
	if (myId < numElements)
	{
		// load shared mem from global mem
		sdata[tid] = dataIn[myId];
		__syncthreads(); // make sure entire block is loaded!

		// do reduction in shared mem
		// Basically same as normal version
		for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
		{
			if (tid < s)
			{
				sdata[tid] += sdata[tid + s];
			}
			__syncthreads(); // make sure all adds at one stage are done!
		}

		// only thread 0 writes result for this block back to global mem
		if (tid == 0)
		{
			dataOut[blockIdx.x] = sdata[0];
		}
	}
}

void opt_reduction_Kernel_Wrapper(dim3 gridSize, dim3 blockSize, int numElements, float *dataIn, float *dataOut)
{
	float *rs; // tmp
	cudaMalloc(&rs, gridSize.x * sizeof(float));
	// Calculate shared memory size
	opt_reduction_Kernel<<<gridSize, blockSize, blockSize.x * sizeof(float)>>>(numElements, dataIn, rs);
	cudaDeviceSynchronize();
	// Only start 1 block (exc)
	opt_reduction_Kernel<<<1, gridSize, gridSize.x * sizeof(float)>>>(numElements, rs, dataOut);

	cudaFree(rs);

}
