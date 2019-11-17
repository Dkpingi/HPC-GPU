/*************************************************************************************************
 *
 *        Computer Engineering Group, Heidelberg University - GPU Computing Exercise 03
 *
 *                           Group : TBD
 *
 *                            File : main.cu
 *
 *                         Purpose : Memory Operations Benchmark
 *
 *************************************************************************************************/

//
// Kernels
//

__global__ void 
globalMemCoalescedKernel(int *src_array,int *copy_array, int size)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int m = floorf(id/size);
    for (int i = 1; i < m + 1; i++)
	{
		if (i*id < size)
		{
		copy_array[i*id] = src_array[i*id];
		}
		syncthreads();
	}
    
}

void 
globalMemCoalescedKernel_Wrapper(dim3 gridDim, dim3 blockDim,int *src_array,int *copy_array, int size) {
	globalMemCoalescedKernel<<< gridDim, blockDim, 0 /*Shared Memory Size*/ >>>(src_array,copy_array,size);
}

__global__ void 
globalMemStrideKernel(/*TODO Parameters*/)
{
    /*TODO Kernel Code*/
}

void 
globalMemStrideKernel_Wrapper(dim3 gridDim, dim3 blockDim /*TODO Parameters*/) {
	globalMemStrideKernel<<< gridDim, blockDim, 0 /*Shared Memory Size*/ >>>( /*TODO Parameters*/);
}

__global__ void 
globalMemOffsetKernel(/*TODO Parameters*/)
{
    /*TODO Kernel Code*/
}

void 
globalMemOffsetKernel_Wrapper(dim3 gridDim, dim3 blockDim /*TODO Parameters*/) {
	globalMemOffsetKernel<<< gridDim, blockDim, 0 /*Shared Memory Size*/ >>>( /*TODO Parameters*/);
}

