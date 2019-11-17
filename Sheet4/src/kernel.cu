/******************************************************************************
 *
 *Computer Engineering Group, Heidelberg University - GPU Computing Exercise 04
 *
 *                  Group : TBD
 *
 *                   File : kernel.cu
 *
 *                Purpose : Memory Operations Benchmark
 *
 ******************************************************************************/

#include <stdio.h>
//
// Test Kernel
//

__global__ void 
globalMem2SharedMem
(int a_size, float* g_array)
{
	int id = threadIdx.x;
	extern __shared__ float s_array[];
	int size = floorf(a_size/sizeof(float));
	for (int i = 1;i*blockDim.x <= size; i++)
	{
		//printf("upper:%d\n",id);
		s_array[id] = g_array[id];
		id  += blockDim.x;
		__syncthreads();
	}
	if (id < size)
	{
		//printf("lower:%d\n",id);
		s_array[id] = g_array[id];
	}
	if (threadIdx.x == 0)
		g_array[0] = 1.0;
}

void globalMem2SharedMem_Wrapper(dim3 gridSize, dim3 blockSize, int shmSize, float* g_array) {
	globalMem2SharedMem<<< gridSize, blockSize, shmSize >>>(shmSize, g_array);
}

__global__ void 
SharedMem2globalMem
(int a_size,float* g_array)
{
	int id = threadIdx.x;
	extern __shared__ float s_array[];
	int size = floorf(a_size/sizeof(float));
	for (int i = 1;i*blockDim.x <= size; i++)
	{
		//printf("upper:%d\n",id);
		g_array[id] = s_array[id];
		id  += blockDim.x;
		__syncthreads();
	}
	if (id < size)
	{
		//printf("lower:%d\n",id);
		g_array[id] = s_array[id];
	}
	//if (threadIdx.x == 0) //PROBABLY NOT NEEDED
		//g_array[0] = 1.0;
}
void SharedMem2globalMem_Wrapper(dim3 gridSize, dim3 blockSize, int shmSize,float* g_array) {
	SharedMem2globalMem<<< gridSize, blockSize, shmSize >>>(shmSize,g_array);
}

__global__ void 
SharedMem2Registers
//(/*TODO Parameters*/)
(int a_size, float* g_array)
{
	int id = threadIdx.x;
	float reg;  // single float should always be stored in register I think
	extern __shared__ float s_array[];
	int size = floorf(a_size/sizeof(float));
	for (int i = 1;i*blockDim.x <= size; i++)
	{
		//printf("upper:%d\n",id);
		reg = s_array[id];
		id  += blockDim.x;
		__syncthreads();
	}
	if (id < size)
	{
		//printf("lower:%d\n",id);
		reg = s_array[id];
	}
	if (threadIdx.x == 0)
		g_array[0] = 1.0;
}
void SharedMem2Registers_Wrapper(dim3 gridSize, dim3 blockSize, int shmSize, float* g_array) {
	SharedMem2Registers<<< gridSize, blockSize, shmSize >>>(shmSize,g_array);
}

__global__ void 
Registers2SharedMem
//(/*TODO Parameters*/)
(int a_size, float* g_array)
{
	int id = threadIdx.x;
	float reg = 4.0;  // single float should always be stored in register I think
	extern __shared__ float s_array[];
	int size = floorf(a_size/sizeof(float));
	for (int i = 1;i*blockDim.x <= size; i++)
	{
		//printf("upper:%d\n",id);
		s_array[id] = reg;
		id  += blockDim.x;
		__syncthreads();
	}
	if (id < size)
	{
		//printf("lower:%d\n",id);
		s_array[id] = reg;
	}
	if (threadIdx.x == 0)
		g_array[0] = 1.0;
}
void Registers2SharedMem_Wrapper(dim3 gridSize, dim3 blockSize, int shmSize, float* g_array) {
	Registers2SharedMem<<< gridSize, blockSize, shmSize >>>(shmSize, g_array);
}

__global__ void 
bankConflictsRead
//(/*TODO Parameters*/)
(int stride,long int* dClocks, float* g_array)
{
	long long start,end;
	int id = stride*threadIdx.x;
	float reg = 4.0;  // single float should always be stored in register I think
	extern __shared__ float s_array[];

	start = clock64();
	for(int i = 0;i<100;i++)
	{
		__syncthreads();
		reg = s_array[id];
	}
	end = clock64();
	long long time = end - start;
	*dClocks += time;
	if (threadIdx.x == 0)
		g_array[0] = reg;
}

void bankConflictsRead_Wrapper(dim3 gridSize, dim3 blockSize, int shmSize, int stride,long int* dClocks,float* g_array) {
	bankConflictsRead<<< gridSize, blockSize, shmSize >>>(stride,dClocks,g_array);
	cudaDeviceSynchronize();
}
