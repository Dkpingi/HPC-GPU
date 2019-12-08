/**************************************************************************************************
 *
 *       Computer Engineering Group, Heidelberg University - GPU Computing Exercise 05
 *
 *                                 Group : TODO
 *
 *                                  File : main.cu
 *
 *                               Purpose : Naive Matrix Multiplication
 *
 *************************************************************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <chCommandLine.h>
#include <chTimer.hpp>
#include <cuda_runtime.h>
#include <cstdio>
#include "mmult_cpu.h"

const static int DEFAULT_MATRIX_WIDTH  = 1024;
const static int DEFAULT_BLOCK_DIM     =   32;

//
// Function Prototypes
//
void printHelp(char * /*programName*/);

//
// matMul_Kernel
//
__global__ void
matMul_Kernel(int matrixSize, float* matrixA, float* matrixB, float* matrixC)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < matrixSize && 
        j < matrixSize) {
        *(matrixC + i*matrixSize + j) = 0;
	for(int k=0;k<matrixSize;k++)
	{
		*(matrixC + i*matrixSize + j) += *(matrixA + i*matrixSize + k)*(*(matrixB + k*matrixSize + j));
	}
    }
}

//
// Shared matMul_Kernel
//
__global__ void
shMatMul_Kernel(int matrixSize, float* matrixA, float* matrixB, float* matrixC)
{
    extern __shared__ float sh_Mem[];
    float* sh_MatrixA = sh_Mem;
    float* sh_MatrixB = sh_MatrixA + blockDim.x*blockDim.y;
    //float* sh_MatrixC = sh_MatrixB + sizeof(float)*blockDim.x*blockDim.y;

    int ig = blockIdx.x * blockDim.x + threadIdx.x;
    int jg = blockIdx.y * blockDim.y + threadIdx.y;
    
    int is = threadIdx.x;
    int js = threadIdx.y;
    //printf("Kernel: %d, Before If\n",ig);
  //  if (ig < matrixSize && 
  //      jg < matrixSize) 
    {
        float c  = 0.0f;
	for(int k = 0; k < gridDim.x;k++)
	{
	    //printf("Kernel: %d, before global load\n",ig);
	    //printf("Kernel: %d, %d, %d, %d \n",ig,ig*matrixSize + k*blockDim.x + js,jg + matrixSize*(k*blockDim.x + is),is*matrixSize + js);
	    if (k*blockDim.x + js < matrixSize)
	    {
                *(sh_MatrixA + is*blockDim.x + js) = *(matrixA + ig*matrixSize + k*blockDim.x + js);
	    //printf("Kernel: %d, after global load\n",ig);
	    }
	    else
	    {
		*(sh_MatrixA + is*blockDim.x + js) = 0.0f;
	    }
	    if (k*blockDim.x + is < matrixSize)
	    {
                *(sh_MatrixB + is*blockDim.x + js) = *(matrixB + jg + matrixSize*(k*blockDim.x + is));
	    //printf("Kernel: %d, after global load\n",ig);
	    }
	    else
	    {
		*(sh_MatrixB + is*blockDim.x + js) = 0.0f;
	    }
	    __syncthreads();
            
	
	    for(int l = 0; l< blockDim.x;l++)
	    {
			c += *(sh_MatrixA + is*blockDim.x + l)*(*(sh_MatrixB + l*blockDim.y + js));
            }
	    __syncthreads();
	}
        if (ig < matrixSize && 
            jg < matrixSize) 
	{
	 //printf("Kernel: %d, Before C\n",ig);
	    *(matrixC + ig*matrixSize + jg) = c;
	}
    }
    //printf("Kernel: %d, done\n",ig);
    
}


// might need later
/**
__global__ void
shMatMul_Kernel(int matrixSize, float* matrixA, float* matrixB, float* matrixC)
{
    extern __shared__ float sh_Mem[];
    float* sh_MatrixA = sh_Mem;
    float* sh_MatrixB = sh_MatrixA + blockDim.x*blockDim.y;
    //float* sh_MatrixC = sh_MatrixB + sizeof(float)*blockDim.x*blockDim.y;

    int ig = blockIdx.x * blockDim.x + threadIdx.x;
    int jg = blockIdx.y * blockDim.y + threadIdx.y;
    
    int is = threadIdx.x;
    int js = threadIdx.y;
    //printf("Kernel: %d, Before If\n",ig);
  //  if (ig < matrixSize && 
  //      jg < matrixSize) 
    {
        float c  = 0.0f;
	for(int k = 0; k < gridDim.x;k++)
	{
	    //printf("Kernel: %d, before global load\n",ig);
	    //printf("Kernel: %d, %d, %d, %d \n",ig,ig*matrixSize + k*blockDim.x + js,jg + matrixSize*(k*blockDim.x + is),is*matrixSize + js);
	    if ((k*blockDim.x + js < matrixSize) && (k*blockDim.x + is < matrixSize))
	    {
            *(sh_MatrixA + is*blockDim.x + js) = *(matrixA + ig*matrixSize + k*blockDim.x + js);
            *(sh_MatrixB + is*blockDim.x + js) = *(matrixB + jg + matrixSize*(k*blockDim.x + is));
	    //printf("Kernel: %d, after global load\n",ig);
	    }
	    else
	    {
		*(sh_MatrixA + is*blockDim.x + js) = 0.0f;
		*(sh_MatrixB + is*blockDim.x + js) = 0.0f;
	    }
	    __syncthreads();
            
	
	    for(int l = 0; l< blockDim.x;l++)
	    {
			c += *(sh_MatrixA + is*blockDim.x + l)*(*(sh_MatrixB + l*blockDim.y + js));
            }
	    __syncthreads();
	}
        if (ig < matrixSize && 
            jg < matrixSize) 
	{
	 //printf("Kernel: %d, Before C\n",ig);
	*(matrixC + ig*matrixSize + jg) = c;
	}
    }
    //printf("Kernel: %d, done\n",ig);
    
}
**/
//
// Main
//
int
main(int argc, char * argv[])
{
    //
    // Show Help
    //
    bool showHelp = chCommandLineGetBool("h", argc, argv);
    if (!showHelp) {
        showHelp = chCommandLineGetBool("help", argc, argv);
    }

    if (showHelp) {
        printHelp(argv[0]);
        exit(0);
    }

    std::cout << "***" << std::endl
              << "*** Starting ..." << std::endl
              << "***" << std::endl;

    ChTimer memCpyH2DTimer, memCpyD2HTimer;
    ChTimer kernelTimer, CPUTimer;

    //
    // Allocate Memory
    //
    int matrixWidth = 0;
    chCommandLineGet<int>(&matrixWidth, "s", argc, argv);
    chCommandLineGet<int>(&matrixWidth, "size", argc, argv);
    matrixWidth = matrixWidth != 0 ? matrixWidth : DEFAULT_MATRIX_WIDTH;

    int matrixSize = matrixWidth * matrixWidth;

    //
    // Host Memory
    //
    bool pinnedMemory = chCommandLineGetBool("p", argc, argv);
    if (!pinnedMemory) {
        pinnedMemory = chCommandLineGetBool("pinned-memory",argc,argv);
    }

    float* h_matrixA = NULL;
    float* h_matrixB = NULL;
    float* h_matrixC = NULL;
    if (!pinnedMemory) {
        // Pageable
        h_matrixA = static_cast<float*>(malloc(
                        static_cast<size_t>(matrixSize * sizeof(*h_matrixA))));
        h_matrixB = static_cast<float*>(malloc(
                        static_cast<size_t>(matrixSize * sizeof(*h_matrixB))));
        h_matrixC = static_cast<float*>(calloc(
                        static_cast<size_t>(matrixSize), sizeof *h_matrixC));

    } else {
        // Pinned
        cudaMallocHost(&h_matrixA, static_cast<size_t>(matrixSize * sizeof(*h_matrixA)));
        cudaMallocHost(&h_matrixB, static_cast<size_t>(matrixSize * sizeof(*h_matrixB)));
        cudaMallocHost(&h_matrixC, static_cast<size_t>(matrixSize * sizeof(*h_matrixC)));
        memset ( h_matrixC, 0, matrixSize * sizeof(*h_matrixC) );
    }

    //
    // Device Memory
    //
    float* d_matrixA = NULL;
    float* d_matrixB = NULL;
    float* d_matrixC = NULL;
    cudaMalloc(&d_matrixA, static_cast<size_t>(matrixSize * sizeof(*d_matrixA)));
    cudaMalloc(&d_matrixB, static_cast<size_t>(matrixSize * sizeof(*d_matrixB)));
    cudaMalloc(&d_matrixC, static_cast<size_t>(matrixSize * sizeof(*d_matrixC)));

    //
    // Check Pointers
    //
    if (h_matrixA == NULL || h_matrixB == NULL || h_matrixC == NULL ||
        d_matrixA == NULL || d_matrixB == NULL || d_matrixC == NULL )
    {
        std::cout << "\033[31m***" << std::endl
                  << "*** Error - Allocation of Memory failed!!!" << std::endl
                  << "***\033[0m" << std::endl;
        exit(-1);
    }

    //
    // Init Matrices
    //
    for (int i = 0; i < matrixSize; i++) {
        int x = i % matrixWidth;
        int y = i / matrixWidth;
        h_matrixA[i] = static_cast<float>(x * y);
        h_matrixB[i] = static_cast<float>(x + y);
    }

    //
    // Copy Data to the Device
    //
    memCpyH2DTimer.start();

    cudaMemcpy(d_matrixA, h_matrixA, static_cast<size_t>(matrixSize * sizeof(*d_matrixA)), 
            cudaMemcpyHostToDevice);
    cudaMemcpy(d_matrixB, h_matrixB, static_cast<size_t>(matrixSize * sizeof(*d_matrixB)), 
            cudaMemcpyHostToDevice);

    memCpyH2DTimer.stop();

    //
    // Get Kernel Launch Parameters
    //
    int blockSize = 0,
        gridSize = 0;

    // Block Dimension / Threads per Block
    chCommandLineGet<int>(&blockSize,"t", argc, argv);
    chCommandLineGet<int>(&blockSize,"threads-per-block", argc, argv);
    blockSize = blockSize != 0 ? blockSize : DEFAULT_BLOCK_DIM;

    if (blockSize > 32) {
        std::cout << "\033[31m***" << std::endl
                  << "*** Error - The number of threads per block is too big" << std::endl
                  << "***\033[0m" << std::endl;
        exit(-1);
    }

    gridSize = ceil(static_cast<float>(matrixWidth) / static_cast<float>(blockSize));

    dim3 grid_dim = dim3(gridSize, gridSize, 1);
    dim3 block_dim = dim3(blockSize, blockSize, 1);

    std::cout << "***" << std::endl
              << "*** Grid Dim:  " << grid_dim.x << "x" << grid_dim.y << "x" << grid_dim.z 
                      << std::endl
              << "*** Block Dim: " << block_dim.x << "x" << block_dim.y << "x" << block_dim.z 
                      << std::endl
              << "***" << std::endl;

    // TODO Calc shared mem size
    int sharedMemSize = 2*pow(blockSize,2)*sizeof(float);
    kernelTimer.start();

    //
    // Launch Kernel
    //
    if (!chCommandLineGetBool("shared", argc, argv)) {
        matMul_Kernel<<<grid_dim, block_dim>>>(matrixWidth, d_matrixA, d_matrixB, d_matrixC);
    } else {
        shMatMul_Kernel<<<grid_dim, block_dim, sharedMemSize>>>(matrixWidth, d_matrixA, d_matrixB, d_matrixC);
    }

    //
    // Synchronize
    //
    cudaDeviceSynchronize();

    //
    // Check for Errors
    //
    cudaError_t cudaError = cudaGetLastError();
    if ( cudaError != cudaSuccess ) {
        std::cout << "\033[31m***" << std::endl
                  << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
                    << std::endl
                  << "***\033[0m" << std::endl;

        return -1;
    }

    kernelTimer.stop();

    //
    // Copy Back Data
    //
    memCpyD2HTimer.start();

    cudaMemcpy(h_matrixC, d_matrixC, static_cast<size_t>(matrixSize * sizeof(*d_matrixC)), 
            cudaMemcpyDeviceToHost);

    memCpyD2HTimer.stop();

    //
    // Check Result
    //
    bool dontCheckResult = chCommandLineGetBool("c", argc, argv);
    if (!dontCheckResult) {
        dontCheckResult = chCommandLineGetBool("no-check", argc, argv);
    }

    if (!dontCheckResult) {
        float* h_matrixD = static_cast<float*>(
                calloc(static_cast<size_t>(matrixSize), sizeof(*h_matrixD)));

	CPUTimer.start();
        MatrixMulOnHostBlocked(h_matrixA, h_matrixB, h_matrixD, 
                static_cast<long>(matrixWidth), 32);
	CPUTimer.stop();
	std::cout<< "CPU took " << 1e3*CPUTimer.getTime() << "ms." << std::endl;
        bool resultOk = MatrixCompare(h_matrixC, h_matrixD, 
                static_cast<long>(matrixWidth));

        if (!resultOk) {
            std::cout << "\033[31m***" << std::endl
                      << "*** Error - The two matrices are different!!!" << std::endl
                      << "***\033[0m" << std::endl;
	    printOutMatrix(h_matrixA, matrixWidth);
	    printOutMatrix(h_matrixB, matrixWidth);
	    printOutMatrix(h_matrixC, matrixWidth);
	    printOutMatrix(h_matrixD, matrixWidth);

            exit(-1);
        }

        free(h_matrixD);
    }

    //
    // Print Meassurement Results
    //
    std::cout << "***" << std::endl
              << "*** Results:" << std::endl
              << "***    Matrix Size: " << matrixSize << std::endl
              << "***    Time to Copy to Device: " << 1e3 * memCpyH2DTimer.getTime()
                << " ms" << std::endl
              << "***    Copy Bandwidth: " 
                << 1e-9 * memCpyH2DTimer.getBandwidth(2 * matrixSize * sizeof(*h_matrixA))
                << " GB/s" << std::endl
              << "***    Time to Copy from Device: " << 1e3 * memCpyD2HTimer.getTime()
                << " ms" << std::endl
              << "***    Copy Bandwidth: " 
                << 1e-9 * memCpyD2HTimer.getBandwidth(matrixSize * sizeof(*h_matrixA))
                << " GB/s" << std::endl
              << "***    Time for Matrix Multiplication: " << 1e3 * kernelTimer.getTime()
                  << " ms" << std::endl
              << "***" << std::endl;
    std::fstream fin,fout;
    fout.open("matMul.txt", std::ios::in | std::ios::out | std::ios::app);
    fin.open("matMul.txt", std::ios::in);
    if (!fout)
    {
    	std::cerr << "file open failed:" << std::endl; 
    }
    if(fin.peek() == std::fstream::traits_type::eof())
    {
	fout << "BoolShMemory" << " " 
               << "MatrixWidth" << " " << "BlockWidth" << " " 
               << "GridWidth" << " " << "H2DTime" << " " 
               << "CalcTime" << " " << "D2HTime" << " "
	       << "\n";
    } 
    fout << chCommandLineGetBool("shared", argc, argv) << " "
	   << matrixSize << " " << blockSize << " " 
           << gridSize << " " << 1e3*memCpyH2DTimer.getTime() << " " 
	   << 1e3*kernelTimer.getTime() << " "  << 1e3*memCpyD2HTimer.getTime()
           << "\n";
    fin.close();
    fout.close();
    if (chCommandLineGetBool("print-matrix", argc, argv) 
       && matrixWidth <= 16) {
        printOutMatrix(h_matrixC, matrixWidth);
    }

    // Free Memory
    if (!pinnedMemory) {
        free(h_matrixA);
        free(h_matrixB);
        free(h_matrixC);
    } else {
        cudaFreeHost(h_matrixA);
        cudaFreeHost(h_matrixB);
        cudaFreeHost(h_matrixC);
    }
    cudaFree(d_matrixA);
    cudaFree(d_matrixB);
    cudaFree(d_matrixC);

    return 0;
}

void
printHelp(char * programName)
{
    std::cout << "Help:" << std::endl
              << "  Usage: " << std::endl
              << "  " << programName << " [-p] [-s <matrix_size>] [-t <threads_per_block>]" 
                << std::endl
              << "                 [-g <blocks_per_grid] [-c] [--print-matrix]" 
                << std::endl
              << "" << std::endl
              << "  -p|--pinned-memory" << std::endl
              << "  Use pinned Memory instead of pageable memory" << std::endl
              << "" << std::endl
              << "  -s <matrix_size>|--size <matix_size>" << std::endl
              << "  The width of the Matrix" << std::endl
              << "" << std::endl
              << "  -t <threads_per_block>|--threads-per-block <threads_per_block>" 
                << std::endl
              << "  The number of threads per block" << std::endl
              << "" << std::endl
              << "  -c|--no-checking" << std::endl
              << "  Do not check the result of the matrix multiplication" << std::endl
              << "" << std::endl
              << "  --print-matrix" << std::endl
              << "  Print the output matrix (only recommended for small matrices)" << std::endl
              << std::endl;
}
