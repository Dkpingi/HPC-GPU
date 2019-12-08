/**************************************************************************************************
 *
 *       Computer Engineering Group, Heidelberg University - GPU Computing ExercthreadXe 05
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

const static int DEFAULT_MATRIX_WIDTH = 1024;
const static int DEFAULT_BLOCK_DIM = 32;

//
// Function Prototypes
//
void printHelp(char * /*programName*/);

/**
* Naive Kernel implementation without shared memory: 5.2
* C = A*B
**/
__global__ void
matMul_Kernel(int matrixSize, float* matrixA, float* matrixB, float* matrixC)
{
    // Jeder Thread berechnet einen Value der Ergebnis Matrix
    // Column
    int elementIdx = blockIdx.x * blockDim.x + threadIdx.x;
    // Row
    int elementIdy = blockIdx.y * blockDim.y + threadIdx.y;

    int elementId = elementIdy * matrixSize + elementIdx;
    // The result value the thread calculates
    float c_res = 0;
    // As usual: Check index does not exceed matrix dimensions
    if (elementIdx < matrixSize && elementIdy < matrixSize) {
        for(int i = 0; i < matrixSize; ++i){
            // Calculate the required index
            c_res += matrixA[elementIdy * matrixSize + i] * matrixB[i * matrixSize + elementIdx];
        }
        // We can write c_res to the result matrix position of the thread
        // Remember: One thread calcs one result of the matrix
        matrixC[elementIdy * matrixSize + elementIdx] = c_res;
    }
}

//
// Shared matMul_Kernel
//
__global__ void
shMatMul_Kernel(int matrixSize, float *matrixA, float *matrixB, float *matrixC)
{

    // HEIGHT = WIDTH
    int TILEWIDTH = blockDim.x;
    int TILEHEIGHT = blockDim.y;
    extern __shared__ float sh_Mem[];
    float *sh_MatrixA = &(sh_Mem[0]);
    float *sh_MatrixB = &(sh_Mem[TILEWIDTH * TILEHEIGHT]);
    // Why would we need shared mem for C?
    // float *sh_MatrixC = &(sh_Mem[2 /*TODO Calc offset*/]);

    int blockX = blockIdx.x;
    int blockY = blockIdx.y;
    int threadX = threadIdx.x;
    int threadY = threadIdx.y;
    int row = blockY * TILEHEIGHT + threadY;
    int col = blockX * TILEWIDTH + threadX;

    float Cvalue = 0.0;
    // Iterate the tiles, loop the submatrices
    for (int m = 0; m < gridDim.x; ++m)
    {
        // Load A and B from device memory to shared memory
        // Each thread loads one element of each sub-matrix
        // (row, col) = row * width + col
        // A.stride * BLOCK_SIZE * row + BLOCK_SIZE * col

        // If the matrix width is not a multiple of tile dim, we have to zero-out
        if (m * TILEWIDTH + threadX < matrixSize && row < matrixSize)
        {
            sh_MatrixA[threadY * TILEWIDTH + threadX] = matrixA[row * matrixSize + (m * TILEWIDTH + threadX)];
        }
        else
        {
            sh_MatrixA[threadY * TILEWIDTH + threadX] = 0.0;
        }
        if (m * TILEWIDTH + threadY < matrixSize && col < matrixSize)
        {
            sh_MatrixB[threadY * TILEWIDTH + threadX] = matrixB[col + (m * TILEWIDTH + threadY) * matrixSize];
        }
        else
        {
            sh_MatrixB[threadY * TILEWIDTH + threadX] = 0.0;
        }
        // Wait for operation to complete
        __syncthreads();

        // Multiply the sub matrices
        for (int e = 0; e < TILEWIDTH; ++e)
        {
            Cvalue += sh_MatrixA[threadY * TILEWIDTH + e] * sh_MatrixB[e * TILEWIDTH + threadX];
        }

        // Lets wait for all to complete before continuing
        __syncthreads();
    }
    if (row < matrixSize &&
        col < matrixSize)
    {
        // Write real result back
        matrixC[row * matrixSize + col] = Cvalue;
    }
}

// Main
int main(int argc, char *argv[])
{
    //
    // Show Help
    //
    bool showHelp = chCommandLineGetBool("h", argc, argv);
    if (!showHelp)
    {
        showHelp = chCommandLineGetBool("help", argc, argv);
    }

    if (showHelp)
    {
        printHelp(argv[0]);
        exit(0);
    }
    std::cout << "**************************************************************************\n";
    std::cout << "***" << std::endl
              << "*** Starting ..." << std::endl
              << "***" << std::endl;

    ChTimer memCpyH2DTimer, memCpyD2HTimer;
    ChTimer kernelTimer, CPUTimer;

    // More would be better, but it just takes too long
    const int ITERATIONS = 10;
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
    if (!pinnedMemory)
    {
        pinnedMemory = chCommandLineGetBool("pinned-memory", argc, argv);
    }

    float *h_matrixA = NULL;
    float *h_matrixB = NULL;
    float *h_matrixC = NULL;
    if (!pinnedMemory)
    {
        // Pageable
        h_matrixA = static_cast<float *>(malloc(
            static_cast<size_t>(matrixSize * sizeof(*h_matrixA))));
        h_matrixB = static_cast<float *>(malloc(
            static_cast<size_t>(matrixSize * sizeof(*h_matrixB))));
        h_matrixC = static_cast<float *>(calloc(
            static_cast<size_t>(matrixSize), sizeof *h_matrixC));
    }
    else
    {
        // Pinned
        cudaMallocHost(&h_matrixA, static_cast<size_t>(matrixSize * sizeof(*h_matrixA)));
        cudaMallocHost(&h_matrixB, static_cast<size_t>(matrixSize * sizeof(*h_matrixB)));
        cudaMallocHost(&h_matrixC, static_cast<size_t>(matrixSize * sizeof(*h_matrixC)));
        memset(h_matrixC, 0, matrixSize * sizeof(*h_matrixC));
    }

    //
    // Device Memory
    //
    float *d_matrixA = NULL;
    float *d_matrixB = NULL;
    float *d_matrixC = NULL;
    cudaMalloc(&d_matrixA, static_cast<size_t>(matrixSize * sizeof(*d_matrixA)));
    cudaMalloc(&d_matrixB, static_cast<size_t>(matrixSize * sizeof(*d_matrixB)));
    cudaMalloc(&d_matrixC, static_cast<size_t>(matrixSize * sizeof(*d_matrixC)));

    //
    // Check Pointers
    //
    if (h_matrixA == NULL || h_matrixB == NULL || h_matrixC == NULL ||
        d_matrixA == NULL || d_matrixB == NULL || d_matrixC == NULL)
    {
        std::cout << "\033[31m***" << std::endl
                  << "*** Error - Allocation of Memory failed!!!" << std::endl
                  << "***\033[0m" << std::endl;
        exit(-1);
    }

    //
    // Init Matrices
    //
    for (int i = 0; i < matrixSize; i++)
    {
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
    chCommandLineGet<int>(&blockSize, "t", argc, argv);
    chCommandLineGet<int>(&blockSize, "threads-per-block", argc, argv);
    blockSize = blockSize != 0 ? blockSize : DEFAULT_BLOCK_DIM;

    if (blockSize > 32)
    {
        std::cout << "\033[31m***" << std::endl
                  << "*** Error - The number of threads per block is too high" << std::endl
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
    int sharedMemSize = 2 * pow(blockSize, 2) * sizeof(float);
    kernelTimer.start();

    //
    // Launch Kernel for i ITERATIONS
    //
    for (int i = 0; i < ITERATIONS; ++i)
    {
        if (!chCommandLineGetBool("shared", argc, argv))
        {
            matMul_Kernel<<<grid_dim, block_dim>>>(matrixWidth, d_matrixA, d_matrixB, d_matrixC);
        }
        else
        {
            shMatMul_Kernel<<<grid_dim, block_dim, sharedMemSize>>>(matrixWidth, d_matrixA, d_matrixB, d_matrixC);
        }

        //
        // Synchronize
        //
        cudaDeviceSynchronize();
    }
    //
    // Check for Errors
    //
    cudaError_t cudaError = cudaGetLastError();
    if (cudaError != cudaSuccess)
    {
        std::cout << "\033[31m***" << std::endl
                  << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
                  << std::endl
                  << "***\033[0m" << std::endl;

        return -1;
    }

    kernelTimer.stop();

    if (!chCommandLineGetBool("shared", argc, argv))
    {
        std::cout << "Using matMul_Kernel\n";
    }
    else
    {
        std::cout << "Using shMatMul_Kernel\n";
    }
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
    if (!dontCheckResult)
    {
        dontCheckResult = chCommandLineGetBool("no-check", argc, argv);
    }

    float cpu_mat_time = 0;
    if (!dontCheckResult)
    {
        std::cout << "Now checking results...\n";
        float *h_matrixD = static_cast<float *>(
            calloc(static_cast<size_t>(matrixSize), sizeof(*h_matrixD)));

        CPUTimer.start();
        MatrixMulOnHostBlocked(h_matrixA, h_matrixB, h_matrixD,
                               static_cast<long>(matrixWidth), 32);
        CPUTimer.stop();
        std::cout << "CPU took " << 1e3 * CPUTimer.getTime() << "ms." << std::endl;
        cpu_mat_time = 1e3 * CPUTimer.getTime();
        bool resultOk = MatrixCompare(h_matrixC, h_matrixD,
                                      static_cast<long>(matrixWidth));

        if (!resultOk)
        {
            std::cout << "\033[31m***" << std::endl
                      << "*** Error - The two matrices are different!!!" << std::endl
                      << "***\033[0m" << std::endl;
            printOutMatrix(h_matrixA, matrixWidth);
            printOutMatrix(h_matrixB, matrixWidth);
            printOutMatrix(h_matrixC, matrixWidth);
            printOutMatrix(h_matrixD, matrixWidth);

            exit(-1);
        }
        else
        {
            std::cout << "The two matrices are the same\n";
        }

        free(h_matrixD);
    }
    else
    {
        std::cout << "Result CPU Check disabled\n";
    }

    //
    // Print Meassurement Results
    //
    std::cout << "***" << std::endl
              << "*** Results:" << std::endl
              << "*** Matrix width:" << matrixWidth << std::endl
              << "***    Matrix Size (width*width): " << matrixSize << std::endl
              << "***    Grid size: " << gridSize << std::endl
              << "***    Block size : " << blockSize << std::endl
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
              << "***    Time for Matrix Multiplication: " << 1e3 * kernelTimer.getTime() / ITERATIONS
              << " ms" << std::endl
              << "***" << std::endl;
    // Print csv style for plotting
    std::cout << "CSV\n matrixWidth, blockSize, H2D_ms, D2H_ms, matrixMult_ms_complCalc, matrixMult_ms, ITERATIONS, CPU Time\n";
    std::cout << matrixWidth << ", " << blockSize << ", " << 1e3 * memCpyH2DTimer.getTime() << ", " << 1e3 * memCpyD2HTimer.getTime() << ", " << 1e3 * kernelTimer.getTime() << ", " << 1e3 * kernelTimer.getTime() / ITERATIONS << ", " << ITERATIONS << ", " << cpu_mat_time << std::endl;
    std::cout << "CSV\n";
    std::cout << "**************************************************************************\n";
    if (chCommandLineGetBool("print-matrix", argc, argv) && matrixWidth <= 16)
    {
        printOutMatrix(h_matrixC, matrixWidth);
    }

    // Free Memory
    if (!pinnedMemory)
    {
        free(h_matrixA);
        free(h_matrixB);
        free(h_matrixC);
    }
    else
    {
        cudaFreeHost(h_matrixA);
        cudaFreeHost(h_matrixB);
        cudaFreeHost(h_matrixC);
    }
    cudaFree(d_matrixA);
    cudaFree(d_matrixB);
    cudaFree(d_matrixC);

    return 0;
}

void printHelp(char *programName)
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
