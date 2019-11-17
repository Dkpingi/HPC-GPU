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

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chCommandLine.h>
#include <chTimer.hpp>
#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>

const static int DEFAULT_MEM_SIZE       = 10*1024*1024; // 10 MB
const static int DEFAULT_NUM_ITERATIONS =         1000;
const static int DEFAULT_BLOCK_DIM      =          128;
const static int DEFAULT_GRID_DIM       =           16;

//
// Function Prototypes
//
void printHelp(char *);

//
// Kernel Wrappers
//


extern void globalMemCoalescedKernel_Wrapper(dim3 gridDim, dim3 blockDim,int *h_array,int *d_array, int size);
extern void globalMemStrideKernel_Wrapper(dim3 gridDim, dim3 blockDim /*TODO Parameters*/);
extern void globalMemOffsetKernel_Wrapper(dim3 gridDim, dim3 blockDim /*TODO Parameters*/);

//
// Main
//
int
main ( int argc, char * argv[] )
{
    // Show Help
    bool optShowHelp = chCommandLineGetBool("h", argc, argv);
    if ( !optShowHelp )
        optShowHelp = chCommandLineGetBool("help", argc, argv);

    if ( optShowHelp ) {
        printHelp ( argv[0] );
        exit (0);
    }

    std::cout << "***" << std::endl
              << "*** Starting ..." << std::endl
              << "***" << std::endl;

    ChTimer memCpyH2DTimer, memCpyD2HTimer, memCpyD2DTimer;
    ChTimer kernelTimer;

    //
    // Get kernel launch parameters and configuration
    //
    int optNumIterations = 0,
        optBlockSize = 0,
        optGridSize = 0;

    // Number of Iterations
    chCommandLineGet<int> ( &optNumIterations,"i", argc, argv );
    chCommandLineGet<int> ( &optNumIterations,"iterations", argc, argv );
    optNumIterations = ( optNumIterations != 0 ) ? optNumIterations : DEFAULT_NUM_ITERATIONS;

    // Block Dimension / Threads per Block
    chCommandLineGet <int> ( &optBlockSize,"t", argc, argv );
    chCommandLineGet <int> ( &optBlockSize,"threads-per-block", argc, argv );
    optBlockSize = optBlockSize != 0 ? optBlockSize : DEFAULT_BLOCK_DIM;

    if ( optBlockSize > 1024 ) {
        std::cout << "\033[31m***" << std::cout
                  << "*** Error - The number of threads per block is too big" 
                    << std::endl
                  << "***\033[0m" << std::endl;

        exit(-1);
    }

    // Grid Dimension
    chCommandLineGet <int> ( &optGridSize,"g", argc, argv );
    chCommandLineGet <int> ( &optGridSize,"grid-dim", argc, argv );
    optGridSize = optGridSize != 0 ? optGridSize : DEFAULT_GRID_DIM;

    // Sync after each Kernel Launch
    bool optSynchronizeKernelLaunch = chCommandLineGetBool ( "y", argc, argv );
    if (!optSynchronizeKernelLaunch)
        optSynchronizeKernelLaunch = chCommandLineGetBool ( "synchronize-kernel",   argc, argv );

    dim3 grid_dim = dim3 ( optGridSize );
    dim3 block_dim = dim3 ( optBlockSize );


    int optStride = 1; //default stride for global-stride test
    chCommandLineGet <int> ( &optStride,"stride", argc, argv );

    int optOffset = 0; //default offset for global-stride test
    chCommandLineGet <int> ( &optOffset,"offset", argc, argv );

    // Allocate Memory (take optStride resp. optOffset into account)
    int optMemorySize = 0;

    if ( chCommandLineGetBool("global-stride", argc, argv) ) {
        // determine memory size from kernel launch parameters and stride
        optMemorySize = optGridSize * optBlockSize * optStride * sizeof optMemorySize;
        std::cout << "*** Ignoring size parameter; using kernel launch parameters and stride, size=" << optMemorySize << std::endl;
    } else if ( chCommandLineGetBool("global-offset", argc, argv) ) {
        // determine memory size from kernel launch parameters and stride
        optMemorySize = optGridSize * optBlockSize + optOffset * sizeof optMemorySize;
        std::cout << "*** Ignoring size parameter; using kernel launch parameters and offset, size=" << optMemorySize << std::endl; 
    } else {
        // determine memory size from parameters
        chCommandLineGet <int> ( &optMemorySize, "s", argc, argv );
        chCommandLineGet <int> ( &optMemorySize, "size", argc, argv );
        optMemorySize = optMemorySize != 0 ? optMemorySize : DEFAULT_MEM_SIZE;
    }

    //
    // Host Memory
    //
    int* h_memoryA = NULL;
    int* h_memoryB = NULL;
    bool optUsePinnedMemory = chCommandLineGetBool ( "p", argc, argv );
    if ( !optUsePinnedMemory )
        optUsePinnedMemory = chCommandLineGetBool ( "pinned-memory",argc,argv );

    if ( !optUsePinnedMemory ) { // Pageable
        std::cout << "***" << " Using pageable memory" << std::endl;
        h_memoryA = static_cast <int*> ( malloc ( static_cast <size_t> ( optMemorySize ) ) );
        h_memoryB = static_cast <int*> ( malloc ( static_cast <size_t> ( optMemorySize ) ) );
    } else { // Pinned
        std::cout << "***" << " Using pinned memory" << std::endl;
        // Allocation of pinned host memory
        cudaMallocHost((void**) &h_memoryA,optMemorySize);
	cudaMallocHost((void**) &h_memoryB,optMemorySize);
    }

    //
    // Device Memory
    //
    int* d_memoryA = NULL;
    int* d_memoryB = NULL;
    // Allocation of device memory
    //cudaMalloc((void**) &d_memoryA,optMemorySize);
cudaError_t rc = cudaMalloc((void **) &d_memoryA, optMemorySize);

if (rc != cudaSuccess)
    printf("Could not allocate memory: %d", rc);

    cudaMalloc((void**) &d_memoryB,optMemorySize);
    printf("%d %d %d %d",!h_memoryA,!h_memoryB,!d_memoryA,!d_memoryB);

    if ( !h_memoryA || !h_memoryB || !d_memoryA || !d_memoryB ) {
        std::cout << "\033[31m***" << std::endl
                  << "*** Error - Memory allocation failed" << std::endl
                  << "***\033[0m" << std::endl;

        exit (-1);
    }

    //
    // Copy
    //
    int optMemCpyIterations = 0;
    chCommandLineGet <int> ( &optMemCpyIterations, "im", argc, argv );
    chCommandLineGet <int> ( &optMemCpyIterations, "memory-copy-iterations", argc, argv );
    optMemCpyIterations = optMemCpyIterations != 0 ? optMemCpyIterations : 1;

    // Host To Device
    memCpyH2DTimer.start ();
    for ( int i = 0; i < optMemCpyIterations; i ++ ) {
        // H2D copy
    }
    memCpyH2DTimer.stop ();

    // Device To Device
    memCpyD2DTimer.start ();
    for ( int i = 0; i < optMemCpyIterations; i ++ ) {
        // D2D copy
    }
    memCpyD2DTimer.stop ();

    // Device To Host
    memCpyD2HTimer.start ();
    for ( int i = 0; i < optMemCpyIterations; i ++ ) {
        // D2H copy
    }
    memCpyD2HTimer.stop ();

    //
    // Check for Errors
    //
    cudaError_t cudaError = cudaGetLastError ();
    if ( cudaError != cudaSuccess ) {
        std::cout << "\033[31m***" << std::endl
                  << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
                  << std::endl
                  << "***\033[0m" << std::endl;

        return -1;
    }

    //
    // Global Memory Tests
    //
    kernelTimer.start();
    for ( int i = 0; i < optNumIterations; i++ ) {
        //
        // Launch Kernel
        //
        if ( chCommandLineGetBool ( "global-coalesced", argc, argv ) ) {

            globalMemCoalescedKernel_Wrapper(grid_dim, block_dim,d_memoryA,d_memoryB,optMemorySize);
        } else if ( chCommandLineGetBool ( "global-stride", argc, argv ) ) {
            globalMemStrideKernel_Wrapper(grid_dim, block_dim /*TODO Parameters*/);
        } else if ( chCommandLineGetBool ( "global-offset", argc, argv ) ) {
            globalMemOffsetKernel_Wrapper(grid_dim, block_dim /*TODO Parameters*/);
        } else {
            break;
        }

        if ( optSynchronizeKernelLaunch ) { // Synchronize after each kernel launch
            cudaDeviceSynchronize ();

            //
            // Check for Errors
            //
            cudaError_t cudaError = cudaGetLastError ();
            if ( cudaError != cudaSuccess ) {
                std::cout << "\033[31m***" << std::endl
                          << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
                              << std::endl
                          << "***\033[0m" << std::endl;

                return -1;
            }
        }
    }

    // Mandatory synchronize after all kernel launches
    cudaDeviceSynchronize();
    kernelTimer.stop();

    //
    // Check for Errors
    //
    cudaError = cudaGetLastError();
    if ( cudaError != cudaSuccess ) {
        std::cout << "\033[31m***" << std::endl
                  << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
                  << std::endl
                  << "***\033[0m" << std::endl;

        return -1;
    }

    // Print Measurement Results
    std::cout   << "Results for cudaMemcpy:" << std::endl
                << "Size: " << std::setw(10) << optMemorySize << "B";
                //<< "***     Time to Copy (H2D): " << 1e6 * memCpyH2DTimer.getTime() << " µs" << std::endl
    std::cout.precision(2);
    std::cout   << ", H2D: " << std::fixed << std::setw(6)
                << 1e-9 * memCpyH2DTimer.getBandwidth ( optMemorySize, optMemCpyIterations ) << " GB/s"
                //<< "***     Time to Copy (D2H): " << 1e6 * memCpyD2HTimer.getTime() << " µs" << std::endl
                << ", D2H: " << std::fixed << std::setw(6) 
                << 1e-9 * memCpyD2HTimer.getBandwidth ( optMemorySize, optMemCpyIterations ) << " GB/s" 
                //<< "***     Time to Copy (D2D): " << 1e6 * memCpyD2DTimer.getTime() << " µs" << std::endl
                << ", D2D: " << std::fixed << std::setw(6)
                << 1e-9 * memCpyD2DTimer.getBandwidth ( optMemorySize, optMemCpyIterations ) << " GB/s" 
                //<< "***     Kernel (Start-Up) Time: " 
                //<< 1e6 * kernelTimer.getTime(optNumIterations) 
                //<< " µs" << std::endl
                << std::endl;

    if ( chCommandLineGetBool ( "global-coalesced", argc, argv ) ) {
        std::cout << "Coalesced copy of global memory, size=" << std::setw(10) << optMemorySize << ", gDim=" << std::setw(5) << grid_dim.x << ", bDim=" << std::setw(5) << block_dim.x;
        //std::cout << ", time=" << kernelTimer.getTime(optNumIterations) << 
        std::cout.precision ( 2 );
        std::cout << ", bw=" << std::fixed << std::setw(6) << optMemorySize / kernelTimer.getTime(optNumIterations) / (1E09) << "GB/s" << std::endl;
    }

    if ( chCommandLineGetBool ( "global-stride", argc, argv ) ) {
        std::cout << "Strided(" << std::setw(3) << optStride << ") copy of global memory, size=" << std::setw(10) << optMemorySize << ", gDim=" << std::setw(5) << grid_dim.x << ", bDim=" << std::setw(5) << block_dim.x;
        //std::cout << ", time=" << kernelTimer.getTime(optNumIterations) << 
        std::cout.precision ( 2 );
        std::cout << ", bw=" << std::fixed << std::setw ( 6 ) << ( grid_dim.x * block_dim.x ) * sizeof (int) / kernelTimer.getTime (optNumIterations ) / ( 1E09 ) << "GB/s" << std::endl;
    }

    if ( chCommandLineGetBool ( "global-offset", argc, argv ) ) {
        std::cout << "Offset(" << std::setw(3) << optOffset << ") copy of global memory, size=" << std::setw(10) << optMemorySize << ", gDim=" << std::setw(5) << grid_dim.x << ", bDim=" << std::setw(5) << block_dim.x;
        //std::cout << ", time=" << kernelTimer.getTime(optNumIterations) << 
        std::cout.precision ( 2 );
        std::cout << ", bw=" << std::fixed << std::setw ( 6 ) << ( grid_dim.x * block_dim.x ) * sizeof (int) / kernelTimer.getTime ( optNumIterations ) / ( 1E09 ) << "GB/s" << std::endl;
    }

    return 0;
}

void
printHelp(char * programName)
{
    std::cout
        << "Usage: " << std::endl
        << "  " << programName << " [-p] [-s <memory_size>] [-i <num_iterations>]" << std::endl
        << "                [-t <threads_per_block>] [-g <blocks_per_grid]" << std::endl
        << "                [-m <memory-copy-iterations>] [-y] [-stride <stride>] [-offset <offset>]" << std::endl
        << "  --global-{coalesced|stride|offset}" << std::endl
        << "    Run kernel analyzing global memory performance" << std::endl
        << "  -p|--pinned-memory" << std::endl
        << "    Use pinned Memory instead of pageable memory" << std::endl
        << "  -y|--synchronize-kernel" << std::endl
        << "    Synchronize device after each kernel launch" << std::endl
        << "  -s <memory_size>|--size <memory_size>" << std::endl
        << "    The amount of memory to allcate" << std::endl
        << "  -t <threads_per_block>|--threads-per-block <threads_per_block>" << std::endl
        << "    The number of threads per block" << std::endl
        << "  -g <blocks_per_grid>|--grid-dim <blocks_per_grid>" << std::endl
        << "     The number of blocks per grid" << std::endl
        << "  -i <num_iterations>|--iterations <num_iterations>" << std::endl
        << "     The number of iterations to launch the kernel" << std::endl
        << "  --im <memory-copy-iterations>|--memory-iterations <memory-copy-iterations>" << std::endl
        << "     The number of times the memory will be copied. Use this to get more stable results." << std::endl
        << "  --stride <stride>" << std::endl
        << "     Stride parameter for global-stride test. Not that size parameter is ignored then." << std::endl
        << "  --offset <offset>" << std::endl
        << "     Offset parameter for global-offset test. Not that size parameter is ignored then." << std::endl
        << "" << std::endl;
}
