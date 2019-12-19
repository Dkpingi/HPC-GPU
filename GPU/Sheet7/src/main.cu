/******************************************************************************
 *
 *           XXXII Heidelberg Physics Graduate Days - GPU Computing
 *
 *                 Gruppe : TODO
 *
 *                   File : main.cu
 *
 *                Purpose : n-Body Computation
 *
 ******************************************************************************/

#include <cmath>
#include <ctime>
#include <iostream>
#include <cstdlib>
#include <chCommandLine.h>
#include <chTimer.hpp>
#include <cstdio>
#include <iomanip>
#include <fstream>

const static int DEFAULT_NUM_ELEMENTS   = 1024;
const static int DEFAULT_NUM_ITERATIONS =    5;
const static int DEFAULT_BLOCK_DIM      =  128;

const static float TIMESTEP =      1e-6; // s
const static float GAMMA    = 6.673e-11; // (Nm^2)/(kg^2)

//
// Structures
//
// Use a SOA (Structure of Arrays)
//
struct Body_t {
	float4* posMass;  /* x = x */
	                  /* y = y */
	                  /* z = z */
	                  /* w = Mass */
	float3* velocity; /* x = v_x*/
	                  /* y = v_y */
	                  /* z= v_z */
	
	Body_t(): posMass(NULL), velocity(NULL) {}
	};

//
// Function Prototypes
//
void printHelp(char *);
void printElement(Body_t, int, int);

//
// Device Functions
//

//
// Calculate the Distance of two points
//
__device__ float
getDistance(float4 a, float4 b)
{
	// TODO: Calculate distance of two particles
	return norm3df(a.x - b.x,a.y - b.y,a.z - b.z);
}

//
// Calculate the forces between two bodies
//
__device__ void
bodyBodyInteraction(float4 bodyA, float4 bodyB, float3& force)
{
	float distance = getDistance(bodyA, bodyB);

	if (distance==0)
	{
		return;
	}
	else
	{	
		float invDistanceCubed = 1/(distance*distance*distance);
		force.x = GAMMA*bodyA.w*bodyB.w*(bodyB.x - bodyA.x)*invDistanceCubed;
		force.y = GAMMA*bodyA.w*bodyB.w*(bodyB.y - bodyA.y)*invDistanceCubed;
		force.z = GAMMA*bodyA.w*bodyB.w*(bodyB.z - bodyA.z)*invDistanceCubed;
	}
	// TODO: Calc Force
}

__device__ void
bodyBodyInteraction_opt(float4 bodyA, float4 bodyB, float3& accel)
{
	float dx = (bodyB.x - bodyA.x);
	float dy = (bodyB.y - bodyA.y);
	float dz = (bodyB.z - bodyA.z);
	float InvDistance = rsqrtf(dx*dx + dy*dy + dz*dz);
	float CinvDistanceCubed = bodyB.w*(InvDistance*InvDistance*InvDistance);
	accel.x += dx*CinvDistanceCubed;
	accel.y += dy*CinvDistanceCubed;
	accel.z += dz*CinvDistanceCubed;
	// TODO: Calc Force
}

//
// Calculate the new velocity of one particle
//
__device__ void
calculateSpeed(float mass, float3& currentSpeed, float3 force)
{	
	currentSpeed.x += (force.x/mass)*TIMESTEP;
	currentSpeed.y += (force.y/mass)*TIMESTEP;
	currentSpeed.z += (force.z/mass)*TIMESTEP;
	
	// TODO: Calculate the new velocity of a particle
}

__device__ void
calculateSpeed_opt(float3& currentSpeed, float3 accel)
{
	currentSpeed.x += accel.x*TIMESTEP;
	currentSpeed.y += accel.y*TIMESTEP;
	currentSpeed.z += accel.z*TIMESTEP;
	
	// TODO: Calculate the new velocity of a particle
}

//
// n-Body Kernel for the speed calculation
//
__global__ void
simpleNbody_Kernel(int numElements, float4* bodyPos, float3* bodySpeed)
{
	int elementId = blockIdx.x * blockDim.x + threadIdx.x;
	
	for(int i = elementId; i<numElements; i+= gridDim.x*blockDim.x)
	{
		float4 elementPosMass = bodyPos[i];
		float3 elementSpeed = bodySpeed[i];
		float3 elementForce = make_float3(0,0,0);
		for (int j = 0; j < numElements; j++) {
			float3 tempForce = make_float3(0,0,0);
			bodyBodyInteraction(elementPosMass, bodyPos[j], tempForce);
			elementForce.x += tempForce.x;
			elementForce.y += tempForce.y;
			elementForce.z += tempForce.z;
		}

		calculateSpeed(elementPosMass.w, elementSpeed, elementForce);

		bodySpeed[i] = elementSpeed;
	}
}

__global__ void
sharedNbody_Kernel(const int numElements, const float4* bodyPos, float3* bodySpeed)
{
	extern __shared__ float4 posMass[];
	const int elementId = blockIdx.x * blockDim.x + threadIdx.x;
	const int tid = threadIdx.x;
	for(int i = elementId; i<numElements; i += blockDim.x*gridDim.x)
	{
		const float4 elementbodyPos = bodyPos[i];
		float3 elementAccel = make_float3(0,0,0);
		#pragma unroll 16 // slight performance increase
		for(int j = 0;j<numElements;j+= blockDim.x)
		{
			posMass[tid] = bodyPos[tid + j];
			posMass[tid].w *= GAMMA;
			__syncthreads();
			if(abs(j-i) < blockDim.x)
			{
				#pragma unroll 32 // significant performance increase 
				for(int k = 0; k<blockDim.x; k++)
				{
					if(i != j + k)
					{
						bodyBodyInteraction_opt(elementbodyPos, posMass[k], elementAccel);
					}
				}
			}
			else
			{
				#pragma unroll 32
				for(int k = 0; k<blockDim.x; k++)
				{
					bodyBodyInteraction_opt(elementbodyPos, posMass[k], elementAccel);
				}

			}
			__syncthreads();
		}
		calculateSpeed_opt(bodySpeed[i], elementAccel);
	}		
}

//
// n-Body Kernel to update the position
// Neended to prevent write-after-read-hazards
//
__global__ void
updatePosition_Kernel(int numElements, float4* bodyPos, float3* bodySpeed)
{
	int elementId = blockIdx.x * blockDim.x + threadIdx.x;

	for(int i = elementId; i<numElements; i+= gridDim.x*blockDim.x)
	{
		bodyPos[i].x += bodySpeed[i].x*TIMESTEP;
		bodyPos[i].y += bodySpeed[i].y*TIMESTEP;
		bodyPos[i].z += bodySpeed[i].z*TIMESTEP;
	}
}

//
// Main
//
int
main(int argc, char * argv[])
{
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
	ChTimer kernelTimer;

	//
	// Allocate Memory
	//
	int numElements = 0;
	chCommandLineGet<int>(&numElements, "s", argc, argv);
	chCommandLineGet<int>(&numElements, "size", argc, argv);
	numElements = numElements != 0 ?
			numElements : DEFAULT_NUM_ELEMENTS;
	//
	// Host Memory
	//
	bool pinnedMemory = chCommandLineGetBool("p", argc, argv);
	if (!pinnedMemory) {
		pinnedMemory = chCommandLineGetBool("pinned-memory",argc,argv);
	}

	Body_t h_particles;
	if (!pinnedMemory) {
		// Pageable
		h_particles.posMass = static_cast<float4*>
				(malloc(static_cast<size_t>
				(numElements * sizeof(*(h_particles.posMass)))));
		h_particles.velocity = static_cast<float3*>
				(malloc(static_cast<size_t>
				(numElements * sizeof(*(h_particles.velocity)))));
	} else {
		// Pinned
		cudaMallocHost(&(h_particles.posMass), 
				static_cast<size_t>
				(numElements * sizeof(*(h_particles.posMass))));
		cudaMallocHost(&(h_particles.velocity), 
				static_cast<size_t>
				(numElements * sizeof(*(h_particles.velocity))));
	}

	// Init Particles
//	srand(static_cast<unsigned>(time(0)));
	srand(0); // Always the same random numbers
	for (int i = 0; i < numElements; i++) {
		h_particles.posMass[i].x = 1e-8*static_cast<float>(rand()); // Modify the random values to
		h_particles.posMass[i].y = 1e-8*static_cast<float>(rand()); // increase the position changes
		h_particles.posMass[i].z = 1e-8*static_cast<float>(rand()); // and the velocity
		h_particles.posMass[i].w =  1e4*static_cast<float>(rand());
		h_particles.velocity[i].x = 0.0f;
		h_particles.velocity[i].y = 0.0f;
		h_particles.velocity[i].z = 0.0f;
	}
	
	//printElement(h_particles, 0, 0);

	// Device Memory
	Body_t d_particles;
	cudaMalloc(&(d_particles.posMass), 
			static_cast<size_t>(numElements * sizeof(*(d_particles.posMass))));
	cudaMalloc(&(d_particles.velocity), 
			static_cast<size_t>(numElements * sizeof(*(d_particles.velocity))));

	if (h_particles.posMass == NULL || h_particles.velocity == NULL ||
		d_particles.posMass == NULL || d_particles.velocity == NULL) {
		std::cout << "\033[31m***" << std::endl
		          << "*** Error - Memory allocation failed" << std::endl
		          << "***\033[0m" << std::endl;

		exit(-1);
	}

	//
	// Copy Data to the Device
	//
	memCpyH2DTimer.start();

	cudaMemcpy(d_particles.posMass, h_particles.posMass, 
			static_cast<size_t>(numElements * sizeof(float4)), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(d_particles.velocity, h_particles.velocity, 
			static_cast<size_t>(numElements * sizeof(float3)), 
			cudaMemcpyHostToDevice);

	memCpyH2DTimer.stop();

	//
	// Get Kernel Launch Parameters
	//
	int blockSize = 0,
	    gridSize = 0,
	    numIterations = 0;

	// Number of Iterations	
	chCommandLineGet<int>(&numIterations,"i", argc, argv);
	chCommandLineGet<int>(&numIterations,"num-iterations", argc, argv);
	numIterations = numIterations != 0 ? numIterations : DEFAULT_NUM_ITERATIONS;

	// Block Dimension / Threads per Block
	chCommandLineGet<int>(&blockSize,"t", argc, argv);
	chCommandLineGet<int>(&blockSize,"threads-per-block", argc, argv);
	blockSize = blockSize != 0 ? blockSize : DEFAULT_BLOCK_DIM;

	if (blockSize > 1024) {
		std::cout << "\033[31m***" << std::endl
		          << "*** Error - The number of threads per block is too big" << std::endl
		          << "***\033[0m" << std::endl;

		exit(-1);
	}

	gridSize = ceil(static_cast<float>(numElements) / static_cast<float>(blockSize));

	dim3 grid_dim = dim3(gridSize);
	dim3 block_dim = dim3(blockSize);

	std::cout << "***" << std::endl;
	std::cout << "*** Grid: " << gridSize << std::endl;
	std::cout << "*** Block: " << blockSize << std::endl;
	std::cout << "***" << std::endl;


	bool optimized = chCommandLineGetBool("opt", argc, argv);
	
	if(!optimized)
	{
	
		kernelTimer.start();

		for (int i = 0; i < numIterations; i ++) {
			simpleNbody_Kernel<<<grid_dim, block_dim>>>(numElements, d_particles.posMass, 
					d_particles.velocity);
			updatePosition_Kernel<<<grid_dim, block_dim>>>(numElements, d_particles.posMass,
					d_particles.velocity);
	
			//cudaMemcpy(h_particles.posMass, d_particles.posMass, sizeof(float4), cudaMemcpyDeviceToHost);
			//cudaMemcpy(h_particles.velocity, d_particles.velocity, sizeof(float3), cudaMemcpyDeviceToHost);
			//printElement(h_particles, 0, i+1);
		}
		// Synchronize
		cudaDeviceSynchronize();

		// Check for Errors
		cudaError_t cudaError = cudaGetLastError();
		if ( cudaError != cudaSuccess ) {
			std::cout << "\033[31m***" << std::endl
			          << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
			          << std::endl
			          << "***\033[0m" << std::endl;

			return -1;
		}

		kernelTimer.stop();

	}
	else
	{
		kernelTimer.start();

		for (int i = 0; i < numIterations; i ++) {
			sharedNbody_Kernel<<<grid_dim, block_dim,block_dim.x*4*sizeof(float)>>>(numElements, 						d_particles.posMass, 
					d_particles.velocity);
			updatePosition_Kernel<<<grid_dim, block_dim>>>(numElements, d_particles.posMass,
					d_particles.velocity);

			//cudaMemcpy(h_particles.posMass, d_particles.posMass, sizeof(float4), cudaMemcpyDeviceToHost);
			//cudaMemcpy(h_particles.velocity, d_particles.velocity, sizeof(float3), cudaMemcpyDeviceToHost);
			//printElement(h_particles, 0, i+1);
		}
		//printElement(h_particles, 0, numIterations);
		// Synchronize
		cudaDeviceSynchronize();

		// Check for Errors
		cudaError_t cudaError = cudaGetLastError();
		if ( cudaError != cudaSuccess ) {
			std::cout << "\033[31m***" << std::endl
			          << "***ERROR*** " << cudaError << " - " << cudaGetErrorString(cudaError)
			          << std::endl
			          << "***\033[0m" << std::endl;
	
			return -1;
		}

		kernelTimer.stop();

	}
	//
	// Copy Back Data
	//
	memCpyD2HTimer.start();
	
	cudaMemcpy(h_particles.posMass, d_particles.posMass, 
			static_cast<size_t>(numElements * sizeof(*(h_particles.posMass))), 
			cudaMemcpyDeviceToHost);
	cudaMemcpy(h_particles.velocity, d_particles.velocity, 
			static_cast<size_t>(numElements * sizeof(*(h_particles.velocity))), 
			cudaMemcpyDeviceToHost);

	memCpyD2HTimer.stop();

	printElement(h_particles, 0, numIterations);
	// Free Memory
	if (!pinnedMemory) {
		free(h_particles.posMass);
		free(h_particles.velocity);
	} else {
		cudaFreeHost(h_particles.posMass);
		cudaFreeHost(h_particles.velocity);
	}

	cudaFree(d_particles.posMass);
	cudaFree(d_particles.velocity);
	
	// Print Meassurement Results
    std::cout << "***" << std::endl
	      << "*** Mode: opt = " << optimized << std::endl
              << "*** Results:" << std::endl
              << "***    Num Elements: " << numElements << std::endl
              << "***    Time to Copy to Device: " << 1e3 * memCpyH2DTimer.getTime()
                << " ms" << std::endl
              << "***    Copy Bandwidth: " 
                << 1e-9 * memCpyH2DTimer.getBandwidth(numElements * sizeof(h_particles))
                << " GB/s" << std::endl
              << "***    Time to Copy from Device: " << 1e3 * memCpyD2HTimer.getTime()
                << " ms" << std::endl
              << "***    Copy Bandwidth: " 
                << 1e-9 * memCpyD2HTimer.getBandwidth(numElements * sizeof(h_particles))
                << " GB/s" << std::endl
              << "***    Time for n-Body Computation: " << 1e3 * kernelTimer.getTime()
                << " ms" << std::endl
		<< "***    Body-body interactions/s [G]: " <<  1e-9*numElements*(numElements-1)*numIterations/kernelTimer.getTime()
                << std::endl
              << "***" << std::endl;
	
	std::fstream fin,fout;
	const char* fname;
	if(!optimized)
	{
		fname = "nBody.txt";
	}
	else
	{
		fname = "nBody_opt.txt";
	}
	fout.open(fname, std::ios::in | std::ios::out | std::ios::app);
	fin.open(fname, std::ios::in);
	if (!fout)
	{
	std::cerr << "file open failed:" << std::endl; 
	}
	if(fin.peek() == std::fstream::traits_type::eof())
	{
		fout << "numElements" << " " << "blockDim" << " " << "bodyBodyInt[G]" << "\n";
	} 
		fout << numElements << " " << block_dim.x << " " << 1e-9*numElements*(numElements-1)*numIterations/kernelTimer.getTime() << "\n"; 
		
	fin.close();
	fout.close();
	return 0;
}

void
printHelp(char * argv)
{
    std::cout << "Help:" << std::endl
              << "  Usage: " << std::endl
              << "  " << argv << " [-p] [-s <num-elements>] [-t <threads_per_block>]"
                  << std::endl
              << "" << std::endl
              << "  -p|--pinned-memory" << std::endl
              << "    Use pinned Memory instead of pageable memory" << std::endl
              << "" << std::endl
              << "  -s <num-elements>|--size <num-elements>" << std::endl
              << "    Number of elements (particles)" << std::endl
              << "" << std::endl
              << "  -t <threads_per_block>|--threads-per-block <threads_per_block>" 
                  << std::endl
              << "    The number of threads per block" << std::endl
              << "" << std::endl;
}

//
// Print one element
//
void
printElement(Body_t particles, int elementId, int iteration)
{
    float4 posMass = particles.posMass[elementId];
    float3 velocity = particles.velocity[elementId];

    std::cout << "***" << std::endl
              << "*** Printing Element " << elementId << " in iteration " << iteration << std::endl
              << "***" << std::endl
              << "*** Position: <" 
                  << std::setw(11) << std::setprecision(9) << posMass.x << "|"
                  << std::setw(11) << std::setprecision(9) << posMass.y << "|"
                  << std::setw(11) << std::setprecision(9) << posMass.z << "> [m]" << std::endl
              << "*** velocity: <" 
                  << std::setw(11) << std::setprecision(9) << velocity.x << "|"
                  << std::setw(11) << std::setprecision(9) << velocity.y << "|"
                  << std::setw(11) << std::setprecision(9) << velocity.z << "> [m/s]" << std::endl
              << "*** Mass: " 
                  << std::setw(11) << std::setprecision(9) << posMass.w << " kg"<< std::endl
              << "***" << std::endl;
}
