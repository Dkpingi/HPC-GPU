

//#define bool int
//#define false 0
//#define true -1

#ifdef __cplusplus
extern "C" {
#endif

long calcTime(timeval* /*tp1*/, timeval* /*tp2*/);
void printTime( char* /*label*/, timeval* /*tp1*/, timeval* /*tp2*/, long /*bytes*/, long /*mflops*/);
bool MatrixCompare(float* /*P*/, float* /*Q*/, long /*matWidth*/);
void MatrixMulOnHost(float* /*M*/, float* /*N*/, float* /*P*/, long /*matWidth*/);
long findMin(long /*a*/, long /*b*/);
void MatrixMulOnHostBlocked(float* /*M*/, float* /*N*/, float* /*P*/, long /*matWidth*/, long /*blockSize*/);
void printOutMatrix(float* /*matrix*/, int /*width*/); 

#ifdef __cplusplus
}
#endif
