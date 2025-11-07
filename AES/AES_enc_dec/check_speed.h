#ifndef __CHECK_SPEED_H__
#define __CHECK_SPEED_H__


#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <time.h>


//general
static clock_t START_CLK;  //secs*(CLOCKS_PER_SEC)
static clock_t END_CLK;
static clock_t DELTA_CLK;

static time_t  START_TIME; //secs
static time_t  END_TIME;
static time_t  DELTA_TIME;

double OP_TIME;

#ifdef _MSC_VER //for Windows
#include <Windows.h>
//Checking Speed
static __int64 start_point, end_point, freq;

#define INNER_CHECK_TIME_START {\
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);\
	QueryPerformanceCounter((LARGE_INTEGER*)&start_point);\
}

#define INNER_CHECK_TIME_END {\
	QueryPerformanceCounter((LARGE_INTEGER*)&end_point);\
	OP_TIME = (double)((double)(end_point - start_point) / freq);\
}
#endif //_MSC_VER 


#ifdef __GNUC__ //for gcc
#include <sys/time.h>
//Checking Speed
static struct timeval start_point, end_point;

#define INNER_CHECK_TIME_START {\
	gettimeofday(&start_point, NULL);\
}

#define INNER_CHECK_TIME_END {\
	gettimeofday(&end_point, NULL);\
	OP_TIME = ((double)(end_point.tv_sec) + (double)(end_point.tv_usec) / 1000000.0 - (double)(start_point.tv_sec) - (double)(start_point.tv_usec) / 1000000.0);\
}

#endif //__GNUC__


#define CHECK_TIME_START {\
	START_CLK  = clock();\
	START_TIME = time(NULL);\
	INNER_CHECK_TIME_START\
}

#define CHECK_TIME_END {\
	END_CLK  = clock();\
	END_TIME = time(NULL);\
	INNER_CHECK_TIME_END\
	DELTA_CLK = END_CLK - START_CLK;\
	DELTA_TIME = END_TIME - START_TIME;\
}

#define PRINT_DELTA_TIME(bits) {\
	printf("%u clocks : ",(unsigned int)DELTA_CLK);\
	printf("%u",(unsigned int)DELTA_TIME);\
	printf("(%0.4f)secs", OP_TIME);\
	if((bits) == 0)\
	{\
		printf("\n");\
	}\
	else\
	{\
		printf("[%0.4f Mbps]\n", ((bits)/OP_TIME)/(1000*1000));\
	}\
}


#define NUM_CHECK_TIME    (10ULL)
#define NUM_BLOCKS_OF_AES (500000ULL)


#ifdef __cplusplus
}
#endif /*extern "C"*/
#endif /*__CHECK_SPEED_H__*/