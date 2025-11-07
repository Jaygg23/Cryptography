#include <stdint.h>
#include <string.h>

#include "aes.h"
#include "check_speed.h"


extern char SUBMISSION_INFO[256];


#define AES_TV_CHECHER(kbsz)\
{\
	AES_STATE_t T;\
	int i;\
	printf("##AES%d##\n", kbsz);\
	for (i = 0; i < NUM_AES##kbsz##_TVS; i++)\
	{\
		printf("  > %d tv : ", i);\
		AES##kbsz##_enc(T, aes##kbsz##_tvs[i].P, aes##kbsz##_tvs[i].K##kbsz);\
		if (memcmp(T, aes##kbsz##_tvs[i].C, sizeof(AES_STATE_t)) != 0)\
		{\
			printf("An Error Occurs while Checking the Encryption Function\n");\
			goto err;\
		}\
		AES##kbsz##_dec(T, aes##kbsz##_tvs[i].C, aes##kbsz##_tvs[i].K##kbsz);\
		if (memcmp(T, aes##kbsz##_tvs[i].P, sizeof(AES_STATE_t)) != 0)\
		{\
			printf("An Error Occurs while Checking the Decrpytion Function\n");\
			goto err;\
		}\
		printf("Pass\n");\
	}\
	printf("  >> Good job!! Passed All the Test Vectors for AES%d\n", kbsz);\
}\


#define AES_PF_CHECHER(kbsz)\
{\
	AES_STATE_t T, P, C;\
	AES##kbsz##_KEY_t K##kbsz;\
	double throughput = 0.;\
	int i, j;\
	printf("##AES%d##\n", kbsz);\
	srand((unsigned int)time(NULL));\
	for (i = 0; i < NUM_CHECK_TIME; i++)\
	{\
		for (j = 0; j < 16; j++)\
		{\
			T[j] = (uint8_t)(rand() & 0xff);\
		}\
		for (j = 0; j < (kbsz / 8); j++)\
		{\
			K##kbsz[j] = (uint8_t)(rand() & 0xff);\
		}\
		CHECK_TIME_START;\
		for (j = 0; j < NUM_BLOCKS_OF_AES; j++)\
		{\
			AES##kbsz##_enc(C, T, K##kbsz);\
			AES##kbsz##_dec(P, C, K##kbsz);\
			if (memcmp(T, P, sizeof(AES_STATE_t)) != 0)\
			{\
				printf("An Error Occurs while Checking the Performance!\n");\
				goto err;\
			}\
			else\
			{\
				memcpy(T, C, sizeof(AES_STATE_t));\
			}\
		}\
		CHECK_TIME_END;\
		printf("  > %d its : ", i);\
		PRINT_DELTA_TIME(2 * NUM_BLOCKS_OF_AES * sizeof(AES_STATE_t) * 8);\
		throughput += OP_TIME;\
	}\
	printf("  >> AES%d : On average, %0.4f Mbps\n", kbsz, ((2 * NUM_BLOCKS_OF_AES * 128) / (throughput / NUM_CHECK_TIME)) / (1000 * 1000));\
}\

int main(void)
{


	printf("Submission Info : %s\n", SUBMISSION_INFO);

	//Step 1 : Check Test Vectors
	printf("Step 1 : Check Test Vectors\n");
	AES_TV_CHECHER(128);
	AES_TV_CHECHER(192);
	AES_TV_CHECHER(256);


	//Step 2 : Check Performance
	printf("Step 2 : Check Performance\n");
	AES_PF_CHECHER(128);
	AES_PF_CHECHER(192);
	AES_PF_CHECHER(256);



	
err:

	return 0;
}