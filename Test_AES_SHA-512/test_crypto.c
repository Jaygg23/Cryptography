#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "crypto_api.h"

// 외부 함수 선언 (aes.c, sha512.c에서 제공)
int test_aes(void);
int test_sha512(void);

int main(void) {
    printf("=======================================\n");
    printf("       CRYPTO FUNCTION TEST START      \n");
    printf("=======================================\n\n");

    printf("[1] AES CTR Mode Test 시작 \n");
    if (test_aes() == 0)
        printf("[AES] 모든 AES CTR Mode 테스트 벡터 통과!\n\n");
    else
        printf("[AES] AES CTR Mode 테스트 실패!\n\n");

    printf("[2] SHA-512 Test 시작...\n");
    if (test_sha512() == 0)
        printf("[SHA512] 모든 SHA-512 테스트 벡터 통과!\n\n");
    else
        printf("[SHA512] SHA-512 테스트 실패!\n\n");

    printf("=======================================\n");
    printf("        모든 CRYPTO TEST 완료!         \n");
    printf("=======================================\n");

    return 0;
}
