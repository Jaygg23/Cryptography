#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "randombytes.h"
#include "sign.h"

#define MLEN 59 // 메시지 길이
#define CTXLEN 14 // 컨텍스트 길이
#define NTESTS 10000 // 테스트 반복 횟수

int main(void)
{
  size_t i, j;
  int ret; // 함수 반환값 저장
  size_t mlen, smlen; // 검증 후 복구된 메시지 길이가 mlen에 채워짐, 서명 결과(sm)의 총 길이
  uint8_t b; // 위변조에서 바이트 단위 변경값
  uint8_t ctx[CTXLEN] = {0}; // 컨텍스트 바이트열 버퍼
  uint8_t m[MLEN + CRYPTO_BYTES]; // 원본 메시지 버퍼
  uint8_t m2[MLEN + CRYPTO_BYTES]; // 검증에 성공했을 때 복구된 메시지를 받는 버퍼
  uint8_t sm[MLEN + CRYPTO_BYTES]; // 서명된 메시지 저장 버퍼
  uint8_t pk[CRYPTO_PUBLICKEYBYTES]; // 공개키 버퍼
  uint8_t sk[CRYPTO_SECRETKEYBYTES]; // 개인키 버퍼

  snprintf((char*)ctx,CTXLEN,"test_dilitium"); //test_dilitium이라는 문자열을 컨텍스트 값으로 해서 ctx 버퍼에 채움

  for(i = 0; i < NTESTS; ++i) { // NTESTS만큼 테스트 반복
    randombytes(m, MLEN); // MLEN 길이의 랜덤 메시지를 m에 저장

    crypto_sign_keypair(pk, sk); // 키 생성
    crypto_sign(sm, &smlen, m, MLEN, ctx, CTXLEN, sk); // 서명
    ret = crypto_sign_open(m2, &mlen, sm, smlen, ctx, CTXLEN, pk); // 검증

    if(ret) {
      fprintf(stderr, "Verification failed\n");
      return -1;
    }
    if(smlen != MLEN + CRYPTO_BYTES) {
      fprintf(stderr, "Signed message lengths wrong\n");
      return -1;
    }
    if(mlen != MLEN) {
      fprintf(stderr, "Message lengths wrong\n");
      return -1;
    }
    for(j = 0; j < MLEN; ++j) {
      if(m2[j] != m[j]) {
        fprintf(stderr, "Messages don't match\n");
        return -1;
      }
    }

    randombytes((uint8_t *)&j, sizeof(j));
    do {
      randombytes(&b, 1);
    } while(!b);
    sm[j % (MLEN + CRYPTO_BYTES)] += b;
    ret = crypto_sign_open(m2, &mlen, sm, smlen, ctx, CTXLEN, pk);
    if(!ret) {
      fprintf(stderr, "Trivial forgeries possible\n");
      return -1;
    }
  }

  printf("CRYPTO_PUBLICKEYBYTES = %d\n", CRYPTO_PUBLICKEYBYTES);
  printf("CRYPTO_SECRETKEYBYTES = %d\n", CRYPTO_SECRETKEYBYTES);
  printf("CRYPTO_BYTES = %d\n", CRYPTO_BYTES);

  return 0;
}
