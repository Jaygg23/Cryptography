#ifndef SHA512_H
#define SHA512_H

#include "crypto_api.h"

#ifdef __cplusplus
extern "C" {
#endif

	// SHA512 관련 함수 선언
	CRYPTO_STATUS sha512_init(SHA512_CTX* ctx);
	CRYPTO_STATUS sha512_update(SHA512_CTX* ctx, const uint8_t* data, size_t len);
	CRYPTO_STATUS sha512_final(SHA512_CTX* ctx, uint8_t* hash);

	// 테스트 함수 (sha512.c 내부 구현)
	int test_sha512(void);

#ifdef __cplusplus
}
#endif

#endif // SHA512_H
