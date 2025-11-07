#ifndef AES_H
#define AES_H

#include "crypto_api.h"

#ifdef __cplusplus
extern "C" {
#endif

	// AES 관련 함수 선언
	CRYPTO_STATUS AES_set_key(AES_CTX* ctx, const uint8_t* key, int key_bits);
	CRYPTO_STATUS AES_encrypt_block(const AES_CTX* ctx, const uint8_t in[AES_BLOCK_SIZE], uint8_t out[AES_BLOCK_SIZE]);
	CRYPTO_STATUS AES_decrypt_block(const AES_CTX* ctx, const uint8_t in[AES_BLOCK_SIZE], uint8_t out[AES_BLOCK_SIZE]);
	CRYPTO_STATUS AES_CTR_crypt(const AES_CTX* ctx, const uint8_t* in, size_t length, uint8_t* out, uint8_t nonce_counter[AES_BLOCK_SIZE]);

	// 테스트 함수 (aes.c 내부 구현)
	int test_aes(void);

#ifdef __cplusplus
}
#endif

#endif // AES_H
