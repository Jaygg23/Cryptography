#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "crypto_api.h"
#include "aes.h"

/*****************************************************
 * AES (Advanced Encryption Standard) 내부 헬퍼 함수들
 * 이 함수들은 AES 암호화의 핵심 구성 요소입니다.
 *****************************************************/

 /**
  * @brief xtimes: GF(2^8) 상에서의 곱셈 연산. MixColumns에서 사용됩니다.
  * * AES는 '갈루아 필드(Galois Field)'라는 특수한 수학적 공간에서 연산을 수행합니다.
  * 이 함수는 입력 바이트에 2를 곱하는 연산을 GF(2^8) 규칙에 따라 수행합니다.
  * 일반적인 곱셈과 달리, 최상위 비트(MSB)가 1이면 XOR 연산(0x1b)이 추가로 발생합니다.
  * 이는 곱셈 결과가 1바이트(8비트)를 넘어가지 않도록 보장하기 위함입니다.
  */
#define xtimes(input) (((input) << 1) ^ (((input) >> 7) * 0x1b))

  /**
   * @brief s_box: SubBytes 연산에 사용되는 치환 테이블.
   * * S-Box는 AES의 비선형성을 제공하는 핵심 요소입니다. 입력 바이트를 미리 정해진 다른 바이트로
   * 완전히 대체하여, 암호문이 원래 평문과 통계적 관계를 갖기 어렵게 만듭니다. (혼돈 효과)
   * 이 테이블 값은 특정 수학적 계산(역원 계산 + 아핀 변환)을 통해 설계되었습니다.
   */
static const uint8_t s_box[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
};

/**
 * @brief Rcon: 라운드 상수(Round Constant). 키 스케줄링에서 사용됩니다.
 *
 * 키 스케줄링 과정에서 각 라운드 키가 이전 라운드 키와 달라지도록 보장하는 역할(라운드 대칭성 파괴)을 합니다.
 * 매 라운드마다 이 값을 XOR하여 라운드 간의 독립성을 높입니다.
 */
static const uint8_t Rcon[11] = { 0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36 };

// --- AES의 4가지 기본 연산 ---
// AES는 16바이트(128비트) 데이터를 4x4 행렬(state)로 간주하고,
// 이 4가지 연산을 여러 라운드에 걸쳐 반복적으로 적용합니다.

/**
 * @brief SubBytes: state의 각 바이트를 S-Box를 이용해 치환합니다.
 * @param state 4x4 크기의 데이터 블록
 */
static void SubBytes(uint8_t state[4][4]) {
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) state[i][j] = s_box[state[i][j]];
}

/**
 * @brief ShiftRows: state의 각 행을 정해진 규칙에 따라 왼쪽으로 이동시킵니다.
 * - 0행: 이동 없음
 * - 1행: 1칸 왼쪽 이동
 * - 2행: 2칸 왼쪽 이동
 * - 3행: 3칸 왼쪽 이동
 * 이 연산은 열(column) 단위로 묶여 있던 데이터들을 여러 열에 걸쳐 섞어주는 역할(확산 효과)을 합니다.
 * @param state 4x4 크기의 데이터 블록
 */
static void ShiftRows(uint8_t state[4][4]) {
    uint8_t temp;
    // 1행
    temp = state[1][0]; state[1][0] = state[1][1]; state[1][1] = state[1][2]; state[1][2] = state[1][3]; state[1][3] = temp;
    // 2행
    temp = state[2][0]; state[2][0] = state[2][2]; state[2][2] = temp;
    temp = state[2][1]; state[2][1] = state[2][3]; state[2][3] = temp;
    // 3행
    temp = state[3][3]; state[3][3] = state[3][2]; state[3][2] = state[3][1]; state[3][1] = state[3][0]; state[3][0] = temp;
}

/**
 * @brief MixColumns: state의 각 열을 GF(2^8) 상에서 특정 행렬과 곱셈 연산을 수행합니다.
 * ShiftRows가 행 단위로 데이터를 섞었다면, MixColumns는 열 단위로 데이터를 섞습니다.
 * 이 두 연산을 통해 16바이트 블록 전체에 걸쳐 데이터가 골고루 섞이게 됩니다. (확산 효과 극대화)
 * @param state 4x4 크기의 데이터 블록
 */
static void MixColumns(uint8_t state[4][4]) {
    uint8_t col[4], res[4];
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) col[i] = state[i][j];
        res[0] = xtimes(col[0]) ^ (xtimes(col[1]) ^ col[1]) ^ col[2] ^ col[3];
        res[1] = col[0] ^ xtimes(col[1]) ^ (xtimes(col[2]) ^ col[2]) ^ col[3];
        res[2] = col[0] ^ col[1] ^ xtimes(col[2]) ^ (xtimes(col[3]) ^ col[3]);
        res[3] = (xtimes(col[0]) ^ col[0]) ^ col[1] ^ col[2] ^ xtimes(col[3]);
        for (int i = 0; i < 4; i++) state[i][j] = res[i];
    }
}

/**
 * @brief AddRoundKey: state와 현재 라운드 키를 XOR 연산합니다.
 * @param state 4x4 크기의 데이터 블록
 * @param round_key 현재 라운드에서 사용할 16바이트 라운드 키
 */
static void AddRoundKey(uint8_t state[4][4], const uint8_t* round_key) {
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) state[j][i] ^= round_key[i * 4 + j];
}


// --- 키 스케줄링용 헬퍼 함수 ---
// 사용자가 제공한 마스터 키로부터 각 라운드에서 사용할 라운드 키들을 생성합니다.

static void RotWord(uint8_t* word) { // 4바이트(1워드)를 왼쪽으로 1바이트씩 순환
    uint8_t temp = word[0];
    word[0] = word[1]; word[1] = word[2]; word[2] = word[3]; word[3] = temp;
}

static void SubWord(uint8_t* word) { // 4바이트(1워드)의 각 바이트를 S-Box로 치환
    word[0] = s_box[word[0]]; word[1] = s_box[word[1]];
    word[2] = s_box[word[2]]; word[3] = s_box[word[3]];
}

/*****************************************************
 * 키 스케줄링(Key Expansion) 함수들
 * 마스터 키로부터 모든 라운드 키를 생성하여 AES_CTX에 저장합니다.
 * 키 길이에 따라 생성되는 라운드 키의 개수와 생성 방식이 약간씩 다릅니다.
 *****************************************************/
static void KeySchedule128(const uint8_t* key, AES_CTX* ctx) {
    uint8_t temp[4]; // 임시 4바이트 워드
    uint8_t* w = ctx->round_keys; // 생성된 라운드 키들이 저장될 배열
    memcpy(w, key, 16); // 첫 16바이트는 마스터 키 그대로 사용

    // 나머지 라운드 키 생성 (4워드씩 생성)
    for (int i = 4; i < 4 * (ctx->Nr + 1); i++) {
        memcpy(temp, &w[(i - 1) * 4], 4); // 이전 워드를 temp에 복사
        if (i % ctx->Nk == 0) { // 특정 조건(Nk의 배수)일 때
            RotWord(temp); // 워드 회전
            SubWord(temp); // S-Box 치환
            temp[0] ^= Rcon[i / ctx->Nk]; // 라운드 상수와 XOR
        }
        // 새 워드 = (i-Nk)번째 워드 XOR temp
        w[i * 4 + 0] = w[(i - ctx->Nk) * 4 + 0] ^ temp[0];
        w[i * 4 + 1] = w[(i - ctx->Nk) * 4 + 1] ^ temp[1];
        w[i * 4 + 2] = w[(i - ctx->Nk) * 4 + 2] ^ temp[2];
        w[i * 4 + 3] = w[(i - ctx->Nk) * 4 + 3] ^ temp[3];
    }
}

static void KeySchedule192(const uint8_t* key, AES_CTX* ctx) {
    uint8_t temp[4];
    uint8_t* w = ctx->round_keys;
    memcpy(w, key, 24); // 마스터 키 24바이트 복사

    for (int i = 6; i < 4 * (ctx->Nr + 1); i++) {
        memcpy(temp, &w[(i - 1) * 4], 4);
        if (i % ctx->Nk == 0) {
            RotWord(temp);
            SubWord(temp);
            temp[0] ^= Rcon[i / ctx->Nk];
        }
        w[i * 4 + 0] = w[(i - ctx->Nk) * 4 + 0] ^ temp[0];
        w[i * 4 + 1] = w[(i - ctx->Nk) * 4 + 1] ^ temp[1];
        w[i * 4 + 2] = w[(i - ctx->Nk) * 4 + 2] ^ temp[2];
        w[i * 4 + 3] = w[(i - ctx->Nk) * 4 + 3] ^ temp[3];
    }
}

static void KeySchedule256(const uint8_t* key, AES_CTX* ctx) {
    uint8_t temp[4];
    uint8_t* w = ctx->round_keys;
    memcpy(w, key, 32); // 마스터 키 32바이트 복사

    for (int i = 8; i < 4 * (ctx->Nr + 1); i++) {
        memcpy(temp, &w[(i - 1) * 4], 4);
        if (i % ctx->Nk == 0) {
            RotWord(temp);
            SubWord(temp);
            temp[0] ^= Rcon[i / ctx->Nk];
        }
        else if (i % ctx->Nk == 4) { // 256비트 키 스케줄에만 있는 추가 규칙
            SubWord(temp);
        }
        w[i * 4 + 0] = w[(i - ctx->Nk) * 4 + 0] ^ temp[0];
        w[i * 4 + 1] = w[(i - ctx->Nk) * 4 + 1] ^ temp[1];
        w[i * 4 + 2] = w[(i - ctx->Nk) * 4 + 2] ^ temp[2];
        w[i * 4 + 3] = w[(i - ctx->Nk) * 4 + 3] ^ temp[3];
    }
}


/*****************************************************
 * Crypto API 함수 구현
 * 헤더 파일(crypto_api.h)에 선언된 함수들을 실제로 구현하는 부분입니다.
 *****************************************************/

 /**
  * @brief AES_set_key: 사용자가 제공한 마스터 키로 AES 컨텍스트를 초기화하고,
  * 모든 라운드 키를 미리 생성합니다.
  * @param ctx AES 상태를 저장할 컨텍스트 구조체 포인터
  * @param key 마스터 키
  * @param key_bits 키의 비트 길이 (128, 192, 256)
  * @return 성공 시 CRYPTO_SUCCESS, 실패 시 오류 코드
  */
CRYPTO_STATUS AES_set_key(AES_CTX* ctx, const uint8_t* key, int key_bits) {
    if (!ctx || !key) return CRYPTO_ERR_NULL_CONTEXT;

    ctx->key_bits = key_bits;
    ctx->Nk = key_bits / 32; // Nk: 키 길이를 32비트 워드 단위로 나타낸 값

    switch (key_bits) {
    case 128:
        ctx->Nr = AES_ROUND_128; // Nr: 라운드 수
        KeySchedule128(key, ctx);
        break;
    case 192:
        ctx->Nr = AES_ROUND_192;
        KeySchedule192(key, ctx);
        break;
    case 256:
        ctx->Nr = AES_ROUND_256;
        KeySchedule256(key, ctx);
        break;
    default:
        return CRYPTO_ERR_INVALID_ARGUMENT; // 지원하지 않는 키 길이
    }
    return CRYPTO_SUCCESS;
}

/**
 * @brief AES_encrypt_block: 16바이트 평문 블록 하나를 암호화합니다.
 * @param ctx 초기화된 AES 컨텍스트
 * @param in 16바이트 평문 블록
 * @param out 암호화된 결과가 저장될 16바이트 버퍼
 * @return 성공 시 CRYPTO_SUCCESS
 */
CRYPTO_STATUS AES_encrypt_block(const AES_CTX* ctx, const uint8_t in[AES_BLOCK_SIZE], uint8_t out[AES_BLOCK_SIZE]) {
    if (!ctx || !in || !out) return CRYPTO_ERR_NULL_CONTEXT;

    // 1. 입력 평문을 4x4 state 행렬로 변환
    uint8_t state[4][4];
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) state[j][i] = in[i * 4 + j];

    // 2. 초기 라운드: AddRoundKey
    AddRoundKey(state, ctx->round_keys);

    // 3. 메인 라운드: (Nr - 1)번 반복
    for (int r = 1; r < ctx->Nr; r++) {
        SubBytes(state);
        ShiftRows(state);
        MixColumns(state);
        AddRoundKey(state, ctx->round_keys + r * 16); // r번째 라운드 키 사용
    }

    // 4. 마지막 라운드: MixColumns 제외
    SubBytes(state);
    ShiftRows(state);
    AddRoundKey(state, ctx->round_keys + ctx->Nr * 16);

    // 5. state 행렬을 출력 버퍼로 변환
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) out[i * 4 + j] = state[j][i];
    return CRYPTO_SUCCESS;
}

// AES_decrypt_block 은 이번 과제에서 필요 없으므로 구현하지 않음
CRYPTO_STATUS AES_decrypt_block(const AES_CTX* ctx, const uint8_t in[AES_BLOCK_SIZE], uint8_t out[AES_BLOCK_SIZE]) {
    return CRYPTO_ERR_INTERNAL_FAILURE;
}


/**
 * @brief AES_CTR_crypt: AES 카운터(CTR) 모드로 암호화 또는 복호화를 수행합니다.
 * * CTR 모드는 블록 암호인 AES를 스트림 암호처럼 사용할 수 있게 해줍니다.
 * Nonce(재사용 금지 값)와 Counter를 합쳐 암호화한 결과를 '키스트림'으로 만들고,
 * 이 키스트림을 평문(또는 암호문)과 XOR하여 암호문(또는 평문)을 생성합니다.
 * 암호화와 복호화 과정이 동일하다는 장점이 있습니다.
 * * @param ctx 초기화된 AES 컨텍스트
 * @param in 입력 데이터 (평문 또는 암호문)
 * @param length 입력 데이터의 길이 (바이트)
 * @param out 출력 데이터가 저장될 버퍼
 * @param nonce_counter 16바이트 Nonce+Counter 블록. 함수 호출 후 자동으로 1씩 증가합니다.
 * @return 성공 시 CRYPTO_SUCCESS
 */
CRYPTO_STATUS AES_CTR_crypt(const AES_CTX* ctx, const uint8_t* in, size_t length, uint8_t* out, uint8_t nonce_counter[AES_BLOCK_SIZE]) {
    if (!ctx) return CRYPTO_ERR_NULL_CONTEXT;
    if ((!in || !out) && length > 0) return CRYPTO_ERR_INVALID_INPUT;
    if (!nonce_counter) return CRYPTO_ERR_INVALID_INPUT;

    uint8_t counter_block[AES_BLOCK_SIZE]; // 현재 카운터 값 (암호화 대상)
    uint8_t keystream_block[AES_BLOCK_SIZE]; // 카운터를 암호화한 결과 (키스트림)
    size_t offset = 0;

    memcpy(counter_block, nonce_counter, AES_BLOCK_SIZE);

    // 데이터를 16바이트 블록 단위로 처리
    while (length > 0) {
        // 1. 현재 카운터 블록을 AES로 암호화하여 키스트림 생성
        AES_encrypt_block(ctx, counter_block, keystream_block);

        // 처리할 데이터 길이 결정 (마지막 블록은 16바이트보다 작을 수 있음)
        size_t block_len = (length < AES_BLOCK_SIZE) ? length : AES_BLOCK_SIZE;

        // 2. 입력 데이터와 키스트림을 XOR하여 출력 생성
        for (size_t i = 0; i < block_len; i++) {
            out[offset + i] = in[offset + i] ^ keystream_block[i];
        }

        length -= block_len;
        offset += block_len;

        // 3. 다음 블록을 위해 카운터 값을 1 증가 (big-endian 방식)
        for (int i = AES_BLOCK_SIZE - 1; i >= 0; i--) {
            if (++counter_block[i] != 0) break; // 자리올림이 없으면 중단
        }
    }

    // 최종적으로 증가된 카운터 값을 원래 nonce_counter 배열에 업데이트
    memcpy(nonce_counter, counter_block, AES_BLOCK_SIZE);
    return CRYPTO_SUCCESS;
}

/*****************************************************
 * NIST SP 800-38A CTR 모드 테스트
 * 이 main 함수는 구현된 AES-CTR 코드가 표준 테스트 값으로
 * 정확하게 동작하는지 검증하는 역할을 합니다.
 *****************************************************/

 // 헬퍼 함수: 데이터를 16진수 문자열로 예쁘게 출력
void print_hex(const char* label, const unsigned char* data, int len) {
    printf("%-14s: ", label);
    for (int i = 0; i < len; i++) {
        printf("%02x", data[i]);
    }
    printf("\n");
}

// 헬퍼 함수: 두 데이터 배열이 일치하는지 비교
int compare_hex(const unsigned char* d1, const unsigned char* d2, int len) {
    return memcmp(d1, d2, len) == 0;
}

int test_aes() {
    printf("=======================================\n");
    printf("  NIST SP 800-38A CTR Mode Test Vector \n");
    printf("=======================================\n");

    AES_CTX ctx;

    // --- AES-128 CTR Test ---
    printf("--- AES-128 CTR Test ---\n");
    uint8_t key128[] = { 0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c };
    uint8_t pt128[] = { 0x6b, 0xc1, 0xbe, 0xe2, 0x2e, 0x40, 0x9f, 0x96, 0xe9, 0x3d, 0x7e, 0x11, 0x73, 0x93, 0x17, 0x2a };
    uint8_t iv128[] = { 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff };
    uint8_t expected_ct128[] = { 0x87, 0x4d, 0x61, 0x91, 0xb6, 0x20, 0xe3, 0x26, 0x1b, 0xef, 0x68, 0x64, 0x99, 0x0d, 0xb6, 0xce };
    uint8_t ct128[AES_BLOCK_SIZE];

    AES_set_key(&ctx, key128, 128); // 1. 키 설정
    AES_CTR_crypt(&ctx, pt128, sizeof(pt128), ct128, iv128); // 2. 암호화 수행

    print_hex("Key", key128, sizeof(key128));
    print_hex("Plaintext", pt128, sizeof(pt128));
    print_hex("IV", (uint8_t[]) { 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff }, 16);
    print_hex("Result", ct128, sizeof(ct128));
    print_hex("Expected", expected_ct128, sizeof(expected_ct128));
    printf("Verification  : %s\n\n", compare_hex(ct128, expected_ct128, sizeof(ct128)) ? "SUCCESS" : "FAILURE");

    // --- AES-192 CTR Test ---
    printf("--- AES-192 CTR Test ---\n");
    uint8_t key192[] = { 0x8e, 0x73, 0xb0, 0xf7, 0xda, 0x0e, 0x64, 0x52, 0xc8, 0x10, 0xf3, 0x2b, 0x80, 0x90, 0x79, 0xe5, 0x62, 0xf8, 0xea, 0xd2, 0x52, 0x2c, 0x6b, 0x7b };
    uint8_t pt192[] = { 0x6b, 0xc1, 0xbe, 0xe2, 0x2e, 0x40, 0x9f, 0x96, 0xe9, 0x3d, 0x7e, 0x11, 0x73, 0x93, 0x17, 0x2a };
    uint8_t iv192[] = { 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff };
    uint8_t expected_ct192[] = { 0x1a, 0xbc, 0x93, 0x24, 0x17, 0x52, 0x1c, 0xa2, 0x4f, 0x2b, 0x04, 0x59, 0xfe, 0x7e, 0x6e, 0x0b };
    uint8_t ct192[AES_BLOCK_SIZE];

    AES_set_key(&ctx, key192, 192);
    AES_CTR_crypt(&ctx, pt192, sizeof(pt192), ct192, iv192);

    print_hex("Key", key192, sizeof(key192));
    print_hex("Plaintext", pt192, sizeof(pt192));
    print_hex("IV", (uint8_t[]) { 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff }, 16);
    print_hex("Result", ct192, sizeof(ct192));
    print_hex("Expected", expected_ct192, sizeof(expected_ct192));
    printf("Verification  : %s\n\n", compare_hex(ct192, expected_ct192, sizeof(ct192)) ? "SUCCESS" : "FAILURE");

    // --- AES-256 CTR Test ---
    printf("--- AES-256 CTR Test ---\n");
    uint8_t key256[] = { 0x60, 0x3d, 0xeb, 0x10, 0x15, 0xca, 0x71, 0xbe, 0x2b, 0x73, 0xae, 0xf0, 0x85, 0x7d, 0x77, 0x81, 0x1f, 0x35, 0x2c, 0x07, 0x3b, 0x61, 0x08, 0xd7, 0x2d, 0x98, 0x10, 0xa3, 0x09, 0x14, 0xdf, 0xf4 };
    uint8_t pt256[] = { 0x6b, 0xc1, 0xbe, 0xe2, 0x2e, 0x40, 0x9f, 0x96, 0xe9, 0x3d, 0x7e, 0x11, 0x73, 0x93, 0x17, 0x2a };
    uint8_t iv256[] = { 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff };
    uint8_t expected_ct256[] = { 0x60, 0x1e, 0xc3, 0x13, 0x77, 0x57, 0x89, 0xa5, 0xb7, 0xa7, 0xf5, 0x04, 0xbb, 0xf3, 0xd2, 0x28 };
    uint8_t ct256[AES_BLOCK_SIZE];

    AES_set_key(&ctx, key256, 256);
    AES_CTR_crypt(&ctx, pt256, sizeof(pt256), ct256, iv256);

    print_hex("Key", key256, sizeof(key256));
    print_hex("Plaintext", pt256, sizeof(pt256));
    print_hex("IV", (uint8_t[]) { 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff }, 16);
    print_hex("Result", ct256, sizeof(ct256));
    print_hex("Expected", expected_ct256, sizeof(expected_ct256));
    printf("Verification  : %s\n\n", compare_hex(ct256, expected_ct256, sizeof(ct256)) ? "SUCCESS" : "FAILURE");

    return 0;
}
