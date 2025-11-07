#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <stdint.h>
#include <stddef.h>
#include "crypto_api.h"
#include "sha512.h"

#define SHA512_BLOCK_SIZE 128 // 128 바이트 = 1024비트
#define SHA512_DIGEST_LENGTH 64

static const uint64_t K[80] = {
    0x428a2f98d728ae22,0x7137449123ef65cd,0xb5c0fbcfec4d3b2f,0xe9b5dba58189dbbc,
    0x3956c25bf348b538,0x59f111f1b605d019,0x923f82a4af194f9b,0xab1c5ed5da6d8118,
    0xd807aa98a3030242,0x12835b0145706fbe,0x243185be4ee4b28c,0x550c7dc3d5ffb4e2,
    0x72be5d74f27b896f,0x80deb1fe3b1696b1,0x9bdc06a725c71235,0xc19bf174cf692694,
    0xe49b69c19ef14ad2,0xefbe4786384f25e3,0x0fc19dc68b8cd5b5,0x240ca1cc77ac9c65,
    0x2de92c6f592b0275,0x4a7484aa6ea6e483,0x5cb0a9dcbd41fbd4,0x76f988da831153b5,
    0x983e5152ee66dfab,0xa831c66d2db43210,0xb00327c898fb213f,0xbf597fc7beef0ee4,
    0xc6e00bf33da88fc2,0xd5a79147930aa725,0x06ca6351e003826f,0x142929670a0e6e70,
    0x27b70a8546d22ffc,0x2e1b21385c26c926,0x4d2c6dfc5ac42aed,0x53380d139d95b3df,
    0x650a73548baf63de,0x766a0abb3c77b2a8,0x81c2c92e47edaee6,0x92722c851482353b,
    0xa2bfe8a14cf10364,0xa81a664bbc423001,0xc24b8b70d0f89791,0xc76c51a30654be30,
    0xd192e819d6ef5218,0xd69906245565a910,0xf40e35855771202a,0x106aa07032bbd1b8,
    0x19a4c116b8d2d0c8,0x1e376c085141ab53,0x2748774cdf8eeb99,0x34b0bcb5e19b48a8,
    0x391c0cb3c5c95a63,0x4ed8aa4ae3418acb,0x5b9cca4f7763e373,0x682e6ff3d6b2b8a3,
    0x748f82ee5defb2fc,0x78a5636f43172f60,0x84c87814a1f0ab72,0x8cc702081a6439ec,
    0x90befffa23631e28,0xa4506cebde82bde9,0xbef9a3f7b2c67915,0xc67178f2e372532b,
    0xca273eceea26619c,0xd186b8c721c0c207,0xeada7dd6cde0eb1e,0xf57d4f7fee6ed178,
    0x06f067aa72176fba,0x0a637dc5a2c898a6,0x113f9804bef90dae,0x1b710b35131c471b,
    0x28db77f523047d84,0x32caab7b40c72493,0x3c9ebe0a15c9bebc,0x431d67c49c100d4c,
    0x4cc5d4becb3e42b6,0x597f299cfc657e2a,0x5fcb6fab3ad6faec,0x6c44198c4a475817
};

#define ROTR(x,n) (((x) >> (n)) | ((x) << (64-(n))))
#define CH(x,y,z) (((x)&(y)) ^ (~(x)&(z)))
#define MAJ(x,y,z) (((x)&(y)) ^ ((x)&(z)) ^ ((y)&(z)))
#define EP0(x) (ROTR(x,28)^ROTR(x,34)^ROTR(x,39))
#define EP1(x) (ROTR(x,14)^ROTR(x,18)^ROTR(x,41))
#define SIG0(x) (ROTR(x,1)^ROTR(x,8)^((x)>>7))
#define SIG1(x) (ROTR(x,19)^ROTR(x,61)^((x)>>6))

static void transform(SHA512_CTX* ctx, const uint8_t data[SHA512_BLOCK_SIZE]) {
    uint64_t W[80];
    register uint64_t a, b, c, d, e, f, g, h;
    register uint64_t t1, t2;

    // --- 1. 초기 16개 워드 (Big-endian → uint64_t 변환)
    const uint64_t* p = (const uint64_t*)data;
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    for (int i = 0; i < 16; i++) {
        uint64_t v = p[i];
        W[i] = ((v & 0x00000000000000FF) << 56) |
            ((v & 0x000000000000FF00) << 40) |
            ((v & 0x0000000000FF0000) << 24) |
            ((v & 0x00000000FF000000) << 8) |
            ((v & 0x000000FF00000000) >> 8) |
            ((v & 0x0000FF0000000000) >> 24) |
            ((v & 0x00FF000000000000) >> 40) |
            ((v & 0xFF00000000000000) >> 56);
    }
#else
    memcpy(W, data, 128);
#endif

    // --- 2. 나머지 64개 워드 계산 (Message Schedule)
    for (int i = 16; i < 80; i++) {
        uint64_t s0 = SIG0(W[i - 15]);
        uint64_t s1 = SIG1(W[i - 2]);
        W[i] = W[i - 16] + s0 + W[i - 7] + s1;
    }

    // --- 3. 초기 해시 상태 로드
    a = ctx->state[0];
    b = ctx->state[1];
    c = ctx->state[2];
    d = ctx->state[3];
    e = ctx->state[4];
    f = ctx->state[5];
    g = ctx->state[6];
    h = ctx->state[7];

    // --- 4. 메인 압축 루프 (80라운드)
    for (int i = 0; i < 80; i++) {
        t1 = h + EP1(e) + CH(e, f, g) + K[i] + W[i];
        t2 = EP0(a) + MAJ(a, b, c);

        // 레지스터 회전 (의미 그대로 유지)
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }

    // --- 5. 중간 해시 상태 업데이트
    ctx->state[0] += a;
    ctx->state[1] += b;
    ctx->state[2] += c;
    ctx->state[3] += d;
    ctx->state[4] += e;
    ctx->state[5] += f;
    ctx->state[6] += g;
    ctx->state[7] += h;
}

CRYPTO_STATUS sha512_init(SHA512_CTX* ctx) { // 초기 해시값 H(0) 설정
    ctx->state[0] = 0x6a09e667f3bcc908; ctx->state[1] = 0xbb67ae8584caa73b;
    ctx->state[2] = 0x3c6ef372fe94f82b; ctx->state[3] = 0xa54ff53a5f1d36f1;
    ctx->state[4] = 0x510e527fade682d1; ctx->state[5] = 0x9b05688c2b3e6c1f;
    ctx->state[6] = 0x1f83d9abfb41bd6b; ctx->state[7] = 0x5be0cd19137e2179;
    ctx->datalen = 0; ctx->bitlen_high = ctx->bitlen_low = 0;

    return CRYPTO_SUCCESS;
}

static void add_bitlen(SHA512_CTX* ctx, size_t len) {
    ctx->bitlen_low += (uint64_t)len * 8;
    if (ctx->bitlen_low < (uint64_t)len * 8) ctx->bitlen_high++;
}

CRYPTO_STATUS sha512_update(SHA512_CTX* ctx, const uint8_t* data, size_t len) {
    size_t i = 0;

    // 남은 데이터가 버퍼에 이미 일부 있는 경우 먼저 채움
    if (ctx->datalen > 0) {
        size_t fill = SHA512_BLOCK_SIZE - ctx->datalen;
        if (len < fill) {
            memcpy(ctx->buffer + ctx->datalen, data, len);
            ctx->datalen += len;
            return;
        }
        memcpy(ctx->buffer + ctx->datalen, data, fill);

        // 버퍼가 SHA512_BLOCK_SIZE(128 바이트 = 1024비트)만큼 차면 transform()을 호출해 한 블록을 처리
        transform(ctx, ctx->buffer);
        add_bitlen(ctx, SHA512_BLOCK_SIZE);
        data += fill;
        len -= fill;
        ctx->datalen = 0;
    }

    // 이제 남은 입력 데이터를 128바이트(=SHA512_BLOCK_SIZE) 단위로 처리
    while (len >= SHA512_BLOCK_SIZE) {
        transform(ctx, data);
        add_bitlen(ctx, SHA512_BLOCK_SIZE);
        data += SHA512_BLOCK_SIZE;
        len -= SHA512_BLOCK_SIZE;
    }

    // 마지막으로 남은 (<128바이트) 부분을 버퍼에 복사
    if (len > 0) {
        memcpy(ctx->buffer, data, len);
        ctx->datalen = len;
    }
    return CRYPTO_SUCCESS;
}

CRYPTO_STATUS sha512_final(SHA512_CTX* ctx, uint8_t* hash) {
    // i = 현재 버퍼에 남아 있는 데이터 바이트 수
    size_t i = ctx->datalen;

    // 남은 메시지 바이트 수를 비트 길이에 더함 (패딩 바이트는 포함하지 않음)
    if (i > 0) add_bitlen(ctx, i);

    // 메시지 끝에 '1' 비트를 추가 (0x80)하고 필요한 만큼 0으로 패딩할 준비
    ctx->buffer[i++] = 0x80;

    // 남은 공간이 16바이트(길이 필드)보다 적으면, 현재 블록을 0으로 채워서 처리하고 새 블록에서 이어서 작업
    if (i > 112) {
        // 남은 부분을 0으로 패딩
        if (i < SHA512_BLOCK_SIZE) memset(ctx->buffer + i, 0, SHA512_BLOCK_SIZE - i);
        transform(ctx, ctx->buffer);
        // 새 블록 시작
        i = 0;
    }

    // 블록의 112바이트(=128-16) 지점까지 0으로 패딩
    if (i < 112) memset(ctx->buffer + i, 0, 112 - i);

    // 마지막 16바이트(128비트)에 메시지 전체 길이를 big-endian으로 기록
    // bitlen_high와 bitlen_low는 각각 64비트이며, 상위 64비트가 먼저 옴
    for (int j = 0; j < 8; ++j) {
        ctx->buffer[112 + j] = (uint8_t)(ctx->bitlen_high >> (56 - 8 * j));
        ctx->buffer[120 + j] = (uint8_t)(ctx->bitlen_low >> (56 - 8 * j));
    }

    // 마지막 블록을 처리
    transform(ctx, ctx->buffer);

    // 내부 상태(state[8])를 big-endian 바이트 배열(64바이트)로 변환하여 해시 결과에 저장
    for (int j = 0; j < 8; ++j) {
        uint64_t v = ctx->state[j];
        hash[j * 8 + 0] = (uint8_t)(v >> 56);
        hash[j * 8 + 1] = (uint8_t)(v >> 48);
        hash[j * 8 + 2] = (uint8_t)(v >> 40);
        hash[j * 8 + 3] = (uint8_t)(v >> 32);
        hash[j * 8 + 4] = (uint8_t)(v >> 24);
        hash[j * 8 + 5] = (uint8_t)(v >> 16);
        hash[j * 8 + 6] = (uint8_t)(v >> 8);
        hash[j * 8 + 7] = (uint8_t)(v >> 0);
    }

    // 민감한 내부 데이터(버퍼, 상태 등)를 0으로 초기화하여 보안상 잔류 데이터 방지
    memset(ctx->buffer, 0, sizeof(ctx->buffer));
    ctx->datalen = 0;
    ctx->bitlen_low = ctx->bitlen_high = 0;
    for (int j = 0; j < 8; ++j) ctx->state[j] = 0;

    return CRYPTO_SUCCESS;
}

// 16진 문자열 → 바이트 배열
static int hex_to_bytes(const char* hex, unsigned char* out) {
    int len = strlen(hex);
    int i;
    for (i = 0; i < len / 2; i++) {
        sscanf(hex + 2 * i, "%2hhx", &out[i]);
    }
    return i;
}

#ifdef _WIN32
#include <windows.h>
double now_sec(void) {
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)freq.QuadPart;
}
#else
#include <sys/time.h>
double now_sec(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}
#endif

// 파일에서 테스트 벡터를 읽고 실행
// 실제 사용할 때는 입력에 대한 해시을 출력하도록 수정
int test_sha512(void) {
    FILE* fp = fopen("SHA512ShortMsg.rsp", "r");
    if (!fp) {
        perror("파일 열기 실패");
        return 1;
    }

    char line[1024];
    int len = 0;
    char msg_hex[1024] = { 0 };
    char md_ref[200] = { 0 };
    unsigned char msg_bytes[512];
    unsigned char digest[64];
    char digest_hex[129];
    int test_count = 0, pass_count = 0;

    double start = now_sec(); // 시간 측정 시작

    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "Len = ", 6) == 0) {
            len = atoi(line + 6);
        }
        else if (strncmp(line, "Msg = ", 6) == 0) {
            sscanf(line + 6, "%s", msg_hex);
        }
        else if (strncmp(line, "MD = ", 5) == 0) {
            sscanf(line + 5, "%s", md_ref);

            int msg_len = len / 8; // Len은 비트 단위
            hex_to_bytes(msg_hex, msg_bytes);

            // 해시 계산
            SHA512_CTX ctx;
            sha512_init(&ctx);
            sha512_update(&ctx, msg_bytes, msg_len);
            sha512_final(&ctx, digest);

            // 결과를 16진 문자열로 변환
            for (int i = 0; i < 64; i++)
                sprintf(digest_hex + i * 2, "%02x", digest[i]);
            digest_hex[128] = '\0';

            // 비교
            test_count++;
            if (strcmp(digest_hex, md_ref) == 0) {
                pass_count++;
            }
            else {
                printf("Test %d 실패\n", test_count);
                printf("Len=%d\nMsg=%s\n예상=%s\n결과=%s\n\n",
                    len, msg_hex, md_ref, digest_hex);
            }
        }
    }

    double end = now_sec();  // 시간 측정 종료
    double elapsed = end - start;

    fclose(fp);
    printf("테스트 완료: %d 중 %d 통과 \n", test_count, pass_count);
    printf("총 실행 시간: %.6f초 (%.3f ms)\n", elapsed, elapsed * 1000.0);

    if (test_count == pass_count)
        printf("모든 SHA-512 테스트 벡터 통과!\n");

    return 0;
}
