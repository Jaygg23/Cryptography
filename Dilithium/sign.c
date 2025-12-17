#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"

/*************************************************
* Name:        crypto_sign_keypair
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
  uint8_t seedbuf[2*SEEDBYTES + CRHBYTES]; // 시드 확장용 임시 버퍼 (rho, rhoprime, K를 뽑음)
  // [2 * 32바이트 + 64바이트] = [2 * 256비트 + 512비트]
  // 한 번의 SHAKE-256 호출로 ρ, ρ′, K를 모두 생성 -> rhoprime 영역을 CRH 입력(64B)용으로 재활용하기 위함
  uint8_t tr[TRBYTES]; // 개인키의 성분 중 하나인 tr로, 해시 H(ρ||t1) 결과를 담는 버퍼
  // tr[64바이트] = tr[512비트]
  // tr은 384비트이지만 SHAKE-256 출력 버퍼를 64바이트 단위로 맞춤
  const uint8_t *rho, *rhoprime, *key; // seedbuf 안을 가리킴
  polyvecl mat[K]; // 길이 L인 다항식 벡터를 K개 가진 배열 -> 행렬 A의 각 행(벡터)
  polyvecl s1, s1hat; // 비밀키 벡터, NTT 표현
  polyveck s2, t1, t0;

  /* Get randomness for rho, rhoprime and key */
  // pseudo 코드 1번에 해당 -> 랜덤 비트 문자열 난수 시드(제타) 생성
  randombytes(seedbuf, SEEDBYTES); // seedbuf 버퍼에 SEEDBYTES(32바이트 = 256비트)의 난수를 채워 넣음
  seedbuf[SEEDBYTES+0] = K; // 난수 자체에 파라미터 정보를 덧붙임
  seedbuf[SEEDBYTES+1] = L; 
  // seedbuf = [32바이트 난수 | 1바이트 K | 1바이트 L] 와 같이 모드 구분 정보를 시드에 덧붙여서 SHAKE-256 입력을 유일하게 만듦

  // pseudo 코드 2번에 해당
  shake256(seedbuf, 2*SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES+2);
  // seedbuf의 앞쪽 SEEDBYTES + 2바이트(난수 시드+파라미터 정보 K, L)를 입력값으로 해시하여, 
  // 그 결과로 2 * SEEDBYTES + CRHBYTES 바이트를 다시 seedbuf에 출력 -> rho: SEEDBYTES, rhoprime: CRHBYTES, key: SEEDBYTES
  rho = seedbuf; 
  rhoprime = rho + SEEDBYTES;
  key = rhoprime + CRHBYTES;
  // seedbuf = | < -------- - rho(SEEDBYTES) ----------->|| < ------rhoprime(CRHBYTES) ------->|| < -- - key(SEEDBYTES) ---> |

  /* Expand matrix */
  // pseudo 코드 4번에 해당
  polyvec_matrix_expand(mat, rho); // rho를 시드로 사용해서 크기 k×l의 행렬 A를 결정적으로 만듦

  /* Sample short vectors s1 and s2 */
  // pseudo 코드 3번에 해당
  polyvecl_uniform_eta(&s1, rhoprime, 0); // 비밀키 벡터 s1 생성
  polyveck_uniform_eta(&s2, rhoprime, L); // 비밀키 벡터 s2 생성
  // polyvecl_uniform_eta(): η(eta) 분포를 따르는 다항식 벡터를 만드는 함수

  /* Matrix-vector multiplication */
  // NTT를 이용한 행렬-벡터 곱(A * s1)
  s1hat = s1; // 비밀 벡터 s1을 s1hat에 복사
  polyvecl_ntt(&s1hat); // s1hat의 각 다항식에 대해 NTT 변환을 수행
  polyvec_matrix_pointwise_montgomery(&t1, mat, &s1hat); // 실제 행렬–벡터 곱셈
  // mat: k x l 크기의 공개 행렬 A, s1hat: 길이 l의 비밀 벡터 s1, t1: 길이 k의 결과 벡터
  // s1hat이 NTT 변환되어 있으므로, 곱셈은 각 계수별로 수행 -> 다항식 곱셈이 단순한 계수 곱으로 바뀜 (Montgomery 곱셈은 모듈러 연산을 빠르게 하기 위한 기법)
  polyveck_reduce(&t1); // NTT 계산과 몽고메리 곱셈 후, 각 계수가 Q로 나눈 나머지로 정리되지 않았을 수 있으므로 reduce()는 모든 계수를 모듈러 Q 범위 안으로 조정
  polyveck_invntt_tomont(&t1); // 역NTT를 수행하여 t1을 시간 영역으로 돌림

  /* Add error vector s2 */
  // pseudo 5번에 해당
  polyveck_add(&t1, &t1, &s2); // As1 + s2 mod q

  /* Extract t1 and write public key */
  // pseudo 코드 6번에 해당
  polyveck_caddq(&t1); // 모든 계수를 [0,q) 범위로 정규화 -> t = t mod q
  polyveck_power2round(&t1, &t0, &t1); // t를 상위비트(t1)와 하위비트(t0)로 분리

  // pseudo 코드 8번에 해당
  pack_pk(pk, rho, &t1); // 공개키를 직렬화하는 단계

  /* Compute H(rho, t1) and write secret key */
  // pseudo 코드 7번에 해당
  shake256(tr, TRBYTES, pk, CRYPTO_PUBLICKEYBYTES); // 공개키 pk 전체에 대해 SHAKE256 해시를 계산하여 그 결과 tr를 저장
  // CRYPTO_PUBLICKEYBYTES = (SEEDBYTES + K*POLYT1_PACKEDBYTES) = (32 + 4*320) 
  // rho = 32바이트, POLYT1_PACKEDBYTES = 320: 다항식 하나를 압축한 바이트 수 => 총 1312바이트
  
  // pseudo 코드 8번에 해당
  pack_sk(sk, rho, tr, key, &t0, &s1, &s2);

  return 0;
}

/*************************************************
* Name:        crypto_sign_signature_internal
*
* Description: Computes signature. Internal API.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *pre:   pointer to prefix string
*              - size_t prelen:  length of prefix string
*              - uint8_t *rnd:   pointer to random seed
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign_signature_internal(uint8_t *sig,
                                   size_t *siglen,
                                   const uint8_t *m,
                                   size_t mlen,
                                   const uint8_t *pre,
                                   size_t prelen,
                                   const uint8_t rnd[RNDBYTES],
                                   const uint8_t *sk)
{
  unsigned int n;
  uint8_t seedbuf[2*SEEDBYTES + TRBYTES + 2*CRHBYTES]; // 여러 시드/해시값을 한 번에 담는 큰 버퍼(2*32 + 48 + 2*64)
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0; // y 벡터 샘플링할 때, 같은 rhoprime로 여러 번 y를 만들기 위해 쓰는 카운터
  polyvecl mat[K], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp; // 챌린지  다항식 c
  keccak_state state;

  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + TRBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;
  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

  /* Compute mu = CRH(tr, pre, msg) */
  // 10. CRH에 tr || M을 시드로 받아, mu 생성
  shake256_init(&state); // Keccak 스폰지 상태(1600비트)를 0으로 초기화하고, 내부 포인터 pos를 0으로 함
  shake256_absorb(&state, tr, TRBYTES); // 입력을 차례대로 스폰지에 흡수 (tr || pre || m)
  shake256_absorb(&state, pre, prelen); 
  shake256_absorb(&state, m, mlen);
  // 중간에 rate가 꽉 차면 Keccak - f[1600] 순열을 한 번 돌리고 계속 이어서 흡수
  shake256_finalize(&state); // 더 이상 absorb를 하지 않는다는 표시 + 패딩 비트를 붙이고 스폰지를 squeeze 모드로 전환
  shake256_squeeze(mu, CRHBYTES, &state); // 스폰지에서 원하는 길이만큼 출력을 뽑아냄. CRHBYTES만큼 뽑아서 mu에 저장

  /* Compute rhoprime = CRH(key, rnd, mu) */
  // 12. rhoprime을 CRH( K || mu )로 생성
  shake256_init(&state);
  shake256_absorb(&state, key, SEEDBYTES);
  shake256_absorb(&state, rnd, RNDBYTES); // #ifdef DILITHIUM_RANDOMIZED_SIGNING가 켜져 있으면 randombytes(rnd, 32);로 진짜 랜덤을 채우고,
                                          // 꺼져 있으면 memset(rnd, 0, 32);로 전부 0을 채워 넣음
  shake256_absorb(&state, mu, CRHBYTES);
  // 의사코드 상에서는 rhoprime을 K를 사용한 다이제스트와 그냥 랜덤값으로 나눴지만,
  // 실제 구현 상에서는 rhoprime = CRH( K || rnd || mu ) 형태로 rnd 값에 0 또는 실제 랜덤값이 들어가게 되며 케이스가 나뉨
  shake256_finalize(&state);
  shake256_squeeze(rhoprime, CRHBYTES, &state);

  /* Expand matrix and transform vectors */
  // 9. rho를 시드로 행렬 A를 결정론적으로 생성
  polyvec_matrix_expand(mat, rho); 
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);

rej: // 13. (z, h)가 NULL일 경우 루프
  /* Sample intermediate vector y */
  // 14. 중간 벡터 y 샘플링
  polyvecl_uniform_gamma1(&y, rhoprime, nonce++);

  /* Matrix-vector multiplication */
  // 15. NTT 변환을 이용해 w를 구함
  z = y;
  polyvecl_ntt(&z); // y -> NTT(y)
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z); // w1= A * y
  polyveck_reduce(&w1); // mod q
  polyveck_invntt_tomont(&w1); // 기존 도메인으로 변경 -> w (메모리 절약을 위해 w1에 저장 후 덮어쓰기)

  /* Decompose w and call the random oracle */
  // 16. w의 상위 성분 w1 구함
  polyveck_caddq(&w1); // 감산(mod q) 과정에서 생길 수 있는 음수 계수를 q를 더해서 전부 0~(q-1) 범위로 만들어 줌
  polyveck_decompose(&w1, &w0, &w1); //다시 더 봐라!! w를 decompose 함수를 사용해서 상위비트(w1)와 하위비트(w0)로 나눔
  polyveck_pack_w1(sig, &w1); // 서명(signature)에 저장될 w1 값을 바이트형식으로 압축

  // 17. H(mu || w1)를  ~c에 저장
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES); // mu라는 바이트 배열을 CRHBYTES 길이만큼 흡수
  shake256_absorb(&state, sig, K*POLYW1_PACKEDBYTES); // w1 값 흡수
  shake256_finalize(&state); // c = H(mu || w1)
  shake256_squeeze(sig, CTILDEBYTES, &state); // CTILDEBYTES 만큼을 sig 버퍼에 저장
  poly_challenge(&cp, sig); // 바이트열 sig(~c)를 사용해서 희소 다항식 cp(c) 생성
  poly_ntt(&cp);

  /* Compute z, reject if it reveals secret */
  // 19. 서명 후보 z = y + c * s1 계산
  polyvecl_pointwise_poly_montgomery(&z, &cp, &s1); // cp * s1을 z에 저장
  polyvecl_invntt_tomont(&z); // 역ntt 변환
  polyvecl_add(&z, &z, &y); // z = (cp * s1) + y
  polyvecl_reduce(&z);
  if(polyvecl_chknorm(&z, GAMMA1 - BETA)) // 거절 조건 : z의 어떤 계수가 GAMMA1 - BETA보다 크면 거절
    goto rej;

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  // 21. 거절 조건 확인
  // (c * s2)를 반영했을 때도 hint가 안전하게 생성되는지 검사하는 reject test
  polyveck_pointwise_poly_montgomery(&h, &cp, &s2); // c * s2
  polyveck_invntt_tomont(&h);
  polyveck_sub(&w0, &w0, &h); // w = w - (c * s2)
  polyveck_reduce(&w0);
  if(polyveck_chknorm(&w0, GAMMA2 - BETA)) // 벡터 w의 어떤 계수라도 주어진 범위(GAMMA2 - BETA) 이상이면 1(true)을 반환 -> 거절
    goto rej;

  /* Compute hints for w1 */
  polyveck_pointwise_poly_montgomery(&h, &cp, &t0); // h = c * t0
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
  if(polyveck_chknorm(&h, GAMMA2)) // c * t0의 어떤 계수라도 GAMMA2 이상이면 거절
    goto rej;

  // 힌트 h 생성
  polyveck_add(&w0, &w0, &h); // w = w0 + (c * t0) 
  n = polyveck_make_hint(&h, &w0, &w1); // 힌트 h 생성 -> h: 힌트 벡터(0 / 1비트), n: 힌트에서 1의 개수
  if(n > OMEGA) // n이 OMEGA보다 크면 거절
    goto rej;

  /* Write signature */
  // 26. 서명 생성
  pack_sig(sig, sig, &z, &h); // 서명 생성
  *siglen = CRYPTO_BYTES; // 서명 길이
  return 0; // 서명 생성 성공
}

/*************************************************
* Name:        crypto_sign_signature
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *ctx:   pointer to contex string
*              - size_t ctxlen:  length of contex string
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success) or -1 (context string too long)
**************************************************/
int crypto_sign_signature(uint8_t *sig, // 서명 결과가 저장될 버퍼
                          size_t *siglen, // 서명 길이
                          const uint8_t *m, // 메시지
                          size_t mlen, // 메시지 길이
                          const uint8_t *ctx, 
                          size_t ctxlen,
                          const uint8_t *sk) // 비밀키
{
  size_t i;
  uint8_t pre[257]; // 서명용 prefix 버퍼 (0x00, ctxlen, ctx 자체)
  uint8_t rnd[RNDBYTES]; // 서명 시 사용하는 랜덤 바이트 버퍼 (RNDBYTES: 32바이트)

  if(ctxlen > 255)
    return -1;

  /* Prepare pre = (0, ctxlen, ctx) */
  // 메시지 앞에 붙여줄 prefix 
  pre[0] = 0; // 같은 해시 함수를 여러 곳에서 쓸 때 서명용 입력/다른 용도임을 구분하기 위한 도메인 태그 (0x00: 서명용 입력)
  pre[1] = ctxlen; // "dilithium"과 같이 쓰이는 짧은 ASCII 문자열 -> 이 서명은 이 용도에서만 유효함을 나타냄
  for(i = 0; i < ctxlen; i++)
    pre[2 + i] = ctx[i]; // ctx를 prefix(pre[]) 버퍼 안에 복사

 // 서명 생성 시 랜덤성을 쓸지 결정론적 모드를 쓸지 결정
#ifdef DILITHIUM_RANDOMIZED_SIGNING // rhoprime을 랜덤으로 생성하는 randomized signing 경우
  randombytes(rnd, RNDBYTES);
#else
  for(i=0;i<RNDBYTES;i++)
    rnd[i] = 0;
#endif

  crypto_sign_signature_internal(sig,siglen,m,mlen,pre,2+ctxlen,rnd,sk); // 실제 서명 생성 알고리즘
  return 0;
}

/*************************************************
* Name:        crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *ctx: pointer to context string
*              - size_t ctxlen: length of context string
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success) or -1 (context string too long)
**************************************************/
int crypto_sign(uint8_t *sm, // 결과를 담을 버퍼 (signed message)
                size_t *smlen, // 결과 전체 길이를 저장할 포인터
                const uint8_t *m, // 서명할 메시지 (입력 메시지)
                size_t mlen, // 메시지 길이
                const uint8_t *ctx,
                size_t ctxlen,
                const uint8_t *sk) // 비밀키
{
  int ret;
  size_t i;

  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i]; // 메시지 m을 sm의 버퍼 뒤쪽(CRYPTO_BYTES 이후 위치)에 역순으로 복사 (CRYPTO_BYTES는 서명의 길이)
    // 서명을 붙이기 전에 메시지를 안전하게 옮겨두기 위해 뒤에서부터 채움
  ret = crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, ctx, ctxlen, sk); // 서명 생성
  // sm + CRYPTO_BYTES: 메시지가 저장된 위치
  // sm + CRYPTO_BYTES 위치에 있는 메시지를 읽어서 서명을 생성한 뒤 그 서명을 sm[0..CRYPTO_BYTES-1]에 저장

  *smlen += mlen; // 서명 길이 + 메시지 길이
  return ret; // crypto_sign_signature() 함수의 반환값을 그대로 전달 -> 서명 생성이 성공했는지를 호출자에게 알려주는 용도
}

/*************************************************
* Name:        crypto_sign_verify_internal
*
* Description: Verifies signature. Internal API.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pre: pointer to prefix string
*              - size_t prelen: length of prefix string
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_verify_internal(const uint8_t *sig,
                                size_t siglen,
                                const uint8_t *m,
                                size_t mlen,
                                const uint8_t *pre,
                                size_t prelen,
                                const uint8_t *pk)
{
  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[CTILDEBYTES];
  uint8_t c2[CTILDEBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;
  keccak_state state;

  if(siglen != CRYPTO_BYTES) // 서명 길이가 다를 경우 거절
    return -1;

  unpack_pk(rho, &t1, pk);
  if(unpack_sig(c, &z, &h, sig)) // 서명 sig를 풀어서 형식이 이상하면 거절
    return -1;
  if(polyvecl_chknorm(&z, GAMMA1 - BETA)) // z가 너무 크면 실패
    return -1;

  /* Compute CRH(H(rho, t1), pre, msg) */
  shake256(mu, TRBYTES, pk, CRYPTO_PUBLICKEYBYTES); // tr 생성 (변수명은 mu이지만 tr이 들어감): tr = CRH(pk)
  shake256_init(&state);
  shake256_absorb(&state, mu, TRBYTES);
  shake256_absorb(&state, pre, prelen);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state); // 28. mu = CRH(tr || pre || m)

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  poly_challenge(&cp, c); // 29. SampleInBall (해시 결과인 ~c(c)를 희소 다항식 c(cp)로 만듦)
  polyvec_matrix_expand(mat, rho); // 27. 행렬 A 생성

  polyvecl_ntt(&z); // z를 NTT 도메인으로 바꾼 후 
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z); // A * z

  poly_ntt(&cp);
  polyveck_shiftl(&t1); // t1에 2^d를 곱함 (계수에 2^d배, 왼쪽 시프트)
  polyveck_ntt(&t1);
  polyveck_pointwise_poly_montgomery(&t1, &cp, &t1); // t1 = c * t1 * 2^d

  polyveck_sub(&w1, &w1, &t1); // Az - (c * t1 * 2^d)
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);

  /* Reconstruct w1 */
  polyveck_caddq(&w1);
  polyveck_use_hint(&w1, &w1, &h); // 30. UseHint로 서명 생성과 동일한 w1 재생성
  polyveck_pack_w1(buf, &w1);

  /* Call random oracle and verify challenge */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES); // 서명 생성과 동일하게 CRH(mu || w1) 연산
  shake256_absorb(&state, buf, K*POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(c2, CTILDEBYTES, &state); // 새로운 챌린지 다항식 c2 생성
  for(i = 0; i < CTILDEBYTES; ++i)
    if(c[i] != c2[i]) // 서명에서 unpack한 챌린지 c와 동일한지 확인
      return -1; // 다르면 서명 위조

  return 0;
}

/*************************************************
* Name:        crypto_sign_verify
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *ctx: pointer to context string
*              - size_t ctxlen: length of context string
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_verify(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *ctx,
                       size_t ctxlen,
                       const uint8_t *pk)
{
  size_t i;
  uint8_t pre[257];

  if(ctxlen > 255)
    return -1;

  // pre 포맷: pre = 0x00 || ctxlen || ctx
  pre[0] = 0;
  pre[1] = ctxlen;
  for(i = 0; i < ctxlen; i++)
    pre[2 + i] = ctx[i]; 

  return crypto_sign_verify_internal(sig,siglen,m,mlen,pre,2+ctxlen,pk);
}

/*************************************************
* Name:        crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *ctx: pointer to context tring
*              - size_t ctxlen: length of context string
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *ctx,
                     size_t ctxlen,
                     const uint8_t *pk)
{
  size_t i;

  if(smlen < CRYPTO_BYTES) // 서명은 최소 CRYPTO_BYTES 바이트
    goto badsig;

  *mlen = smlen - CRYPTO_BYTES; // sm = [ signature | message ]
  if(crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, ctx, ctxlen, pk)) // 서명 검증
    goto badsig;
  else {
    /* All good, copy msg, return 0 */
    for(i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i]; // 검증 성공 시, 서명 뒤의 메시지 부분을 출력 버퍼 m으로 복사
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = 0;
  for(i = 0; i < smlen; ++i) 
    m[i] = 0; // 검증 실패 시, 출력 버퍼를 0으로 덮음 (쓰레기 데이터가 남지 않도록 하는 보안 습관)

  return -1;
}