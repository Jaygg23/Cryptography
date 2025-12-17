#include <stdint.h>
#include "params.h"
#include "rounding.h"

/*************************************************
* Name:        power2round
*
* Description: For finite field element a, compute a0, a1 such that
*              a mod^+ Q = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
*              Assumes a to be standard representative.
*
* Arguments:   - int32_t a: input element
*              - int32_t *a0: pointer to output element a0
*
* Returns a1.
**************************************************/
int32_t power2round(int32_t *a0, int32_t a)  {
  int32_t a1;

  a1 = (a + (1 << (D-1)) - 1) >> D; // a1 = (a+2^{D-1}=1)/2^D의 floor
  *a0 = a - (a1 << D); // a0 = a-a1*2^D
  // a를 2^D로 나눈 몫과 나머지로 분리
  // + (1 << (D-1)) - 1의 목적: a1이 a / 2^D 를 단순한 내림(floor)으로 가져가면 a0의 절댓값 범위가 커질 수 있는데,
  // 적당한 바이어스를 주어 a0가 작은 절대값 범위(대칭적) 내에 들어가게끔 조정하기 위함 => a0는 절댓값 기준 작은 값이 되어 비밀키 저장에 유리
  return a1;
}

/*************************************************
* Name:        decompose
*
* Description: For finite field element a, compute high and low bits a0, a1 such
*              that a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
*              if a1 = (Q-1)/ALPHA where we set a1 = 0 and
*              -ALPHA/2 <= a0 = a mod^+ Q - Q < 0. Assumes a to be standard
*              representative.
*
* Arguments:   - int32_t a: input element
*              - int32_t *a0: pointer to output element a0
*
* Returns a1.
**************************************************/
int32_t decompose(int32_t *a0, int32_t a) {
  int32_t a1;

  a1  = (a + 127) >> 7; // a / (2*GAMMA2) 비슷한 값을 만들기 위한 고정소수점(fixed-point) 전처리 단계
#if GAMMA2 == (Q-1)/32 // Dilithium 파라미터에 따라 GAMMA2 값이 달라짐
  a1  = (a1*1025 + (1 << 21)) >> 22; // a1에 어떤 비율을 곱해서 다시 나눠서, 원하는 비율의 근사값을 만듦
  a1 &= 15; // a1 값을 0~15로 제한 (하위 4비트만 남김)
#elif GAMMA2 == (Q-1)/88
  a1  = (a1*11275 + (1 << 23)) >> 24;
  a1 ^= ((43 - a1) >> 31) & a1; // 43 - a1이 음수면 a1 ^= a1 -> 0
  // >>31이 0xFFFFFFFF
  // a1 > 43이면 0으로 만듦 (a1을 0~43 범위 안으로 함)
#endif

  *a0  = a - a1*2*GAMMA2; // a에서 high part를 뺀 나머지 -> 하위 비트
  *a0 -= (((Q-1)/2 - *a0) >> 31) & Q; // a0를 (-Q/2, Q/2] 같은 대칭 범위로 옮기기
  return a1; // 상위 비트
}

/*************************************************
* Name:        make_hint
*
* Description: Compute hint bit indicating whether the low bits of the
*              input element overflow into the high bits.
*
* Arguments:   - int32_t a0: low bits of input element
*              - int32_t a1: high bits of input element
*
* Returns 1 if overflow.
**************************************************/
unsigned int make_hint(int32_t a0, int32_t a1) {
  if(a0 > GAMMA2 || a0 < -GAMMA2 || (a0 == -GAMMA2 && a1 != 0))
    return 1;

  return 0;
}

/*************************************************
* Name:        use_hint
*
* Description: Correct high bits according to hint.
*
* Arguments:   - int32_t a: input element
*              - unsigned int hint: hint bit
*
* Returns corrected high bits.
**************************************************/
int32_t use_hint(int32_t a, unsigned int hint) {
  int32_t a0, a1;

  a1 = decompose(&a0, a); // 계수 a를 decompose 해서 a1, a0로 나눔
  if(hint == 0) // 힌트가 0이면 보정 없이 그대로 리턴
    return a1;

#if GAMMA2 == (Q-1)/32
  if(a0 > 0)
    return (a1 + 1) & 15;
  else
    return (a1 - 1) & 15;
#elif GAMMA2 == (Q-1)/88 // a1 범위가 0~43이므로 43에서 +1이면 0으로, 0에서 -1이면 43으로 래핑 처리
  if(a0 > 0)
    return (a1 == 43) ?  0 : a1 + 1;
  else
    return (a1 ==  0) ? 43 : a1 - 1;
#endif
}
