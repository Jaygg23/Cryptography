#include <stdint.h>
#include "params.h"
#include "reduce.h"

/*************************************************
* Name:        montgomery_reduce
*
* Description: For finite field element a with -2^{31}Q <= a <= Q*2^31,
*              compute r \equiv a*2^{-32} (mod Q) such that -Q < r < Q.
*
* Arguments:   - int64_t: finite field element a
*
* Returns r.
**************************************************/
int32_t montgomery_reduce(int64_t a) { // 나눗셈 없이 mod Q 연산을 빠르게 계산하는 알고리즘(나눗셈 대신 곱셈과 시프트만으로)
  int32_t t;

  t = (int64_t)(int32_t)a*QINV;
  t = (a - (int64_t)t*Q) >> 32;
  return t;
}

/*************************************************
* Name:        reduce32
*
* Description: For finite field element a with a <= 2^{31} - 2^{22} - 1,
*              compute r \equiv a (mod Q) such that -6283008 <= r <= 6283008.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t reduce32(int32_t a) {
  int32_t t;

  t = (a + (1 << 22)) >> 23;
  t = a - t*Q;
  return t;
}

/*************************************************
* Name:        caddq
*
* Description: Add Q if input coefficient is negative.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t caddq(int32_t a) { // 조건부로 +Q 를 더해서 음수를 양수 대표로 바꾸기
// a가 음수일 때 Q를 더해 a를 [0, Q-1] 범위의 양수 대표로 만듦
  a += (a >> 31) & Q; 
  // (a >> 31) : 32비트 산술 시프트(부호 비트 유지)
  //	만약 a가 음수이면(a >> 31)는 - 1(0xFFFFFFFF)
  //	a가 0 이상이면(a >> 31)는 0
  // (a >> 31) & Q : 위 결과와 Q를 AND 처리
  //	음수이면 0xFFFFFFFF & Q = Q
  //	비음수이면 0 & Q = 0
  // a가 음수일 경우 a+=Q, 음수가 아닐 경우 변화 없음
  return a;
}

/*************************************************
* Name:        freeze
*
* Description: For finite field element a, compute standard
*              representative r = a mod^+ Q.
*
* Arguments:   - int32_t: finite field element a
*
* Returns r.
**************************************************/
int32_t freeze(int32_t a) {
  a = reduce32(a);
  a = caddq(a);
  return a;
}
