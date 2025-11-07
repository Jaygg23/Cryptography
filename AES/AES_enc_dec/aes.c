#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include <stdint.h>
#include "aes.h"

//학번_이름
char SUBMISSION_INFO[256] = "2023KU0203_송윤우";

/*****************************************************함수 구현********************************************************/

//#define xtimes(f) ((((f) >> 7 & 0x01) == 1) ? ((f) << 1) ^ 0x1b : (f) << 1) // if문 사용시 시간 더 걸림

#define xtimes(input) (((input)<<1)^(((input)>>7)*0x1b)) // if문 대신 함수로 구현

const unsigned char s_box[256] = {
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

const unsigned char Inv_s_box[256] =
{
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d
};

const unsigned char RC[10] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36 };

void AddRoundKey(unsigned char string[], unsigned char key[])
{
    string[0] ^= key[0];
    string[1] ^= key[1];
    string[2] ^= key[2];
    string[3] ^= key[3];
    string[4] ^= key[4];
    string[5] ^= key[5];
    string[6] ^= key[6];
    string[7] ^= key[7];
    string[8] ^= key[8];
    string[9] ^= key[9];
    string[10] ^= key[10];
    string[11] ^= key[11];
    string[12] ^= key[12];
    string[13] ^= key[13];
    string[14] ^= key[14];
    string[15] ^= key[15];
}

void SubBytes(unsigned char string[])
{
    string[0] = s_box[string[0]];
    string[1] = s_box[string[1]];
    string[2] = s_box[string[2]];
    string[3] = s_box[string[3]];
    string[4] = s_box[string[4]];
    string[5] = s_box[string[5]];
    string[6] = s_box[string[6]];
    string[7] = s_box[string[7]];
    string[8] = s_box[string[8]];
    string[9] = s_box[string[9]];
    string[10] = s_box[string[10]];
    string[11] = s_box[string[11]];
    string[12] = s_box[string[12]];
    string[13] = s_box[string[13]];
    string[14] = s_box[string[14]];
    string[15] = s_box[string[15]];
}

void InvSubBytes(unsigned char string[]) {
    string[0] = Inv_s_box[string[0]];
    string[1] = Inv_s_box[string[1]];
    string[2] = Inv_s_box[string[2]];
    string[3] = Inv_s_box[string[3]];
    string[4] = Inv_s_box[string[4]];
    string[5] = Inv_s_box[string[5]];
    string[6] = Inv_s_box[string[6]];
    string[7] = Inv_s_box[string[7]];
    string[8] = Inv_s_box[string[8]];
    string[9] = Inv_s_box[string[9]];
    string[10] = Inv_s_box[string[10]];
    string[11] = Inv_s_box[string[11]];
    string[12] = Inv_s_box[string[12]];
    string[13] = Inv_s_box[string[13]];
    string[14] = Inv_s_box[string[14]];
    string[15] = Inv_s_box[string[15]];
}

void ShiftRows(unsigned char string[])
{
    unsigned char buf;

    buf = string[1];
    string[1] = string[5];
    string[5] = string[9];
    string[9] = string[13];
    string[13] = buf;

    buf = string[2];
    string[2] = string[10];
    string[10] = buf;
    buf = string[6];
    string[6] = string[14];
    string[14] = buf;

    buf = string[15];
    string[15] = string[11];
    string[11] = string[7];
    string[7] = string[3];
    string[3] = buf;
}

void InvShiftRows(unsigned char string[])
{
    unsigned char buf;

    buf = string[13];
    string[13] = string[9];
    string[9] = string[5];
    string[5] = string[1];
    string[1] = buf;

    buf = string[14];
    string[14] = string[6];
    string[6] = buf;
    buf = string[10];
    string[10] = string[2];
    string[2] = buf;

    buf = string[3];
    string[3] = string[7];
    string[7] = string[11];
    string[11] = string[15];
    string[15] = buf;
}

void MixColumns(unsigned char string[])
{
    unsigned char buf[4] = { 0x00, };

    buf[0] = string[0];
    buf[1] = string[1];
    buf[2] = string[2];
    buf[3] = string[3];
    string[0] = xtimes((buf[0] ^ buf[1])) ^ buf[1] ^ buf[2] ^ buf[3];
    string[1] = buf[0] ^ xtimes((buf[1] ^ buf[2])) ^ buf[2] ^ buf[3];
    string[2] = buf[0] ^ buf[1] ^ xtimes((buf[2] ^ buf[3])) ^ buf[3];
    string[3] = xtimes((buf[0] ^ buf[3])) ^ buf[0] ^ buf[1] ^ buf[2];

    buf[0] = string[4];
    buf[1] = string[5];
    buf[2] = string[6];
    buf[3] = string[7];
    string[4] = xtimes((buf[0] ^ buf[1])) ^ buf[1] ^ buf[2] ^ buf[3];
    string[5] = buf[0] ^ xtimes((buf[1] ^ buf[2])) ^ buf[2] ^ buf[3];
    string[6] = buf[0] ^ buf[1] ^ xtimes((buf[2] ^ buf[3])) ^ buf[3];
    string[7] = xtimes((buf[0] ^ buf[3])) ^ buf[0] ^ buf[1] ^ buf[2];

    buf[0] = string[8];
    buf[1] = string[9];
    buf[2] = string[10];
    buf[3] = string[11];
    string[8] = xtimes((buf[0] ^ buf[1])) ^ buf[1] ^ buf[2] ^ buf[3];
    string[9] = buf[0] ^ xtimes((buf[1] ^ buf[2])) ^ buf[2] ^ buf[3];
    string[10] = buf[0] ^ buf[1] ^ xtimes((buf[2] ^ buf[3])) ^ buf[3];
    string[11] = xtimes((buf[0] ^ buf[3])) ^ buf[0] ^ buf[1] ^ buf[2];

    buf[0] = string[12];
    buf[1] = string[13];
    buf[2] = string[14];
    buf[3] = string[15];
    string[12] = xtimes((buf[0] ^ buf[1])) ^ buf[1] ^ buf[2] ^ buf[3];
    string[13] = buf[0] ^ xtimes((buf[1] ^ buf[2])) ^ buf[2] ^ buf[3];
    string[14] = buf[0] ^ buf[1] ^ xtimes((buf[2] ^ buf[3])) ^ buf[3];
    string[15] = xtimes((buf[0] ^ buf[3])) ^ buf[0] ^ buf[1] ^ buf[2];
}

//InvMixColumns의 연산에 사용
unsigned char GF_mul(unsigned char a, unsigned b)
{
    unsigned char h = 0x00;
    int coef; //계수
    for (int i = 7; i >= 0; i--)
    {
        coef = (a >> i) & 0x01;
        h = xtimes(h);
        if (coef == 1)
        {
            h ^= b;
        }
    }
    return h;
}

void InvMixColumns(unsigned char string[]) {
    unsigned char buf[4] = { 0x00, };

    buf[0] = string[0];
    buf[1] = string[1];
    buf[2] = string[2];
    buf[3] = string[3];
    string[0] = GF_mul(0x0e, buf[0]) ^ GF_mul(0x0b, buf[1]) ^ GF_mul(0x0d, buf[2]) ^ GF_mul(0x09, buf[3]);
    string[1] = GF_mul(0x09, buf[0]) ^ GF_mul(0x0e, buf[1]) ^ GF_mul(0x0b, buf[2]) ^ GF_mul(0x0d, buf[3]);
    string[2] = GF_mul(0x0d, buf[0]) ^ GF_mul(0x09, buf[1]) ^ GF_mul(0x0e, buf[2]) ^ GF_mul(0x0b, buf[3]);
    string[3] = GF_mul(0x0b, buf[0]) ^ GF_mul(0x0d, buf[1]) ^ GF_mul(0x09, buf[2]) ^ GF_mul(0x0e, buf[3]);

    buf[0] = string[4];
    buf[1] = string[5];
    buf[2] = string[6];
    buf[3] = string[7];
    string[4] = GF_mul(0x0e, buf[0]) ^ GF_mul(0x0b, buf[1]) ^ GF_mul(0x0d, buf[2]) ^ GF_mul(0x09, buf[3]);
    string[5] = GF_mul(0x09, buf[0]) ^ GF_mul(0x0e, buf[1]) ^ GF_mul(0x0b, buf[2]) ^ GF_mul(0x0d, buf[3]);
    string[6] = GF_mul(0x0d, buf[0]) ^ GF_mul(0x09, buf[1]) ^ GF_mul(0x0e, buf[2]) ^ GF_mul(0x0b, buf[3]);
    string[7] = GF_mul(0x0b, buf[0]) ^ GF_mul(0x0d, buf[1]) ^ GF_mul(0x09, buf[2]) ^ GF_mul(0x0e, buf[3]);

    buf[0] = string[8];
    buf[1] = string[9];
    buf[2] = string[10];
    buf[3] = string[11];
    string[8] = GF_mul(0x0e, buf[0]) ^ GF_mul(0x0b, buf[1]) ^ GF_mul(0x0d, buf[2]) ^ GF_mul(0x09, buf[3]);
    string[9] = GF_mul(0x09, buf[0]) ^ GF_mul(0x0e, buf[1]) ^ GF_mul(0x0b, buf[2]) ^ GF_mul(0x0d, buf[3]);
    string[10] = GF_mul(0x0d, buf[0]) ^ GF_mul(0x09, buf[1]) ^ GF_mul(0x0e, buf[2]) ^ GF_mul(0x0b, buf[3]);
    string[11] = GF_mul(0x0b, buf[0]) ^ GF_mul(0x0d, buf[1]) ^ GF_mul(0x09, buf[2]) ^ GF_mul(0x0e, buf[3]);

    buf[0] = string[12];
    buf[1] = string[13];
    buf[2] = string[14];
    buf[3] = string[15];
    string[12] = GF_mul(0x0e, buf[0]) ^ GF_mul(0x0b, buf[1]) ^ GF_mul(0x0d, buf[2]) ^ GF_mul(0x09, buf[3]);
    string[13] = GF_mul(0x09, buf[0]) ^ GF_mul(0x0e, buf[1]) ^ GF_mul(0x0b, buf[2]) ^ GF_mul(0x0d, buf[3]);
    string[14] = GF_mul(0x0d, buf[0]) ^ GF_mul(0x09, buf[1]) ^ GF_mul(0x0e, buf[2]) ^ GF_mul(0x0b, buf[3]);
    string[15] = GF_mul(0x0b, buf[0]) ^ GF_mul(0x0d, buf[1]) ^ GF_mul(0x09, buf[2]) ^ GF_mul(0x0e, buf[3]);
}

/*****************************************************AES128********************************************************/

//관련 데이터 타입 정의 내역
//128-bit block
//typedef uint8_t AES_STATE_t[16];
//128-bit masterkey
//typedef uint8_t AES128_KEY_t[16];

void KeySchedule128(unsigned char key[][16])
{
    int cnt_i;

    for (cnt_i = 0; cnt_i < 10; cnt_i++)
    {
        key[cnt_i + 1][0] = key[cnt_i][0] ^ s_box[key[cnt_i][13]] ^ RC[cnt_i];
        key[cnt_i + 1][1] = key[cnt_i][1] ^ s_box[key[cnt_i][14]];
        key[cnt_i + 1][2] = key[cnt_i][2] ^ s_box[key[cnt_i][15]];
        key[cnt_i + 1][3] = key[cnt_i][3] ^ s_box[key[cnt_i][12]];

        key[cnt_i + 1][4] = key[cnt_i + 1][0] ^ key[cnt_i][4];
        key[cnt_i + 1][5] = key[cnt_i + 1][1] ^ key[cnt_i][5];
        key[cnt_i + 1][6] = key[cnt_i + 1][2] ^ key[cnt_i][6];
        key[cnt_i + 1][7] = key[cnt_i + 1][3] ^ key[cnt_i][7];

        key[cnt_i + 1][8] = key[cnt_i + 1][4] ^ key[cnt_i][8];
        key[cnt_i + 1][9] = key[cnt_i + 1][5] ^ key[cnt_i][9];
        key[cnt_i + 1][10] = key[cnt_i + 1][6] ^ key[cnt_i][10];
        key[cnt_i + 1][11] = key[cnt_i + 1][7] ^ key[cnt_i][11];

        key[cnt_i + 1][12] = key[cnt_i + 1][8] ^ key[cnt_i][12];
        key[cnt_i + 1][13] = key[cnt_i + 1][9] ^ key[cnt_i][13];
        key[cnt_i + 1][14] = key[cnt_i + 1][10] ^ key[cnt_i][14];
        key[cnt_i + 1][15] = key[cnt_i + 1][11] ^ key[cnt_i][15];
    }
}

void AES128_enc(AES_STATE_t C, AES_STATE_t P, AES128_KEY_t K128)
{
    unsigned char RoundKey[11][16];
    for (int i = 0; i < 16; i++)
    {
        RoundKey[0][i] = K128[i];
        C[i] = P[i];
    }

    KeySchedule128(RoundKey);

    AddRoundKey(C, RoundKey[0]);

    //1~9 라운드
    for (int cnt_i = 0; cnt_i < 9; cnt_i++)
    {
        SubBytes(C);
        ShiftRows(C);
        MixColumns(C);
        AddRoundKey(C, RoundKey[cnt_i + 1]);
    }

    //마지막 10 라운드
    SubBytes(C);
    ShiftRows(C);
    AddRoundKey(C, RoundKey[10]);
}


void AES128_dec(AES_STATE_t P, AES_STATE_t C, AES128_KEY_t K128)
{
    unsigned char RoundKey[11][16];
    for (int i = 0; i < 16; i++)
    {
        RoundKey[0][i] = K128[i];
        P[i] = C[i];
    }

    KeySchedule128(RoundKey);

    AddRoundKey(P, RoundKey[10]);

    //1~9 라운드
    for (int cnt_i = 9; cnt_i > 0; cnt_i--)
    {
        InvSubBytes(P);
        InvShiftRows(P);
        AddRoundKey(P, RoundKey[cnt_i]);
        InvMixColumns(P);
    }

    //마지막 10 라운드
    InvShiftRows(P);
    InvSubBytes(P);
    AddRoundKey(P, RoundKey[0]);
}

/*****************************************************AES192********************************************************/

//관련 데이터 타입 정의 내역
//128-bit block
//typedef uint8_t AES_STATE_t[16];
//192-bit masterkey
//typedef uint8_t AES192_KEY_t[24];

void KeySchedule192(unsigned char buf[][16], unsigned char key[24]) {
    /*[0]*/
    buf[0][0] = key[0];
    buf[0][1] = key[1];
    buf[0][2] = key[2];
    buf[0][3] = key[3];

    buf[0][4] = key[4];
    buf[0][5] = key[5];
    buf[0][6] = key[6];
    buf[0][7] = key[7];

    buf[0][8] = key[8];
    buf[0][9] = key[9];
    buf[0][10] = key[10];
    buf[0][11] = key[11];

    buf[0][12] = key[12];
    buf[0][13] = key[13];
    buf[0][14] = key[14];
    buf[0][15] = key[15];
    /*[1]*/
    buf[1][0] = key[16];
    buf[1][1] = key[17];
    buf[1][2] = key[18];
    buf[1][3] = key[19];

    buf[1][4] = key[20];
    buf[1][5] = key[21];
    buf[1][6] = key[22];
    buf[1][7] = key[23];

    buf[1][8] = key[0] ^ s_box[key[21]] ^ RC[0];
    buf[1][9] = key[1] ^ s_box[key[22]];
    buf[1][10] = key[2] ^ s_box[key[23]];
    buf[1][11] = key[3] ^ s_box[key[20]];

    buf[1][12] = buf[1][8] ^ key[4];
    buf[1][13] = buf[1][9] ^ key[5];
    buf[1][14] = buf[1][10] ^ key[6];
    buf[1][15] = buf[1][11] ^ key[7];
    /*[2]*/
    buf[2][0] = buf[1][12] ^ key[8];
    buf[2][1] = buf[1][13] ^ key[9];
    buf[2][2] = buf[1][14] ^ key[10];
    buf[2][3] = buf[1][15] ^ key[11];

    buf[2][4] = buf[2][0] ^ key[12];
    buf[2][5] = buf[2][1] ^ key[13];
    buf[2][6] = buf[2][2] ^ key[14];
    buf[2][7] = buf[2][3] ^ key[15];

    buf[2][8] = buf[2][4] ^ buf[1][0];
    buf[2][9] = buf[2][5] ^ buf[1][1];
    buf[2][10] = buf[2][6] ^ buf[1][2];
    buf[2][11] = buf[2][7] ^ buf[1][3];

    buf[2][12] = buf[2][8] ^ buf[1][4];
    buf[2][13] = buf[2][9] ^ buf[1][5];
    buf[2][14] = buf[2][10] ^ buf[1][6];
    buf[2][15] = buf[2][11] ^ buf[1][7];
    ////////////////////////////////////////
    /*[3]*/
    buf[3][0] = buf[1][8] ^ s_box[buf[2][13]] ^ RC[1];
    buf[3][1] = buf[1][9] ^ s_box[buf[2][14]];
    buf[3][2] = buf[1][10] ^ s_box[buf[2][15]];
    buf[3][3] = buf[1][11] ^ s_box[buf[2][12]];

    buf[3][4] = buf[3][0] ^ buf[1][12];
    buf[3][5] = buf[3][1] ^ buf[1][13];
    buf[3][6] = buf[3][2] ^ buf[1][14];
    buf[3][7] = buf[3][3] ^ buf[1][15];

    buf[3][8] = buf[3][4] ^ buf[2][0];
    buf[3][9] = buf[3][5] ^ buf[2][1];
    buf[3][10] = buf[3][6] ^ buf[2][2];
    buf[3][11] = buf[3][7] ^ buf[2][3];

    buf[3][12] = buf[3][8] ^ buf[2][4];
    buf[3][13] = buf[3][9] ^ buf[2][5];
    buf[3][14] = buf[3][10] ^ buf[2][6];
    buf[3][15] = buf[3][11] ^ buf[2][7];
    /*[4]*/
    buf[4][0] = buf[3][12] ^ buf[2][8];
    buf[4][1] = buf[3][13] ^ buf[2][9];
    buf[4][2] = buf[3][14] ^ buf[2][10];
    buf[4][3] = buf[3][15] ^ buf[2][11];

    buf[4][4] = buf[4][0] ^ buf[2][12];
    buf[4][5] = buf[4][1] ^ buf[2][13];
    buf[4][6] = buf[4][2] ^ buf[2][14];
    buf[4][7] = buf[4][3] ^ buf[2][15];
    ////////////////////////////////////////
    buf[4][8] = buf[3][0] ^ s_box[buf[4][5]] ^ RC[2];
    buf[4][9] = buf[3][1] ^ s_box[buf[4][6]];
    buf[4][10] = buf[3][2] ^ s_box[buf[4][7]];
    buf[4][11] = buf[3][3] ^ s_box[buf[4][4]];

    buf[4][12] = buf[4][8] ^ buf[3][4];
    buf[4][13] = buf[4][9] ^ buf[3][5];
    buf[4][14] = buf[4][10] ^ buf[3][6];
    buf[4][15] = buf[4][11] ^ buf[3][7];
    /*[5]*/
    buf[5][0] = buf[4][12] ^ buf[3][8];
    buf[5][1] = buf[4][13] ^ buf[3][9];
    buf[5][2] = buf[4][14] ^ buf[3][10];
    buf[5][3] = buf[4][15] ^ buf[3][11];

    buf[5][4] = buf[5][0] ^ buf[3][12];
    buf[5][5] = buf[5][1] ^ buf[3][13];
    buf[5][6] = buf[5][2] ^ buf[3][14];
    buf[5][7] = buf[5][3] ^ buf[3][15];

    buf[5][8] = buf[5][4] ^ buf[4][0];
    buf[5][9] = buf[5][5] ^ buf[4][1];
    buf[5][10] = buf[5][6] ^ buf[4][2];
    buf[5][11] = buf[5][7] ^ buf[4][3];

    buf[5][12] = buf[5][8] ^ buf[4][4];
    buf[5][13] = buf[5][9] ^ buf[4][5];
    buf[5][14] = buf[5][10] ^ buf[4][6];
    buf[5][15] = buf[5][11] ^ buf[4][7];
    ////////////////////////////////////////
    /*[6]*/
    buf[6][0] = buf[4][8] ^ s_box[buf[5][13]] ^ RC[3];
    buf[6][1] = buf[4][9] ^ s_box[buf[5][14]];
    buf[6][2] = buf[4][10] ^ s_box[buf[5][15]];
    buf[6][3] = buf[4][11] ^ s_box[buf[5][12]];

    buf[6][4] = buf[6][0] ^ buf[4][12];
    buf[6][5] = buf[6][1] ^ buf[4][13];
    buf[6][6] = buf[6][2] ^ buf[4][14];
    buf[6][7] = buf[6][3] ^ buf[4][15];

    buf[6][8] = buf[6][4] ^ buf[5][0];
    buf[6][9] = buf[6][5] ^ buf[5][1];
    buf[6][10] = buf[6][6] ^ buf[5][2];
    buf[6][11] = buf[6][7] ^ buf[5][3];

    buf[6][12] = buf[6][8] ^ buf[5][4];
    buf[6][13] = buf[6][9] ^ buf[5][5];
    buf[6][14] = buf[6][10] ^ buf[5][6];
    buf[6][15] = buf[6][11] ^ buf[5][7];
    /*[7]*/
    buf[7][0] = buf[6][12] ^ buf[5][8];
    buf[7][1] = buf[6][13] ^ buf[5][9];
    buf[7][2] = buf[6][14] ^ buf[5][10];
    buf[7][3] = buf[6][15] ^ buf[5][11];

    buf[7][4] = buf[7][0] ^ buf[5][12];
    buf[7][5] = buf[7][1] ^ buf[5][13];
    buf[7][6] = buf[7][2] ^ buf[5][14];
    buf[7][7] = buf[7][3] ^ buf[5][15];
    ////////////////////////////////////////
    buf[7][8] = buf[6][0] ^ s_box[buf[7][5]] ^ RC[4];
    buf[7][9] = buf[6][1] ^ s_box[buf[7][6]];
    buf[7][10] = buf[6][2] ^ s_box[buf[7][7]];
    buf[7][11] = buf[6][3] ^ s_box[buf[7][4]];

    buf[7][12] = buf[7][8] ^ buf[6][4];
    buf[7][13] = buf[7][9] ^ buf[6][5];
    buf[7][14] = buf[7][10] ^ buf[6][6];
    buf[7][15] = buf[7][11] ^ buf[6][7];
    /*[8]*/
    buf[8][0] = buf[7][12] ^ buf[6][8];
    buf[8][1] = buf[7][13] ^ buf[6][9];
    buf[8][2] = buf[7][14] ^ buf[6][10];
    buf[8][3] = buf[7][15] ^ buf[6][11];

    buf[8][4] = buf[8][0] ^ buf[6][12];
    buf[8][5] = buf[8][1] ^ buf[6][13];
    buf[8][6] = buf[8][2] ^ buf[6][14];
    buf[8][7] = buf[8][3] ^ buf[6][15];

    buf[8][8] = buf[8][4] ^ buf[7][0];
    buf[8][9] = buf[8][5] ^ buf[7][1];
    buf[8][10] = buf[8][6] ^ buf[7][2];
    buf[8][11] = buf[8][7] ^ buf[7][3];

    buf[8][12] = buf[8][8] ^ buf[7][4];
    buf[8][13] = buf[8][9] ^ buf[7][5];
    buf[8][14] = buf[8][10] ^ buf[7][6];
    buf[8][15] = buf[8][11] ^ buf[7][7];
    ////////////////////////////////////////
    /*[9]*/
    buf[9][0] = buf[7][8] ^ s_box[buf[8][13]] ^ RC[5];
    buf[9][1] = buf[7][9] ^ s_box[buf[8][14]];
    buf[9][2] = buf[7][10] ^ s_box[buf[8][15]];
    buf[9][3] = buf[7][11] ^ s_box[buf[8][12]];

    buf[9][4] = buf[9][0] ^ buf[7][12];
    buf[9][5] = buf[9][1] ^ buf[7][13];
    buf[9][6] = buf[9][2] ^ buf[7][14];
    buf[9][7] = buf[9][3] ^ buf[7][15];

    buf[9][8] = buf[9][4] ^ buf[8][0];
    buf[9][9] = buf[9][5] ^ buf[8][1];
    buf[9][10] = buf[9][6] ^ buf[8][2];
    buf[9][11] = buf[9][7] ^ buf[8][3];

    buf[9][12] = buf[9][8] ^ buf[8][4];
    buf[9][13] = buf[9][9] ^ buf[8][5];
    buf[9][14] = buf[9][10] ^ buf[8][6];
    buf[9][15] = buf[9][11] ^ buf[8][7];
    /*[10]*/
    buf[10][0] = buf[9][12] ^ buf[8][8];
    buf[10][1] = buf[9][13] ^ buf[8][9];
    buf[10][2] = buf[9][14] ^ buf[8][10];
    buf[10][3] = buf[9][15] ^ buf[8][11];

    buf[10][4] = buf[10][0] ^ buf[8][12];
    buf[10][5] = buf[10][1] ^ buf[8][13];
    buf[10][6] = buf[10][2] ^ buf[8][14];
    buf[10][7] = buf[10][3] ^ buf[8][15];
    ////////////////////////////////////////
    buf[10][8] = buf[9][0] ^ s_box[buf[10][5]] ^ RC[6];
    buf[10][9] = buf[9][1] ^ s_box[buf[10][6]];
    buf[10][10] = buf[9][2] ^ s_box[buf[10][7]];
    buf[10][11] = buf[9][3] ^ s_box[buf[10][4]];

    buf[10][12] = buf[10][8] ^ buf[9][4];
    buf[10][13] = buf[10][9] ^ buf[9][5];
    buf[10][14] = buf[10][10] ^ buf[9][6];
    buf[10][15] = buf[10][11] ^ buf[9][7];
    /*[11]*/
    buf[11][0] = buf[10][12] ^ buf[9][8];
    buf[11][1] = buf[10][13] ^ buf[9][9];
    buf[11][2] = buf[10][14] ^ buf[9][10];
    buf[11][3] = buf[10][15] ^ buf[9][11];

    buf[11][4] = buf[11][0] ^ buf[9][12];
    buf[11][5] = buf[11][1] ^ buf[9][13];
    buf[11][6] = buf[11][2] ^ buf[9][14];
    buf[11][7] = buf[11][3] ^ buf[9][15];

    buf[11][8] = buf[11][4] ^ buf[10][0];
    buf[11][9] = buf[11][5] ^ buf[10][1];
    buf[11][10] = buf[11][6] ^ buf[10][2];
    buf[11][11] = buf[11][7] ^ buf[10][3];

    buf[11][12] = buf[11][8] ^ buf[10][4];
    buf[11][13] = buf[11][9] ^ buf[10][5];
    buf[11][14] = buf[11][10] ^ buf[10][6];
    buf[11][15] = buf[11][11] ^ buf[10][7];
    ////////////////////////////////////////
    /*[12]*/
    buf[12][0] = buf[10][8] ^ s_box[buf[11][13]] ^ RC[7];
    buf[12][1] = buf[10][9] ^ s_box[buf[11][14]];
    buf[12][2] = buf[10][10] ^ s_box[buf[11][15]];
    buf[12][3] = buf[10][11] ^ s_box[buf[11][12]];

    buf[12][4] = buf[12][0] ^ buf[10][12];
    buf[12][5] = buf[12][1] ^ buf[10][13];
    buf[12][6] = buf[12][2] ^ buf[10][14];
    buf[12][7] = buf[12][3] ^ buf[10][15];

    buf[12][8] = buf[12][4] ^ buf[11][0];
    buf[12][9] = buf[12][5] ^ buf[11][1];
    buf[12][10] = buf[12][6] ^ buf[11][2];
    buf[12][11] = buf[12][7] ^ buf[11][3];

    buf[12][12] = buf[12][8] ^ buf[11][4];
    buf[12][13] = buf[12][9] ^ buf[11][5];
    buf[12][14] = buf[12][10] ^ buf[11][6];
    buf[12][15] = buf[12][11] ^ buf[11][7];
}

void AES192_enc(AES_STATE_t C, AES_STATE_t P, AES192_KEY_t K192)
{
    unsigned char RoundKey_0[24];
    unsigned char buf[13][16]; //라운드 키를 저장할 배열

    for (int i = 0; i < 24; i++)
    {
        RoundKey_0[i] = K192[i]; //MasterKey를 라운드 키 0으로 저장
    }

    for (int i = 0; i < 16; i++) {
        C[i] = P[i];
    }

    KeySchedule192(buf, RoundKey_0);

    AddRoundKey(C, buf[0]);

    //1~11 라운드
    for (int cnt_i = 0; cnt_i < 11; cnt_i++)
    {
        SubBytes(C);
        ShiftRows(C);
        MixColumns(C);
        AddRoundKey(C, buf[cnt_i + 1]);
    }

    //마지막 12 라운드
    SubBytes(C);
    ShiftRows(C);
    AddRoundKey(C, buf[12]);
}


void AES192_dec(AES_STATE_t P, AES_STATE_t C, AES192_KEY_t K192)
{
    unsigned char RoundKey_0[24];
    unsigned char buf[13][16];

    for (int i = 0; i < 24; i++)
    {
        RoundKey_0[i] = K192[i];
    }

    for (int i = 0; i < 16; i++) {
        P[i] = C[i];
    }

    KeySchedule192(buf, RoundKey_0);

    AddRoundKey(P, buf[12]);

    //1~11 라운드
    for (int cnt_i = 11; cnt_i > 0; cnt_i--)
    {
        InvSubBytes(P);
        InvShiftRows(P);
        AddRoundKey(P, buf[cnt_i]);
        InvMixColumns(P);
    }

    //마지막 12 라운드
    InvShiftRows(P);
    InvSubBytes(P);
    AddRoundKey(P, buf[0]);
}

/*****************************************************AES256********************************************************/

//관련 데이터 타입 정의 내역
//128-bit block
//typedef uint8_t AES_STATE_t[16];
//256-bit masterkey
//typedef uint8_t AES256_KEY_t[32];

void KeySchedule256(unsigned char key[][16])
{
    int cnt_i;

    for (cnt_i = 1; cnt_i < 7; cnt_i++) {
        /*[짝]*/
        key[cnt_i * 2][0] = key[cnt_i * 2 - 2][0] ^ s_box[key[cnt_i * 2 - 1][13]] ^ RC[cnt_i - 1];
        key[cnt_i * 2][1] = key[cnt_i * 2 - 2][1] ^ s_box[key[cnt_i * 2 - 1][14]];
        key[cnt_i * 2][2] = key[cnt_i * 2 - 2][2] ^ s_box[key[cnt_i * 2 - 1][15]];
        key[cnt_i * 2][3] = key[cnt_i * 2 - 2][3] ^ s_box[key[cnt_i * 2 - 1][12]];

        key[cnt_i * 2][4] = key[cnt_i * 2][0] ^ key[cnt_i * 2 - 2][4];
        key[cnt_i * 2][5] = key[cnt_i * 2][1] ^ key[cnt_i * 2 - 2][5];
        key[cnt_i * 2][6] = key[cnt_i * 2][2] ^ key[cnt_i * 2 - 2][6];
        key[cnt_i * 2][7] = key[cnt_i * 2][3] ^ key[cnt_i * 2 - 2][7];

        key[cnt_i * 2][8] = key[cnt_i * 2][4] ^ key[cnt_i * 2 - 2][8];
        key[cnt_i * 2][9] = key[cnt_i * 2][5] ^ key[cnt_i * 2 - 2][9];
        key[cnt_i * 2][10] = key[cnt_i * 2][6] ^ key[cnt_i * 2 - 2][10];
        key[cnt_i * 2][11] = key[cnt_i * 2][7] ^ key[cnt_i * 2 - 2][11];

        key[cnt_i * 2][12] = key[cnt_i * 2][8] ^ key[cnt_i * 2 - 2][12];
        key[cnt_i * 2][13] = key[cnt_i * 2][9] ^ key[cnt_i * 2 - 2][13];
        key[cnt_i * 2][14] = key[cnt_i * 2][10] ^ key[cnt_i * 2 - 2][14];
        key[cnt_i * 2][15] = key[cnt_i * 2][11] ^ key[cnt_i * 2 - 2][15];

        /*[홀]*/
        key[cnt_i * 2 + 1][0] = s_box[key[cnt_i * 2][12]] ^ key[cnt_i * 2 - 1][0];
        key[cnt_i * 2 + 1][1] = s_box[key[cnt_i * 2][13]] ^ key[cnt_i * 2 - 1][1];
        key[cnt_i * 2 + 1][2] = s_box[key[cnt_i * 2][14]] ^ key[cnt_i * 2 - 1][2];
        key[cnt_i * 2 + 1][3] = s_box[key[cnt_i * 2][15]] ^ key[cnt_i * 2 - 1][3];

        key[cnt_i * 2 + 1][4] = key[cnt_i * 2 + 1][0] ^ key[cnt_i * 2 - 1][4];
        key[cnt_i * 2 + 1][5] = key[cnt_i * 2 + 1][1] ^ key[cnt_i * 2 - 1][5];
        key[cnt_i * 2 + 1][6] = key[cnt_i * 2 + 1][2] ^ key[cnt_i * 2 - 1][6];
        key[cnt_i * 2 + 1][7] = key[cnt_i * 2 + 1][3] ^ key[cnt_i * 2 - 1][7];

        key[cnt_i * 2 + 1][8] = key[cnt_i * 2 + 1][4] ^ key[cnt_i * 2 - 1][8];
        key[cnt_i * 2 + 1][9] = key[cnt_i * 2 + 1][5] ^ key[cnt_i * 2 - 1][9];
        key[cnt_i * 2 + 1][10] = key[cnt_i * 2 + 1][6] ^ key[cnt_i * 2 - 1][10];
        key[cnt_i * 2 + 1][11] = key[cnt_i * 2 + 1][7] ^ key[cnt_i * 2 - 1][11];

        key[cnt_i * 2 + 1][12] = key[cnt_i * 2 + 1][8] ^ key[cnt_i * 2 - 1][12];
        key[cnt_i * 2 + 1][13] = key[cnt_i * 2 + 1][9] ^ key[cnt_i * 2 - 1][13];
        key[cnt_i * 2 + 1][14] = key[cnt_i * 2 + 1][10] ^ key[cnt_i * 2 - 1][14];
        key[cnt_i * 2 + 1][15] = key[cnt_i * 2 + 1][11] ^ key[cnt_i * 2 - 1][15];
    }
    /*[14]*/
    key[14][0] = key[12][0] ^ s_box[key[13][13]] ^ RC[6];
    key[14][1] = key[12][1] ^ s_box[key[13][14]];
    key[14][2] = key[12][2] ^ s_box[key[13][15]];
    key[14][3] = key[12][3] ^ s_box[key[13][12]];

    key[14][4] = key[14][0] ^ key[12][4];
    key[14][5] = key[14][1] ^ key[12][5];
    key[14][6] = key[14][2] ^ key[12][6];
    key[14][7] = key[14][3] ^ key[12][7];

    key[14][8] = key[14][4] ^ key[12][8];
    key[14][9] = key[14][5] ^ key[12][9];
    key[14][10] = key[14][6] ^ key[12][10];
    key[14][11] = key[14][7] ^ key[12][11];

    key[14][12] = key[14][8] ^ key[12][12];
    key[14][13] = key[14][9] ^ key[12][13];
    key[14][14] = key[14][10] ^ key[12][14];
    key[14][15] = key[14][11] ^ key[12][15];
}

void AES256_enc(AES_STATE_t C, AES_STATE_t P, AES256_KEY_t K256)
{
    unsigned char RoundKey[15][16];
    for (int i = 0; i < 16; i++)
    {
        RoundKey[0][i] = K256[i];
        RoundKey[1][i] = K256[i + 16];
        C[i] = P[i];
    }

    KeySchedule256(RoundKey);

    AddRoundKey(C, RoundKey[0]);

    //1~13 라운드
    for (int cnt_i = 0; cnt_i < 13; cnt_i++)
    {
        SubBytes(C);
        ShiftRows(C);
        MixColumns(C);
        AddRoundKey(C, RoundKey[cnt_i + 1]);
    }

    //마지막 14 라운드
    SubBytes(C);
    ShiftRows(C);
    AddRoundKey(C, RoundKey[14]);
}


void AES256_dec(AES_STATE_t P, AES_STATE_t C, AES256_KEY_t K256)
{
    unsigned char RoundKey[15][16];
    for (int i = 0; i < 16; i++)
    {
        RoundKey[0][i] = K256[i];
        RoundKey[1][i] = K256[i + 16];
        P[i] = C[i];
    }

    KeySchedule256(RoundKey);

    AddRoundKey(P, RoundKey[14]);

    //1~13 라운드
    for (int cnt_i = 13; cnt_i > 0; cnt_i--)
    {
        InvSubBytes(P);
        InvShiftRows(P);
        AddRoundKey(P, RoundKey[cnt_i]);
        InvMixColumns(P);
    }

    //마지막 14 라운드
    InvShiftRows(P);
    InvSubBytes(P);
    AddRoundKey(P, RoundKey[0]);
}
