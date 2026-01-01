#include <stdio.h>    
#include "test.h"    
#include <stdbool.h>
#include <stdint.h>

#define N 50000
#define D 1000
static inline uint32_t xorshift32(uint32_t* s) {
    uint32_t x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    return x;
}
void make_collision_fail_sample_bias(int* out, uint32_t seed) {
    const int HEAVY = 0;        // 가장 자주 나오는 심볼
    const double HEAVY_PCT = 1;   // 35%를 0으로 고정 -> p_max=0.35

    for (int i = 0; i < N; i++) {
        uint32_t r = xorshift32(&seed) % 100;
        if ((int)r < HEAVY_PCT) {
            out[i] = HEAVY;  // 자주 반복 -> 첫 충돌 매우 빨라짐
        }
        else {
            out[i] = (int)(xorshift32(&seed) % 16);
        }
    }
}
void make_sample(int out[N], uint32_t seed)
{
    uint32_t st = seed ? seed : 1u;

    // ====== 5개 파라미터(함수 내부에서 고정) ======
    const int LAG_C = 128;   // 근접 지연(충돌 망가뜨리기)
    const int BURST_LEN = 16;   // burst 길이
    const int BURST_RATE = 5;    // burst 시작 확률(%)

    const int LAG_M = 1024; // 원거리 지연
    const int COPY_PCT_M = 20;   // 원거리 복사 확률(%)
    // =============================================

    // 1) 전체 i.i.d. 랜덤(0~15)
    for (int i = 0; i < N; i++) {
        out[i] = (int)(xorshift32(&st) & 0x0F);
    }

    int start = D + 1;
    if (start < LAG_C) start = LAG_C;

    int burst_left = 0;

    for (int i = start; i < N; i++) {
        // (A) burst 시작 여부 결정
        if (burst_left == 0) {
            int rr = (int)(xorshift32(&st) % 100u);
            if (rr < BURST_RATE && i - LAG_C >= 0) {
                burst_left = BURST_LEN;
            }
        }

        // (B) burst 중이면 근접 복사(=Collision에 강력한 근접 중복 유발)
        if (burst_left > 0) {
            out[i] = out[i - LAG_C];
            burst_left--;
            continue;
        }

        // (C) burst가 아니면, 기존 장거리 복사는 "가볍게"만
        if (i - LAG_M >= 0) {
            int r = (int)(xorshift32(&st) % 100u);
            if (r < COPY_PCT_M) {
                out[i] = out[i - LAG_M];
            }
        }
        // else: 랜덤 유지
    }
}


int main(void)
{
    int s[N];
    unsigned char cState1[8] = { 1,0,1,0,1,0,1,0 };  // 초기 상태 (0만 아니면 됨)

    for (int i = 0; i < 50000; i++)
    {
        s[i] = 0;
        for (int n = 0; n < 4; n++)
        {
            // 기존 구조 유지: 최상위 비트를 출력으로 사용
            s[i] |= (cState1[7] << (4 - 1 - n));

            // 8차 다항식(예): P(x)=x^8 + x^6 + x^5 + x^4 + 1
            // taps(0-index): 7,5,4,3
            unsigned char temp = (cState1[7] ^ cState1[5] ^ cState1[4] ^ cState1[3]) & 1;

            // 시프트 구조 그대로, 범위만 7로 축소
            for (int j = 7; j > 0; j--)
            {
                cState1[j] = cState1[j - 1];
            }
            cState1[0] = temp;
        }
    }
    int idx = 0;  // b[] 인덱스
    int b[N] = { 0 };
    unsigned char cState[8] = { 1,0,1,0,1,0,1,0 };  // 초기상태 (전부 0만 아니면 됨)

    for (int i = 0; i < N; i++)
    {
        // 출력 비트: MSB
        b[i] = cState[7];

        // 8차 다항식: P(x)=x^8 + x^6 + x^5 + x^4 + 1
        // taps(0-based): 7,5,4,3
        unsigned char temp =
            (cState[7] ^ cState[5] ^ cState[4] ^ cState[3]) & 1;

        // 시프트 구조 유지: 7<-6<-...<-1<-0, 새 비트는 0에
        for (int j = 7; j > 0; j--)
            cState[j] = cState[j - 1];

        cState[0] = temp;
    }
    static int out[N];
    //make_sample(out,123456789u);
    make_collision_fail_sample_bias(out, 123456789u);
    //123456789u
    //56789u1234
    double H = -1;
    int k = 16;
    printf("\n");
    double P_value = -1;
    autocorrelation(b, N, &P_value);
    printf("p=%f\n", P_value);
    //Collision(out, N, k, &H); //ttak
    Collision(s, N, k, &H);
    printf("H = %f\n", H);
    //Compression(out, N, k, &H);
    Compression(s, N, k, &H);
    printf("H = %f\n", H);
}