/* mconf.h : Cephes 기본 설정 */
#ifndef MCON1F_H
#define MCONF1_H
#include <stdio.h>
#include <math.h>
#include <float.h>

/* Cephes 오류 처리 (간단 버전) */
static void mtherr(const char* name, int code) {
    // 원래는 다양한 에러 코드 처리, 여기선 간단히 출력만
    fprintf(stderr, "Cephes error in %s (code %d)\n", name, code);
}

#endif
#pragma once
