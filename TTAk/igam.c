#pragma once
#include <math.h>
#include "mconf1.h"

/* 함수 원형 */
double igam(double a, double x);
double igamc(double a, double x);

static const double MACHEP = 1.11022302462515654042E-16;
static const double BIG = 4.503599627370496e15;
static const double BIGINV = 2.22044604925031308085e-16;

/* 하위 불완전 감마 함수 */
double igam(double a, double x) {
    double ans, ax, c, r;

    if ((x <= 0) || (a <= 0)) {
        return 0.0;
    }

    if ((x > 1.0) && (x > a)) {
        return 1.0 - igamc(a, x);  // 여기서 igamc 호출
    }

    ax = a * log(x) - x - lgamma(a);
    if (ax < -709.78271289338399) {
        return 0.0;
    }
    ax = exp(ax);

    r = a;
    c = 1.0;
    ans = 1.0;

    do {
        r += 1.0;
        c *= x / r;
        ans += c;
    } while (c / ans > MACHEP);

    return ans * ax / a;
}

/* 상위 불완전 감마 함수 */
double igamc(double a, double x) {
    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    if ((x <= 0) || (a <= 0)) {
        return 1.0;
    }

    if ((x < 1.0) || (x < a)) {
        return 1.0 - igam(a, x);  // 여기서 igam 호출
    }

    ax = a * log(x) - x - lgamma(a);
    if (ax < -709.78271289338399) {
        return 0.0;
    }
    ax = exp(ax);

    /* continued fraction */
    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1 / qkm1;

    do {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        yc = y * c;
        pk = pkm1 * z - pkm2 * yc;
        qk = qkm1 * z - qkm2 * yc;
        if (qk != 0) {
            r = pk / qk;
            t = fabs((ans - r) / r);
            ans = r;
        }
        else {
            t = 1.0;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if (fabs(pk) > BIG) {
            pkm2 *= BIGINV;
            pkm1 *= BIGINV;
            qkm2 *= BIGINV;
            qkm1 *= BIGINV;
        }
    } while (t > MACHEP);

    return ans * ax;
}

