/* F(z,t,u) 구현 */
#include <stdio.h>   
#include "test.h"  
#include <math.h>   
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
double igam(double a, double x);
double igamc(double a, double x);


double factorial(int n)
{
    if (n < 0) return 0;   // 음수 입력 방어

    double result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}
static double CF(double z, int k)
{
    const double inv = 1.0 / z;
    double ig = igamc(k + 1.0, inv);
    double f = factorial(k);
    double p = pow(inv, -(double)(k + 1));
    double e = exp(inv);

    /*printf("z=%.6e inv=%.6e\n", z, inv);
    printf("ig=%.6e  \n", ig);
    printf("f =%.6e  \n", f);
    printf("p =%.6e  \n", p);
    printf("e =%.6e  \n", e);
    */

    double re = (f * ig) * p * e; //nist
    //double re = ( ig) * p * e;
    return re;
}
static double Nist_CF(double q, int k)
{

    //return q + (2 * q * q) + (2 * q * q * q);

    return
        (factorial(16) * q * q * q * q * q * q * q * q * q * q * q * q * q * q * q * q * q) +
        (factorial(16) / factorial(1)) * q * q * q * q * q * q * q * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(2)) * q * q * q * q * q * q * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(3)) * q * q * q * q * q * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(4)) * q * q * q * q * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(5)) * q * q * q * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(6)) * q * q * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(7)) * q * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(8)) * q * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(9)) * q * q * q * q * q * q * q * q +
        (factorial(16) / factorial(10)) * q * q * q * q * q * q * q +
        (factorial(16) / factorial(11)) * q * q * q * q * q * q +
        (factorial(16) / factorial(12)) * q * q * q * q * q +
        (factorial(16) / factorial(13)) * q * q * q * q +
        (factorial(16) / factorial(14)) * q * q * q +
        (factorial(16) / factorial(15)) * q * q +
        (factorial(16) / factorial(16)) * q;

}
static double P_RHS(double p, int k, double target)
{
    if (k <= 1) return NAN;
    if (!(p > 0.0 && p < 1.0)) {

        return NAN;
    }

    double q = (1.0 - p) / (double)(k - 1);


    double pinv = 1.0 / p;
    double qinv = 1.0 / q;

    double A = p * (qinv * qinv) * (1.0 + (1.0 / (double)k) * (pinv - qinv));
    double B = p * qinv * (1.0 / (double)k) * (pinv - qinv);

    double rhs = A * Nist_CF(q, k) - B;
    return (rhs - target);
}
bool bisection_root1(double lo, double hi, double eps, int max_iter,
    int k, double target,
    double* root_out)
{
    /* 출력 포인터 유효성 및 입력값 검사 */
    if (!root_out || !isfinite(lo) || !isfinite(hi) ||
        !isfinite(eps) || eps <= 0.0)
        return false;

    /* 반복 횟수는 양수여야 함 */
    if (max_iter <= 0) return false;

    /* lo <= hi 가 되도록 정렬 */
    if (lo > hi) {
        double t = lo;
        lo = hi;
        hi = t;
    }

    /* 구간 양 끝에서 함수값 계산 */
    double flo = P_RHS(lo, k, target);
    double fhi = P_RHS(hi, k, target);

    /* 함수값이 NaN/Inf이면 수치적으로 실패 */
    if (!isfinite(flo) || !isfinite(fhi)) return false;

    /* lo 자체가 이미 근에 충분히 가까운 경우 */
    if (fabs(flo) <= eps) {
        *root_out = lo;
        return true;
    }

    /* hi 자체가 이미 근에 충분히 가까운 경우 */
    if (fabs(fhi) <= eps) {
        *root_out = hi;
        return true;
    }

    /*
     * 중간값 정리 조건 검사
     * [lo, hi] 구간에 근이 있으려면
     * f(lo)와 f(hi)는 서로 부호가 달라야 한다
     */
    if ((flo > 0.0 && fhi > 0.0) ||
        (flo < 0.0 && fhi < 0.0))
    {
        printf("중간값 정리x\n");
        return false;
    }

    /* 이분법 반복 */
    for (int i = 0; i < max_iter; i++) {

        /* 구간 중간점 */
        double mid = 0.5 * (lo + hi);

        /* 중간점에서 함수값 계산 */
        double fmid = P_RHS(mid, k, target);
        if (!isfinite(fmid)) return false;

        /*
         * 종료 조건:
         * 1) 함수값이 eps 이내로 0에 근접
         * 2) 구간 길이가 충분히 작아짐
         */
        if (fabs(fmid) <= eps || 0.5 * (hi - lo) <= eps) {
            *root_out = mid;
            return true;
        }

        /*
         * 부호 변화를 이용해 다음 구간 선택
         * f(lo)와 f(mid)의 부호가 다르면 근은 [lo, mid]
         * 아니면 [mid, hi]
         */
        if ((flo < 0.0 && fmid > 0.0) ||
            (flo > 0.0 && fmid < 0.0)) {
            hi = mid;
            fhi = fmid;
        }
        else {
            lo = mid;
            flo = fmid;
        }
    }

    /*
     * 최대 반복 횟수에 도달했을 경우
     * 현재 구간의 중간값을 근의 근사값으로 반환
     */
    *root_out = 0.5 * (lo + hi);
    return true;
}

void Collision(int s[], int L, int k, double* H)//50000,k=2^l(4)=16
{
    int v = 0;
    int index = 0; //
    int* t = (int*)malloc(L * sizeof(int));
    while (index < L) {
        int temp_j = -1;

        for (int j = index + 1; j < L; j++) {

            int found = 0;

            for (int i = index; i < j; i++) {
                if (s[i] == s[j]) {
                    //충돌발생
                    found = 1;
                    break;
                }
            }
            if (found) {
                temp_j = j;
                break;
            }
        }

        if (temp_j == -1) break;

        t[v] = temp_j - index + 1;
        v++;
        index = temp_j + 1;
    }
    if (v > 0) v--;
    int sum1 = 0;
    for (int i = 0; i < v; i++)
    {
        sum1 += t[i];
    }
    double X = (double)sum1 / (double)v;
    double sum2 = 0;
    for (int i = 0; i < v; i++)
    {
        sum2 += (t[i] - X) * (t[i] - X);
    }
    double sigma = sqrt(sum2 / (v));
    double X_prime = X - 2.576 * sigma / sqrt(v);

    printf("v(충돌횟수)=%d\n", v);
    printf("X(표본평균)=%f\n", X);
    printf("sigma(편차)=%f\n", sigma);
    printf("X_prime(X의 99% 신뢰구간 하한)=%f\n", X_prime);
    double  p = 0.3;

    //printf("(0.732을 식에 대입한 값)-X_prime= %f\n", P_RHS(0.732, k, X_prime));
    //nist문서에는 근이 0.5~1 사이 라고 언급 이유x
    //0.99까지는 이진탐색O, 0.999부터는 X
    //ttak:1e-4, 1-(1e-4), 1e-4
    //1e-4, 1- (1e-4)
    // 상한이 절대1이되면 안됨 1/(1-p)가 정의가 안됨

    //printf("(0.99을 식에 대입한 값)-X_prime= %f\n", P_RHS(0.99, k, X_prime));

    //printf("(0.99을 식에 대입한 값)-X_prime= %f\n", P_RHS(0.99, k, X_prime));

    printf("(1/k을 식에 대입한 값)-X_prime= %f\n", P_RHS(1.0/k, k, X_prime));
    printf("(0.9999을 식에 대입한 값)-X_prime= %f\n", P_RHS(0.9999, k, X_prime));
    // (1e-5) =0.00001
     // (1e-7) =0.0000001
    bool ok = bisection_root1(1.0/k, 0.9999, 0.0001, 200,
        k, X_prime,
        &p);
    *H = -((log(p) / log(2.0)));




    if (ok)
    {
        printf("이진탐색 성공!\n");
        printf("p=%f\n", p);
        *H = -log(p) / log(2.0);
        printf("(p을 식에 대입한 값)-X_prime= %f\n", P_RHS(p, k, X_prime));
    }
    else
    {
        printf("이진탐색 실패!\n");
        *H = -log(k) / log(2.0);
    }

  

    free(t);

}
double findMax(double arr[], int size) {
    double max = arr[0]; // 첫 번째 원소로 초기화
    for (int i = 1; i < size; i++) {
        if (arr[i] > max) {
            max = arr[i]; // 더 큰 값 발견 시 갱신
        }
    }
    return max;
}
static double F(double z, int t, int u) {
    double one_minus = 1.0 - z;
    if (u < t) {
        /* z^2 * (1-z)^(u-1) */
        return (z * z) * pow(one_minus, (u - 1));
    }
    else {
        /* u == t: z * (1-z)^(t-1) */
        return z * pow(one_minus, (t - 1));
    }
}

/* G(z) 구현: d, L을 받아 v=L-d로 계산 */
static double G(double z, int L, int d) {
    int v = L - d;
    double sum = 0.0;

    for (int t = d + 1; t <= L; t++) {
        for (int u = 1; u <= t; u++) {
            double fu = F(z, t, u);
            sum += (log(u) / log(2)) * fu;
        }
    }

    return sum / (double)v;
}
static double G_RHS(double p, int k, int  L, int d, double target)
{
    double q = (1 - p) / (k - 1);
    double RHS = G(p, L, d) + (k - 1) * G(q, L, d);
    return RHS - target;
}
bool bisection_root2(double lo, double hi, double eps, int max_iter,
    int k, int L, int d, double target,
    double* root_out)
{
    if (!root_out || !isfinite(lo) || !isfinite(hi) ||
        !isfinite(eps) || eps <= 0.0 || max_iter <= 0)
        return false;

    if (lo > hi) { double t = lo; lo = hi; hi = t; }

    double flo = G_RHS(lo, k, L, d, target);
    double fhi = G_RHS(hi, k, L, d, target);
    if (!isfinite(flo) || !isfinite(fhi)) return false;
    printf("flo = %f", flo);
    printf("fhi = %f", fhi);

    /* 단조 감소 함수에서 근 존재 조건 */
    if (flo < 0.0 || fhi > 0.0) {
        printf("근 구간 아님 (단조감소 가정 위배)\n");
        return false;
    }

    for (int i = 0; i < max_iter; i++) {
        double mid = 0.5 * (lo + hi);
        double fmid = G_RHS(mid, k, L, d, target);
        if (!isfinite(fmid)) return false;

        if (fabs(fmid) <= eps || 0.5 * (hi - lo) <= eps) {
            *root_out = mid;
            return true;
        }

        /* 단조 감소 성질 이용 */
        if (fmid > 0.0) {
            lo = mid;   // p가 아직 작음
        }
        else {
            hi = mid;   // p가 너무 큼
        }
    }

    *root_out = 0.5 * (lo + hi);
    return true;
}
void Compression(int s[], int L, int k, double* H)
{
    int* dict = (int*)malloc(sizeof(int) * k);// 1)

    int d = 1000;
    for (int i = 0; i < d; i++) {
        int si = s[i];
        if (0 <= si && si < k) dict[si] = i;
    }
    int v = L - d;
    int* D = (int*)malloc(sizeof(int) * v);
    if (!D) { *H = -2; free(dict); return; }
    for (int i = d; i < L; i++) //3.2
    {
        int si = s[i];
        int idx = i - d;

        if (dict[si] != 0)// 3.2.1)
            //샘플이 등장한 경우
        {
            D[idx] = (i - dict[si]);
            //이전 등장 위치간 간격
            //반복이 많이 일어나면 작아짐
            dict[si] = i;
            //마지막 위치 갱신
        }
        else//3.2.2) 
         //첫등장
        {
            dict[si] = i;
            D[idx] = i;
        }
    }

    /* 4) b = floor(log2(k-1)) + 1 */
    double b = floor(log((double)(k - 1)) / log(2.0)) + 1.0;
    printf("b(b = floor(log2(k-1)) + 1) = %f\n", b);
    /* 5) Xbar = (1/v) * sum log2(D[i]) */
    double sum_log = 0.0;
    double sum_log2 = 0.0;
    for (int i = 0; i < v; i++) {

        double ld = log(D[i]) / log(2.0); /* log2 */
        sum_log += ld;
        sum_log2 += ld * ld;
    }
    double X = sum_log / (double)v;
    printf("X(표본평균) = %f\n", X);
    double c = 0.7
        - (0.8 / b)
        + ((4.0 + (32.0 / b)) * pow((double)v, -3.0 / b)) / 15.0;
    printf("c= %f\n", c);
    double var = (sum_log2 / (double)v) - (X * X);
    if (var < 0.0) var = 0.0;
    double sigma = c * sqrt(var);
    double X_prime = X - (2.576 * sigma) / sqrt(v);
    double p = 0.3;
    printf("sigma(편차) = %f\n", sigma);
    printf("X_prime(X의 99% 신뢰구간 하한) = %f\n", X_prime);
    bool ok = bisection_root2(1.0 / k, 0.5, 0.0001, 3,
        k, L, d, X_prime,
        &p);
    if (ok)
    {
        printf("성공\n");
        printf("근 =%f\n", p);
        printf("(p을 식에 대입한 값)-X_prime= %f\n", G_RHS(p, k, L,d, X_prime));
        *H = -log(p) / log(2);
    }
    else {
        printf("실패\n");
    }

    free(D);
    free(dict);

}
void autocorrelation(int b[], int N, double* P_value)
{
    int* Z = (int*)malloc((N / 4 + 1) * sizeof(int));
    for (int tau = 1; tau <= N / 4; tau++) {
        int sum = 0;
        for (int i = 0; i < N / 4; i++) {
            sum += b[i] ^ b[i + tau];
        }
        Z[tau] = sum;
    }// 2)
    double mu = N / 8.0;
    double sigma = sqrt(N / 16.0);// 3)
    printf("mu(평균)= %f, sigma(편차) =%f\n\n", mu, sigma);
    int tau0 = 1;
    double temp = fabs(Z[1] - mu);
    for (int tau = 2; tau <= N / 4; tau++) {
        double diff = fabs(Z[tau] - mu);
        if (diff > temp) {
            temp = diff;
            tau0 = tau;
        }
    }// 4)
    tau0--;//(인덱스반영)
    printf("tau0(크게만드는tau)= %d\n\n", tau0);
    int c = 0;
    for (int i = N / 2; i < 3 * N / 4; i++) {
        c += b[i] ^ b[i + tau0];
    }// 5)
    printf("c(뒷비트열에서의합)= %d\n\n", c);
    double S_obs = fabs((c - mu) / sigma);
    printf("S_obs= %f\n\n", S_obs);
    *P_value = erfc(S_obs / sqrt(2.0));
}

