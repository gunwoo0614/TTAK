#include <stdio.h>   
#include "test.h"  
#include "mconf1.h" 
#include <stdint.h>
#include <stdbool.h>
#include <fenv.h>

   
double igam(double a, double x);
double igamc(double a, double x);
void bit_string_conversion(int L, int l, int s[], int b[]) //[L*l]
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < l; j++)
        {
            b[(i) * l + j] = s[i] >> (l - j-1) & 1;
        }
    }
}

void monobit(int b[], int N, double* P_value) //int b[n]
{
	int ctr = 0;
	for (int i = 0; i < N; i++)
	{
		if (b[i])
		{
			ctr++;
		}
	}
	printf("ctr(개수)  =%d\n\n", ctr);
	double myu = N / 2.0;
	printf("myu(평균)  =%f\n\n", myu);
	double sigma = sqrt(N / 4.0);
	printf("sigma(편차)  =%f\n\n", sigma);
	double S_obs = fabs((ctr - myu) / sigma);
	printf("S_obs  =%f\n\n", S_obs);
	*P_value = erfc(S_obs / sqrt(2));
}
void poker(int b[], int N, double* P_value)//9/19
{
    int n = N / 4;
    int* d = (int*)malloc(n * sizeof(int));
    for (int j = 0; j < n; j++)
    {
        d[j] = (b[4 * j] << 3) | (b[4 * j + 1] << 2) | (b[4 * j + 2] << 1) | b[4 * j + 3];

    }// 1) 4비트 샘플로 변형
    int E[16] = { 0 }; // 2)
    for (int j = 0; j < n; j++)
    {
        E[d[j]]++;
    }// 3) 0,1...15 의 개수 각각 카운팅
    for (int j = 0; j < 16; j++)
    {
        printf("%d       ", E[j]);
    }
    printf("\n");
    int sum = 0;
    for (int i = 0; i <= 15; i++)
    {
        sum += (E[i] * E[i]);
    }
    double z = (sum * 16) / (N / 4.0) - (N / 4.0); // 4)
    printf("\n z(통계량)=%f\n\n", z);
    *P_value = igamc(15 / 2.0, z / 2.0); //5)
}
void run(int b[], int N, double* P_value1, double* P_value0)
{
    int R[2][7] = { 0 };// 1)
    int run_val = b[0];
    int run_len = 1;
    for (int i = 1; i < N; i++)
    {
        if (b[i] == run_val) {
            run_len++;
        }
        else {
            if (run_len > 6) run_len = 6;
            R[run_val][run_len]++;   // run 기록
            run_val = b[i];
            run_len = 1;
        }
    } // 3)
    if (run_len > 6) run_len = 6; //0또는 1이 계속 반복될 때
    R[run_val][run_len]++;
    for (int i = 1; i <= 6; i++)
    {
        printf("%d  ", R[0][i]);
    }
    printf("\n\n");
    for (int i = 1; i <= 6; i++)
    {
        printf("%d  ", R[1][i]);
    }
    double mu[7] = { 0 }; // 4)
    double sum = 0;
    for (int i = 1; i <= 5; i++)
    {
        mu[i] = (N - i + 3) / (double)(1 << (i + 2));
    } // 4.1)
    mu[6] = ((N - 4) / (double)(1 << 7)) - (1 / (double)(1 << (N + 2)));
    printf("\n mu[6] = %f\n\n", mu[6]); // 문서에서는 단순 mu
    double z1 = 0;
    double z0 = 0;
    for (int i = 1; i <= 6; i++)
    {
        z1 = z1 + ((R[1][i] - mu[i]) * (R[1][i] - mu[i])) / mu[i];
        z0 = z0 + ((R[0][i] - mu[i]) * (R[0][i] - mu[i])) / mu[i];
    }// 5)
    printf("z1(통게량) = %f\n\n", z1);
    printf("z0(통게량) = %f\n\n", z0);
    *P_value1 = igamc(5 / 2.0, z1 / 2.0);
    *P_value0 = igamc(5 / 2.0, z0 / 2.0); // 6)
}
void long_run(int b[], int N, double* P_value)
{
    int MR0 = 0;
    int MR1 = 0;
    int run_val = b[0];
    int run_len = 1;
    for (int i = 1; i < N; i++)
    {
        if (b[i] == run_val)
        {
            run_len++;
        }
        else
        {
            if (run_val == 0) {
                if (run_len > MR0) { //가장 긴 런 반영
                    MR0 = run_len;
                }
            }
            if (run_val == 1) {
                if (run_len > MR1) {
                    MR1 = run_len;
                }
            }
            run_val = b[i];
            run_len = 1;
        }
    }
    if (run_val == 0)//0,1 계속 반복 될때
    {
        if (run_len > MR0)
        {
            MR0 = run_len;
        }
    }
    if (run_val == 1)
    {
        if (run_len > MR1)
        {
            MR1 = run_len;
        }
    }
    int M;
    if (MR0 > MR1)
    {
        M = MR0;
    }
    else
    {
        M = MR1;
    }
    printf("\n M(최대런)= %d\n", M);
    *P_value = (N - M) / (double)(1 << M);
}

int longest_run_of_ones_in_block(int b[], int start)
{
    int max = 0;
    int ctr = 0;
    for (int i = 0; i < 128; i++) {
        if (b[start + i] == 1) {
            ctr++;
            if (ctr > max) max = ctr;
        }
        else ctr = 0;
    }
    return max;
}
int category_index1(int L) {
    if (L <= 1) return 0;   // 
    if (L == 2) return 1;
    if (L == 3) return 2;
    return 3;               
}
int category_index2(int L) {
    if (L <= 4) return 0;   // ≤4
    if (L == 5) return 1;
    if (L == 6) return 2;
    if (L == 7) return 3;
    if (L == 8) return 4;
    return 5;               // ≥9
}
void NIST_long_run1(int b[], double* P_value)
{
    double PI[4] = { 0.2148, 0.3672, 0.2305, 0.1875};
    int N_B = 128 / 8; //블록 개수
    int v[4] = { 0 };
    // 각 블록의 최장 런 길이 카테고리별로 집계
    for (int i = 0; i < N_B; i++)
    {
        int start = i * 8;
        int L = longest_run_of_ones_in_block(b, start);
        int idx = category_index1(L);
        v[idx]++;
    }
    for (int i = 0; i < 4; i++)
    {
        printf("v[%d] = %d\n\n ", i, v[i]);
    }
    double z = 0.0;
    for (int i = 0; i < 4; i++) {
        double Ei = PI[i] * (double)N_B;   // 기대값
        double diff = (double)v[i] - Ei;
        z += (diff * diff) / Ei;
    }
    printf("z(통계량) = %f\n\n", z);
    *P_value = igamc(3.0 / 2.0, z / 2.0);
}
void NIST_long_run2(int b[], double* P_value)
{
    double PI[6] = { 0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124 };
    int N_B = 50000 / 128; //블록 개수
    int v[6] = { 0 };
    // 각 블록의 최장 런 길이 카테고리별로 집계
    for (int i = 0; i < N_B; i++)
    {
        int start = i * 128;
        int L = longest_run_of_ones_in_block(b, start);
        int idx = category_index2(L);
        v[idx]++;
    }
    for (int i = 0; i < 6; i++)
    {
        printf("v[%d] = %d\n\n ",i  ,v[i]);
    }
    double z = 0.0;
    for (int i = 0; i < 6; i++) {
        double Ei = PI[i] * (50000 / 128);   // 기대값
        double diff = (double)v[i] - Ei;
        z += (diff * diff) / Ei;
    }
    printf("z(통계량) = %f\n\n", z);
    *P_value = igamc(5 / 2.0, z / 2.0);
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


void repetition_count(int s[], int L, int C, double* TF)
{ 
    *TF = 0.0;          
    int t = s[0];      
    int count = 1;     

    for (int i = 1; i < L; i++)
    {
        if (s[i] == t)
        {
            count++;
            if (count == C)
            {
               
                *TF = 1.0; 
                return;
            }
        }
        else
        {
            t = s[i];     
            count = 1;    
        }
    }
    printf("count=%d\n", count);
}
void adaptive_proportion(int s[], int W, int C, double* TF)
{
    *TF = 0.0;         
    int t = s[0];     
    int count = 1;     

    for (int i = 1; i < W; i++)
    {
        if (s[i] == t)
        {
            count++;
         
        }
    }
    if (count == C)
    {
        *TF = 1.0; 
        return;
    }
    printf("count=%d\n", count);
}

void The_Most_Common_Value(int s[], int L, int k, double* H)//50000,k=2^l(4)=16
{
    int* cnt = (int*)calloc(k, sizeof(int));
    for (int i = 0; i < L; i++) {
        int v = s[i];
        if (v >= 0 && v < k) {
            cnt[v]++;
        }
    } // 빈도수 세기
    int max = 0;
    for (int i = 0; i < k; i++) {
        if (cnt[i] > max) max = cnt[i];
    }// 1) max(p)찾기

    double p_hat = max / (double)L;
    printf("p_hat=%f\n", p_hat);
    double temp = p_hat + 2.576 * sqrt(p_hat * (1.0 - p_hat) / (double)L);
    // P(-2.576<= Z <=2.576) =0.99
    double p_u = (temp > 1.0) ? 1.0 : temp;
    printf("p_u=%f\n", p_u);
    *H = -(log(p_u) / log(2.0)); // 3)
}







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
        return false;

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

    //printf("(0.096을 식에 대입한 값)-X_prime= %f\n", P_RHS(0.096, k, X_prime));
    // (1e-5) =0.00001
     // (1e-7) =0.0000001
    bool ok = bisection_root1(1.0/k, 1-(1e-7), (1e-7), 200,
        k, X_prime, 
        &p);
    *H = -((log(p) / log(2.0)));


 
   
    if (ok)
    {
        printf("이진탐색 성공!\n");
        printf("p=%f\n", p);
        printf("(p을 식에 대입한 값)-X_prime= %f\n", P_RHS(p, k, X_prime));
    }
    *H = -log(p) / log(2.0);

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
void Markov(int s[], int L, int k, double* H)//50000,k=2^l(4)=16
{
    int d = 128;// 1)
    double alpha = pow(0.99, k * k);
    printf("alpha=%f\n", alpha);
    if (pow(0.99, k * k) >= pow(0.99, d))
    {
        alpha = pow(0.99, d);
    }// 2)
    int* O = (int*)calloc(k, sizeof(int));
    double* P = (double*)malloc(k * sizeof(double)); //3)
    for (int i = 0; i < L; i++)//문서에서는 k
    {
        int v = s[i];
        if (v >= 0 && v < k)
        {
            O[v]++;
        }
    }// 4)
    double eps = sqrt(-log(1 - alpha) / (log(2) * 2 * L));// 5)
    printf("eps=%f\n", eps);
    printf("\n");
    for (int i = 0; i < k; i++)
    {
        P[i] = 1;
        if (1 > (O[i] / (double)L) + eps)
        {
            P[i] = (O[i] / (double)L) + eps;
            printf("P[%d]=%f   ", i, P[i]);
        }
    }// 6)
    printf("\n");
    O[s[L - 1]] = O[s[L - 1]] - 1;// 7)
    printf("O[s[L-1]]=%d\n", O[s[L - 1]]);
    int** R = (int**)calloc(k, sizeof(int*));
    for (int i = 0; i < k; i++) {
        R[i] = (int*)calloc(k, sizeof(int));
    }
    double** B = (double**)malloc(k * sizeof(double*)); // 8)
    for (int i = 0; i < k; i++) {
        B[i] = (double*)malloc(k * sizeof(double));
    } // 8)
    double* E = (double*)malloc(k * sizeof(double));// 8)
    for (int t = 0; t < L - 1; t++) {
        int i = s[t];     // 현재 값
        int j = s[t + 1];   // 다음 값
        if (i < k && j < k) {
            R[i][j]++;
        }
    }// 9)
    printf("\n");
    for (int i = 0; i < k; i++)
    {
        E[i] = sqrt((-log(1 - alpha) / log(2)) / (double)(2 * O[i]));
        printf("E[%d] = %f   ", i, E[i]);
    }// 10)
    printf("\n");
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            if (O[i] == 0) B[i][j] = 0;
            else
            {
                B[i][j] = 1 < (R[i][j] / (double)O[i]) + E[i] ? 1 : (R[i][j] / (double)O[i]) + E[i];
            }
        }
    }// 11)
    double* M = (double*)malloc(k * sizeof(double));
    double* p_prime = (double*)malloc(k * sizeof(double));
    for (int m = 0; m < d-1 ; m++)
    {
        for (int n = 0; n < k; n++)
        {
            for (int i = 0; i < k; i++)
            {
                p_prime[i] = P[i] * B[i][n];
            }
            M[n] = findMax(p_prime, k);
        }
        for (int i = 0; i < k; i++)
        {
            P[i] = M[i];
        }
    }
    double p_max = findMax(P, k);
    printf("\n-log_2(p_max)=%f\n", -log(p_max) / log(2));
    *H = (-1.0 / d) * log(p_max) / log(2);
}



/* F(z,t,u) 구현 */
static double F(double z, int t, int u) {
     double one_minus = 1.0 - z;

    if (u < t) {
        /* z^2 * (1-z)^(u-1) */
        return (z * z) * pow(one_minus,(u - 1));
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

    for (int t = d+1 ; t <= L; t++) {
        for (int u = 1; u <= t; u++) {
          double fu = F(z, t, u);
            sum += (log(u)/log(2)) * fu;
        }
    }

    return sum / (double)v;
}
static double G_RHS(double p, int k,int  L, int d, double target)
{
    double q = (1 - p) / (k - 1);
    double RHS = G(p, L, d) + (k - 1) * G(q, L, d);
    return RHS-target;
}
bool bisection_root2(double lo, double hi, double eps, int max_iter,
    int k, int L, int d, double target,
    double* root_out)
{
    if (!root_out || !isfinite(lo) || !isfinite(hi) || !isfinite(eps) || eps <= 0.0)
        return false;
    if (max_iter <= 0) return false;

    if (lo > hi) { double t = lo; lo = hi; hi = t; }

    double flo = G_RHS(lo, k, L,d, target);
    double fhi = G_RHS(hi, k, L,d, target);
    if (!isfinite(flo) || !isfinite(fhi)) return false;


    if (fabs(flo) <= eps) { *root_out = lo; return true; }
    if (fabs(fhi) <= eps) { *root_out = hi; return true; }

    //중간값 정리에의해 0보다 작아야 근이 존재
    if (flo * fhi > 0.0) {
        printf("중간값정리만족x");
        return false;
    }

    for (int i = 0; i < max_iter; i++) {
        double mid = 0.5 * (lo + hi);
        double fmid = G_RHS(mid, k, L,d, target);
        if (!isfinite(fmid)) return false;

        if (fabs(fmid) <= eps || 0.5 * (hi - lo) <= eps) {
            *root_out = mid;
            return true;
        }

        if (flo * fmid < 0.0) {
            hi = mid;
            fhi = fmid;
        }
        else {
            lo = mid;
            flo = fmid;
        }
    }
    *root_out = 0.5 * (lo + hi);
    return true;
}
void Compression(int s[], int L, int k, double* H)
{
    int* dict = (int*)malloc(sizeof(int) * k);// 1)
    if (!dict) { *H = NAN; return; }
    int d = 1000;
    for (int i = 0; i < d; i++) {
        int si = s[i];
        if (0 <= si && si < k) dict[si] = i;
    }
    int v = L - d;
    int* D = (int*)malloc(sizeof(int) * v);
    if (!D) { *H = NAN; free(dict); return; }
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
    bool ok = bisection_root2(1.0/16.0,1 , 1e-4, 1,
        k,L,d, X_prime,
        &p);
    if (ok)
    {
        printf("성공\n");
        printf("근 =%f\n", p);
        *H = -log(p) / log(2);
    }
    else {
        printf("실패\n");
    }

    printf("sigma(편차) = %f\n", sigma);
    printf("X_prime(X의 99% 신뢰구간 하한) = %f\n", X_prime);
    printf("X_prime(X의 99% 신뢰구간 하한) = %f\n", X_prime);

    free(D);
    free(dict);

}






double g_func(int i)
{
    if (i <= 1) {

        return 0.0;
    }
    double sum = 0.0;
    for (int z = 1; z < i; z++)
        sum += (1.0 / z);

    return sum / log(2.0);
}

void entropy_test(int b[], int N, int T, int Q, int K, double* H)

{
    int* b1 = (int*)calloc(N + 1, sizeof(int));
    if (!b1) return;
    for (int i = 1; i <= N; i++) {
        b1[i] = b[i - 1] & 1;  // bit로 강제
    }
    int total_words = Q + K;
    int* w = (int*)calloc(total_words+1, sizeof(int));  
    int* D = (int*)calloc(K+   1, sizeof(int));               
    if (!w || !D) { free(w); free(D); if (H) *H = NAN; return; }

    // w[1]부터 채우기: i=1..total_words
    // 문서(MSB)
    for (int i = 1; i <=total_words; i++) {
        int word = 0;
        for (int j = 0; j < T; j++) {
            word = (word << 1) | (b1[(i-1) * T +1+ j]);
        }
        w[i] = word;
    }  // 1)

    /* LSB
    for (int i = 1; i <= total_words; i++) {
        int word = 0;
        for (int j = 0; j < T; j++) {
            word |= (b1[(i - 1) * T + 1 + j] & 1) << j;
        }
        w[i] = word;
    }*/


    // m은 1-based 기준: m = Q+1 .. Q+K
    for (int m = Q+1 ; m <= Q + K; m++) {
        int found = 0;
        int distance = 0;

        // i는 1..m-1 (가까운 과거부터 찾는다면 i=1부터)
        for (int i = 1; i <= m-1 ; i++) {
            if (w[m] == w[m - i]) {
                found = 1;
                distance = i;
                break;
            }
        }
        int idx = m - Q; 
        if (found)
            D[idx] = distance;   // 3.2)
        else
            D[idx] = m;          // 3.1) 
      
    }
    for (int i = 1; i <= K; i++)
    {
        printf("%d  ", D[i]);
    }
    double sum = 0.0;
    int cnt = 0;
    for (int j = 1; j <=K; j++) {
        //if (D[j] == 1) cnt++;
        sum += g_func(D[j]);
    }  
 
    //printf("cnt=%d\n", cnt);
    *H = sum / K;
    printf("sum:%f\n", sum);

    free(w);
    free(D);
    free(b1);
}

void Mutual_Information(int b[], int N, int T, int Q, int K, double* H) //T, Q=4, K=50
{ 
    *H = 0;
    // N = (Q + K) * T, 216 = (50+4)*4
    int L = N / T; //54
    int* w = (int*)malloc(L * sizeof(int));
    for (int i = 0; i < L; i++) {
        int temp = 0;
        for (int j = 0; j < T; j++) {
            temp = (temp << 1) | b[i * T + j];
        }
        w[i] = temp;
    }// 1)
    int* E = (int*)calloc(pow(2, T), sizeof(int));
    for (int i = 0; i < Q; i++)
    {
        E[w[i]]++;
    }// 3)
    printf("단계 3\n");
    for (int i = 0; i < 1 << T; i++)
    {
        printf("%d   ", E[i]);
    }
    int i = Q;
    int count = 0;
    while (i < K + Q)
    {
        count++;
        E[w[i]]++;
        *H += ((log(i)) / ((log(2)) * E[w[i]]));
        i++;
    }
    
    printf("\n단계 5\n");
    for (int i = 0; i < 1 << T; i++)
    {
        printf("%d   ", E[i]);
    }
    printf("\n");
    printf("\n수행이후의 H= %f\n", *H);
    *H /= count;
}

double double_correlation(uint8_t col1[256], uint8_t col2[256]) {
    double mean1 = 0, mean2 = 0;
    for (int i = 0; i < 256; i++) {
        mean1 += col1[i];
        mean2 += col2[i];
    }
    mean1 /= 256;
    mean2 /= 256;

    double num = 0, den1 = 0, den2 = 0;
    for (int i = 0; i < 256; i++) {
        double d1 = col1[i] - mean1;
        double d2 = col2[i] - mean2;
        num += d1 * d2;
        den1 += d1 * d1;
        den2 += d2 * d2;
    }
    if (den1 == 0 || den2 == 0) return 0.0;
    return num / sqrt(den1 * den2);
}

// 2.3.2 Shannon 엔트로피 계산
double shannon_entropy(const uint8_t col[256]) {
    int freq[256] = { 0 };

    // 빈도 세기
    for (int i = 0; i < 256; i++) {
        freq[col[i]]++;
    }

    // H = - sum p log2 p
    double H = 0.0;
    for (int t = 0; t < 256; t++) {
        if (freq[t] > 0) {
            double prob = (double)freq[t] / 256.0;
            H -=prob * log2(prob);
        }
    }
    return H;
}


// 2.3.3 Min-엔트로피 계산
double min_entropy(uint8_t col[256]) {
    int freq[256] = { 0 };
    int maxf = 0;//최대 빈도
    for (int i = 0; i < 256; i++) {
        int v = col[i];
        if (++freq[v] > maxf) maxf = freq[v];
    }
    double pmax = (double)maxf / 256;
    return -log2(pmax);
}

 
void Byte_Correlation(uint8_t X[256][8], int min_distub, int use_shannon, double* e_f)
{
    int valid_cols[8] = { 0 };
    uint8_t Xprime[256][8] = { 0 };
    int k = 0; // 살아남은 열 개수
    for (int j = 0; j < 8; j++) {
        int seen[256] = { 0 };
        int distinct = 0;
        for (int i = 0; i < 256; i++) {
            if (!seen[X[i][j]]) {
                seen[X[i][j]] = 1;
                distinct++;
            }

        }// 1.1)
        if (distinct >= min_distub) {
            for (int i = 0; i < 256; i++) {
                Xprime[i][k] = X[i][j];
            }
            k++;
        }// 1.2)
        printf("X_%d일때 취하는 값의 종류%d\n", j + 1, distinct);
    }

    if (k == 0) {
        *e_f = 0.0;
        return;
    }
    // k열로 압축된 Xprime 이후…

    // (2) Pearson 상관계수 먼저
    double C[8][8] = { 0.0 };

    for (int i = 0; i < k; i++) {
        uint8_t ci[256];
        for (int r = 0; r < 256; r++) {
            ci[r] = Xprime[r][i];
        }

        for (int j = 0; j < k; j++) {
            if (i == j) {
                C[i][j] = 0.0;
            }
            else {
                uint8_t cj[256];
                for (int r = 0; r < 256; r++) {
                    cj[r] = Xprime[r][j];
                }
                C[i][j] = double_correlation(ci, cj);
            }
        }
    }
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            printf("c[%d][%d]=%f    ", i + 1, j + 1, C[i][j]);
        }
        printf("\n");
    }
    double e[8];

    for (int j = 0; j < k; j++) {
        uint8_t col[256];
        for (int r = 0; r < 256; r++) {
            col[r] = Xprime[r][j];
        }
        e[j] = shannon_entropy(col);
    }
    printf("\n");
    printf("Shannon 엔트로피(단계3.1)\n");
    for (int i = 0; i < k; i++)
    {
        printf("e[%d] =%f  ", i + 1, e[i]);

    }
    printf("\n");
    printf("재산정 Shannon 엔트로피(단계4.2)\n");
    double e1[8];
    for (int i = 0; i < k; i++) {
        double sum_abs = 0.0;

        for (int j = 0; j < k; j++) {
            if (i == j) continue;
            sum_abs += fabs(C[i][j]);          // |C_{i,j}| 누적
        }

        double Ci_avg = sum_abs / (double)k;   // (1/k) * C'_i  (그림 식 그대로면 /k)
        e1[i] = e[i] * (1.0 - Ci_avg);         // e'_i = e_i * (1 - (1/k)C'_i)
    }
    for (int i = 0; i < k; i++)
    {
        printf("e`[%d] =%f  ", i + 1, e1[i]);
    }
    double sum1 = 0;
    for (int i = 0; i < k; i++)
    {
        sum1 += e1[i];
    }
    printf("\n");
    printf("최종 Shannon 엔트로피(단계5)\n");
    printf("e_f=%f", sum1);











    e[0] = 0.0;
    e[1] = 0.0;
    e[2] = 0.0;
    for (int j = 0; j < k; j++) {
        uint8_t col[256];
        for (int r = 0; r < 256; r++) {
            col[r] = Xprime[r][j];
        }
        e[j] = min_entropy(col);
    }
    printf("\n");

    printf("\n");
    printf("\n");
    printf("\n");
    printf("최소 엔트로피(단계3.2)\n");
    for (int i = 0; i < k; i++)
    {
        printf("e[%d] =%f  ", i + 1, e[i]);
    }
  
    // (4) 보정 e'_i = e_i * (1 - 평균 |C_ij|)
    double e2[8];
    for (int i = 0; i < k; i++) {
        double sum_abs = 0.0;

        for (int j = 0; j < k; j++) {
            if (i == j) continue;              
            sum_abs += fabs(C[i][j]);          // |C_{i,j}| 누적
        }

        double Ci_avg = sum_abs / (double)k;   // (1/k) * C'_i  (그림 식 그대로면 /k)
        e2[i] = e[i] * (1.0 - Ci_avg);         // e'_i = e_i * (1 - (1/k)C'_i)
    }
    printf("\n");
    printf("재산정 Shannon 엔트로피(단계4.2)\n");
    for (int i = 0; i < k; i++)
    {
        printf("e`[%d] =%f  ", i + 1, e2[i]);
    }
    double sum2 = 0;
    for (int i = 0; i < k; i++)
    {
        sum2 += e2[i];
    }
    printf("\n");
    printf("최종 최소 엔트로피(단계5)\n");
    printf("e_f=%f", sum2);

}