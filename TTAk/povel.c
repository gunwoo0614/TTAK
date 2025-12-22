/* polevl.c : 다항식 평가 */
double polevl(double x, const double coef[], int N) {
    double ans;
    int i;
    ans = coef[0];
    for (i = 1; i <= N; i++)
        ans = ans * x + coef[i];
    return ans;
}

double p1evl(double x, const double coef[], int N) {
    double ans;
    int i;
    ans = x + coef[0];
    for (i = 1; i < N; i++)
        ans = ans * x + coef[i];
    return ans;
}
