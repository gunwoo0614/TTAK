/* gamma.c : 감마 함수와 로그 감마 근사 */
#include "mconf1.h"

/* log gamma 함수 */
double lgam(double x) {
    return lgamma(x); // C 표준 math.h 함수 이용
}

/* gamma 함수 */
double gammafn(double x) {
    return tgamma(x); // C 표준 math.h 함수 이용
}
