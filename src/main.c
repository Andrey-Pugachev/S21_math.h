#include <stdio.h>
#include <math.h>
#include "s21_math.h"

int main() {
    double x;
    scanf("%lf", &x);
    printf("%Lf    s21\n", s21_acos(x));
    printf("%f    math\n", acos(x));
    return 0;
}
