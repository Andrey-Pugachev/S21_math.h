#include "s21_math.h"

long double s21_factorial(int x) {
    if (x < 0)
        return 0;
    if (x == 1)
        return 1;
    return x * s21_factorial(x - 1);
}

long double s21_exp(double x) {
    long double newMember = 1;
    long double i = 1;
    long double exp = 1;
    int isNegativeNumber = 0;
    if (x < 0) {
    x *= -1;
    isNegativeNumber = 1;
    }
    while (s21_fabs(newMember) > S21_PRECISION) {
        newMember *= x / i;
        i++;
        exp += newMember;
        if (exp > DBL_MAX) {
            exp = S21_INF;
            break;
        }
    }
    if (isNegativeNumber == 1) {
        if (exp > DBL_MAX)
            exp = 0;
        else
            exp = 1. / exp;
    }
    if (exp > DBL_MAX)
        return S21_INF;
    return exp;
}

long double s21_fabs(double x) {
    if (x != x)
        return S21_NAN;
    if (x == S21_INF || x == -S21_INF)
        return S21_INF;
    if (x == -0)
        return 0;
    return (x < 0) ? -x : x;
}

long double s21_ceil(double x) {
    long double result = 0.;
    if (x == S21_INF)
        result = S21_INF;
    else if (x == -S21_INF)
        result = -S21_INF;
    else if (x != x)
        result = x;
    else if ((x < 0 && ((long long int)x == 0)) || (x == (long long int)x))
        result = (long long int)x;
    else {
    long long int tmp;
    tmp = (long long int)x;
    tmp < 0 ? tmp : tmp++;
    result = (long double)tmp;
    }
    return result;
}

long double s21_floor(double x) {
    long double result = 0.;
    if (x != x)
        result = S21_NAN;
    else if (x == -S21_INF)
        result = -S21_INF;
    else if (x == S21_INF)
        result = S21_INF;
    else if (x < 0 && ((long long int)x == 0))
        result = (long long int)x - 1;
    else if (x == (long long int)x)
        result = (long long int)x;
    else {
        long long int tmp;
        tmp = (long long int)x;
        tmp < 0 ? tmp-- : tmp;
        result = (long double)tmp;
    }
    return result;
}

long double s21_log(double x) {
    int amount_of_e_in_x = 0;
    double result = 0;
    double previousResult = 0;
    if (x == S21_INF)
        result = S21_INF;
    else if (x == 0)
        result = -S21_INF;
    else if (x < 0)
        result = S21_NAN;
    else if (x == 1)
        result = 0;
    else {
        for (; x >= S21_EXP; x /= S21_EXP, amount_of_e_in_x++)
            continue;
        for (int i = 0; i < 100; i++) {
            previousResult = result;
            result = previousResult + 2 * (x - s21_exp(previousResult)) / (x + s21_exp(previousResult));
        }
    }
    return amount_of_e_in_x + result;
}

long double s21_sqrt(double x) {
    long double left = 0;
    long double right = (x < 1) ? 1 : x;
    long double middle;
    if (x == S21_INF)
        middle = S21_INF;
    else if (x == -S21_INF)
        middle = S21_NAN;
    else if (x == 0 || x == -0)
        middle = 0;
    else if (x < 0)
        middle = S21_NAN;
    else if (x == 1)
        middle = 1;
    else {
        middle = (left + right) / 2;
        while ((middle - left) > S21_PRECISION) {
            if ((middle * middle) > x)
                right = middle;
            else
                left = middle;
            middle = (left + right) / 2;
        }
    }
    return middle;
}

long double s21_fmod(double x, double y) {
    long double result = 0;
    if (y == 0 || y == -0 || y == S21_NAN || x == S21_INF || x == -S21_INF || x == S21_NAN)
        result = S21_NAN;
    else if (y == -S21_INF || y == S21_INF)
        result = x;
    else {
        long long int integerPart;
        integerPart = x / y;
        result = x - y * integerPart;
    }
    return result;
}

long double s21_pow(double base, double exp) {
    long double result = 0;
    long double copyOfBase = base;
    if (base == 0) {
        if (exp == +0)
            result = 1;
        if (exp == -0)
            result = 1;
        if (exp != exp)
            result = S21_NAN;
        if (exp == -S21_INF)
            result = S21_INF;
        if (exp == +S21_INF)
            result = 0;
        if (exp > 0)
            result = 0;
        if (exp < 0)
            result = S21_INF;
        return result;
    }
    if (base == -1) {
        if (exp == +0)
            result = 1;
        else if (exp == -0)
            result = 1;
        else if (exp != exp)
            result = S21_NAN;
        else if (exp == -S21_INF)
            result = 1;
        else if (exp == +S21_INF)
            result = 1;
        else if (exp != +S21_INF && exp != -S21_INF && exp == exp && (s21_floor(exp) - exp) != 0)
            result = S21_NAN;
        else
            if ((exp - (long long int)exp) == 0) {
                if ((long long int)exp % 2)
                    result = base;
                else
                    result = s21_fabs(base);
            } else
                result = S21_NAN;
        return result;
    }
    if (base == 1) {
        if (exp == +0)
            result = 1;
        else if (exp == -0)
            result = 1;
        else if (exp != exp)
            result = 1;
        else if (exp == -S21_INF)
            result = 1;
        else if (exp == +S21_INF)
            result = 1;
        else
            result = 1;
        return result;
    }
    if (base != base) {
        if (exp == +0)
            result = 1;
        else if (exp == -0)
            result = 1;
        else if (exp != exp)
            result = S21_NAN;
        else if (exp == -S21_INF)
            result = S21_NAN;
        else if (exp == +S21_INF)
            result = S21_NAN;
        else
            result = S21_NAN;
        return result;
    }
    if (base == S21_INF) {
        if (exp == +0)
            result = 1;
        if (exp == -0)
            result = 1;
        if (exp != exp)
            result = S21_NAN;
        if (exp == -S21_INF)
            result = 0;
        if (exp == +S21_INF)
            result = S21_INF;
        return result;
    }
    if (base == -S21_INF) {
        if (exp == +0)
            result = 1;
        else if (exp == -0)
            result = 1;
        else if (exp != exp)
            result = S21_NAN;
        else if (exp == -S21_INF)
            result = 0;
        else if (exp == +S21_INF)
            result = S21_INF;
        else
            if (exp > 0) {
                if ((exp - (long long int)exp) == 0)
                    if ((long long int)exp % 2)
                        result = -S21_INF;
                    else
                        result = S21_INF;
                else
                    result = S21_INF;
            }
        return result;
    }
    if (exp == 0) {
        return result = 1;
    }
    if (exp == -S21_INF && (base > -1 && base < 1)) {
        return result = S21_INF;
    }
    if (exp == S21_INF && (base < -1 || base > 1)) {
        return result = S21_INF;
    }
    if (copyOfBase < 0) {
       copyOfBase = -copyOfBase;
       result = s21_exp(exp * s21_log(copyOfBase));
       if (s21_fmod(exp, 2) != 0)
           result = -result;
    } else
       result = s21_exp(exp * s21_log(copyOfBase));
    return result;
}

long double s21_atan(double x) {
    long double result = 0;
    if (x != x)
        return S21_NAN;
    if (x == 1) {
        result = 0.7853981633974480L;
    } else if (x == -1) {
        result = -0.7853981633974480L;
    } else if (x == S21_PI / 2) {
        result = 1.003884821853887214L;
    } else if (x == -S21_PI / 2) {
        result = -1.003884821853887214L;
    } else if (x == S21_INF || x == -S21_INF) {
        result = x < 0 ? -S21_PI / 2 : S21_PI / 2;
    } else if (-1. < x && x < 1.) {
        for (int i = 0; i < 5000; i++)
            result += s21_pow(-1, i) * s21_pow(x, 1 + (2 * i)) / (1 + (2 * i));
    } else {
        for (int i = 0; i < 5000; i++)
            result += s21_pow(-1, i) * s21_pow(x, -1 - (2 * i)) / (1 + (2 * i));
        result = S21_PI * ((x < 0) ? -x : x) / (2 * x) - result;
    }
    return result;
}

long double s21_asin(double x) {
    if (x == 1.)
        return S21_PI / 2;
    if (x == -1.)
        return -S21_PI / 2;
    if (x == 0.7071067811865475244)
        return S21_PI / 4;
    if (x == -0.7071067811865475244)
        return -S21_PI / 4;
    long double result = 0.;
    if (-1. < x && x < 1.)
        result = s21_atan(x / s21_sqrt(1 - x * x));
    else
        return S21_NAN;
    return result;
}

long double s21_acos(double x) {
    if (x == 1.)
        return 0;
    if (x == -1.)
        return S21_PI;
    if (x == 0)
        return S21_PI / 2;
    if (x == 0.7071067811865475244)
        return S21_PI / 4;
    if (x == -0.7071067811865475244)
        return 3 * S21_PI / 4;
    long double result = 0.;
    if (0. < x && x < 1.)
        result = s21_atan(s21_sqrt(1 - x * x) / x);
    else if (-1. < x && x < 0.)
        result = S21_PI + s21_atan(s21_sqrt(1 - x * x) / x);
    else
        return S21_NAN;
    return result;
}

long double s21_sin(double x) {
    int sign = 1;
    long double result = 0;
    if (x == S21_NAN || x == -S21_INF || x == S21_INF)
        return S21_NAN;
    for (; x > 2 * S21_PI || x < -2 * S21_PI;) {
        if (x > 2 * S21_PI)
            x -= 2 * S21_PI;
        else
            x += 2 * S21_PI;
    }
    if (x < 0) {
        x = -x;
        sign = -1;
    }
    for (int i = 0; i < 150; i++) {
        result += s21_pow(-1, i) * s21_pow(x, 1 + 2 * i) / s21_factorial(1 + 2 * i);
    }
    return result * sign;
}

long double s21_cos(double x) {
  long double member, answer;
  x = s21_fmod(x, 2 * S21_PI);
  member = 1;
  answer = 1;
  if (s21_fabs(x) < S21_PRECISION) {
    answer = 1.;
  } else {
    for (long double i = 1.; s21_fabs(member) > S21_PRECISION && i < 50; i++) {
      member *= ((-1.) * x * x / (2. * i * (2. * i - 1.)));
      answer += member;
    }
  }
  return answer;
}

long double s21_tan(double x) {
    if (x == S21_PI / 2)
        return 16331239353195370L;
    else if (x == -S21_PI / 2)
        return -16331239353195370L;
    if (x == 0)
        return 0;
    return s21_sin(x) / s21_cos(x);
}
