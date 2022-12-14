#include <stdio.h>
#include <float.h>
#include <math.h>

#define S21_PRECISION 1e-9 //Максимальная точность для функций нана.
#define S21_INF 1.0 / 0.0
#define S21_NAN 0.0 / 0.0
#define S21_EXP 2.7182818284590452353602874713526624
#define S21_PI 3.14159265358979323846264338327950288

//Факториал числа  -  это произведение всех неотрицательных чисел до него и включая его (0! = 1, 1! = 1, 2! = 2, 3! = 6, ...).
long double s21_factorial(int x);
//Возведенное значение e, в заданную степень.
long double s21_exp(double x);
//Вычисляет абсолютное значение числа с плавающей точкой, то есть модуля.
long double s21_fabs(double x);
// //Возвращает ближайшее целое число, не меньшее заданного значения, то есть округляет в большую сторону (1.2 = 2, -1.8 = -1).
long double s21_ceil(double x);
// //Возвращает ближайшее целое число, не превышающее заданное значение, то есть округляет в меньшую сторону (1.8 = 1, -1.2 = -2).
long double s21_floor(double x);
// //Вычисляет натуральный логарифм то есть функцыя возвращает число в которое нужно возвести константу e, что бы получить число x.
long double s21_log(double x);
// //Вычисляет квадратный корень (реализацыя методом дихотомии).
long double s21_sqrt(double x);
// //Остаток операции деления с плавающей точкой (аналог %).
long double s21_fmod(double x, double y);
// //Возводит число в заданную степень
long double s21_pow(double base, double exp);
// //Вычисляет аргтангенс (получает на вход угол заданный в радианах).
long double s21_atan(double x);
// //Вычисляет арксинус (получает на вход угол заданный в радианах).
// //long double s21_asin(double x);
// //Вычисляет арккосинус (получает на вход угол заданный в радианах).
// //long double s21_acos(double x);
// //Вычисляет синус (получает на вход угол заданный в радианах).
long double s21_sin(double x);
// //Вычисляет косинус (получает на вход угол заданный в радианах).
long double s21_cos(double x);
// //Вычисляет тангенс (получает на вход угол заданный в радианах).
// long double s21_tan(double x);

int main() {
    double x;
    scanf("%lf", &x);
    printf("%Lf    s21\n", s21_cos(x));
    printf("%f    math\n", cos(x));
    return 0;
}

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
    if (exp > DBL_MAX) // ????????????????
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

long double s21_ceil(double x) { //как обработать -0?
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
    int amount_of_e_in_x = 0; //Переменная для хранения целой части той степени при возведдении в которую мы получим x
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
    long double right = (x < 1) ? 1 : x; // так как при значениях x меньше еденицы функцыя ведет себя неоднозначно, например корень из 0.1 = 0.316... и это больше чем 0.1, то мы вводим это условие.
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
    if (base == 0) { //Отработка краевых значений при основании 0
        if (exp == +0)
            result = 1;
        if (exp == -0)
            result = 1;
        if (exp != exp) //Единственный вариант сравнить NaN это проверить не равно ли оно само себе?
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
    if (base == -1) { // Отработка краевых значений при основании -1
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
    if (base == 1) { // Отработка краевых значений при основании 1
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
    if (base != base) { //Отработка краевых значений при основании nan
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
    if (base == S21_INF) { //Отработка краевых значений при основании inf
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
    if (base == -S21_INF) { //Отработка краевых значений при основании -inf
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
    if (exp == 0) { //Отработка краевых значений при степени 0
        return result = 1;
    }
    if (exp == -S21_INF && (base > -1 && base < 1)) { //Отработка краевых значений при степени -inf
        return result = S21_INF;
    }
    if (exp == S21_INF && (base < -1 || base > 1)) { //Отработка краевых значений при степени inf
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

// //=============================================================================================================================================================================================

long double s21_atan(double x) {
    /*
    |x| < 1    <=>    -1 < x < 1
    |x| > 1    <=>    x < -1 && x > 1
    В полной окружности 360 градусов    -    6.28319... радиан <=> 2pi.
    */
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

// long double s21_asin(double x) {

//     return;
// }

// long double s21_acos(double x) {

//     return;
// }

long double s21_sin(double x) {
    long double result = 0;
    if (x == S21_NAN || x == -S21_INF || x == S21_INF)
        return S21_NAN;
    for (; x > 2 * S21_PI || x < -2 * S21_PI;) {
        if (x > 2 * S21_PI)
            x -= 2 * S21_PI;
        else
            x += 2 * S21_PI;
    }
    for (int i = 0; i < 150; i++) { //колличество итерацый задаём не более 150 иначе функцыя ломается
        result += s21_pow(-1, i) * s21_pow(x, 1 + 2 * i) / s21_factorial(1 + 2 * i);
    }
    return result;
}

long double s21_cos(double x) {
    long double result = 0;
    for (; x > 2 * S21_PI || x < -2 * S21_PI;) {
        if (x > 2 * S21_PI)
            x -= 2 * S21_PI;
        else
            x += 2 * S21_PI;
    }
    if (x < 0)
    x = -x;
    for (int i = 0; i < 10; i++) {
        result += s21_pow(-1, i) * s21_pow(x, 2 * i) / s21_factorial(2 * i);
    }
    return result;
}

// long double s21_tan(double x) {
//     return s21_sin(x) / s21_cos(x);
// }



/*
1. Что такое INF и NaN?
2. Как проверяете свои функции?
3. Правильно ли я прописал константы (какие еще нужно добавить)?
4. Как побитово выделить целую часть?
5. Как пишутся автотесты?
6. Необходимо пофиксить pow(), что бы нормально работали тригонометрические функции.
*/
