main.o : main.c s21_math.h
	gcc -Wall -Werror -Wextra -c main.c -o main.o
s21_math.o : s21_math.c s21_math.h
	gcc -Wall -Werror -Wextra -c s21_math.c -o s21_math.o
s21_math : main.o s21_math.o
	gcc main.o s21_math.o -o s21_math
clean :
	rm s21_math *.o
