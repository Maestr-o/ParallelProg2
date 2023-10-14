#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

bool is_prime(long num); // проверка на простоту
void generate_combinations_util(int x, int n, int current, int* combination, int index,
    long** result, int* num_combinations); // рекурсивная функция генерации комбинации простых чисел
long** generate_combinations(int x, int n, int* num_combinations); // обертывает вызов рекурсивной функции
void free_combinations(long** combinations, int num_combinations); // освобождает память
long max(long *a, long sz); // поиск максимального значения
long seq_calc(long n, long x, long exp); // последовательное вычисление
void run();

int main() {
    run();
    return 0;
}

void run() {
    long x; // количество чисел в комбинации
    long n; // максимальное значение числа
    long exp; // степень
    long result = 0;

    printf("Enter N: ");
    scanf("%ld", &n);
    printf("Enter quantity of prime numbers: ");
    scanf("%ld", &x);
    printf("Enter exponent: ");
    scanf("%ld", &exp);
    
    double start, end; // последовательное вычисление
    start = omp_get_wtime();
    result = seq_calc(n, x, exp);
    end = omp_get_wtime();
    printf("Seq: time = %.3lf sec, result = %ld", (end - start), result);
}

bool is_prime(long num) {
    if (num < 2) return false;
    for (long i = 2; i * 2 <= num; i++) {
        if (num % i == 0) return false;
    }
    return true;
}

void generate_combinations_util(int x, int n, int current, int* combination, int index, long** result, int* num_combinations) {
    if (index == x) {
        for (int i = 0; i < x; i++) {
            result[*num_combinations][i] = combination[i];
        }
        (*num_combinations)++;
        return;
    }

    for (int i = current; i <= n; i++) {
        if (is_prime(i)) {
            combination[index] = i;
            generate_combinations_util(x, n, i, combination, index + 1, result, num_combinations);
        }
    }
}

long** generate_combinations(int x, int n, int* num_combinations) {
    *num_combinations = 0;
    int max_combinations = 1;
    for (int i = 0; i < x; i++) {
        max_combinations *= (n - 1);
    }

    long** result = (long**)malloc(max_combinations * sizeof(long*));
    for (int i = 0; i < max_combinations; i++) {
        result[i] = (long*)malloc(x * sizeof(long));
    }

    int* combination = (int*)malloc(x * sizeof(int));
    generate_combinations_util(x, n, 2, combination, 0, result, num_combinations);
    free(combination);

    return result;
}

void free_combinations(long** combinations, int num_combinations) {
    for (int i = 0; i < num_combinations; i++) {
        free(combinations[i]);
    }
    free(combinations);
}

long max(long *a, long sz) {
    long i = 1, m = a[0];
    while (i < sz) {
        if (a[i] > m)
            m = a[i];
        i++;
    }
    free(a);
    return m;
}

long seq_calc(long n, long x, long exp) {
    int num_combinations, ans;
    long** combinations = generate_combinations(x, n, &num_combinations);
    long* results = calloc(sizeof(long), num_combinations);
    for (int i = 0; i < num_combinations; i++) {
        for (unsigned int j = 0; j <= sizeof(combinations[i]) / sizeof(long); j++) {
            results[i] += pow(combinations[i][j], exp);
        }
        if (results[i] >= n)
            results[i] = 0;
    }
    ans = max(results, num_combinations);
    free_combinations(combinations, num_combinations);
    free(results);
    return ans;
}
