#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

bool is_prime(long num); // проверка на простоту
void generate_combinations_util(long x, long n, long current, long* combination, long index,
    long** result, long* num_combinations); // рекурсивная функция генерации комбинации простых чисел
long** generate_combinations(long x, long n, long* num_combinations); // обертывает вызов рекурсивной функции
void free_combinations(long** combinations, long num_combinations); // освобождает память
long max(long *a, long sz); // поиск максимального значения
long calc(long n, long x, long exp); // вычисление
void run_seq(long n, long x, long exp); // запуск последовательного вычисления
void run_par(long n, long x, long exp); // запуск параллельного вычисления

int main(int argc, char *argv[]) {
    #ifdef parallel
    MPI_Init(&argc, &argv);
    #endif
    long x; // количество чисел в комбинации
    long n; // максимальное значение числа
    long exp; // степень
    
    printf("Enter N: ");
    if (!scanf("%ld", &n)) {
        printf("Input error\n");
        return 0;
    }
    printf("Enter quantity of prime numbers: ");
    if (!scanf("%ld", &x)) {
        printf("Input error\n");
        return 0;
    }
    printf("Enter exponent: ");
    if (!scanf("%ld", &exp)) {
        printf("Input error\n");
        return 0;
    }
    #ifdef seq
    run_seq(n, x, exp);
    #endif
    #ifdef parallel
    run_par(n, x, exp);
    MPI_Finalize();
    #endif
    return 0;
}

void run_seq(long n, long x, long exp) {
    long result = 0;
    double start, end;
    start = omp_get_wtime();
    result = calc(n, x, exp);
    end = omp_get_wtime();
    printf("Seq: time = %.3lf sec, result = %ld\n", (end - start), result);
}

void run_par(long n, long x, long exp) {
    int rank, size;
    long result = 0;
    double start, end;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&exp, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    start = omp_get_wtime();
    result = calc(n, x, exp);
    end = omp_get_wtime();
    printf("Parallel: time = %.3lf sec, result = %ld\n", (end - start), result);
}

bool is_prime(long num) {
    if (num < 2) return false;
    long i;
    for (i = 2; i * 2 <= num; i++) {
        if (num % i == 0) return false;
    }
    return true;
}

void generate_combinations_util(long x, long n, long current, long* combination, long index, long** result, long* num_combinations) {
    long i;
    if (index == x) {
        for (i = 0; i < x; i++) {
            result[*num_combinations][i] = combination[i];
        }
        (*num_combinations)++;
        return;
    }
    for (i = current; i <= n; i++) {
        if (is_prime(i)) {
            combination[index] = i;
            generate_combinations_util(x, n, i, combination, index + 1, result, num_combinations);
        }
    }
}

long** generate_combinations(long x, long n, long* num_combinations) {
    *num_combinations = 0;
    long max_combinations = 1, i;
    for (i = 0; i < x; i++) {
        max_combinations *= (n - 1);
    }

    long** result = (long**)malloc(max_combinations * sizeof(long*));
    for (i = 0; i < max_combinations; i++) {
        result[i] = (long*)malloc(x * sizeof(long));
    }

    long* combination = (long*)malloc(x * sizeof(long));
    generate_combinations_util(x, n, 2, combination, 0, result, num_combinations);
    free(combination);

    return result;
}

void free_combinations(long** combinations, long num_combinations) {
    long i;
    for (i = 0; i < num_combinations; i++) {
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
    return m;
}

long calc(long n, long x, long exp) {
    long num_combinations, ans;
    long** combinations = generate_combinations(x, n, &num_combinations);
    long* results = calloc(sizeof(long), num_combinations);
    long i;
    for (i = 0; i < num_combinations; i++) {
        unsigned long j;
        for (j = 0; j <= sizeof(combinations[i]) / sizeof(long); j++) {
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
