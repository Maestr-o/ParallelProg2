#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#define MAX_THREADS 8

bool is_prime(long num); // проверка на простоту
void generate_combinations_util(long x, long n, long current, long *combination, long index,
                                long **result, long *num_combinations);             // рекурсивная функция генерации комбинации простых чисел
long **generate_combinations(long x, long start, long end, long *num_combinations); // обертка для рекурсивной функции генерации
void free_combinations(long **combinations, long num_combinations);                 // освобождение памяти
long max(long *a, long sz);                                                         // поиск максимального значения
void run_seq(long n, long x, long exp);                                             // запуск последовательного вычисления
void run_par(long n, long x, long exp);                                             // запуск параллельного вычисления
long calc_seq(long n, long x, long exp);                                            // последовательное вычисление
long calc_par(long n, long x, long exp);                                            // параллельное вычисление

int main(int argc, char *argv[])
{
    long x;   // количество чисел в комбинации
    long n;   // максимальное значение числа
    long exp; // степень

#ifdef parallel
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0 && size > MAX_THREADS)
    {
        printf("Limit of threads = %d\n", MAX_THREADS);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    if (rank == 0)
    {
#endif
        printf("Enter N: ");
        if (!scanf("%ld", &n))
        {
            printf("Input error\n");
#ifdef parallel
            MPI_Finalize();
#endif
            return 0;
        }
        printf("Enter quantity of prime numbers: ");
        if (!scanf("%ld", &x))
        {
            printf("Input error\n");
#ifdef parallel
            MPI_Finalize();
#endif
            return 0;
        }
        printf("Enter exponent: ");
        if (!scanf("%ld", &exp))
        {
            printf("Input error\n");
#ifdef parallel
            MPI_Finalize();
#endif
            return 0;
        }
#ifdef parallel
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&exp, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    run_par(n, x, exp);
    MPI_Finalize();
#endif
#ifdef seq
    run_seq(n, x, exp);
#endif
    return 0;
}

void run_seq(long n, long x, long exp)
{
    long result = 0;
    double start, end;
    start = omp_get_wtime();
    result = calc_seq(n, x, exp);
    end = omp_get_wtime();
    printf("Seq: time = %.3lf sec, result = %ld\n", (end - start), result);
}

void run_par(long n, long x, long exp)
{
    int rank, size;
    long result = 0;
    double start = 0.0, end = 0.0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == 0)
        start = omp_get_wtime();
    result = calc_par(n, x, exp);
    if (rank == 0)
    {
        end = omp_get_wtime();
        printf("Parallel: time = %.3lf sec, result = %ld\n", (end - start), result);
    }
}

bool is_prime(long num)
{
    if (num < 2)
        return false;
    long i;
    for (i = 2; i * 2 <= num; i++)
    {
        if (num % i == 0)
            return false;
    }
    return true;
}

long **generate_combinations(long x, long start, long end, long *num_combinations)
{
    *num_combinations = 0;
    long max_combinations = 1, i;
    for (i = 0; i < x; i++)
        max_combinations *= (end - start - 1);

    long **result = (long **)malloc(max_combinations * sizeof(long *));
    if (result == NULL)
    {
        fprintf(stderr, "Error allocate memory\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < max_combinations; i++)
        result[i] = (long *)malloc(x * sizeof(long));

    long *combination = (long *)malloc(x * sizeof(long));
    generate_combinations_util(x, end - start, start, combination, start, result, num_combinations);
    free(combination);

    return result;
}

void generate_combinations_util(long x, long n, long current, long *combination, long index, long **result, long *num_combinations)
{
    long i;
    if (index == x)
    {
        for (i = 0; i < x; i++)
            result[*num_combinations][i] = combination[i];
        (*num_combinations)++;
        return;
    }
    for (i = current; i <= n; i++)
    {
        if (is_prime(i))
        {
            combination[index] = i;
            generate_combinations_util(x, n, i, combination, index + 1, result, num_combinations);
        }
    }
}

void free_combinations(long **combinations, long num_combinations)
{
    long i;
    for (i = 0; i < num_combinations; i++)
        free(combinations[i]);
    free(combinations);
}

long max(long *a, long sz)
{
    long i = 1, m = a[0];
    while (i < sz)
    {
        if (a[i] > m)
            m = a[i];
        i++;
    }
    return m;
}

long calc_seq(long n, long x, long exp)
{
    long num_combinations, ans;
    long **combinations = generate_combinations(x, 0, n, &num_combinations);
    long *results = calloc(sizeof(long), num_combinations);
    long i;
    for (i = 0; i < num_combinations; i++)
    {
        unsigned long j;
        for (j = 0; j <= sizeof(combinations[i]) / sizeof(long); j++)
        {
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

long calc_par(long n, long x, long exp)
{
    int rank, size;
    long result;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunk_size = n / size;
    int start_index = rank * chunk_size;
    int end_index = start_index + chunk_size;

    long num_combinations;
    long **combinations = generate_combinations(x, start_index, end_index, &num_combinations);
    long *results = calloc(sizeof(long), num_combinations);
    long i;
    for (i = 0; i < num_combinations; i++)
    {
        unsigned long j;
        for (j = 0; j <= sizeof(combinations[i]) / sizeof(long); j++)
            results[i] += pow(combinations[i][j], exp);
        if (results[i] >= n)
            results[i] = 0;
    }
    result = max(results, num_combinations);
    free_combinations(combinations, num_combinations);
    free(results);
    return result;
}
