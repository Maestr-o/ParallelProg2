#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void run();
bool is_prime(unsigned long long int num); // проверка числа на простоту
unsigned long long int* generate_primes(unsigned long long int n); // генерация массива простых чисел
unsigned long long int seq_calc(unsigned long long int num_primes, unsigned long long int exp); // последовательно вычисление

int main() {
    run();
    return 0;
}

void run() {
    unsigned long long int N, num_primes, exp, result = 0;
    printf("Enter N: ");
    scanf("%I64u", &N);
    printf("Enter quantity of prime numbers: ");
    scanf("%I64u", &num_primes);
    printf("Enter exponent: ");
    scanf("%I64u", &exp);

    double start, end; // последовательное вычисление
    start = omp_get_wtime();
    result = seq_calc(num_primes, exp);
    end = omp_get_wtime();
    printf("Seq: time = %.3lf sec, result = ", (end - start));
    if (result < N)
        printf("%I64u\n", result);
    else
        printf("No\n");

    
}

bool is_prime(unsigned long long int num) {
    if (num < 2) return false;
    for (unsigned long long int i = 2; i * i <= num; i++) {
        if (num % i == 0) return false;
    }
    return true;
}

unsigned long long int* generate_primes(unsigned long long int n) {
    unsigned long long int* primes = malloc(n * sizeof(unsigned long long int));
    unsigned long long int count = 0;
    unsigned long long int current_num = 2;
    while (count < n) {
        if (is_prime(current_num)) {
            primes[count] = current_num;
            count++;
        }
        current_num++;
    }
    return primes;
}

unsigned long long int seq_calc(unsigned long long int num_primes, unsigned long long int exp) {
    unsigned long long int result;
    unsigned long long int *primes = generate_primes(num_primes);
    for (unsigned long long int i = 0; i < num_primes; i++) {
        result += pow(primes[i], exp);
    }
    free(primes);
    return result;
}