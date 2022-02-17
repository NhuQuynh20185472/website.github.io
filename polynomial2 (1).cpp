#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#define INF 999.0

struct PolynomialFunction {
    double coeff[100];
    int deg;
} func;
typedef struct PolynomialFunction polynomial;

double roots[10][10];

// Calculate value of f at x
double calc(polynomial func, double x) {
    double value = 0;
    int i;
    for (i = 0; i <= func.deg; i++) {
        value += func.coeff[i] * pow(x, i);
    }
    return value;
}

// Calculate value of df at x
polynomial diff(polynomial func) {
    polynomial dfunc;
    int i;
    dfunc.deg = func.deg - 1;
    for (i = 0; i <= dfunc.deg; i++) {
        dfunc.coeff[i] = func.coeff[i + 1] * (i + 1);
    }
    return dfunc;
}

// Merge two sorted array
void merge(double arr1[], double arr2[], double arr3[]) {
    int i, j, k, m, n;
    m = (int) arr1[0];
    n = (int) arr2[0];
    i = j = k = 1;
    arr3[0] = (double) m + n;
    for (i = 1; i <= m && j <= n;) {
        if (arr1[i] < arr2[j]) {
            if ((arr3[k - 1] != arr1[i]) || (k == 1)) {
                arr3[k] = arr1[i];
                k++;
            } else {
                arr3[0] -= 1;
            }
            i++;
        } else {
            if ((arr3[k - 1] != arr2[j]) || (k == 1)) {
                arr3[k] = arr2[j];
                k++;
            } else {
                arr3[0] -= 1;
            }
            j++;
        }
    }
    while (i <= m) {
        arr3[k] = arr1[i];
        k++;
        i++;
    }
    while (j <= n) {
        arr3[k] = arr2[j];
        k++;
        j++;
    }
}

//Sign of a number
double sign(double x) {
    if (x >= 0) return 1;
    else return -1;
}

// Find range that contains all roots
void rootsRange(polynomial func, double a, double b, double range[]) {
    int i;
    double max_coeff = fabs(func.coeff[0]);
    for (i = 1; i <= func.deg; i++) {
        if (max_coeff < fabs(func.coeff[i])) {
            max_coeff = fabs(func.coeff[i]);
        }
    }
    range[0] = -1 - max_coeff / fabs(func.coeff[func.deg]);
    if (a > range[0]) {
        range[0] = a;
    }
    range[1] = 1 + max_coeff / fabs(func.coeff[func.deg]);
    if (b < range[1]) {
        range[1] = b;
    }
}

//Newton Raphson method 
double solve(polynomial func, double a, double b, double e) {
    double s, x, m1, m2, err;
    polynomial dfunc = diff(func), ddfunc = diff(dfunc);
    int k = 0;
    if (fabs(calc(dfunc, a)) > fabs(calc(dfunc, b))) {
        m1 = calc(dfunc, b);
    } else {
        m1 = calc(dfunc, a);
    }
    if (fabs(calc(ddfunc, a)) > fabs(calc(ddfunc, b))) {
        m2 = calc(ddfunc, a);
    } else {
        m2 = calc(ddfunc, b);
    }
    err = sqrt(fabs(2.0 * m1 * e / m2));
    printf("\nError = %lf", err);
    // Choose Fourier point
    if (sign(calc(func, a)) == sign(calc(ddfunc, a))) x = a;
    else x = b;
    do {
        s = x;
        x = x - calc(func, x) / calc(dfunc, x);
        printf("\n    Iteration %d:", k + 1);
        printf("\n    x = %lf", x);
        k++;
    } while (fabs(x - s) >= err);
    printf("\n   Total iterations: %d", k);
    printf("\n   x = %lf", x);
    return (x);
}

void nrsolve(polynomial func, double a, double b, double e) {
    if (func.deg == 1) {
        roots[1][0] = 1;
        roots[1][1] = -func.coeff[0] / func.coeff[1];
        int i;
        printf("\nExecute current function: ");
        for (i = 0; i <= func.deg; i++) {
            printf(" %lfx^%d ", func.coeff[i], i);
        }
    } else
    if (func.deg == 0) {
        roots[0][0] = 0;
        int i;
        printf("\nExecute current function: ");
        for (i = 0; i <= func.deg; i++) {
            printf(" %lfx^%d ", func.coeff[i], i);
        }
    } else {
        nrsolve(diff(func), a, b, e);
        int i, index = 1, n;
        printf("\n");
        printf("\nExecute current function: ");
        for (i = 0; i <= func.deg; i++) {
            printf(" %lfx^%d ", func.coeff[i], i);
        }
        n = func.deg;
        double ab[3] = {2, a, b};
        double sp[(int) roots[n - 1][0] + (int) roots[n - 2][0] + 1];
        double point[(int) sp[0] + 3];
        merge(roots[n - 1], roots[n - 2], sp);
        merge(sp, ab, point);
        printf("\nRange: ");
        for (i = 1; i <= (int) point[0]; i++) {
            printf(" %lf ", point[i]);
        }
        for (i = 1; i < (int) point[0]; i++) {
            if (sign(calc(func, point[i] + e / 2)) != sign(calc(func, point[i + 1] - e / 2))) {
                printf("\nDivision range from %lf to %lf", point[i] + e / 2, point[i + 1] - e / 2);
                roots[n][index++] = solve(func, point[i] + e / 2, point[i + 1] - e / 2, e);
            }
            // speacial point is root ???
            if (fabs(calc(func, point[i])) <= fabs(calc(diff(func), point[i]) * e)) {
                roots[n][index++] = point[i];
            }
        }
        if (fabs(calc(func, point[(int) point[0]])) <= fabs(calc(diff(func), point[(int) point[0]]) * e)) {
            roots[n][index++] = point[(int) point[0]];
        }
        roots[n][0] = index - 1;
        printf("\nRoot Set: ");
        for (i = 1; i <= (int) roots[n][0]; i++) {
            printf("\n   x%d = %lf ", i, roots[n][i]);
        }
    }
}

main() {
    int i;
    double D[2], e;
    FILE * f;
    f = fopen("input.txt", "r");
    while (true) {
        fscanf(f, "%lf", & func.coeff[i]);
        if (feof(f)) break;
        else i++;
    }
    func.deg = i;
    fclose(f);
    printf("Function: ");
    for (i = 0; i <= func.deg; i++) {
        printf(" %lfx^%d ", func.coeff[i], i);
    }
    printf("\nInput error: ");
    scanf("%lf", & e);
    rootsRange(func, -INF, INF, D);
    printf("\nRange [%lf, %lf]", D[0], D[1]);
    nrsolve(func, D[0], D[1], e);
    getch();
}
