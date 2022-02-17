#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<string.h>
#define INF 999.0

// Polynomial Function
struct PolynomialFunction {
    double coeff[100];
    int deg;
} func;
typedef struct PolynomialFunction polynomial;

// Roots of f and its derivative
double roots[50][50];

// Error of roots
double err;

// Number of iterations
int n;

// Display digits
int digits;

// File
FILE * f;

// String display
char disp[51];

// Display decimals (maximum 18 digits)
char * display(double var);

// Merge two sorted array
void merge(double arr1[], double arr2[], double arr3[]);

// Calculate value of f at x
double calc(polynomial func, double x);

// Derivative of f
polynomial diff(polynomial func);

// Sign of a number
int sign(double x);

// Find range that contains all roots
void rootsRange(polynomial func, double a, double b, double range[]);

// Use bisection method to narrow division range 
void bisection(polynomial func, double a, double b, double distance, double range[]);

// Newton Raphson method: two consecutive approximate according to formula
double nrsolve1(polynomial func, double a, double b, double e);

// Newton Raphson method: two consecutive approximate according to constant error
double nrsolve2(polynomial func, double a, double b, double e);

// Newton Raphson method: iterate according to a constant number of interations
double nriter(polynomial func, double a, double b, int n);

// solve equation f(x) = 0, first step is find the division range then use Newton Raphson method
void solve(polynomial func, double a, double b, double e, int choice);

main() {
    int i = 0, choice;
    double range[2];
    bool continueloop = true;
    while (true) {
        // Input file
        f = fopen("input.txt", "r");
        while (true) {
            fscanf(f, "%lf", & func.coeff[i]);
            if (feof(f)) break;
            else i++;
        }
        func.deg = i;
        fclose(f);
        printf("\n1. Newton Raphson method: two consecutive approximate according to formula");
        printf("\n2. Newton Raphson method: two consecutive approximate according to constant error");
        printf("\n3. Newton Raphson method: iterate according to a constant number of interations");
        printf("\nMake your choice: ");
        scanf("%d", & choice);
        switch (choice) {
        case 1:
            f = fopen("output.txt", "w");
            fprintf(f, "\nNewton Raphson method: two consecutive approximate according to formula");
            printf("\nInput number of decimal digits display (maximum 18): ");
            scanf("%d", & digits);
            if ((digits > 18) || (digits <= 0)) {
                digits = 18;
            }
            rootsRange(func, -INF, INF, range);
            fprintf(f, "\nRange from %s ", display(range[0]));
            fprintf(f, "to %s", display(range[1]));
            printf("\nRange from %s ", display(range[0]));
            printf("to %s", display(range[1]));
            err = 0;
            while (err == 0) {
                printf("\nInput error (none zero): ");
                scanf("%lf", & err);
            }
            solve(func, range[0], range[1], err, 1);
            fclose(f);
            break;
        case 2:
            f = fopen("output.txt", "w");
            fprintf(f, "\nNewton Raphson method: two consecutive approximate according to constant error");
            printf("\nInput number of decimal digits display (maximum 18): ");
            scanf("%d", & digits);
            if ((digits > 18) || (digits <= 0)) {
                digits = 18;
            }
            rootsRange(func, -INF, INF, range);
            fprintf(f, "\nRange from %s ", display(range[0]));
            fprintf(f, "to %s", display(range[1]));
            printf("\nRange from %s ", display(range[0]));
            printf("to %s", display(range[1]));
            err = 0;
            while (err == 0) {
                printf("\nInput error (none zero): ");
                scanf("%lf", & err);
            }
            solve(func, range[0], range[1], err, 2);
            fclose(f);
            break;
        case 3:
            f = fopen("output.txt", "w");
            fprintf(f, "\nNewton Raphson method: iterate according to a constant number of interations");
            printf("\nInput number of decimal digits display (maximum 18): ");
            scanf("%d", & digits);
            if ((digits > 18) || (digits <= 0)) {
                digits = 18;
            }
            rootsRange(func, -INF, INF, range);
            fprintf(f, "\nRange from %s ", display(range[0]));
            fprintf(f, "to %s", display(range[1]));
            printf("\nRange from %s ", display(range[0]));
            printf("to %s", display(range[1]));
            err = 1e-6;
            n = 0;
            while (n <= 0) {
                printf("\nInput number of iterations (positive): ");
                scanf("%d", & n);
            }
            solve(func, range[0], range[1], err, 3);
            fclose(f);
            break;
        default:
            continueloop = false;
            break;
        }
        printf("\n\nDONE!\n");
        if (continueloop == false) {
            break;
        }
    }
    printf("\nExit!");
    getch();
}

// Display decimals (maximum 18 digits)
char * display(double var) {
    int i;
    long long intpart = (long long) var, decpart = (long long)(fabs(var * pow(10, digits) - (double) intpart * pow(10, digits)));
    char istring[25], dstring[25];
    lltoa(intpart, istring, 10);
    lltoa(decpart, dstring, 10);
    if ((var < 0) && (intpart == 0)) {
        strcpy(disp, "-");
        strcat(disp, istring);
    } else {
        strcpy(disp, istring);
    }
    strcat(disp, ".");
    if (strlen(dstring) < digits) {
        for (i = 0; i < digits - strlen(dstring); i++) {
            strcat(disp, "0");
        }
    }
    strcat(disp, dstring);
    return disp;
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

// Calculate value of f at x
double calc(polynomial func, double x) {
    double value = 0;
    int i;
    for (i = 0; i <= func.deg; i++) {
        value += func.coeff[i] * pow(x, i);
    }
    return value;
}

// Derivative of f
polynomial diff(polynomial func) {
    polynomial dfunc;
    int i;
    dfunc.deg = func.deg - 1;
    for (i = 0; i <= dfunc.deg; i++) {
        dfunc.coeff[i] = func.coeff[i + 1] * (i + 1);
    }
    return dfunc;
}

// Sign of a number
int sign(double x) {
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

// Use bisection method to narrow division range 
void bisection(polynomial func, double a, double b, double distance, double range[]) {
    int i;
    double temp;
    fprintf(f, "\n\t\tBisection:");
    printf("\n\t\tBisection:");
    while (b - a > distance) {
        temp = (b + a) / 2;
        if (sign(calc(func, temp)) == sign(calc(func, b))) {
            b = temp;
        } else {
            a = temp;
        }
        fprintf(f, "\n\t\t\t%s ", display(a));
        fprintf(f, "%s", display(b));
        printf("\n\t\t\t%s ", display(a));
        printf("%s", display(b));
    }
    range[0] = a;
    range[1] = b;
}

// Newton Raphson method: two consecutive approximate according to formula
double nrsolve1(polynomial func, double a, double b, double e) {
    double s, x, m1, m2, err;
    polynomial dfunc = diff(func), ddfunc = diff(dfunc);
    int k = 0;
    // m1 = min|df(x)|
    if (fabs(calc(dfunc, a)) > fabs(calc(dfunc, b))) {
        m1 = calc(dfunc, b);
    } else {
        m1 = calc(dfunc, a);
    }
    // m2 = max|ddf(x)|
    if (fabs(calc(ddfunc, a)) > fabs(calc(ddfunc, b))) {
        m2 = calc(ddfunc, a);
    } else {
        m2 = calc(ddfunc, b);
    }
    err = sqrt(fabs(2.0 * m1 * e / m2));
    fprintf(f, "\n\t\tError = %s", display(err));
    printf("\n\t\tError = %s", display(err));
    // Choose Fourier point
    if (sign(calc(func, a)) == sign(calc(ddfunc, a))) {
        x = a;
    } else {
        x = b;
    }
    do {
        s = x;
        x = x - calc(func, x) / calc(dfunc, x);
        fprintf(f, "\n\t\t\tIteration %d: x = %s", k + 1, display(x));
        printf("\n\t\t\tIteration %d: x = %s", k + 1, display(x));
        k++;
    } while (fabs(x - s) >= err);
    fprintf(f, "\n\t\tTotal iterations: %d", k);
    printf("\n\t\tTotal iterations: %d", k);
    fprintf(f, "\n\t\tx = %s", display(x));
    printf("\n\t\tx = %s", display(x));
    return x;
}

// Newton Raphson method: two consecutive approximate according to constant error
double nrsolve2(polynomial func, double a, double b, double e) {
    double s, x, err;
    polynomial dfunc = diff(func), ddfunc = diff(dfunc);
    int k = 0;
    err = e;
    fprintf(f, "\n\t\tError = %s", display(err));
    printf("\n\t\tError = %s", display(err));
    // Choose Fourier point
    if (sign(calc(func, a)) == sign(calc(ddfunc, a))) {
        x = a;
    } else {
        x = b;
    }
    do {
        s = x;
        x = x - calc(func, x) / calc(dfunc, x);
        fprintf(f, "\n\t\t\tIteration %d: x = %s", k + 1, display(x));
        printf("\n\t\t\tIteration %d: x = %s", k + 1, display(x));
        k++;
    } while (fabs(x - s) >= err);
    fprintf(f, "\n\t\tTotal iterations: %d", k);
    printf("\n\t\tTotal iterations: %d", k);
    fprintf(f, "\n\t\tx = %s", display(x));
    printf("\n\t\tx = %s", display(x));
    return x;
}

// Newton Raphson method: iterate according to a constant number of interations
double nriter(polynomial func, double a, double b, int n) {
    polynomial dfunc = diff(func), ddfunc = diff(dfunc);
    int k = 0;
    double s, x;
    // Choose Fourier point
    if (sign(calc(func, a)) == sign(calc(ddfunc, a))) {
        x = a;
    } else {
        x = b;
    }
    do {
        s = x;
        x = x - calc(func, x) / calc(dfunc, x);
        fprintf(f, "\n\t\t\tIteration %d: x = %s", k + 1, display(x));
        printf("\n\t\t\tIteration %d: x = %s", k + 1, display(x));
        k++;
    } while (k < n);
    fprintf(f, "\n\t\tAfter %d iterations: x = %s", k, display(x));
    printf("\n\t\tAfter %d iterations: x = %s", k, display(x));
    return x;
}

// solve equation f(x) = 0, 
// first step is find the division range
// then use Newton Raphson method
void solve(polynomial func, double a, double b, double e, int choice) {
    if (func.deg == 0) {
        roots[0][0] = 0;
    } else
    if (func.deg == 1) {
        roots[1][0] = 1;
        roots[1][1] = -func.coeff[0] / func.coeff[1];
        int i;
        fprintf(f, "\n\nCurrent executing function: ");
        printf("\n\nCurrent executing function: ");
        for (i = 0; i <= func.deg; i++) {
            fprintf(f, " %sx^%d ", display(func.coeff[i]), i);
            printf(" %sx^%d ", display(func.coeff[i]), i);
        }
        fprintf(f, "\n\nRoot Set: ");
        printf("\n\nRoot Set: ");
        for (i = 1; i <= (int) roots[func.deg][0]; i++) {
            fprintf(f, "\n\tx%d = %s", i, display(roots[func.deg][i]));
            printf("\n\tx%d = %s", i, display(roots[func.deg][i]));
        }
    } else {
        solve(diff(func), a, b, e, choice);
        int i, index = 1;
        double ab[3] = {2, a, b};
        // Narrowed range after bisection 
        double range[2];
        // Special points is extrema and inflection points
        double special[(int) roots[func.deg - 1][0] + (int) roots[func.deg - 2][0] + 1];
        merge(roots[func.deg - 1], roots[func.deg - 2], special);
        // Range to be consider division
        double point[(int) special[0] + 3];
        merge(special, ab, point);

        fprintf(f, "\n\nCurrent executing function: ");
        printf("\n\nCurrent executing function: ");
        for (i = 0; i <= func.deg; i++) {
            fprintf(f, " %sx^%d ", display(func.coeff[i]), i);
            printf(" %sx^%d ", display(func.coeff[i]), i);
        }
        for (i = 1; i < (int) point[0]; i++) {
            if (sign(calc(func, point[i] + e / 2)) != sign(calc(func, point[i + 1] - e / 2))) {
                fprintf(f, "\n\n\tDivision range from %s ", display(point[i] + e / 2));
                fprintf(f, "to %s", display(point[i + 1] - e / 2));
                printf("\n\n\tDivision range from %s ", display(point[i] + e / 2));
                printf("to %s", display(point[i + 1] - e / 2));
                bisection(func, point[i] + e / 2, point[i + 1] - e / 2, 0.5, range);
                fprintf(f, "\n\t\tBisection range method from %s ", display(range[0]));
                fprintf(f, "to %s", display(range[1]));
                printf("\n\t\tBisection range method from %s ", display(range[0]));
                printf("to %s", display(range[1]));
                if (choice == 1) {
                    roots[func.deg][index++] = nrsolve1(func, range[0], range[1], e);
                } else
                if (choice == 2) {
                    roots[func.deg][index++] = nrsolve2(func, range[0], range[1], e);
                } else
                if (choice == 3) {
                    roots[func.deg][index++] = nriter(func, range[0], range[1], n);
                }
            }
            // Is special point a root ?
            if (fabs(calc(func, point[i])) <= fabs(calc(diff(func), point[i]) * e)) {
                roots[func.deg][index++] = point[i];
            }
        }
        if (fabs(calc(func, point[(int) point[0]])) <= fabs(calc(diff(func), point[(int) point[0]]) * e)) {
            roots[func.deg][index++] = point[(int) point[0]];
        }
        roots[func.deg][0] = index - 1;
        fprintf(f, "\n\nRoots Set: ");
        printf("\n\nRoots Set: ");
        if (roots[func.deg][0] == 0) {
            fprintf(f, "\nNo root found");
            printf("\nNo root found");
        } else {
            for (i = 1; i <= (int) roots[func.deg][0]; i++) {
                fprintf(f, "\n\tx%d = %s", i, display(roots[func.deg][i]));
                printf("\n\tx%d = %s", i, display(roots[func.deg][i]));
            }
        }
    }
}
