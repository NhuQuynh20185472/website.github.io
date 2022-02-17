#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#define INF 999.0

// Polynomial Function
struct PolynomialFunction {
    double coeff[100];
    int deg;
} func;
typedef struct PolynomialFunction polynomial;

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

double min(polynomial func, double a, double b) {
    double x0 = a, e = 1e-3, xmin = a, alpha = (b - a) / 10000;
    polynomial dfunc = diff(func);
    do {
        x0 = x0 + e;
        while ((fabs(calc(dfunc, x0)) > e) && (x0 <= b)) {
            x0 = x0 + alpha * fabs(calc(dfunc, x0));
        }
        if (x0 > b) break;
        else
        if (calc(func, x0) < calc(func, xmin)) {
            xmin = x0;
        }
    } while (true);
    if (calc(func, xmin) < calc(func, b)) return xmin;
    else return b;
}

double max(polynomial func, double a, double b) {
    double x0 = a, e = 1e-3, xmax = a, alpha = (b - a) / 10000;
    polynomial dfunc = diff(func);
    do {
        x0 = x0 + e;
        while ((fabs(calc(dfunc, x0)) > e) && (x0 <= b)) {
            x0 = x0 + alpha * fabs(calc(dfunc, x0));
        }
        if (x0 > b) break;
        else
        if (calc(func, x0) > calc(func, xmax)) {
            xmax = x0;
        }
    } while (true);
    if (calc(func, xmax) > calc(func, b)) return xmax;
    else return b;
}

//Newton Raphson method
double sol(polynomial func, double a, double b, double e) {
    double s, x0, m1, m2, c;
    polynomial dfunc = diff(func), ddfunc = diff(dfunc);
    int i = 1;
    if (sign(calc(func, a)) == sign(calc(ddfunc, a))) {
        x0 = a;
    } else {
        x0 = b;
    }
    if (calc(dfunc, a) < 0) {
        m1 = calc(dfunc, max(dfunc, a, b));
    }
    else {
        m1 = calc(dfunc, min(dfunc, a, b));
	}
    if (calc(ddfunc, a) < 0) {
        m2 = calc(ddfunc, min(ddfunc, a, b));
    }
    else {
        m2 = calc(ddfunc, max(ddfunc, a, b));
	}
    c = sqrt(fabs(2.0 * m1 * e / m2));
	printf("\nError = %lf", c);
    do {
        s = x0;
        x0 = x0 - calc(func, x0) / calc(dfunc, x0);
        printf("\nLan thu %d:\n", i);
        printf("x = %8.15lf", x0);
        i++;
    } while (fabs(x0 - s) >= c);
    printf("\nVay so lan lap: %d \n", i - 1);
    printf("x = %lf", x0);
    return (x0);
}

main() {
    FILE * f;
    int i;
    double D[2], err, a, b;
    f = fopen("input.txt", "r");
    while (true) {
        fscanf(f, "%lf", &func.coeff[i]);
        if (feof(f)) break;
        else i++;
    }
    func.deg = i;
    fclose(f);
    printf("Ham so: ");
	for (i = 0; i <= func.deg; i++) {
		printf(" %lfx^%d ", func.coeff[i], i);
	}
	err = 0;
	while (err == 0) {
		printf("\nNhap sai so khac 0: ");
		scanf("%lf", &err);
	}
    rootsRange(func, -INF, INF, D);
    printf("\nKhoang chua nghiem [%lf, %lf]", D[0], D[1]);
    int n = func.deg * func.deg * func.deg;
    double range[n];
    printf("\n");
    for (i = 0; i <= n; i++) {
        range[i] = D[0] + i * (D[1] - D[0]) / n;
//        printf("%lf ", range[i]);
    }
    for (i = 0; i < n; i++) {
        if (sign(calc(func, range[i])) != sign(calc(func, range[i + 1]))) {
            a = min(func, range[i], range[i + 1]);
            b = max(func, range[i], range[i + 1]);
            if (sign(calc(diff(diff(func)), a)) == sign(calc(diff(diff(func)), b))) {
                if (a > b) {
                    printf("\nGiai phuong trinh tren [%lf, %lf]", b, a);
                    sol(func, b, a, err);
                } else {
                    printf("\nGiai phuong trinh tren [%lf, %lf]", a, b);
                    sol(func, a, b, err);
                }
            }
        }
    }
    getch();
}
