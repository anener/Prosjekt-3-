#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10

using namespace std;

double gammln(double);
void gaulag(double *, double *, int, double);
void gaussLaguerre(int);
double intFunction_laguerre(double, double, double, double, double, double);
