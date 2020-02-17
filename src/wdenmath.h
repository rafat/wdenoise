/*
Copyright (c) 2019, Rafat Hussain
*/
#ifndef WDENMATH_H_
#define WDENMATH_H_

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "../header/wavelib.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define PIVAL 3.14159265358979323846264338327950288

#define XINFVAL 1.79e+308

#define XNINFVAL 2.2251e-308

#ifndef DBL_MAX_EXP
#define DBL_MAX_EXP   +1024
#endif

#ifndef DBL_MIN_EXP
#define DBL_MIN_EXP   -1021
#endif

typedef struct custom_funcprior_set custom_funcprior;

struct custom_funcprior_set{
	void(*funcprior) (double *x,int N, double *z, double *w, double *f, void *params);// M Functions in N variables
	void *params;
};

#define FUNCPRIOR_EVAL(F,x,N,z,w,f) (*((F)->funcprior))(x,N,z,w,(f),(F)->params)

typedef struct custom_function_set custom_function;

struct custom_function_set{
	double(*funcpt) (double *x, int N, void *params);// Function in N variables
	void *params;
};

typedef struct custom_gradient_set custom_gradient;

struct custom_gradient_set{
	void(*funcgrad) (double *x, int N, double *g, void *params);
	void *params;
};

#define FUNCPT_EVAL(F,x,N) (*((F)->funcpt))(x,N,(F)->params)

#define FUNCGRAD_EVAL(F,x,N,g) (*((F)->funcgrad))(x,N,(g),(F)->params)

int compare_double(const void* a, const void* b);

int compare_double_desc(const void* a, const void* b);

void range(double *x, int N, double *max, double *min);

double erf__(double x);

double erfc1__(int ind, double x);
/*
double erf(double x);

double erfc(double x);
*/
double erfcx(double x);

double erfinv(double x);

double erfcinv(double x);

double dnorm(double x, double mu, double sigma);

double pnorm(double x, double mu, double sigma);

double qnorm(double p, double mu, double sigma);

void dnorm_vec(double *x,int N, double mu, double sigma,double *oup);

void pnorm_vec(double *x,int N, double mu, double sigma, double *oup);

void qnorm_vec(double *p,int N, double mu, double sigma,double *oup);

void beta_cauchy(double *x, int N, double *beta);

void beta_laplace(double *x, int N, double s, double a, double *beta);

void vecbinsolv(double *zf, int nz, custom_funcprior *func, double *tlo, double *thi, int nits, double *z, double *w, double *tsol);

void postmean_cauchy(double *x, int N, double w, double *oup);

void postmean_laplace(double *x, int N, double s, double ax, double w, double *oup);

void postmean(double *x, int N, double s, const char* prior, double a, double w, double *oup);

void cauchy_medzero(double *x, int N,  double *z, double *w, double *res,void *params);

void cauchy_threshzero(double *z, int N, double *w, double *y);

void laplace_threshzero(double *x, int N, double s, double a, double *w, double *z);

void postmed_cauchy(double *x, int nx, double wx, double *oup);

void postmed_laplace(double *x, int N, double s, double ax, double w, double *oup);

void postmed(double *x, int N, double s, const char* prior, double a, double w, double *oup);

void isotone(double *y, int N, double *w, int increasing, double *z);

void wfromt(double *tt, int N, double s, const char* prior, double a, double *w);

double wfromx(double *x, int N, double s, const char* prior, double a, int uthresh);

int wmonfromx(double *x, int N, const char* prior, double a, double tol, int maxits,double *w);

void beta_cauchy_func(double *x, int N, double *z, double *w, double *beta,void *params);

void beta_laplace_func(double *x, int N, double *z, double *w, double *beta, void *params);

void cauchy_threshzero_func(double *x, int N, double *z, double *w, double *beta, void *params);

void laplace_threshzero_func(double *x, int N, double *z, double *w, double *beta, void *params);

void tfromw(double *w, int N, double s, const char* prior, int bayesfac, double a, double *tt);

void threshld(double *x, int N, double td, const char* thresh,double *z);

void wpost_laplace(double *w, int N, double *x, double s, double a, double *z);

double negloglik_laplace(double *xx, int N, void *params);

double negloglik_laplace2(double *xx, int M, void *params) ;

void wandafromx(double  *x, int N,double s,int universalthresh,double *a, double *w);

void tfromx(double *x, int N, double s, const char* prior, int bayesfac,double a, int uthresh,double *tt);

int unique_sort(double *x, int N);

double zetafromx(double *xd, int nx, double *cs,const char* prior, double a, double *pilo );

double macheps();

void mdisplay(double *A, int row, int col);

int mvalue(int N);

int grad_cd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
	double eps3, double *f);

double signx(double x);

int grad_calc2(custom_function *funcpt, double *x, int N, double *dx, double eps3, double *f);

int grad_calc(custom_function *funcpt, double *x, int N, double *dx, double eps2, double *f);

int grad_fd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
	double eps2, double *f);

int grad_rextp(custom_function *funcpt, double *x, int N, int r, int v, double *a, double *h, double *xp, double *xm, double meps, double *f);

double l2norm(double *vec, int N);

void minfun2d(custom_function *funcpt, custom_gradient *funcgrad, double *x, int NN, int MM, double *xl, double *xu);


#ifdef __cplusplus
}
#endif

#endif /* WDENMATH_H_ */