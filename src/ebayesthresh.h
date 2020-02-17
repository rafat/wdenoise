#ifndef EBAYESTHRESH_H_
#define EBAYESTHRESH_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "wdenmath.h"

#ifdef __cplusplus
extern "C" {
#endif

double mean(double* vec, int N);

double var(double* vec, int N);

double median(double *x, int N);

double mad(double *x, int N);

double medad(const double *x, int N, int low, int high);

int minindex(double *arr, int N);

void ebayesthresh(double *x, int N, const char* prior, double *a, int bayesfac, double *sdev, int ls,
	const char* thresrule, int universalthresh, int stabadjustment, int isotone,  double *muhat);

void ebayesthresh_wavelet(double *signal, int N, int J, char *wname,char *method, char *ext, const  char *level,
	double vscale, const char* prior, double *a, int bayesfac, const char* thresrule, int isotone, int stabadjustment, int universalthresh, double *denoised);

void denoise_ebayes(double *signal, int N, int J, char *wname, char *method, char *ext, const char *level, double vscale,
	const char* prior, double *a, int bayesfac, const char* thresrule, double *denoised);

#ifdef __cplusplus
}
#endif



#endif /* EBAYESTHRESH_H_ */
