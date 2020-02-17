#ifndef WDENOISE_H_
#define WDENOISE_H_

#include "wavelib.h"

#define BLOCK_THRESHOLD_COEFFICIENT 4.5052409999999999
#define W_PI 3.14159265358979323846

#ifdef __cplusplus
extern "C" {
#endif

typedef struct wdenoise_set* wdenoise_object;

wdenoise_object wdenoise_init(int length, int J,const char* wname);

struct wdenoise_set{
	int N; //signal length
	int J; // Levels of Wavelet decomposition
    double qval; // FDR shrink q value
	char wname[10]; //Wavelet name
	char wmethod[10]; //Wavelet decomposition method - dwt or swt
	char cmethod[10]; //Cnvolution Method - direct or fft . Available only for modwt.
	// SWT and DWT only use direct method.
	char ext[10]; // Signal Extension - sym or per
	char thresh[20]; // thresholding - soft or hard
	char level[10]; // Noise Estimation level - first or all
	char dmethod[20]; //Denoising Method
	// EBayes Parameters
	char prior[20];// laplace or cauchy
	char ebthresh[20];// One of median, mean, hard or soft. Default "median" calculates postmedian threshold
	int uthresh; // 1 - adjust threshold for standard deviation , 0 - keep weight between [0,1] . Default 1
	double scalefac;// scale factor. Used only if Laplace prior is used. Default value is 0.5
	// Set double *a in ebayesthresh function to NULL if you want it to be estimated.
	int bayesfac; // bayes factor threshold is calculated instead of posterior median function if bayesfac is 1
	int isotone;// use weighted least squares monotone regression to calculate EBayes weights if isotone is 1. Default 0
	int stabadjustment;// Set to 1 . To-Do list : A more flexible implementation that allows more values of std dev
};

void visushrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void visushrink2(double *data,int rows,int cols, int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void ebayesthresh_wavelet2(double *data, int rows,int cols, int J, char *wname,char *method, char *ext, const  char *level,
	double vscale, const char* prior, double *a, int bayesfac, const char* thresrule, int isotone, int stabadjustment, int universalthresh, double *denoised);

void minimaxshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void minimaxshrink2(double *data,int rows,int cols, int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void sureshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void sureshrink2(double *data,int rows,int cols, int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void modwtshrink(double *signal, int N, int J, const char *wname, const char *cmethod, const char *ext, const char *thresh, double *denoised);

void modwtshrink2(double *data, int rows,int cols, int J, const char *wname, const char *thresh, double *denoised);

void fdrshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,double qval,const char *thresh,const char *level,double *denoised);

void neighcoeffshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void neighblockshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void blockthreshshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised);

void wdenoise(wdenoise_object obj, double *signal,double *denoised);

void setWDenoiseMethod(wdenoise_object obj, const char *dmethod);

void setWDenoiseWTMethod(wdenoise_object obj, const char *wmethod);

void setWDenoiseWTExtension(wdenoise_object obj, const char *extension);

void setWDenoiseParameters(wdenoise_object obj, const char *thresh,const char *level);

void setWDenoiseEbayesParameters(wdenoise_object obj, const char *prior,const char *ebthresh,int *bayesfac, int *isotone,int *universalthresh);

void setWDenoiseEbayesScale(wdenoise_object obj, double *scalefac);

void setWDenoiseQval(wdenoise_object obj,double qval);

void wdenoise_free(wdenoise_object object);

#ifdef __cplusplus
}
#endif



#endif /* WDENOISE_H_ */
