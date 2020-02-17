#ifndef TEST2D_H_
#define TEST2D_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "../header/wdenoise.h"
#include "mt19937ar.h"
/*
#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#endif
*/
#ifdef __cplusplus
extern "C" {
#endif

double rmse(int N, double *x, double *y);

double corrcoef(int N, double *x, double *y);

void getImageMSEandPSNR(unsigned char *img, unsigned char *outimg,int rows, int cols, int channels, double *mse, double *psnr);

double absmax(double *array, int N);

void rescale(double *x,int N);


#ifdef __cplusplus
}
#endif


#endif /* TEST2D*/
