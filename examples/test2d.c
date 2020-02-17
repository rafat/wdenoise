#include "test2d.h"

double rmse(int N, double *x, double *y) {
	double rms;
	int i;

	rms = 0.0;

	for (i = 0; i < N; ++i) {
		rms += (x[i] - y[i]) * (x[i] - y[i]);
	}

	rms = sqrt(rms / (double)N);

	return rms;
}

double corrcoef(int N, double *x, double *y) {
	double cc, xm, ym, tx, ty, num, den1, den2;
	int i;
	xm = ym = 0.0;
	for (i = 0; i < N; ++i) {
		xm += x[i];
		ym += y[i];
	}

	xm = xm / N;
	ym = ym / N;
	num = den1 = den2 = 0.0;

	for (i = 0; i < N; ++i) {
		tx = x[i] - xm;
		ty = y[i] - ym;
		num += (tx*ty);
		den1 += (tx*tx);
		den2 += (ty*ty);
	}

	cc = num / sqrt(den1*den2);

	return cc;
}

void getImageMSEandPSNR(unsigned char *img, unsigned char *outimg,int rows, int cols, int channels, double *mse, double *psnr) {
	int i,K;
	double sum = 0.0,temp,mval;
	unsigned char *pi, *po;

	if (channels == 1 || channels == 2) {
		K = 1;
	} else if (channels == 3 || channels == 4) {
		K = 3;
	} else {
		printf("This example only supports images with 1,2,3 or 4 channels. \n");
		exit(-1);
	}

	mval = -10;

	if (K == 1) {
		pi = img;
		po = outimg;
		for(i = 0; i < rows*cols;++i) {
			temp = (double) (*pi - *po);
			sum += temp*temp;
			if ((double)*pi > mval) {
				mval = (double)*pi;
			}
			pi+=channels;
			po+=channels;
		}
	} else {
		pi = img;
		po = outimg;
		for(i = 0; i < rows*cols*K;++i) {
			temp = (double) (*pi - *po);
			sum += temp*temp;
			if ((double)*pi > mval) {
				mval = (double)*pi;
			}
			temp = (double) (*(pi+1) - *(po+1));
			sum += temp*temp;
			if ((double)*(pi+1) > mval) {
				mval = (double)*(pi+1);
			}
			temp = (double) (*(pi+2) - *(po+2));
			sum += temp*temp;
			if ((double)*(pi+2) > mval) {
				mval = (double)*(pi+2);
			}
			pi+=channels;
			po+=channels;
		}
	}

	*mse = sum/(rows*cols*channels);

	temp = mval*mval/(*mse);
	*psnr = 10 * log(temp)/log(10);
}

double absmax(double *array, int N) {
	double max;
	int i;

	max = 0.0;
	for (i = 0; i < N; ++i) {
		if (fabs(array[i]) >= max) {
			max = fabs(array[i]);
		}
	}

	return max;
}

void rescale(double *x,int N) {
	int i;
	double xmin,xmax,ymin,ymax, mult;
	ymin = 0;
	ymax = 255;
	xmin = 10000;
	xmax = -10000;

	for (i = 0; i < N;++i) {
		if (x[i] < xmin) {
			xmin = x[i];
		}
		if (x[i] > xmax) {
			xmax = x[i];
		}
	}

	mult = (ymax - ymin)/(xmax - xmin);

	for(i = 0; i < N;++i) {
		x[i] = (x[i] - xmin) * mult;
	}
}
