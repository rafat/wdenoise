#include "ebayesthresh.h"

double mean(double* vec, int N) {
	int i;
	double m;
	m = 0.0;

	for (i = 0; i < N; ++i) {
		m+= vec[i];
	}
	m = m / N;
	return m;
}

double var(double* vec, int N) {
	double v,temp,m;
	int i;
	v = 0.0;
	m = mean(vec,N);

	for (i = 0; i < N; ++i) {
		temp = vec[i] - m;
		v+= temp*temp;
	}

	v = v / N;

	return v;

}

double median(double *x, int N) {
	double sigma;

    qsort(x, N, sizeof(double), compare_double);

    if ((N % 2) == 0) {
		sigma = (x[N/2 - 1] + x[N/2] ) / 2.0;
	} else {
		sigma = x[N/2];
	}

	return sigma;
}

double mad(double *x, int N) {
	double sigma;
	int i;

	sigma = median(x,N);

	for(i = 0; i < N;++i) {
		x[i] = (x[i] - sigma) > 0 ? (x[i] - sigma) : -(x[i] - sigma);
	}

	sigma = median(x,N);

	return sigma;
}

double medad(const double *x, int N, int low, int high) {
	double m, constant,temp,center;
	double *x2;
	int rem, N2,i;
	constant = 1.4826;
	rem = N % 2;

	x2 = (double*)calloc(N, sizeof(double));

	if ((low == 1 || high == 1) && rem == 0) {
		if (low == 1 && high == 1) {
			printf("Both low and high cannot be true \n");
			exit(1);
		}

		for (i = 0; i < N; ++i) {
			x2[i] = x[i];
		}
		center = median(x2, N);
		N2 = N / 2 + high - 1;

		for (i = 0; i < N; ++i) {
			temp = x[i] - center;
			x2[i] = fabs(temp);
		}
		qsort(x2, N, sizeof(double), compare_double);
		m = constant * x2[N2];

	}
	else {
		for (i = 0; i < N; ++i) {
			x2[i] = x[i];
		}
		center = median(x2, N);
		for (i = 0; i < N; ++i) {
			temp = x[i] - center;
			x2[i] = fabs(temp);
		}
		m = median(x2, N);
		m *= constant;
	}

	free(x2);
	return m;
}

int minindex(double *arr, int N) {
	double min;
	int index,i;

	min = DBL_MAX;
	index = 0;
	for(i = 0; i < N;++i) {
		if (arr[i] < min) {
			min = arr[i];
			index = i;
		}
	}

	return index;

}

void ebayesthresh(double *x, int N, const char* prior, double *a, int bayesfac, double *sdev, int ls,
	const char* thresrule, int universalthresh, int stabadjustment, int isotone,  double *muhat) {

	int i,stabadjustment_condition,maxits;
	double m_sdev,ax,wx,sx,tol,tt,tcor,sdr;
	double *s;
	char pr[2];
	char th[2];

	s = (double*)calloc(ls, sizeof(double));

	memcpy(&pr, prior, sizeof(pr) - 1);
	pr[1] = '\0';

	memcpy(&th, thresrule, sizeof(th) - 1);
	th[1] = '\0';

	if (ls == 1) {
		if (sdev == NULL) {
			sdr = medad(x, N, 0,0);
		}
		stabadjustment_condition = 1;
	}
	else {
		if (!strcmp(prior, "cauchy")) {
			printf("Standard deviation has to be homogeneous for Cauchy prior.\n");
			exit(1);
		}
		if (ls != N) {
			printf("Standard deviation has to be homogeneous or has the same length as observations. \n");
			exit(1);
		}
		stabadjustment_condition = stabadjustment;
	}

	if (stabadjustment_condition) {
		if (sdev == NULL) {
			m_sdev = sdr;
			s[0] = 1.0;
		}
		else {
			m_sdev = mean(sdev, ls);
			for (i = 0; i < ls; ++i) {
				s[i] = sdev[i] / m_sdev;
			}
		}
			
		for (i = 0; i < N; ++i) {
			x[i] /= m_sdev;
		}
	}
	else {
		if (sdev == NULL) {
			s[0] = sdr;
		}
		else {
			for (i = 0; i < ls; ++i) {
				s[i] = sdev[i];
			}
		}
		
	}
	/*
	For now Implementation works only for scalar sdev
	*/
	sx = s[0];
	//printf("sx %g sdev %g \n", sx,*sdev);
	if (!strcmp(pr, "l") && (a == NULL)) {
	
		wandafromx(x, N, sx, universalthresh, &ax, &wx);
		//printf("estimated a val %g w %g \n",ax,wx);
	}
	else if (isotone) {
		ax = *a;
		tol = 1e-08;
		maxits = 20;
		wmonfromx(x, N, prior, ax, tol, maxits, &wx);
		//printf("estimated a val %g w %g \n",ax,wx);
	}
	else {
		ax = *a;
		wx = wfromx(x, N, sx, prior, ax, universalthresh);
		//printf("estimated a val %g w %g \n",ax,wx);
		//printf("wx %g \n", wx);
	}

	if (strcmp(th, "m")) {
		tfromw(&wx, 1, sx, prior, bayesfac, ax, &tt);

		if (stabadjustment_condition) {
			tcor = tt * m_sdev;
		}
		else {
			tcor = tt;
		}
	}

	//printf("ax %g wx %g \n", ax,wx);
	
	if (!strcmp(thresrule, "median")) {
		postmed(x, N, sx, prior, ax, wx, muhat);
	}
	else if (!strcmp(thresrule, "mean")) {
		postmean(x, N, sx, prior, ax, wx, muhat);
	}
	else if (!strcmp(thresrule, "hard") || !strcmp(thresrule, "soft")) {
		threshld(x, N, tt, thresrule, muhat);
	}
	else {
		muhat = NULL;
	}

	if (stabadjustment_condition) {
		for (i = 0; i < N; ++i) {
			muhat[i] *= m_sdev;
		}
	}

	free(s);
}

void ebayesthresh_wavelet(double *signal, int N, int J, char *wname,char *method, char *ext, const  char *level,
	double vscale, const char* prior, double *a, int bayesfac, const char* thresrule, int isotone, int stabadjustment, int universalthresh, double *denoised) {
	int filt_len, iter, i, dlen, dwt_len, sgn, MaxIter, it,ls;
	double sigma, td, tmp;
	wave_object wave;
	wt_object wt;
	double *dout, *lnoise,*muhat;

	wave = wave_init(wname);

	filt_len = wave->filtlength;

	MaxIter = (int)(log((double)N / ((double)filt_len - 1.0)) / log(2.0));

	if (J > MaxIter) {
		printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
		exit(-1);
	}

	wt = wt_init(wave, method, N, J);
	if (!strcmp(method, "dwt")) {
		setDWTExtension(wt, ext);
		dwt(wt, signal);
	}
	else if (!strcmp(method, "swt")) {
		swt(wt, signal);
	}
	else {
		printf("Acceptable WT methods are - dwt,swt and modwt\n");
		exit(-1);
	}

	lnoise = (double*)malloc(sizeof(double) * J);

	//Set sigma

	iter = wt->length[0];
	dlen = wt->length[J];

	dout = (double*)malloc(sizeof(double) * dlen);
	muhat = (double*)malloc(sizeof(double) * dlen);

	if (!strcmp(level, "first")) {
		for (i = 1; i < J; ++i) {
			iter += wt->length[i];
		}

		for (i = 0; i < dlen; ++i) {
			dout[i] = wt->output[iter + i];
		}

		//sigma = median(dout, dlen) / 0.6745;
		sigma = medad(dout, dlen, 0, 0);
		//printf("sigma %g \n", sigma);
		for (it = 0; it < J; ++it) {
			lnoise[it] = sigma;
		}
	}
	else if (!strcmp(level, "all")){
		for (it = 0; it < J; ++it) {
			dlen = wt->length[it + 1];
			for (i = 0; i < dlen; ++i) {
				dout[i] = wt->output[iter + i];
			}
			//sigma = median(dout, dlen) / 0.6745;
			sigma = medad(dout, dlen, 0, 0);
			lnoise[it] = sigma;
			iter += dlen;
		}

	}
	else {
		sigma = vscale;
		for (it = 0; it < J; ++it) {
			lnoise[it] = sigma;
		}
	}


	dwt_len = wt->outlength;
	iter = wt->length[0];

	ls = 1;
	

	for (it = 0; it < J; ++it) {
		sigma = lnoise[it];
		dlen = wt->length[it + 1];

		for (i = 0; i < dlen; ++i) {
			dout[i] = wt->output[iter + i];
		}
		ebayesthresh(dout, dlen, prior, a, bayesfac, &sigma, ls, thresrule, universalthresh, stabadjustment, isotone, muhat);

		for (i = 0; i < dlen; ++i) {
			wt->output[iter + i] = muhat[i];
		}

		iter += wt->length[it + 1];
	}

	if (!strcmp(method, "dwt")) {
		idwt(wt, denoised);
	}
	else if (!strcmp(method, "swt")) {
		iswt(wt, denoised);
	}

	free(muhat);
	free(dout);
	free(lnoise);
	wave_free(wave);
	wt_free(wt);
}

void denoise_ebayes(double *signal, int N, int J, char *wname, char *method, char *ext, const char *level, double vscale,
	const char* prior, double *a, int bayesfac, const char* thresrule,  double *denoised) {
	double *swtsignal, *swtdenoised;
	
	int i, ret, div, len, N2;
	int isotone;
	int stabadjustment;
	int universalthresh;

	isotone = 0;
	stabadjustment = 1;
	universalthresh = 1; // 1 - adjust threshold for standard deviation , 0 - keep weight between [0,1] 

	if (!strcmp(method, "swt")) {
		div = 1;
		for (i = 0; i < J; ++i) {
			div *= 2;
		}
		ret = N % div;
		if (ret != 0) {
			len = div - ret;
			N2 = N + len;
			swtsignal = (double*)malloc(sizeof(double) * N2);
			swtdenoised = (double*)malloc(sizeof(double) * N2);
			for (i = 0; i < N; ++i) {
				swtsignal[i] = signal[i];
			}
			if (N % 2 == 0) {
				for (i = 0; i < len; ++i) {
					swtsignal[N + i] = signal[i];
				}
			}
			else {
				swtsignal[N] = swtsignal[N - 1];
				for (i = 1; i < len; ++i) {
					swtsignal[N + i] = signal[i - 1];
				}
			}
			
			ebayesthresh_wavelet(swtsignal,N2,J,wname,method,ext,level,vscale,prior,a,bayesfac,thresrule,isotone,stabadjustment,universalthresh,swtdenoised);

			for (i = 0; i < N; ++i) {
				denoised[i] = swtdenoised[i];
			}

			free(swtsignal);
			free(swtdenoised);
		}
		else {
			ebayesthresh_wavelet(signal, N, J, wname, method, ext, level, vscale, prior, a, bayesfac, thresrule, isotone, stabadjustment, universalthresh, denoised);
		}

	}
	else if (!strcmp(method, "dwt")) {
		ebayesthresh_wavelet(signal, N, J, wname, method, ext, level, vscale, prior, a, bayesfac, thresrule, isotone, stabadjustment, universalthresh, denoised);
	}
	else {
		printf("Acceptable WT methods are - dwt and swt\n");
		exit(-1);
	}

}
