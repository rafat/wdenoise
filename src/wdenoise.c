#include "../header/wdenoise.h"
#include "ebayesthresh.h"

wdenoise_object wdenoise_init(int length, int J,const char* wname) {
	wdenoise_object obj = NULL;

	obj = (wdenoise_object)malloc(sizeof(struct wdenoise_set) +sizeof(double)*10);

	obj->N = length;
	obj->J = J;
	obj->qval = 0.05;

	strcpy(obj->wname,wname);

	//Set Default Values
	strcpy(obj->dmethod,"sureshrink");
	strcpy(obj->ext,"sym");
	strcpy(obj->level,"all");
	strcpy(obj->thresh,"soft");
	strcpy(obj->wmethod,"dwt");
	strcpy(obj->cmethod,"direct");

	// EBayes Params

	obj->scalefac = 0.5;
	obj->uthresh = 1;
	obj->bayesfac = 0;
	obj->isotone = 0;
	obj->stabadjustment = 1;

	strcpy(obj->prior,"laplace");
	strcpy(obj->ebthresh,"median");

	return obj;
}

static wt_object getFWT(wave_object wave,double *signal,const char *wname,const char *method,const char* ext,int N,int J) {
	wt_object wt;

	int filt_len,MaxIter;

	filt_len = wave->filtlength;

	MaxIter = (int) (log((double)N / ((double)filt_len - 1.0)) / log(2.0));

	if (J > MaxIter) {
		printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n",MaxIter);
		exit(-1);
	}

	wt = wt_init(wave,method,N,J);
	if(!strcmp(method,"dwt")) {
		setDWTExtension(wt,ext);
		dwt(wt,signal);
	} else if(!strcmp(method,"swt")) {
		swt(wt,signal);
	} else {
		printf("Acceptable WT methods are - dwt,swt and modwt\n");
		exit(-1);
	}

	return wt;
}

static wt2_object getFWT2(wave_object wave,double *data,int rows, int cols,const char *wname,const char *method,const char* ext,int J) {
	wt2_object wt;
	int filt_len,MaxIter, MaxRows, MaxCols;
	double *wavecoeffs;

	filt_len = wave->filtlength;

	MaxRows = (int) (log((double)rows / ((double)filt_len - 1.0)) / log(2.0));
	MaxCols = (int) (log((double)cols / ((double)filt_len - 1.0)) / log(2.0));

	MaxIter = (MaxRows < MaxCols) ? MaxRows : MaxCols;

	if (J > MaxIter) {
		printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n",MaxIter);
		exit(-1);
	}

	wt = wt2_init(wave,method,rows,cols,J);

	return wt;
}

static void getIFWT(wt_object wt,const char *method,double *denoised) {
	if(!strcmp(method,"dwt")) {
		idwt(wt,denoised);
	} else if(!strcmp(method,"swt")) {
		iswt(wt,denoised);
	}
}

static void getIFWT2(wt2_object wt,const char *method,double *wavecoeffs,double *denoised) {
	if(!strcmp(method,"dwt")) {
		idwt2(wt,wavecoeffs,denoised);
	} else if(!strcmp(method,"swt")) {
		iswt2(wt,wavecoeffs,denoised);
	}
}

static void noiseest(wt_object wt,double *dout,int J,const char *level,double *lnoise) {
	int i,iter,dlen,it;
	double sigma;

	iter = wt->length[0];
	dlen = wt->length[J];

	if(!strcmp(level,"first")) {
		for (i = 1; i < J; ++i) {
			iter += wt->length[i];
		}

		for(i = 0; i < dlen;++i) {
			dout[i] = fabs(wt->output[iter+i]);
		}

		sigma = median(dout,dlen) / 0.6745;
		for(it = 0; it < J;++it) {
			lnoise[it] = sigma;
		}
	} else if(!strcmp(level,"all")){
		for(it = 0; it < J;++it) {
			dlen = wt->length[it+1];
			for(i = 0; i < dlen;++i) {
				dout[i] = fabs(wt->output[iter+i]);
			}
			sigma = median(dout,dlen) / 0.6745;
			lnoise[it] = sigma;
			iter += dlen;
		}

	} else {
		printf("Acceptable Noise estimation level values are - first and all \n");
		exit(-1);
	}
}

static void noiseest2(wt2_object wt,double *wavecoeffs,double *dout,int J,const char *level,double *lnoise) {
	int i,aLH,aHL,aHH,vsize,it,dlen;
	double sigmah,sigmav,sigmad;
	// lnoise has dimensions 3 * J lnh,lnv,lnd


	if(!strcmp(level,"first")) {
		dlen = wt->dimensions[2*J-2] * wt->dimensions[2*J-1];

		// Horizontal Component Noise Estimation
		aLH = wt->coeffaccess[3*(J-1)+1];

		for(i = 0; i < dlen;++i) {
			dout[i] = fabs(wavecoeffs[aLH+i]);
		}

		sigmah = median(dout,dlen) / 0.6745;

		// Vertical Component Noise Estimation
		aHL = wt->coeffaccess[3*(J-1) + 2];

		for(i = 0; i < dlen;++i) {
			dout[i] = fabs(wavecoeffs[aHL+i]);
		}

		sigmav = median(dout,dlen) / 0.6745;

		// Diagonal Component Noise Estimation
		aHH = wt->coeffaccess[3*(J-1) + 3];

		for(i = 0; i < dlen;++i) {
			dout[i] = fabs(wavecoeffs[aHH+i]);
		}

		sigmad = median(dout,dlen) / 0.6745;

		for(i = 0; i < J;++i) {
			lnoise[3*i] = sigmah;
			lnoise[3*i+1] = sigmav;
			lnoise[3*i+2] = sigmad;
		}


	} else if(!strcmp(level,"all")) {
		for(it = 0; it < J;++it) {
			dlen = wt->dimensions[2*it] * wt->dimensions[2*it+1];
			// Horizontal Component Noise Estimation
			aLH = wt->coeffaccess[3*it+1];

			for(i = 0; i < dlen;++i) {
				dout[i] = fabs(wavecoeffs[aLH+i]);
			}

			lnoise[3*it] = median(dout,dlen) / 0.6745;


			// Vertical Component Noise Estimation
			aHL = wt->coeffaccess[3*it+2];

			for(i = 0; i < dlen;++i) {
				dout[i] = fabs(wavecoeffs[aHL+i]);
			}

			lnoise[3*it+1] = median(dout,dlen) / 0.6745;

			// Diagonal Component Noise Estimation
			aHH = wt->coeffaccess[3*it+3];

			for(i = 0; i < dlen;++i) {
				dout[i] = fabs(wavecoeffs[aHH+i]);
			}

			lnoise[3*it+2] = median(dout,dlen) / 0.6745;


		}

	} else {
		printf("Acceptable Noise estimation level values are - first and all \n");
		exit(-1);
	}

}

static void noiseest2a(wt2_object wt,double *wavecoeffs,double *dout,int J,const char *level,double *lnoise) {
	int i,aLH,aHL,aHH,vsize,it,dlen;
	double sigma,sigma2,sigmah,sigmav,sigmad;
	double meanh,meanv,meand,xm,xmsum;
	// lnoise has dimensions 3 * J lnh,lnv,lnd


	if(!strcmp(level,"first")) {
		dlen = wt->dimensions[2*J-2] * wt->dimensions[2*J-1];

		// Diagonal Component Median noise estimator
		aHH = wt->coeffaccess[3*(J-1) + 3];

		for(i = 0; i < dlen;++i) {
			dout[i] = fabs(wavecoeffs[aHH+i]);
		}

		sigma = median(dout,dlen) / 0.6745;
		sigma2 = sigma *sigma;

		// Horizontal Component Noise Estimation
		aLH = wt->coeffaccess[3*(J-1)+1];

		meanh = mean(wavecoeffs+aLH,dlen);
		xmsum = 0.0;
		for(i = 0; i < dlen;++i) {
			xm = wavecoeffs[aLH+i] - meanh;
			xmsum += xm*xm;
		}
		xmsum /= (double) dlen;

		xmsum = sqrt(xmsum);

		sigmah = sigma2 / xmsum;

		// Vertical Component Noise Estimation
		aHL = wt->coeffaccess[3*(J-1) + 2];

		meanv = mean(wavecoeffs+aHL,dlen);
		xmsum = 0.0;
		for(i = 0; i < dlen;++i) {
			xm = wavecoeffs[aHL+i] - meanv;
			xmsum += xm*xm;
		}
		xmsum /= (double) dlen;

		xmsum = sqrt(xmsum);

		sigmav = sigma2 / xmsum;

		// Diagonal Component Noise Estimation
		aHH = wt->coeffaccess[3*(J-1) + 3];

		meand = mean(wavecoeffs+aHH,dlen);
		xmsum = 0.0;
		for(i = 0; i < dlen;++i) {
			xm = wavecoeffs[aHH+i] - meand;
			xmsum += xm*xm;
		}
		xmsum /= (double) dlen;

		xmsum = sqrt(xmsum);

		sigmad = sigma2 / xmsum;

		for(i = 0; i < J;++i) {
			lnoise[3*i] = sigmah;
			lnoise[3*i+1] = sigmav;
			lnoise[3*i+2] = sigmad;
		}


	} else if(!strcmp(level,"all")) {
		for(it = 0; it < J;++it) {
			dlen = wt->dimensions[2*it] * wt->dimensions[2*it+1];

			// Diagonal Component Median Noise Estimation
			aHH = wt->coeffaccess[3*it+3];

			for(i = 0; i < dlen;++i) {
				dout[i] = fabs(wavecoeffs[aHH+i]);
			}

			sigma = median(dout,dlen) / 0.6745;
			sigma2 = sigma *sigma;

			// Horizontal Component Noise Estimation
			aLH = wt->coeffaccess[3*it+1];

			meanh = mean(wavecoeffs+aLH,dlen);
			xmsum = 0.0;
			for(i = 0; i < dlen;++i) {
				xm = wavecoeffs[aLH+i] - meanh;
				xmsum += xm*xm;
			}
			xmsum /= (double) dlen;

			xmsum = sqrt(xmsum);

			lnoise[3*it] = sigma2 / xmsum;


			// Vertical Component Noise Estimation
			aHL = wt->coeffaccess[3*it+2];

			meanv = mean(wavecoeffs+aHL,dlen);
			xmsum = 0.0;
			for(i = 0; i < dlen;++i) {
				xm = wavecoeffs[aHL+i] - meanv;
				xmsum += xm*xm;
			}
			xmsum /= (double) dlen;

			xmsum = sqrt(xmsum);

			lnoise[3*it+1] = sigma2 / xmsum;

			// Diagonal Component Noise Estimation
			aHH = wt->coeffaccess[3*it+3];

			meand = mean(wavecoeffs+aHH,dlen);
			xmsum = 0.0;
			for(i = 0; i < dlen;++i) {
				xm = wavecoeffs[aHH+i] - meand;
				xmsum += xm*xm;
			}
			xmsum /= (double) dlen;

			xmsum = sqrt(xmsum);

			lnoise[3*it+2] = sigma2 / xmsum;

		}

	} else {
		printf("Acceptable Noise estimation level values are - first and all \n");
		exit(-1);
	}

}

static void threshrule(double *x, int N, double td, const char* thresh) {
	int i,sgn;

	if (!strcmp(thresh, "hard")) {
		for (i = 0; i < N; ++i) {
			x[i] = fabs(x[i]) >= td ? x[i] : 0;
		}
	}
	else if (!strcmp(thresh, "soft")) {
		for (i = 0; i < N; ++i) {
			if (fabs(x[i]) < td) {
				x[i] = 0.0;
			}
			else {
				sgn = x[i] > 0 ? 1 : -1;
				if (x[i] == 0) {
					sgn = 0;
				}
				x[i] = sgn * (fabs(x[i]) - td);
			}
		}
	}

}

void visushrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int iter,i,dlen,dwt_len,sgn,it;
	double sigma,td,tmp;
	wt_object wt;
	double *dout,*lnoise;
	wave_object wave;
	wave = wave_init(wname);
	wt = getFWT(wave,signal,wname,method,ext,N,J);

	lnoise = (double*)malloc(sizeof(double) * J);

	//Set sigma

	iter = wt->length[0];
	dlen = wt->length[J];

	dout = (double*)malloc(sizeof(double) * dlen);

	noiseest(wt,dout,J,level,lnoise);

	dwt_len = wt->outlength;
	iter = wt->length[0];
	for(it = 0; it < J;++it) {
		sigma = lnoise[it];
		dlen = wt->length[it+1];
		td = sqrt(2.0 * log(dwt_len)) * sigma;

		threshrule(wt->output+iter,dlen,td,thresh);

		iter += wt->length[it+1];
	}

	getIFWT(wt,method,denoised);

	wave_free(wave);
	free(dout);
	free(lnoise);
	wt_free(wt);
}

void visushrink2(double *data,int rows,int cols, int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int i,dlen,dwt_len,sgn,it,aLH,aHL,aHH;
	double tdh,tdv,tdd;
	wt2_object wt;
	double *wavecoeffs,*dout,*lnoise;
	wave_object wave;
	wave = wave_init(wname);

	wt = getFWT2(wave,data,rows,cols,wname,method,ext,J);

	if(!strcmp(method,"dwt")) {
		setDWT2Extension(wt,ext);
		wavecoeffs = dwt2(wt,data);
	} else if(!strcmp(method,"swt")) {
		wavecoeffs = swt2(wt,data);
	} else {
		printf("Acceptable WT methods are - dwt,swt and modwt\n");
		exit(-1);
	}


	lnoise = (double*)malloc(sizeof(double) * J * 3);

	dlen = wt->dimensions[2*J-2] * wt->dimensions[2*J-1];

	dout = (double*)malloc(sizeof(double) * dlen);

	//noiseest2(wt,wavecoeffs,dout,J,level,lnoise);
	noiseest2a(wt,wavecoeffs,dout,J,level,lnoise);

	dwt_len = wt->outlength;

	for(it = 0; it < J;++it) {
		aLH = wt->coeffaccess[3*it+1];
		aHL = wt->coeffaccess[3*it+2];
		aHH = wt->coeffaccess[3*it+3];

		//tdh = sqrt(2.0 * log(dwt_len)) * lnoise[3*it];
		//tdv = sqrt(2.0 * log(dwt_len)) * lnoise[3*it+1];
		//tdd = sqrt(2.0 * log(dwt_len)) * lnoise[3*it+2];

		if (!strcmp(thresh,"soft")) {
			tdh = 1.5 * lnoise[3*it];
			tdv = 1.5 * lnoise[3*it+1];
			tdd = 1.5 * lnoise[3*it+2];
		} else if (!strcmp(thresh,"hard")) {
			tdh = 3.0 * lnoise[3*it];
			tdv = 3.0 * lnoise[3*it+1];
			tdd = 3.0 * lnoise[3*it+2];
		}


		dlen = wt->dimensions[2*it] * wt->dimensions[2*it+1];
		threshrule(wavecoeffs+aLH,dlen,tdh,thresh);
		threshrule(wavecoeffs+aHL,dlen,tdv,thresh);
		threshrule(wavecoeffs+aHH,dlen,tdd,thresh);

	}

	getIFWT2(wt,method,wavecoeffs,denoised);

	wave_free(wave);
	free(wavecoeffs);
	free(dout);
	free(lnoise);
	wt2_free(wt);
}

void ebayesthresh_wavelet2(double *data, int rows,int cols, int J, char *wname,char *method, char *ext, const  char *level,
	double vscale, const char* prior, double *a, int bayesfac, const char* thresrule, int isotone, int stabadjustment, int universalthresh, double *denoised)  {
	int i,dlen,dwt_len,sgn,it,aLH,aHL,aHH,ls;
	double tdh,tdv,tdd,sigma;
	wt2_object wt;
	double *wavecoeffs,*dout,*lnoise,*muhat;
	wave_object wave;
	wave = wave_init(wname);

	wt = getFWT2(wave,data,rows,cols,wname,method,ext,J);
	ls = 1;

	if(!strcmp(method,"dwt")) {
		setDWT2Extension(wt,ext);
		wavecoeffs = dwt2(wt,data);
	} else if(!strcmp(method,"swt")) {
		wavecoeffs = swt2(wt,data);
	} else {
		printf("Acceptable WT methods are - dwt,swt and modwt\n");
		exit(-1);
	}


	lnoise = (double*)malloc(sizeof(double) * J * 3);

	dlen = wt->dimensions[2*J-2] * wt->dimensions[2*J-1];
	muhat = (double*)malloc(sizeof(double) * dlen);

	dout = (double*)malloc(sizeof(double) * dlen);

	//noiseest2(wt,wavecoeffs,dout,J,level,lnoise);
	noiseest2a(wt,wavecoeffs,dout,J,level,lnoise);

	for(it = 0; it < J;++it) {
		aLH = wt->coeffaccess[3*it+1];
		aHL = wt->coeffaccess[3*it+2];
		aHH = wt->coeffaccess[3*it+3];

		dlen = wt->dimensions[2*it] * wt->dimensions[2*it+1];

		for(i = 0; i < dlen;++i) {
			dout[i] = wavecoeffs[aLH+i];
		}

		sigma = lnoise[3*it];
		ebayesthresh(dout, dlen, prior, a, bayesfac, &sigma, ls, thresrule, universalthresh, stabadjustment, isotone, muhat);

		for(i = 0; i < dlen;++i) {
			wavecoeffs[aLH+i] = muhat[i];
		}

		for(i = 0; i < dlen;++i) {
			dout[i] = wavecoeffs[aHL+i];
		}

		sigma = lnoise[3*it+1];
		ebayesthresh(dout, dlen, prior, a, bayesfac, &sigma, ls, thresrule, universalthresh, stabadjustment, isotone, muhat);

		for(i = 0; i < dlen;++i) {
			wavecoeffs[aHL+i] = muhat[i];
		}

		for(i = 0; i < dlen;++i) {
			dout[i] = wavecoeffs[aHH+i];
		}

		sigma = lnoise[3*it+2];
		ebayesthresh(dout, dlen, prior, a, bayesfac, &sigma, ls, thresrule, universalthresh, stabadjustment, isotone, muhat);

		for(i = 0; i < dlen;++i) {
			wavecoeffs[aHH+i] = muhat[i];
		}

	}

	getIFWT2(wt,method,wavecoeffs,denoised);

	wave_free(wave);
	free(wavecoeffs);
	free(dout);
	free(muhat);
	free(lnoise);
	wt2_free(wt);
}

static void getPTSR(wt_object wt,double *ptsr) {
	int i,k,J,iter;
	double sum,max,fval;
	J = wt->J;

	iter = wt->length[0];

	for(i = 0; i < J;++i) {
		sum = 0.0;
		max = 0.0;
		for(k = 0; k < wt->length[i+1];++k) {
			fval = fabs(wt->output[iter+k]);
			sum += fval;
			max = max < fval ? fval : max;
		}

		ptsr[i] = max / sum;

		iter += wt->length[i+1];
	}

}

static int getDLevel(double *ptsr, int J) {
	int level,i;

	double SAFL;

	SAFL = 0.2;
	level = 0;
	for(i = J-1; i >= 0; --i) {
		if (ptsr[i] <= SAFL) {
			level++;
		} else {
			return level;
		}
	}

	return level;
}

void safshrink(double *signal,int N,const char *wname,const char *method,const char *ext,double *denoised) {
	wave_object wave;
	wt_object wt;
	int filt_len, J, level, i;
	double *ptsr;
	double *thlo,*thhi,*mu,*sigma;

	filt_len = wave->filtlength;

	J = (int) (log((double)N / ((double)filt_len - 1.0)) / log(2.0));

	wave = wave_init(wname);
	wt = wt_init(wave,method,N,J);

	if(!strcmp(method,"dwt")) {
		setDWTExtension(wt,ext);
		dwt(wt,signal);
	} else if(!strcmp(method,"swt")) {
		swt(wt,signal);
	} else {
		printf("Acceptable WT methods are - dwt,swt and modwt\n");
		exit(-1);
	}

	ptsr = (double*)calloc(J,sizeof(double));

	// get PTSR

	getPTSR(wt,ptsr);

	// get Decomposition level
	level = getDLevel(ptsr,J);

	// Calculate Thresholds

	thhi = (double*)calloc(level,sizeof(double));
	thlo = (double*)calloc(level,sizeof(double));
	mu = (double*)calloc(level,sizeof(double));
	sigma = (double*)calloc(level,sizeof(double));

	for(i = level; i > 0;--i) {

	}


	if(!strcmp(method,"dwt")) {
		idwt(wt,denoised);
	} else if(!strcmp(method,"swt")) {
		iswt(wt,denoised);
	}

	free(ptsr);
	free(thhi);
	free(thlo);
	free(mu);
	free(sigma);
	wave_free(wave);
	wt_free(wt);
}

void minimaxshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int iter,i,dlen,dwt_len,sgn,it;
	double sigma,td,tmp;
	wt_object wt;
	double *dout,*lnoise;
	wave_object wave;
	wave = wave_init(wname);
	wt = getFWT(wave,signal,wname,method,ext,N,J);

	lnoise = (double*)malloc(sizeof(double) * J);

	//Set sigma

	iter = wt->length[0];
	dlen = wt->length[J];

	dout = (double*)malloc(sizeof(double) * dlen);

	noiseest(wt,dout,J,level,lnoise);

	dwt_len = wt->outlength;
	iter = wt->length[0];
	for(it = 0; it < J;++it) {
		sigma = lnoise[it];
		dlen = wt->length[it+1];
		td = dwt_len >= 32 ? sigma * (0.3936 + 0.1829*(log(dwt_len) / log(2))) : 0.0;

		threshrule(wt->output+iter,dlen,td,thresh);

		iter += wt->length[it+1];
	}

	getIFWT(wt,method,denoised);

	wave_free(wave);
	free(dout);
	free(lnoise);
	wt_free(wt);
}

void minimaxshrink2(double *data,int rows,int cols, int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int i,dlen,dwt_len,sgn,it,aLH,aHL,aHH;
	double tdh,tdv,tdd;
	wt2_object wt;
	double *wavecoeffs,*dout,*lnoise;
	wave_object wave;
	wave = wave_init(wname);

	wt = getFWT2(wave,data,rows,cols,wname,method,ext,J);

	if(!strcmp(method,"dwt")) {
		setDWT2Extension(wt,ext);
		wavecoeffs = dwt2(wt,data);
	} else if(!strcmp(method,"swt")) {
		wavecoeffs = swt2(wt,data);
	} else {
		printf("Acceptable WT methods are - dwt,swt and modwt\n");
		exit(-1);
	}


	lnoise = (double*)malloc(sizeof(double) * J * 3);

	dlen = wt->dimensions[2*J-2] * wt->dimensions[2*J-1];

	dout = (double*)malloc(sizeof(double) * dlen);

	noiseest2(wt,wavecoeffs,dout,J,level,lnoise);

	dwt_len = wt->outlength;

	for(it = 0; it < J;++it) {
		aLH = wt->coeffaccess[3*it+1];
		aHL = wt->coeffaccess[3*it+2];
		aHH = wt->coeffaccess[3*it+3];

		tdh = dwt_len >= 32 ? lnoise[3*it] * (0.3936 + 0.1829*(log(dwt_len) / log(2))) : 0.0;
		tdv = dwt_len >= 32 ? lnoise[3*it+1] * (0.3936 + 0.1829*(log(dwt_len) / log(2))) : 0.0;
		tdd = dwt_len >= 32 ? lnoise[3*it+2] * (0.3936 + 0.1829*(log(dwt_len) / log(2))) : 0.0;

		dlen = wt->dimensions[2*it] * wt->dimensions[2*it+1];
		threshrule(wavecoeffs+aLH,dlen,tdh,thresh);
		threshrule(wavecoeffs+aHL,dlen,tdv,thresh);
		threshrule(wavecoeffs+aHH,dlen,tdd,thresh);

	}

	getIFWT2(wt,method,wavecoeffs,denoised);

	wave_free(wave);
	free(wavecoeffs);
	free(dout);
	free(lnoise);
	wt2_free(wt);
}


void sureshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int i,it,len,dlen,dwt_len,min_index,sgn,iter;
	double sigma,norm,td,tv,te,ct,thr,temp,x_sum;
	wt_object wt;
	double *dout,*risk,*dsum,*lnoise;
	wave_object wave;
	wave = wave_init(wname);
	wt = getFWT(wave,signal,wname,method,ext,N,J);

	len = wt->length[0];
	dlen = wt->length[J];

	dout = (double*)malloc(sizeof(double) * dlen);
	risk = (double*)malloc(sizeof(double) * dlen);
	dsum = (double*)malloc(sizeof(double) * dlen);
	lnoise = (double*)malloc(sizeof(double) * J);

	iter = wt->length[0];

	if(!strcmp(level,"first")) {
		for (i = 1; i < J; ++i) {
			iter += wt->length[i];
		}

		for(i = 0; i < dlen;++i) {
			dout[i] = fabs(wt->output[iter+i]);
		}

		sigma = median(dout,dlen) / 0.6745;
		for(it = 0; it < J;++it) {
			lnoise[it] = sigma;
		}
	} else if(!strcmp(level,"all")){
		for(it = 0; it < J;++it) {
			dlen = wt->length[it+1];
			for(i = 0; i < dlen;++i) {
				dout[i] = fabs(wt->output[iter+i]);
			}
			sigma = median(dout,dlen) / 0.6745;
			lnoise[it] = sigma;
			iter += dlen;
		}

	} else {
		printf("Acceptable Noise estimation level values are - first and all \n");
		exit(-1);
	}


	for(it = 0; it < J;++it) {
		dwt_len = wt->length[it+1];
		sigma = lnoise[it];

		if ( sigma < 0.00000001) {
			td = 0;
		} else {
			tv = sqrt(2.0 * log(dwt_len));
			norm = 0.0;
			for(i = 0; i < dwt_len;++i) {
				norm += (wt->output[len+i] *wt->output[len+i] /(sigma*sigma));
			}
			te =(norm - (double) dwt_len)/(double) dwt_len;
			ct = pow(log((double) dwt_len)/log(2.0),1.5)/sqrt((double) dwt_len);

			if (te < ct) {
				td = tv;
			} else {
				x_sum = 0.0;

				for(i = 0; i < dwt_len;++i) {
					dout[i] = fabs(wt->output[len+i]/sigma);
				}

				qsort(dout, dwt_len, sizeof(double), compare_double);
				for(i = 0; i < dwt_len;++i) {
					dout[i] = (dout[i]*dout[i]);
					x_sum += dout[i];
					dsum[i] = x_sum;
				}

				for(i = 0;i < dwt_len;++i) {
					risk[i] = ((double)dwt_len - 2 * ((double)i + 1) +dsum[i] +
							dout[i]*((double)dwt_len - 1 -(double) i))/(double)dwt_len;
				}
				min_index = minindex(risk,dwt_len);
				thr = sqrt(dout[min_index]);
				td = thr < tv ? thr : tv;
			}
		}

		td = td * sigma;

		threshrule(wt->output+len,dlen,td,thresh);


		len += wt->length[it+1];
	}

	getIFWT(wt,method,denoised);

	wave_free(wave);
	free(dout);
	free(dsum);
	free(risk);
	free(lnoise);
	wt_free(wt);
}

void sureshrink2(double *data,int rows,int cols, int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int i,dlen,dwt_len,sgn,it,aLH,aHL,aHH,min_index;
	double tdh,tdv,tdd,sigmah,sigmav,sigmad,tv,norm,te,ct,x_sum,thr;
	wt2_object wt;
	double *wavecoeffs,*dout,*lnoise,*dsum,*risk;
	wave_object wave;
	wave = wave_init(wname);

	wt = getFWT2(wave,data,rows,cols,wname,method,ext,J);

	if(!strcmp(method,"dwt")) {
		setDWT2Extension(wt,ext);
		wavecoeffs = dwt2(wt,data);
	} else if(!strcmp(method,"swt")) {
		wavecoeffs = swt2(wt,data);
	} else {
		printf("Acceptable WT methods are - dwt,swt and modwt\n");
		exit(-1);
	}


	lnoise = (double*)malloc(sizeof(double) * J * 3);

	dlen = wt->dimensions[2*J-2] * wt->dimensions[2*J-1];

	dout = (double*)malloc(sizeof(double) * dlen);
	dsum = (double*)malloc(sizeof(double) * dlen);
	risk = (double*)malloc(sizeof(double) * dlen);

	noiseest2(wt,wavecoeffs,dout,J,level,lnoise);
	//printf("Noise Estimated\n");

	wt2_summary(wt);

	for(it = 0; it < J;++it) {
		aLH = wt->coeffaccess[3*it+1];
		aHL = wt->coeffaccess[3*it+2];
		aHH = wt->coeffaccess[3*it+3];
		dwt_len = wt->dimensions[2*it] * wt->dimensions[2*it+1];
		sigmah = lnoise[3*it];
		sigmav = lnoise[3*it+1];
		sigmad = lnoise[3*it+2];
		//printf("sigmah dwt_len %d \n",dwt_len);

		if ( sigmah < 0.00000001) {
			tdh = 0;
		} else {
			tv = sqrt(2.0 * log(dwt_len));
			norm = 0.0;
			for(i = 0; i < dwt_len;++i) {
				norm += wavecoeffs[aLH+i] * wavecoeffs[aLH+i] / (sigmah * sigmah);
			}
			te =(norm - (double) dwt_len)/(double) dwt_len;
			ct = pow(log((double) dwt_len)/log(2.0),1.5)/sqrt((double) dwt_len);
			//printf("tv %g te %g ct %g \n\n",tv,te,ct);
			if (te < ct) {
				tdh = tv;
			} else {
				x_sum = 0.0;

				for(i = 0; i < dwt_len;++i) {
					dout[i] = fabs(wavecoeffs[aLH+i]/sigmah);
				}
				//printf("dout \n");
				qsort(dout, dwt_len, sizeof(double), compare_double);
				//printf("qsort\n");
				for(i = 0; i < dwt_len;++i) {
					dout[i] = (dout[i]*dout[i]);
					x_sum += dout[i];
					dsum[i] = x_sum;
				}
				//printf("dsum\n");
				for(i = 0;i < dwt_len;++i) {
					risk[i] = ((double)dwt_len - 2 * ((double)i + 1) +dsum[i] +
							dout[i]*((double)dwt_len - 1 -(double) i))/(double)dwt_len;
				}
				//printf("risk\n");
				min_index = minindex(risk,dwt_len);
				thr = sqrt(dout[min_index]);
				tdh = thr < tv ? thr : tv;
			}
		}

		//printf("sigmav\n");

		if ( sigmav < 0.00000001) {
			tdv = 0;
		} else {
			tv = sqrt(2.0 * log(dwt_len));
			norm = 0.0;
			for(i = 0; i < dwt_len;++i) {
				norm += wavecoeffs[aHL+i] * wavecoeffs[aHL+i] / (sigmav * sigmav);
			}
			te =(norm - (double) dwt_len)/(double) dwt_len;
			ct = pow(log((double) dwt_len)/log(2.0),1.5)/sqrt((double) dwt_len);
			//printf("tv %g te %g ct %g \n\n",tv,te,ct);
			if (te < ct) {
				tdv = tv;
			} else {
				x_sum = 0.0;

				for(i = 0; i < dwt_len;++i) {
					dout[i] = fabs(wavecoeffs[aHL+i]/sigmav);
				}
				qsort(dout, dwt_len, sizeof(double), compare_double);
				for(i = 0; i < dwt_len;++i) {
					dout[i] = (dout[i]*dout[i]);
					x_sum += dout[i];
					dsum[i] = x_sum;
				}
				for(i = 0;i < dwt_len;++i) {
					risk[i] = ((double)dwt_len - 2 * ((double)i + 1) +dsum[i] +
							dout[i]*((double)dwt_len - 1 -(double) i))/(double)dwt_len;
				}
				min_index = minindex(risk,dwt_len);
				thr = sqrt(dout[min_index]);
				tdv = thr < tv ? thr : tv;
			}
		}
		//printf("sigmad\n");

		if ( sigmad < 0.00000001) {
			tdd = 0;
		} else {
			tv = sqrt(2.0 * log(dwt_len));
			norm = 0.0;
			for(i = 0; i < dwt_len;++i) {
				norm += wavecoeffs[aHH+i] * wavecoeffs[aHH+i] / (sigmad * sigmad);
			}
			te =(norm - (double) dwt_len)/(double) dwt_len;
			ct = pow(log((double) dwt_len)/log(2.0),1.5)/sqrt((double) dwt_len);
			if (te < ct) {
				tdd = tv;
			} else {
				x_sum = 0.0;

				for(i = 0; i < dwt_len;++i) {
					dout[i] = fabs(wavecoeffs[aHH+i]/sigmad);
				}
				qsort(dout, dwt_len, sizeof(double), compare_double);
				for(i = 0; i < dwt_len;++i) {
					dout[i] = (dout[i]*dout[i]);
					x_sum += dout[i];
					dsum[i] = x_sum;
				}
				for(i = 0;i < dwt_len;++i) {
					risk[i] = ((double)dwt_len - 2 * ((double)i + 1) +dsum[i] +
							dout[i]*((double)dwt_len - 1 -(double) i))/(double)dwt_len;
				}
				min_index = minindex(risk,dwt_len);
				thr = sqrt(dout[min_index]);
				tdd = thr < tv ? thr : tv;
			}
		}

		tdh = tdh * sigmah;
		tdv = tdv * sigmav;
		tdd = tdd * sigmad;

		//printf("tdh %g tdv %g tdd %g \n\n",tdh,tdv,tdd);
		threshrule(wavecoeffs+aLH,dwt_len,tdh,thresh);
		threshrule(wavecoeffs+aHL,dwt_len,tdv,thresh);
		threshrule(wavecoeffs+aHH,dwt_len,tdd,thresh);
		//printf("Threshrule\n");

	}

	getIFWT2(wt,method,wavecoeffs,denoised);

	wave_free(wave);
	free(wavecoeffs);
	free(dout);
	free(lnoise);
	free(dsum);
	free(risk);
	wt2_free(wt);
}

void modwtshrink(double *signal, int N, int J, const char *wname, const char *cmethod, const char *ext, const char *thresh, double *denoised) {
	int filt_len, iter, i, dlen, sgn, MaxIter, it;
	double sigma, td, tmp, M, llen;
	wave_object wave;
	wt_object wt;
	double *dout, *lnoise;

	wave = wave_init(wname);

	filt_len = wave->filtlength;

	MaxIter = (int)(log((double)N / ((double)filt_len - 1.0)) / log(2.0));

	if (J > MaxIter) {
		printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
		exit(-1);
	}

	wt = wt_init(wave, "modwt", N, J);

	if (!strcmp(ext, "sym") && !strcmp(cmethod,"fft")) {
		setWTConv(wt, "fft");
		setDWTExtension(wt, "sym");
	}
	else if (!strcmp(ext, "sym") && !strcmp(cmethod, "direct")) {
		printf("Symmetric Extension is not available for direct method");
		exit(-1);
	}
	else if (!strcmp(ext, "per") && !strcmp(cmethod, "direct"))  {
		setWTConv(wt, "direct");
		setDWTExtension(wt, "per");
	}
	else if (!strcmp(ext, "per") && !strcmp(cmethod, "fft"))  {
		setWTConv(wt, "fft");
		setDWTExtension(wt, "per");
	}
	else {
		printf("Signal extension can be either per or sym");
		exit(-1);
	}

	modwt(wt, signal);

	lnoise = (double*)malloc(sizeof(double)* J);

	//Set sigma

	iter = wt->length[0];
	dlen = wt->length[J];
	dout = (double*)malloc(sizeof(double)* dlen);

	for (it = 0; it < J; ++it) {
		dlen = wt->length[it + 1];
		for (i = 0; i < dlen; ++i) {
			dout[i] = fabs(wt->output[iter + i]);
		}

		sigma = sqrt(2.0) * median(dout, dlen) / 0.6745;
		lnoise[it] = sigma;
		iter += dlen;
	}

	M = pow(2.0,J);
	llen = log((double)wt->modwtsiglength);
	// Thresholding

	iter = wt->length[0];
	for (it = 0; it < J; ++it) {
		sigma = lnoise[it];
		dlen = wt->length[it + 1];
		td = sqrt(2.0 * llen / M) * sigma;

		threshrule(wt->output+iter,dlen,td,thresh);

		iter += wt->length[it + 1];
		M /= 2.0;
	}

	imodwt(wt, denoised);

	free(dout);
	free(lnoise);
	wave_free(wave);
	wt_free(wt);
}

void modwtshrink2(double *data, int rows,int cols, int J, const char *wname, const char *thresh, double *denoised) {
	int filt_len, i, dlen, sgn, MaxIter, it,MaxRows,MaxCols,aLH,aHL,aHH;
	double sigma, tdh,tdv,tdd, tmp, M, llen;
	wave_object wave;
	wt2_object wt;
	double *dout, *lnoise,*wavecoeffs;

	wave = wave_init(wname);

	filt_len = wave->filtlength;

	MaxRows = (int) (log((double)rows / ((double)filt_len - 1.0)) / log(2.0));
	MaxCols = (int) (log((double)cols / ((double)filt_len - 1.0)) / log(2.0));

	MaxIter = (MaxRows < MaxCols) ? MaxRows : MaxCols;

	if (J > MaxIter) {
		printf("\n Error - The Signal Can only be iterated %d times using this wavelet. Exiting\n", MaxIter);
		exit(-1);
	}

	wt = wt2_init(wave, "modwt", rows,cols, J);

	wavecoeffs = modwt2(wt, data);

	lnoise = (double*)malloc(sizeof(double)* J * 3);

	//Set sigma

	dlen = wt->dimensions[2*J-2] * wt->dimensions[2*J-1];

	dout = (double*)malloc(sizeof(double) * dlen);

	noiseest2(wt,wavecoeffs,dout,J,"all",lnoise);

	for (it = 0; it < J * 3; ++it) {
		lnoise[it] *= sqrt(2.0);
	}

	M = pow(2.0,J);
	llen = log((double)wt->rows * wt->cols);
	// Thresholding


	for (it = 0; it < J; ++it) {
		aLH = wt->coeffaccess[3*it+1];
		aHL = wt->coeffaccess[3*it+2];
		aHH = wt->coeffaccess[3*it+3];
		dlen = wt->dimensions[2*it] * wt->dimensions[2*it+1];

		tdh = sqrt(2.0 * llen / M) * lnoise[3*it];
		tdv = sqrt(2.0 * llen / M) * lnoise[3*it+1];
		tdd = sqrt(2.0 * llen / M) * lnoise[3*it+2];

		threshrule(wavecoeffs+aLH,dlen,tdh,thresh);
		threshrule(wavecoeffs+aHL,dlen,tdv,thresh);
		threshrule(wavecoeffs+aHH,dlen,tdd,thresh);

		M /= 2.0;
	}

	imodwt2(wt,wavecoeffs, denoised);

	free(dout);
	free(lnoise);
	wave_free(wave);
	wt2_free(wt);
}


void fdrshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,double qval,const char *thresh,const char *level,double *denoised) {
	int iter,i,dlen,dwt_len,sgn,it,it2;
	double sigma,mu,td,tmp;
	wt_object wt;
	double *dout,*lnoise,*prob,*qvec;
	wave_object wave;
	wave = wave_init(wname);
	wt = getFWT(wave,signal,wname,method,ext,N,J);

	//Set sigma

	lnoise = (double*)malloc(sizeof(double) * J);

	iter = wt->length[0];
	dlen = wt->length[J];

	dout = (double*)malloc(sizeof(double) * dlen);
	prob = (double*)malloc(sizeof(double) * dlen);
	qvec = (double*)malloc(sizeof(double) * dlen);

	noiseest(wt,dout,J,level,lnoise);


	dwt_len = wt->outlength;
	iter = wt->length[0];
	mu = 0.0;
	for(it = 0; it < J;++it) {
		dlen = wt->length[it+1];
		sigma = lnoise[it];
		for (i = 0; i < dlen;++i) {
			dout[i] = fabs(wt->output[iter + i]);
			prob[i] = ((i+1) * qval)/(2 * dlen);
		}
		qsort(dout, dlen, sizeof(double), compare_double_desc);
		qnorm_vec(prob,dlen,mu,sigma,qvec);

		it2 = 0;

		while (dout[it2] >= -qvec[it2]) {
			it2++;
		}

		if (it2 > 0) {
			td = -qvec[it2-1];
		} else {
			td = -qvec[0];
		}

		threshrule(wt->output+iter,dlen,td,thresh);

		iter += wt->length[it+1];
	}

	getIFWT(wt,method,denoised);

	wave_free(wave);
	free(dout);
	free(lnoise);
	free(prob);
	free(qvec);
	wt_free(wt);
}

void neighcoeffshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int iter,i,dlen,dwt_len,sgn,it,it2,L,j;
	double sigma,td,tmp,b1,SS;
	wt_object wt;
	double *dout,*lnoise;
	wave_object wave;
	wave = wave_init(wname);
	wt = getFWT(wave,signal,wname,method,ext,N,J);


	//Set sigma

	lnoise = (double*)malloc(sizeof(double) * J);

	iter = wt->length[0];
	dlen = wt->length[J];

	dout = (double*)malloc(sizeof(double) * (dlen+2));

	noiseest(wt,dout,J,level,lnoise);

	L = 3;
	dwt_len = wt->outlength;
	iter = wt->length[0];

	for(it = 0; it < J;++it) {
		dlen = wt->length[it+1];
		sigma = lnoise[it];
		td = 2 * sigma *sigma * log((double)dlen);

		dout[0] = wt->output[iter+dlen-1];
		dout[dlen+1] = wt->output[iter];

		memcpy(dout+1,wt->output+iter,sizeof(double)*dlen);

		if (!strcmp(thresh, "hard")) {
			for (i = 0; i < dlen; ++i) {
				SS = 0;
				for(j = i; j < i+L;++j) {
					SS += (dout[j]*dout[j]);
				}
				if (SS <= td) {
					wt->output[iter + i] = 0.0;
				}
			}
		} else if (!strcmp(thresh, "soft")) {
			for (i = 0; i < dlen; ++i) {
				SS = 0;
				for(j = i; j < i+L;++j) {
					SS += (dout[j]*dout[j]);
				}
				if (SS != 0) {
					b1 = 1.0 - td/SS;
					tmp = b1 > 0.0 ? b1 : 0.0;
				} else {
					tmp = 0.0;
				}

				wt->output[iter + i] *= tmp;
			}
		}


		iter += wt->length[it+1];
	}

	getIFWT(wt,method,denoised);

	wave_free(wave);
	free(dout);
	free(lnoise);
	wt_free(wt);
}

void neighblockshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int iter,i,dlen,dwt_len,sgn,it,it2,L,j;
	double sigma,td,b1,ssq;
	int L0,L1,L02,blocks;
	wt_object wt;
	double *dout,*lnoise,*S,*beta;
	wave_object wave;
	wave = wave_init(wname);
	wt = getFWT(wave,signal,wname,method,ext,N,J);

	//Set sigma

	lnoise = (double*)malloc(sizeof(double) * J);

	iter = wt->length[0];
	dlen = wt->length[J];

	dout = (double*)malloc(sizeof(double) * (dlen*2));
	S = (double*)malloc(sizeof(double) * dlen);
	beta = (double*)malloc(sizeof(double) * dlen);

	noiseest(wt,dout,J,level,lnoise);

	dwt_len = wt->outlength;
	iter = wt->length[0];

	for(it = 0; it < J;++it) {
		dlen = wt->length[it+1];
		sigma = lnoise[it];

		L0 = (int) floor(log((double)dlen)/2.0);
		L02 = (int) floor((double)L0/2.0);
		L1 = L02 > 1 ? L02:1;
		L = L0 + 2*L1;
		blocks = (int) floor(dlen/L0);

		td = (double) BLOCK_THRESHOLD_COEFFICIENT * sigma *sigma * L;

		for(i = 0; i < L1;++i) {
			dout[dlen+L1+i] = wt->output[iter+i];
			dout[dlen+i] = wt->output[iter+dlen-L1+i];
		}

		memcpy(dout+L1,wt->output+iter,sizeof(double)*dlen);

		for(i = 0; i < blocks;++i) {
			S[i] = 0;
			for(j = i*L0;j < i*L0+L;++j) {
				S[i] += (dout[j]*dout[j]);
			}

			if (!strcmp(thresh, "soft")) {
				b1 = 1.0 - td/S[i];
				beta[i] = b1 > 0.0 ? b1 : 0.0;
			} else if (!strcmp(thresh, "hard")) {
				beta[i] = S[i] > td ? 1 : 0;
			}

		}

		for(i = 0; i < blocks;++i) {
			for(j = i*L0; j < i*L0+L0;++j) {
				wt->output[iter+j] *= beta[i];
			}
		}

		if ( dlen > blocks*L0) {
			L02 = (int) floor((L-(dlen-blocks*L0))/2.0);
			ssq=0;
			for(i = 0; i < L02;++i) {
				ssq+=(wt->output[iter+dlen-L+L02+i]*wt->output[iter+dlen-L+L02+i] + wt->output[iter+i]*wt->output[iter+i]);
			}

			if (!strcmp(thresh, "soft")) {
				b1 = 1.0 - td/ssq;
				if (b1 < 0.0) {
					b1 = 0.0;
				}
			} else if (!strcmp(thresh, "hard")) {
				b1= ssq > td ? 1 : 0;
			}

			for(i = blocks*L0; i < dlen;++i) {
				wt->output[iter+i] *= b1;
			}
		}


		iter += wt->length[it+1];
	}

	getIFWT(wt,method,denoised);

	wave_free(wave);
	free(dout);
	free(lnoise);
	free(S);
	free(beta);
	wt_free(wt);
}

void blockthreshshrink(double *signal,int N,int J,const char *wname,const char *method,const char *ext,const char *thresh,const char *level,double *denoised) {
	int iter,i,dlen,dwt_len,sgn,it,it2,L,j;
	double sigma,td,b1,ssq;
	int blocks,L02;
	wt_object wt;
	double *lnoise,*S,*beta;
	wave_object wave;
	wave = wave_init(wname);
	wt = getFWT(wave,signal,wname,method,ext,N,J);


	//Set sigma

	lnoise = (double*)malloc(sizeof(double) * J);

	iter = wt->length[0];
	dlen = wt->length[J];


	S = (double*)malloc(sizeof(double) * dlen);
	beta = (double*)malloc(sizeof(double) * dlen);

	noiseest(wt,S,J,level,lnoise);

	dwt_len = wt->outlength;
	iter = wt->length[0];

	for(it = 0; it < J;++it) {
		dlen = wt->length[it+1];
		sigma = lnoise[it];

		L = (int) floor(log((double)dlen));
		blocks = (int) floor(dlen/L);

		td = (double) BLOCK_THRESHOLD_COEFFICIENT * sigma *sigma * L;


		for(i = 0; i < blocks;++i) {
			S[i] = 0;
			for(j = i*L;j < (i+1)*L;++j) {
				S[i] += (wt->output[iter+j]*wt->output[iter+j]);
			}

			if (!strcmp(thresh, "soft")) {
				if (S[i] != 0) {
					b1 = 1.0 - td/S[i];
					beta[i] = b1 > 0.0 ? b1 : 0.0;
				} else {
					beta[i] = 0;
				}

			} else if (!strcmp(thresh, "hard")) {
				beta[i] = S[i] > td ? 1 : 0;
			}

			for(j = i*L;j < (i+1)*L;++j) {
				wt->output[iter+j] *= beta[i];
			}

		}


		if ( dlen > blocks*L) {
			L02 = dlen-blocks*L;
			ssq=0;
			for(i = 0; i < L02;++i) {
				ssq+=(wt->output[iter+blocks*L+i]*wt->output[iter+blocks*L+i] );
			}

			if (!strcmp(thresh, "soft")) {
				if ( ssq != 0) {
					b1 = 1.0 - td/ssq;
					if (b1 < 0.0) {
						b1 = 0.0;
					}
				} else {
					b1 = 0;
				}
			} else if (!strcmp(thresh, "hard")) {
				b1= ssq > td ? 1 : 0;
			}
			for(i = 0; i < L02;++i) {
				wt->output[iter+blocks*L+i] *= b1;
			}
		}


		iter += wt->length[it+1];
	}

	getIFWT(wt,method,denoised);

	wave_free(wave);
	free(lnoise);
	free(S);
	free(beta);
	wt_free(wt);
}


void wdenoise(wdenoise_object obj, double *signal,double *denoised) {
	double vscale;
	if(!strcmp(obj->dmethod,"sureshrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("sureshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		sureshrink(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->thresh,obj->level,denoised);
	} else if(!strcmp(obj->dmethod,"visushrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("visushrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		visushrink(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->thresh,obj->level,denoised);
	} else if(!strcmp(obj->dmethod,"minimaxshrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("minimaxshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		minimaxshrink(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->thresh,obj->level,denoised);
	} else if(!strcmp(obj->dmethod,"modwtshrink")) {
		if (strcmp(obj->wmethod, "modwt")) {
			printf("modwtshrink method only works with modwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		modwtshrink(signal,obj->N,obj->J,obj->wname,obj->cmethod,obj->ext,obj->thresh,denoised);
	} else if(!strcmp(obj->dmethod,"fdrshrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("fdrshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		fdrshrink(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->qval,obj->thresh,obj->level,denoised);
	} else if(!strcmp(obj->dmethod,"neighcoeffshrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("neighcoeffshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		neighcoeffshrink(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->thresh,obj->level,denoised);
	} else if(!strcmp(obj->dmethod,"neighblockshrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("neighblockshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		neighblockshrink(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->thresh,obj->level,denoised);
	} else if(!strcmp(obj->dmethod,"blockthreshshrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("blockthreshshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		blockthreshshrink(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->thresh,obj->level,denoised);
	} else if(!strcmp(obj->dmethod,"ebayesshrink")) {
		if (!strcmp(obj->wmethod, "modwt")) {
			printf("ebayesshrink method only works with swt and dwt. Please use setDenoiseWTMethod to set the correct method\n");
			exit(-1);
		}
		vscale = 0.0;
		if (obj->scalefac == -1) {
			if (!strcmp(obj->prior, "cauchy")) {
				printf("Error : You need to specify a scale value for cauchy prior. The dfault is 0.5");
				exit(-1);
			}
			ebayesthresh_wavelet(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->level,vscale,obj->prior,
		NULL,obj->bayesfac,obj->ebthresh,obj->isotone,obj->stabadjustment,obj->uthresh,denoised);
		} else {
			ebayesthresh_wavelet(signal,obj->N,obj->J,obj->wname,obj->wmethod,obj->ext,obj->level,vscale,obj->prior,
		&obj->scalefac,obj->bayesfac,obj->ebthresh,obj->isotone,obj->stabadjustment,obj->uthresh,denoised);
		}

	} else {
		printf("Acceptable Denoising methods are - sureshrink, visushrink, modwtshrink, fdrshrink, \n");
		printf("neighcoeffshrink, neighblockshrink and blockthreshshrink\n");
		exit(-1);
	}
	//ebayesthresh_wavelet(double *signal, int N, int J, char *wname,char *method, char *ext, const  char *level,
	//double vscale, const char* prior, double *a, int bayesfac, const char* thresrule, int isotone, int stabadjustment, int universalthresh, double *denoised)
}

void setWDenoiseMethod(wdenoise_object obj, const char *dmethod) {
	if (!strcmp(dmethod, "sureshrink")) {
		strcpy(obj->dmethod, "sureshrink");
	}
	else if (!strcmp(dmethod, "visushrink")) {
		strcpy(obj->dmethod, "visushrink");
	}
	else if (!strcmp(dmethod, "minimaxshrink")) {
		strcpy(obj->dmethod, "minimaxshrink");
	}
	 else if (!strcmp(dmethod, "modwtshrink")) {
		strcpy(obj->dmethod, "modwtshrink");
	}
	else if (!strcmp(dmethod, "fdrshrink")) {
		strcpy(obj->dmethod, "fdrshrink");
	}
	else if (!strcmp(dmethod, "neighcoeffshrink")) {
		strcpy(obj->dmethod, "neighcoeffshrink");
	}
	else if (!strcmp(dmethod, "neighblockshrink")) {
		strcpy(obj->dmethod, "neighblockshrink");
	}
	else if (!strcmp(dmethod, "blockthreshshrink")) {
		strcpy(obj->dmethod, "blockthreshshrink");
	}
	else if (!strcmp(dmethod, "ebayesshrink")) {
		strcpy(obj->dmethod, "ebayesshrink");
	}
	else {
		printf("Acceptable Denoising methods are - sureshrink, visushrink, modwtshrink, fdrshrink, \n");
		printf("neighcoeffshrink, neighblockshrink, ebayesshrink and blockthreshshrink\n");
		exit(-1);
	}
}

void setWDenoiseWTMethod(wdenoise_object obj, const char *wmethod) {
	if (!strcmp(wmethod, "dwt")) {
		strcpy(obj->wmethod, "dwt");
	}
	else if (!strcmp(wmethod, "swt")) {
		strcpy(obj->wmethod, "swt");
	}
	 else if (!strcmp(wmethod, "modwt")) {
		strcpy(obj->wmethod, "modwt");
	}
	else {
		printf("Wavelet decomposition method can be one of dwt, modwt or swt.\n");
		exit(-1);
	}
}

void setWDenoiseWTExtension(wdenoise_object obj, const char *extension) {
	if (!strcmp(extension, "sym")) {
		strcpy(obj->ext, "sym");
	}
	else if (!strcmp(extension, "per")) {
		strcpy(obj->ext, "per");
	}
	else {
		printf("Signal extension can be either per or sym");
		exit(-1);
	}
}

void setWDenoiseParameters(wdenoise_object obj, const char *thresh,const char *level) {

	//Set thresholding
	if (!strcmp(thresh, "soft")) {
		strcpy(obj->thresh, "soft");
	}
	else if (!strcmp(thresh, "hard")) {
		strcpy(obj->thresh, "hard");
	}
	else {
		printf("Thresholding Method - soft or hard");
		exit(-1);
	}

	// Set Noise estimation at the first level or at all levels

	if (!strcmp(level, "first")) {
		strcpy(obj->level, "first");
	}
	else if (!strcmp(level, "all")) {
		strcpy(obj->level, "all");
	}
	else {
		printf("Noise Estimation at level - first or all");
		exit(-1);
	}
}

void setWDenoiseEbayesParameters(wdenoise_object obj, const char *prior,const char *ebthresh, int *bayesfac, int *isotone, int *universalthresh) {

	if (prior) {
		if (!strcmp(prior, "laplace") || !strcmp(prior, "cauchy")) {
			strcpy(obj->prior,prior);
		} else {
			printf("Only two priors are acceptable - laplace or cauchy");
			exit(-1);
		}
	}

	if (ebthresh) {
		if (!strcmp(ebthresh, "median") || !strcmp(ebthresh, "mean") || !strcmp(ebthresh, "hard") || !strcmp(ebthresh, "soft")) {
			strcpy(obj->ebthresh,ebthresh);
		} else {
			printf("Only four thresholds are acceptable - median, mean, hard or soft");
			exit(-1);
		}
	}


	if (bayesfac) {
		obj->bayesfac = *bayesfac;
	}

	if (isotone) {
		obj->isotone = *isotone;
	}

	if (universalthresh) {
		obj->uthresh = *universalthresh;
	}

}

void setWDenoiseEbayesScale(wdenoise_object obj, double *scalefac) {

	if (scalefac) {
		obj->scalefac = *scalefac;
	} else {
		obj->scalefac = -1;
		//printf("%g \n",obj->scalefac);
	}

}

void setWDenoiseQval(wdenoise_object obj,double qval) {
	if (qval > 0 && qval <= 0.5) {
		obj->qval = qval;
	} else {
		printf("qval values should be between 0 and 0.5 in this implementation");
		exit(-1);
	}
}

void wdenoise_free(wdenoise_object object) {
	free(object);
}
