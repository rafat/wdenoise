#include "test2d.h"

void denoisebayestests() {
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp,*efp,*cfp,*mfp,*hfp,*dfp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db4";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// visushrink and sureshrink
	char *ext = "per";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "first"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("../data/heavysine1024_noisy.txt", "r");
	i = 0;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}

	while (!feof(ifp)) {
		fscanf(ifp, "%lf \n", &temp[i]);
		i++;
	}

	fclose(ifp);
	N = i;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for(i = 0; i < N;++i) {
		inp[i] = temp[i];
	}

	J = 4;
	

	// Exponential Denoising

	obj = wdenoise_init(N,J,wname);
	
	setWDenoiseMethod(obj,"ebayesshrink");
	
	
	setWDenoiseWTMethod(obj,method);
	setWDenoiseWTExtension(obj,ext);

	setWDenoiseEbayesScale(obj,NULL);
	
	setWDenoiseEbayesParameters(obj,"laplace","median",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	efp = fopen("denoisedexponential.txt", "w");

	for (i = 0; i < N; ++i) {
		sig[i] = oup[i];
		fprintf(efp, "%lf \n", oup[i]);
	}

	fclose(efp);

	wdenoise_free(obj);

	// Cauchy Denoising

	obj = wdenoise_init(N,J,wname);

	setWDenoiseMethod(obj,"ebayesshrink");
	
	
	setWDenoiseWTMethod(obj,method);
	setWDenoiseWTExtension(obj,ext);

	setWDenoiseEbayesParameters(obj,"cauchy","median",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	cfp = fopen("denoisedcauchy.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(cfp, "%lf \n", oup[i]);
	}

	fclose(cfp);

	wdenoise_free(obj);

	// Postmean denoising

	obj = wdenoise_init(N,J,wname);

	setWDenoiseMethod(obj,"ebayesshrink");
	
	
	setWDenoiseWTMethod(obj,method);
	setWDenoiseWTExtension(obj,ext);

	setWDenoiseEbayesParameters(obj,"laplace","mean",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	mfp = fopen("denoisedpostmean.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(mfp, "%lf \n", oup[i]);
	}

	fclose(mfp);

	dfp = fopen("difflaplacemean.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(dfp, "%lf \n", sig[i] - oup[i]);
	}

	fclose(dfp);

	wdenoise_free(obj);

	// Exphard denoising

	obj = wdenoise_init(N,J,wname);

	setWDenoiseMethod(obj,"ebayesshrink");
	
	
	setWDenoiseWTMethod(obj,method);
	setWDenoiseWTExtension(obj,ext);

	setWDenoiseEbayesParameters(obj,"laplace","hard",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	hfp = fopen("denoisedexphard.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(hfp, "%lf \n", oup[i]);
	}

	fclose(hfp);
	
	printf("Denoised Outputs written to files - denoisedexponential.txt,\n");
	printf("denoisedcauchy.txt, denoisedpostmean.txt, denoisedexphard.txt \n");

	wdenoise_free(obj);

	free(sig);
	free(inp);
	free(oup);
}

int main() {
	denoisebayestests();
}
