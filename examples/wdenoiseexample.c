#include "test2d.h"

void wdenoisetest() {
	// gcc -Wall -I../header -L../Bin denoisetest.c -o denoise -lwauxlib -lwavelib -lm
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// the rest
	char *ext = "sym";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("../data/pieceregular1024.txt", "r");
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
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for(i = 0; i < N;++i) {
		sig[i] = temp[i];
	}

	ifp = fopen("../data/PieceRegular10.txt", "r");
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

	for(i = 0; i < N;++i) {
		inp[i] = temp[i];
	}
	obj = wdenoise_init(N,J,wname);
	setWDenoiseMethod(obj,"fdrshrink");
	// modwt works only with modwtshrink method
	setWDenoiseWTMethod(obj,method);// Default is dwt. the other options are swt and modwt
	setWDenoiseWTExtension(obj,ext);// Default is sym. the other option is per
	setWDenoiseParameters(obj,thresh,level);// Default for thresh is soft. Other option is hard
	// Default for level is all. The other option is first

	wdenoise(obj,inp,oup);


	printf("Signal - Noisy Signal Stats \n");
	printf("RMSE %g\n",rmse(N,sig,inp));
	printf("Corr Coeff %g\n",corrcoef(N,sig,inp));

	printf("Signal - DeNoised Signal Stats \n");
	printf("RMSE %g\n",rmse(N,sig,oup));
		printf("Corr Coeff %g\n",corrcoef(N,sig,oup));

	free(sig);
	free(inp);
	wdenoise_free(obj);
	free(oup);
}

int main() {
	wdenoisetest();
}

