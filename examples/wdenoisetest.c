#include "test2d.h"

void visushrinktest() {
    double *sig, *inp, *oup;
	int i, N, J;
	FILE *ifp,*ofp;

	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "soft";
	char *level = "first";
	char *denoise_method = "minimax";

	ifp = fopen("pieceregular1024.txt", "r");
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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		sig[i] = temp[i];
	}

	ifp = fopen("PieceRegular10.txt", "r");
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

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}

	//denoise(inp, N, J, denoise_method, wname, method, ext, thresh, level, oup);

	// Alternative to denoise_object
	// Just use visushrink and sureshrink functions
	visushrink(inp,N,J,wname,method,ext,thresh,level,oup);
	//sureshrink(inp,N,J,wname,method,ext,thresh,level,oup);

	ofp = fopen("denoiseds.txt", "w");
	for (i = 0; i < N; ++i) {
		fprintf(ofp, "%lf \n", oup[i]);
	}
	fclose(ofp);


	printf("Signal - Noisy Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, inp));
	printf("Corr Coeff %g\n", corrcoef(N, sig, inp));

	printf("Signal - DeNoised Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, oup));
	printf("Corr Coeff %g\n", corrcoef(N, sig, oup));

	free(sig);
	free(inp);
	free(oup);
}

void fdrshrinktest() {
    double *sig, *inp, *oup;
	int i, N, J;
	double qval;
	FILE *ifp,*ofp;

	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "hard";
	char *level = "first";
	char *denoise_method = "minimax";
	qval = 0.05;

	ifp = fopen("pieceregular1024.txt", "r");
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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		sig[i] = temp[i];
	}

	ifp = fopen("PieceRegular10.txt", "r");
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

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}

	//denoise(inp, N, J, denoise_method, wname, method, ext, thresh, level, oup);

	// Alternative to denoise_object
	// Just use visushrink and sureshrink functions
	fdrshrink(inp,N,J,wname,method,ext,qval,thresh,level,oup);
	//sureshrink(inp,N,J,wname,method,ext,thresh,level,oup);

	ofp = fopen("denoisedfdr05.txt", "w");
	for (i = 0; i < N; ++i) {
		fprintf(ofp, "%lf \n", oup[i]);
	}
	fclose(ofp);


	printf("Signal - Noisy Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, inp));
	printf("Corr Coeff %g\n", corrcoef(N, sig, inp));

	printf("Signal - DeNoised Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, oup));
	printf("Corr Coeff %g\n", corrcoef(N, sig, oup));

	free(sig);
	free(inp);
	free(oup);
}

void neighcoeffshrinktest() {
    double *sig, *inp, *oup;
	int i, N, J;
	FILE *ifp,*ofp;

	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "hard";
	char *level = "first";

	ifp = fopen("pieceregular1024.txt", "r");
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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		sig[i] = temp[i];
	}

	ifp = fopen("PieceRegular10.txt", "r");
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

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}

	//denoise(inp, N, J, denoise_method, wname, method, ext, thresh, level, oup);

	// Alternative to denoise_object
	// Just use visushrink and sureshrink functions
	neighcoeffshrink(inp,N,J,wname,method,ext,thresh,level,oup);
	//sureshrink(inp,N,J,wname,method,ext,thresh,level,oup);

	ofp = fopen("denoisedneighcoeffhard.txt", "w");
	for (i = 0; i < N; ++i) {
		fprintf(ofp, "%lf \n", oup[i]);
	}
	fclose(ofp);


	printf("Signal - Noisy Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, inp));
	printf("Corr Coeff %g\n", corrcoef(N, sig, inp));

	printf("Signal - DeNoised Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, oup));
	printf("Corr Coeff %g\n", corrcoef(N, sig, oup));

	free(sig);
	free(inp);
	free(oup);
}

void neighblockshrinktest() {
    double *sig, *inp, *oup;
	int i, N, J;
	FILE *ifp,*ofp;

	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "soft";
	char *level = "first";

	ifp = fopen("pieceregular1024.txt", "r");
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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		sig[i] = temp[i];
	}

	ifp = fopen("PieceRegular10.txt", "r");
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

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}

	//denoise(inp, N, J, denoise_method, wname, method, ext, thresh, level, oup);

	// Alternative to denoise_object
	// Just use visushrink and sureshrink functions
	neighblockshrink(inp,N,J,wname,method,ext,thresh,level,oup);
	//sureshrink(inp,N,J,wname,method,ext,thresh,level,oup);

	ofp = fopen("denoisedneighblocksoft.txt", "w");
	for (i = 0; i < N; ++i) {
		fprintf(ofp, "%lf \n", oup[i]);
	}
	fclose(ofp);


	printf("Signal - Noisy Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, inp));
	printf("Corr Coeff %g\n", corrcoef(N, sig, inp));

	printf("Signal - DeNoised Signal Stats \n");
	printf("RMSE %g\n", rmse(N, sig, oup));
	printf("Corr Coeff %g\n", corrcoef(N, sig, oup));

	free(sig);
	free(inp);
	free(oup);
}

void shrinktest() {
    double *inp, *oup;
	int i, N, J;
	double qval;
	FILE *ifp,*ofp;

	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "soft";
	char *level = "first";
	//char *denoise_method = "minimax";
	qval = 0.05;


	N = 1024;
	J = 5;

	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);


	ifp = fopen("heavysine1024_noisy.txt", "r");
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

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}

	//denoise(inp, N, J, denoise_method, wname, method, ext, thresh, level, oup);

	// Alternative to denoise_object
	// Just use visushrink and sureshrink functions
	//fdrshrink(inp,N,J,wname,method,ext,qval,thresh,level,oup);
	blockthreshshrink(inp,N,J,wname,method,ext,thresh,level,oup);

	ofp = fopen("denoisedheavisinenblockthresh.txt", "w");
	for (i = 0; i < N; ++i) {
		fprintf(ofp, "%lf \n", oup[i]);
	}
	fclose(ofp);

	free(inp);
	free(oup);
}

void sorttest() {
	int i;
	int N = 10;
	double a[10] = {-2,-7,4,2,1,89,0,-23,6.4,-7};

	qsort(a, N, sizeof(double), compare_double_desc);

	for(i = 0;i < N;++i) {
		printf("%g ",a[i]);
	}
}

void qnormvectest() {
	int N,i;
	double q, mu, sigma;
	double *p,*qn;
	N = 10;
	q=0.05;

	p = (double*)malloc(sizeof(double)*N);
	qn = (double*)malloc(sizeof(double)*N);

	mu = 0;
	sigma = 2.7;

	for(i = 1; i <= N;++i) {
		p[i-1] = (i * q)/(2 * N); 
	}

	for(i = 0;i < N;++i) {
		printf("%g ",p[i]);
	}
	printf("\n");

	qnorm_vec(p,N,mu,sigma,qn);

	for(i = 0;i < N;++i) {
		printf("%g ",qn[i]);
	}
	printf("\n");

	qnorm_vec(p,N,mu,1.0,qn);

	for(i = 0;i < N;++i) {
		printf("%g ",qn[i]);
	}
	printf("\n");

	free(p);
	free(qn);
}

void neighblock(double *x, int N,double sigma) {
	int L0,L1,L,L02,blocks,i,j;
	double td,b1,ssq;
	double *temp,*S,*beta;

	L0 = (int) floor(log((double)N)/2.0);
	L02 = (int) floor((double)L0/2.0);
	L1 = L02 > 1 ? L02:1;
	L = L0 + 2*L1;
	blocks = (int) floor(N/L0);

	temp = (double*)malloc(sizeof(double)*(N+2*L1));
	S = (double*)calloc(blocks,sizeof(double));
	beta = (double*)calloc(blocks,sizeof(double));

	td = (double) BLOCK_THRESHOLD_COEFFICIENT * sigma *sigma * L;


	for(i = 0; i < L1;++i) {
		temp[N+L1+i] = x[i]; 
		temp[i] = x[N-L1+i];
	}

	for(i = 0; i < N;++i) {
		temp[L1+i] = x[i];
	}

	for(i = 0; i < blocks;++i) {
		for(j = i*L0;j < i*L0+L;++j) {
			S[i] += (temp[j]*temp[j]);
		}
		b1 = 1.0 - td/S[i];
		beta[i] = b1 > 0.0 ? b1 : 0.0;
		printf("SS : %g beta : %g \n",S[i],beta[i]);
	}

	for(i = 0; i < blocks;++i) {
		for(j = i*L0; j < i*L0+L0;++j) {
			x[j] *= beta[i];
		}
	}

	if ( N > blocks*L0) {
		L02 = (int) floor((L-(N-blocks*L0))/2.0);
		ssq=0;
		for(i = 0; i < L02;++i) {
			ssq+=(x[N-L+L02+i]*x[N-L+L02+i] + x[i]*x[i]);
		}
		b1 = 1.0 - td/ssq;
		if (b1 < 0.0) {
			b1 = 0.0;
		}

		for(i = blocks*L0; i < N;++i) {
			x[i] *= b1;
		}
	}

	free(temp);
	free(S);
	free(beta);
}

void neighcoeff(double *x, int N,double sigma) {
	int L,i,j;
	double td,b1;
	double *temp,*S,*beta;

	L = 3;

	td = 2 * sigma *sigma * log((double)N);

	temp = (double*)malloc(sizeof(double)*(N+2));
	S = (double*)calloc(N,sizeof(double));

	temp[0] = x[N-1];
	temp[N+1] = x[0];

	memcpy(temp+1,x,sizeof(double)*N);

	for(i = 0; i < N;++i) {
		for(j = i; j < i+L;++j) {
			S[i] += (temp[j]*temp[j]);
		}
		if (S[i] <= td) {
			x[i] = 0.0;
		}
	}


	free(temp);
	free(S);
}

void makeHeavisine(int N, double *oup) {
	int i;
	double t;// 1/N:1/N:N

	for(i = 0; i < N;++i) {
		t = (double)(i+1)/(double)N;
		oup[i] = 4.*sin(4*W_PI*t) - signx(t - .3) - signx(.72 - t) + 5;
		oup[i] = (0.6/9.0)*oup[i] + 0.2;
	}
}

void neightests() {
	int i, N;
	double sigma;
	double *x;

	N = 100;

	x = (double*)malloc(sizeof(double)*N);

	makeHeavisine(N,x);

	sigma = 1.0;

	neighblock(x,N,sigma);

	printf("\n\n");
	for(i = 0; i < N;++i) {
		printf("%g ",x[i]);
	}

	free(x);
}

int wdenoisetest() {
	// gcc -Wall -I../header -L../Bin denoisetest.c -o denoise -lwauxlib -lwavelib -lm
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// visushrink and sureshrink
	char *ext = "sym";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("pieceregular1024.txt", "r");
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

	ifp = fopen("PieceRegular10.txt", "r");
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
	setWDenoiseMethod(obj,"fdrshrink");// sureshrink is also the default. The other option with dwt and swt is visushrink.
	// modwt works only with modwtshrink method
	setWDenoiseWTMethod(obj,method);// Default is dwt. the other options are swt and modwt
	setWDenoiseWTExtension(obj,ext);// Default is sym. the other option is per
	setWDenoiseParameters(obj,thresh,level);// Default for thresh is soft. Other option is hard
	// Default for level is all. The other option is first

	wdenoise(obj,inp,oup);

	// Alternative to denoise_object
	// Just use visushrink, modwtshrink and sureshrink functions
	//visushrink(inp,N,J,wname,method,ext,thresh,level,oup);
	//sureshrink(inp,N,J,wname,method,ext,thresh,level,oup);
	// modwtshrink(sig,N,J,wname,cmethod,ext,thresh,oup); See modwtdenoisetest.c
	//ofp = fopen("denoiseds.txt", "w");


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
	return 0;
}

void ebayesparamstest() {
	int i,N,J;
	wdenoise_object obj;

	char *wname = "db5";
	double scalefac = 0.8;
	int bayesfac = 1;
	int isotone = 0;
	int uthresh = 1;

	obj = wdenoise_init(N,J,wname);

	// void setWDenoiseEbayesParameters(wdenoise_object obj, const char *prior,const char *ebthresh, double scalefac,
	// int *bayesfac, int *isotone, int *stabadjustment, int *universalthresh);

	setWDenoiseEbayesParameters(obj,"cauchy","mean",&bayesfac,&isotone,&uthresh);

	//setWDenoiseEbayesScale(obj,NULL);

	printf("%s %s %g %d %d %d %d",obj->prior,obj->ebthresh,obj->scalefac,obj->bayesfac,obj->isotone,obj->stabadjustment,obj->uthresh);

	wdenoise_free(obj);
}

int ebayesshrinktest() {
	// gcc -Wall -I../header -L../Bin denoisetest.c -o denoise -lwauxlib -lwavelib -lm
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db5";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// visushrink and sureshrink
	char *ext = "sym";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("pieceregular1024.txt", "r");
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

	ifp = fopen("PieceRegular10.txt", "r");
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
	
	setWDenoiseMethod(obj,"ebayesshrink");// sureshrink is also the default. The other option with dwt and swt is visushrink.
	
	// modwt works only with modwtshrink method
	setWDenoiseWTMethod(obj,method);// Default is dwt. the other options are swt and modwt
	setWDenoiseWTExtension(obj,ext);// Default is sym. the other option is per
	//setWDenoiseParameters(obj,thresh,level);// Default for thresh is soft. Other option is hard
	// Default for level is all. The other option is first

	setWDenoiseEbayesParameters(obj,"laplace","mean",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	// Alternative to denoise_object
	// Just use visushrink, modwtshrink and sureshrink functions
	//visushrink(inp,N,J,wname,method,ext,thresh,level,oup);
	//sureshrink(inp,N,J,wname,method,ext,thresh,level,oup);
	// modwtshrink(sig,N,J,wname,cmethod,ext,thresh,oup); See modwtdenoisetest.c
	//ofp = fopen("denoiseds.txt", "w");


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
	return 0;
}

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

	ifp = fopen("heavysine1024_noisy.txt", "r");
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

	wdenoise_free(obj);

	free(sig);
	free(inp);
	free(oup);
}

void denoiseexponentialtest() {
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp,*efp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db4";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// visushrink and sureshrink
	char *ext = "sym";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("PieceRegular10.txt", "r");//fopen("ebdata.txt", "r");
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
	
	setWDenoiseEbayesParameters(obj,"laplace","median",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	efp = fopen("denoisedexponential.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(efp, "%lf \n", oup[i]);
	}

	fclose(efp);

	wdenoise_free(obj);

	free(sig);
	free(inp);
	free(oup);
}

void denoisecauchytest() {
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp,*efp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db4";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// visushrink and sureshrink
	char *ext = "sym";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("PieceRegular10.txt", "r");//fopen("ebdata.txt", "r");
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
	
	setWDenoiseEbayesParameters(obj,"cauchy","median",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	efp = fopen("denoisedcauchy.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(efp, "%lf \n", oup[i]);
	}

	fclose(efp);

	wdenoise_free(obj);

	free(sig);
	free(inp);
	free(oup);
}

void denoisepostmeantest() {
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp,*efp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db4";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// visushrink and sureshrink
	char *ext = "sym";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("PieceRegular10.txt", "r");//fopen("ebdata.txt", "r");
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
	
	setWDenoiseEbayesParameters(obj,"laplace","mean",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	efp = fopen("denoisedpostmean.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(efp, "%lf \n", oup[i]);
	}

	fclose(efp);

	wdenoise_free(obj);

	free(sig);
	free(inp);
	free(oup);
}

void denoiseexphardtest() {
	double *sig,*inp,*oup;
	int i,N,J;
	FILE *ifp,*efp;

	wdenoise_object obj;
	double temp[2400];

	char *wname = "db4";
	char *method = "dwt";// Available - dwt, swt and modwt. modwt works only with modwtshrink. The other two methods work with
	// visushrink and sureshrink
	char *ext = "sym";// sym and per work with dwt. swt and modwt only use per extension when called through denoise.
	// You can use sy extension if you directly call modwtshrink with cmethod set to fft. See modwtdenoisetest.c file 
	char *thresh = "soft";// soft or hard
	char *level = "all"; // noise estimation at "first" or "all" levels. modwt only has the option of "all" 

	ifp = fopen("PieceRegular10.txt", "r");//fopen("ebdata.txt", "r");
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
	
	setWDenoiseEbayesParameters(obj,"laplace","soft",NULL,NULL,NULL);

	wdenoise(obj,inp,oup);

	efp = fopen("denoisedexphard.txt", "w");

	for (i = 0; i < N; ++i) {
		fprintf(efp, "%lf \n", oup[i]);
	}

	fclose(efp);

	wdenoise_free(obj);

	free(sig);
	free(inp);
	free(oup);
}

static void ebayesthresh_lctests() {
	int i, N;
	FILE *ifp;
	double temp2[200];
	double *x, *y, *z,*z2,*z3;
	const char *prior = "laplace";
	double a, sdev, err;
	int bayesfac;
	int ls;
	const char* thresrule = "median";
	int universalthresh;
	int stabadjustment;
	int isotone;

	double ey[100] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0.461501932251305, 0, -1.18081131503738,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0.87949485647377, 0, 0, 0, 0, 0, 2.18630700818864, 0, 0,
		2.83745453837867, 0, -3.69477351123594, 0, 0, 0, 0, 0, 0, 0,
		3.36467539856769, 0, -4.07317568417785, 0, -0.141700670497741,
		-0.343683094401104, 0, -1.65608219851322, 0, 0, -3.5191992174366,
		0, 0, -2.66613200774895, 1.45289341778632, 1.15387225239285,
		0, 0, 0, 0, 1.59321391647936, 0, 0, 0, -0.925903664068815, 0, 0,
		-3.59828116583749, 2.11687940619439, 0, -2.04737147690392,
		-1.08310726061405, 0, 3.13925967166472, 0 };

	ifp = fopen("e1.csv", "r");
	i = 0;
	N = 100;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}

	for (i = 0; i < N; ++i) {
		fscanf(ifp, "%lf \n", &temp2[i]);
	}

	fclose(ifp);

	x = (double*)malloc(sizeof(double)* N);
	y = (double*)malloc(sizeof(double)* N);
	z = (double*)malloc(sizeof(double)* N);
	z2 = (double*)malloc(sizeof(double)* N);
	z3 = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		x[i] = temp2[i];
	}
	a = 0.5;
	bayesfac = 0;
	ls = 1;
	universalthresh = 1;
	stabadjustment = 1;
	isotone = 0;
	sdev = 1.0;

	ebayesthresh(x, N, prior, &a, bayesfac, &sdev, ls, thresrule, universalthresh, stabadjustment, isotone, y);

	prior = "cauchy";

	ebayesthresh(x, N, prior, &a, bayesfac, &sdev, ls, thresrule, universalthresh, stabadjustment, isotone, z);

	prior = "laplace";

	thresrule = "mean";

	ebayesthresh(x, N, prior, &a, bayesfac, &sdev, ls, thresrule, universalthresh, stabadjustment, isotone, z2);

	thresrule = "hard";

	ebayesthresh(x, N, prior, &a, bayesfac, &sdev, ls, thresrule, universalthresh, stabadjustment, isotone, z3);

	for (i = 0; i < N; ++i) {
		printf(" %g %g %g %g \n", y[i], z[i],z2[i],z3[i]);
	}

	free(x);
	free(y);
	free(z);
	free(z2);
	free(z3);
}

static void stb_read() {
	int width, height, channels;

	unsigned char *img = stbi_load("lena512_20db.bmp",&width,&height,&channels,0);
	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", width, height, channels);

	stbi_image_free(img);
}

static void convert_to_gray() {
	int width, height, channels;
	size_t img_size, gray_image_size;
	int gray_channels;

	unsigned char *img = stbi_load("empire.jpg",&width,&height,&channels,0);
	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", width, height, channels);

	gray_channels = channels == 4 ? 2 : 1;

	img_size = width * height * channels;
	gray_image_size = width * height * gray_channels;

	unsigned char* gray_img = malloc(gray_image_size);

	for (unsigned char *p = img, *pg = gray_img; p != img + img_size; p += channels, pg += gray_channels) {
		*pg = (uint8_t) ((*p + *(p + 1) + *(p + 2)) / 3.0);
		if (channels == 4) {
			*(pg + 1) = *(p + 3);
		}
	}

	stbi_write_jpg("empire_gray.jpg", width, height, gray_channels, gray_img, 100);

	stbi_image_free(img);
	stbi_image_free(gray_img);
}

static void stb_read_write() {
	int width, height, channels;

	unsigned char *img = stbi_load("lenna512.gif",&width,&height,&channels,0);
	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", width, height, channels);

	stbi_write_png("lenna512.png", width, height, channels, img, width*channels);

	stbi_image_free(img);
}

static double max_and_clean(double *x, int N) {
	double max = 0.0;

	for (int i = 0; i < N;++i) {
		if (x[i] < 0.0) {
			x[i] = 0.0;
		}
		if (x[i] > max) {
			max = x[i];
		}

	}

	return max;
}

static void rescale(double *x,int N) {
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

static void rescale2(double *x, int N) {
	int i;

	for(i = 0; i < N;++i) {
		x[i] = x[i] < 0.0 ? 0.0 : (x[i] > 255.0 ? 255.0 : x[i]);
	}
}

static void imagedwttest() {
    wave_object obj;
	wt2_object wt;
	int i, k, J, rows, cols,channels, N, N2;
	size_t img_size,out_img_size;
	double *inpr,*inpg,*inpb, *wavecoeffsr,*wavecoeffsg,*wavecoeffsb;
	double *cLLr,*cLLg,*cLLb,*cLHr,*cLHg,*cLHb;
	double maxr,maxg,maxb;
	int ir, ic;

	unsigned char *img = stbi_load("lenna512.gif",&cols,&rows,&channels,0);
	unsigned char *outimg,*outimg2;
	unsigned char *p;

    N = rows*cols;
	img_size = N * channels;

	char *name = "db2";
	obj = wave_init(name);// Initialize the wavelet
	inpr = (double*)calloc(N, sizeof(double));
	inpg = (double*)calloc(N, sizeof(double));
	inpb = (double*)calloc(N, sizeof(double));
    
	J = 1;

	p = img;

	for (i = 0; i < N; ++i) {
		inpr[i] = (double) *p;
		inpg[i] = (double) *(p+1);
		inpb[i] = (double) *(p+2);
		p+=4;
	}

	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	// Red Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsr = dwt2(wt, inpr);


	cLLr = getWT2Coeffs(wt, wavecoeffsr, J, "A", &ir, &ic);

	cLHr = getWT2Coeffs(wt, wavecoeffsr, J, "H", &ir, &ic);


	wt2_free(wt);


	// green Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsg = dwt2(wt, inpg);

	cLLg = getWT2Coeffs(wt, wavecoeffsg, J, "A", &ir, &ic);

	cLHg = getWT2Coeffs(wt, wavecoeffsg, J, "H", &ir, &ic);

	wt2_free(wt);

	// Blue Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsb = dwt2(wt, inpb);

	cLLb = getWT2Coeffs(wt, wavecoeffsb, J, "A", &ir, &ic);

	cLHb = getWT2Coeffs(wt, wavecoeffsb, J, "H", &ir, &ic);

	wt2_free(wt);

	N2 = ir * ic;

	
	out_img_size = ir * ic * 3;

	outimg = malloc(out_img_size);

	maxr = max_and_clean(cLLr,N2);

	maxg = max_and_clean(cLLg,N2);

	maxb = max_and_clean(cLLb,N2);

	for (p = outimg; p != outimg + out_img_size; p += 3,cLLr++,cLLg++,cLLb++) {
		*p = (uint8_t) (255.0 * *cLLr/maxr);
		*(p+1) = (uint8_t) (255.0 * *cLLg/maxg);
		*(p+2) = (uint8_t) (255.0 * *cLLb/maxb);
	}

	stbi_write_png("lenna512J1.png", ic, ir, 3, outimg, ic*3);

	outimg2 = malloc(out_img_size);

	maxr = max_and_clean(cLHr,N2);

	maxb = max_and_clean(cLHb,N2);

	maxg = max_and_clean(cLHg,N2);

	//amax = absmax(diff, rows*cols);

	for (p = outimg2; p != outimg2 + out_img_size; p += 3,cLHr++,cLHg++,cLHb++) {
		*p = (uint8_t) (255.0 * *cLHr/maxr);
		*(p+1) = (uint8_t) (255.0 * *cLHg/maxg);
		*(p+2) = (uint8_t) (255.0 * *cLHb/maxb);
	}

	stbi_write_png("lenna512J1H.png", ic, ir, 3, outimg2, ic*3);

	wave_free(obj);
	free(inpr);
	free(inpg);
	free(inpb);
	free(wavecoeffsr);
	free(wavecoeffsg);
	free(wavecoeffsb);
	stbi_image_free(img);
	stbi_image_free(outimg);
	stbi_image_free(outimg2);
}

static void imagedwttest2() {
    wave_object obj;
	wt2_object wt;
	int i, k, J, rows, cols,channels, N, N2;
	size_t img_size,out_img_size;
	double *inpr,*inpg,*inpb, *wavecoeffsr,*wavecoeffsg,*wavecoeffsb;
	double *cLLr,*cLLg,*cLLb,*cLHr,*cLHg,*cLHb;
	double maxr,maxg,maxb;
	int ir, ic;

	unsigned char *img = stbi_load("empire.jpg",&cols,&rows,&channels,0);
	unsigned char *outimg,*outimg2;
	unsigned char *p;

    N = rows*cols;
	img_size = N * channels;

	char *name = "db2";
	obj = wave_init(name);// Initialize the wavelet
	inpr = (double*)calloc(N, sizeof(double));
	inpg = (double*)calloc(N, sizeof(double));
	inpb = (double*)calloc(N, sizeof(double));
    
	J = 1;

	p = img;

	for (i = 0; i < N; ++i) {
		inpr[i] = (double) *p;
		inpg[i] = (double) *(p+1);
		inpb[i] = (double) *(p+2);
		p+=channels;
	}

	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	// Red Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsr = dwt2(wt, inpr);


	cLLr = getWT2Coeffs(wt, wavecoeffsr, J, "A", &ir, &ic);

	cLHr = getWT2Coeffs(wt, wavecoeffsr, J, "D", &ir, &ic);


	wt2_free(wt);


	// green Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsg = dwt2(wt, inpg);

	cLLg = getWT2Coeffs(wt, wavecoeffsg, J, "A", &ir, &ic);

	cLHg = getWT2Coeffs(wt, wavecoeffsg, J, "D", &ir, &ic);

	wt2_free(wt);

	// Blue Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsb = dwt2(wt, inpb);

	cLLb = getWT2Coeffs(wt, wavecoeffsb, J, "A", &ir, &ic);

	cLHb = getWT2Coeffs(wt, wavecoeffsb, J, "D", &ir, &ic);

	wt2_free(wt);

	N2 = ir * ic;

	
	out_img_size = ir * ic * 3;

	outimg = malloc(out_img_size);

	maxr = max_and_clean(cLLr,N2);

	maxg = max_and_clean(cLLg,N2);

	maxb = max_and_clean(cLLb,N2);

	for (p = outimg; p != outimg + out_img_size; p += 3,cLLr++,cLLg++,cLLb++) {
		*p = (uint8_t) (255.0 * *cLLr/maxr);
		*(p+1) = (uint8_t) (255.0 * *cLLg/maxg);
		*(p+2) = (uint8_t) (255.0 * *cLLb/maxb);
	}

	stbi_write_png("empireJ1.png", ic, ir, 3, outimg, ic*3);

	outimg2 = malloc(out_img_size);

	maxr = max_and_clean(cLHr,N2);

	maxb = max_and_clean(cLHb,N2);

	maxg = max_and_clean(cLHg,N2);

	//amax = absmax(diff, rows*cols);

	for (p = outimg2; p != outimg2 + out_img_size; p += 3,cLHr++,cLHg++,cLHb++) {
		*p = (uint8_t) (255.0 * *cLHr/maxr);
		*(p+1) = (uint8_t) (255.0 * *cLHg/maxg);
		*(p+2) = (uint8_t) (255.0 * *cLHb/maxb);
	}

	stbi_write_png("empireJ1D.png", ic, ir, 3, outimg2, ic*3);

	wave_free(obj);
	free(inpr);
	free(inpg);
	free(inpb);
	free(wavecoeffsr);
	free(wavecoeffsg);
	free(wavecoeffsb);
	stbi_image_free(img);
	stbi_image_free(outimg);
	stbi_image_free(outimg2);
}

static double absmax(double *array, int N) {
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

static void imrecontest() {
	wave_object obj;
	wt2_object wt;
	int i, k, J, rows, cols,channels, N, N2;
	size_t img_size;
	double *inpr,*inpg,*inpb, *wavecoeffsr,*wavecoeffsg,*wavecoeffsb;
	double *oupr, *oupg, *oupb,*diff;
	double maxr,maxg,maxb;
	int ir, ic;

	unsigned char *img = stbi_load("empire.jpg",&cols,&rows,&channels,0);
	unsigned char *p;

    N = rows*cols;
	img_size = N * channels;

	char *name = "db2";
	obj = wave_init(name);// Initialize the wavelet
	inpr = (double*)calloc(N, sizeof(double));
	inpg = (double*)calloc(N, sizeof(double));
	inpb = (double*)calloc(N, sizeof(double));

	oupr = (double*)calloc(N, sizeof(double));
	oupg = (double*)calloc(N, sizeof(double));
	oupb = (double*)calloc(N, sizeof(double));

	diff = (double*)calloc(N, sizeof(double));
    
	J = 1;

	p = img;

	for (i = 0; i < N; ++i) {
		inpr[i] = (double) *p;
		inpg[i] = (double) *(p+1);
		inpb[i] = (double) *(p+2);
		p+=channels;
	}

	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	// Red Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsr = dwt2(wt, inpr);

	idwt2(wt,wavecoeffsr,oupr);

	for(i = 0; i < N;++i) {
		diff[i] = oupr[i] - inpr[i];
	}

	maxr = absmax(diff,N);

	wt2_free(wt);


	// green Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsg = dwt2(wt, inpg);

	idwt2(wt,wavecoeffsg,oupg);

	for(i = 0; i < N;++i) {
		diff[i] = oupg[i] - inpg[i];
	}

	maxg = absmax(diff,N);

	wt2_free(wt);

	// Blue Element

	wt = wt2_init(obj, "dwt", rows,cols, J);

	wavecoeffsb = dwt2(wt, inpb);

	idwt2(wt,wavecoeffsb,oupb);

	for(i = 0; i < N;++i) {
		diff[i] = oupb[i] - inpb[i];
	}

	maxb = absmax(diff,N);

	wt2_free(wt);

	printf("Recon maxr : %g maxg : %g maxb : %g",maxr,maxg,maxb);

	wave_free(obj);
	free(inpr);
	free(inpg);
	free(inpb);
	free(oupr);
	free(oupg);
	free(oupb);
	free(diff);
	stbi_image_free(img);
}

static void addNoise(double *img,double std,int rows,int cols,long int i,double *noisyimg) {
	int size,j;
	double a,b,pi,nd;

	mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) i);
	size = rows * cols;
	pi = 3.141592653589793115997963468544185161590576171875;
	for(j = 0; j < size;++j) {
		a=mt_genrand_res53();
		b=mt_genrand_res53();
		nd = (double)(std)*sqrt(-2.0*log(a))*cos(2.0*pi*b);
		noisyimg[j] = img[j] + nd;
	}

}

static void getMSEandPSNR(double *img, double *outimg,int rows, int cols, int channels, double *mse, double *psnr) {
	int i;
	double sum = 0.0,temp,mval;

	mval = -10;

	for (i = 0; i < rows*cols*channels;++i) {
		temp = img[i] - outimg[i];
		sum += temp*temp;
		if (img[i] > mval) {
			mval = img[i];
		}
	}

	*mse = sum/(rows*cols*channels);

	temp = mval*mval/(*mse);
	*psnr = 10 * log(temp)/log(10);
}

void addNoisetest() {
	int i, rows, cols, channels,N;
	double sigma,mse,psnr;
	size_t img_size;
	unsigned char *outimg;
	double *inp, *oup;
	unsigned char *img = stbi_load("Boats.png",&cols,&rows,&channels,0);

	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	N = rows*cols;
	img_size = N * channels;
	outimg = malloc(img_size);

	inp = (double*)calloc(N*channels, sizeof(double));
	oup = (double*)calloc(N*channels, sizeof(double));

	// Convert to double

	for (i = 0; i < (int) img_size;++i) {
		inp[i] = (double) img[i];
	}

	printf("Input to Double \n");

	// Add noise

	sigma = 15;

	for(i = 0; i < channels;++i) {
		addNoise(inp+i*N,sigma,rows,cols,(long int)i,oup+i*N);
	}

	printf("Noise added\n");

	rescale2(oup,N*channels);

	for (i = 0; i < (int) img_size;++i) {
		outimg[i] = (uint8_t) oup[i];
	}

	printf("outimg converted to unit8 \n");

	stbi_write_png("Boats_noisy15.png", cols, rows, channels, outimg, cols*channels);

	getMSEandPSNR(inp,oup,rows,cols,channels,&mse,&psnr);

	printf("MSE %g PSNR %g \n",mse,psnr);

	free(inp);
	free(oup);
	stbi_image_free(img);
	stbi_image_free(outimg);
}

void visushrink2test() {
	int rows, cols, channels, N, K, J;
	int i;
	size_t img_size;
	double maxr,maxg,maxb;
	double *inp,*inpr,*inpg,*inpb, *oup, *oupr, *oupg, *oupb;
	unsigned char *p,*pd,*outimg;
	unsigned char *img = stbi_load("lena512_20db.bmp",&cols,&rows,&channels,0);
	char *wname = "sym5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "soft";
	char *level = "all";

	
	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	J = 2;
    N = rows*cols;
	img_size = N * channels;
	outimg = malloc(img_size);

	if (channels == 1 || channels == 2) {
		K = 1;
	} else if (channels == 3 || channels == 4) {
		K = 3;
	} else {
		printf("This example only supports images with 1,2,3 or 4 channels. \n");
		exit(-1);
	}

	if (K == 1) {
		inp = (double*)calloc(N, sizeof(double));
		oup = (double*)calloc(N, sizeof(double));
		p = img;

		for(i = 0; i < N;++i) {
			inp[i] = *p;
			p+=channels;
		}

		visushrink2(inp,rows,cols,J,wname,method,ext,thresh,level,oup);

		p = img;
		for(pd = outimg; pd != outimg +img_size;pd+=channels,p+=channels,oup++) {
			if (*oup > 255) {
				*oup = 255;
			}
			if (*oup < 0) {
				*oup = 0;
			}
			*pd = (uint8_t) *oup;
			if (channels == 2) {
				*(pd+1) = *(p+1);
			}
		}

		stbi_write_bmp("lena512visu.bmp",cols,rows,channels,outimg);

		free(inp);
		free(oup);
	} else {
		inpr = (double*)calloc(N, sizeof(double));
		inpg = (double*)calloc(N, sizeof(double));
		inpb = (double*)calloc(N, sizeof(double));
		oupr = (double*)calloc(N, sizeof(double));
		oupg = (double*)calloc(N, sizeof(double));
		oupb = (double*)calloc(N, sizeof(double));

		p = img;

		for (i = 0; i < N; ++i) {
			inpr[i] = (double) *p;
			inpg[i] = (double) *(p+1);
			inpb[i] = (double) *(p+2);
			p+=channels;
		}

		

		visushrink2(inpr,rows,cols,J,wname,method,ext,thresh,level,oupr);
		
		visushrink2(inpg,rows,cols,J,wname,method,ext,thresh,level,oupg);
		
		visushrink2(inpb,rows,cols,J,wname,method,ext,thresh,level,oupb);
		

		p = img;
		//maxr = max_and_clean(oupr,N);
		//maxg = max_and_clean(oupg,N);
		//maxb = max_and_clean(oupb,N);

		rescale(oupr,N);
		printf("Rescale R \n");
		rescale(oupg,N);
		printf("Rescale G \n");
		rescale(oupb,N);
		printf("Rescale B \n");
		
		i = 0;
		for(pd = outimg; pd != outimg+img_size;pd+=channels,p+=channels) {
			*pd = (uint8_t) oupr[i];
			*(pd+1) = (uint8_t) oupg[i];
			*(pd+2) = (uint8_t) oupb[i];
			if (channels == 4) {
				*(pd+3) = *(p+3);
			}
			i++;
		}

		stbi_write_bmp("lena512visu2.bmp",cols,rows,channels,outimg);

		printf("Image written i %d\n",i);

		free(inpr);
		free(inpg);
		free(inpb);
		printf("Inputs Freed \n");
		free(oupr);
		free(oupg);
		free(oupb);
		printf("Outputs freed \n");
	}

	

	stbi_image_free(img);
	stbi_image_free(outimg);
}

void shrink2test() {
	int rows, cols, channels, N, K, J;
	int i;
	size_t img_size;
	double maxr,maxg,maxb;
	double *inp,*inpr,*inpg,*inpb, *oup, *oupr, *oupg, *oupb;
	unsigned char *p,*pd,*outimg,*diffimg,*pdiff;
	unsigned char *img = stbi_load("building1_10.png",&cols,&rows,&channels,0);
	char *wname = "sym5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "soft";
	char *level = "first";

	
	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	J = 3;
    N = rows*cols;
	img_size = N * channels;
	outimg = malloc(img_size);
	diffimg = malloc(img_size);

	if (channels == 1 || channels == 2) {
		K = 1;
	} else if (channels == 3 || channels == 4) {
		K = 3;
	} else {
		printf("This example only supports images with 1,2,3 or 4 channels. \n");
		exit(-1);
	}

	if (K == 1) {
		inp = (double*)calloc(N, sizeof(double));
		oup = (double*)calloc(N, sizeof(double));
		p = img;

		for(i = 0; i < N;++i) {
			inp[i] = *p;
			p+=channels;
		}

		sureshrink2(inp,rows,cols,J,wname,method,ext,thresh,level,oup);

		//modwtshrink2(inp,rows,cols,J,wname,thresh,oup);

		p = img;
		for(pd = outimg; pd != outimg +img_size;pd+=channels,p+=channels,oup++) {
			if (*oup > 255) {
				*oup = 255;
			}
			if (*oup < 0) {
				*oup = 0;
			}
			*pd = (uint8_t) *oup;
			if (channels == 2) {
				*(pd+1) = *(p+1);
			}
		}

		stbi_write_bmp("lena512visu.bmp",cols,rows,channels,outimg);

		free(inp);
		free(oup);
	} else {
		inpr = (double*)calloc(N, sizeof(double));
		inpg = (double*)calloc(N, sizeof(double));
		inpb = (double*)calloc(N, sizeof(double));
		oupr = (double*)calloc(N, sizeof(double));
		oupg = (double*)calloc(N, sizeof(double));
		oupb = (double*)calloc(N, sizeof(double));

		p = img;

		for (i = 0; i < N; ++i) {
			inpr[i] = (double) *p;
			inpg[i] = (double) *(p+1);
			inpb[i] = (double) *(p+2);
			p+=channels;
		}

		

		visushrink2(inpr,rows,cols,J,wname,method,ext,thresh,level,oupr);
		//modwtshrink2(inpr,rows,cols,J,wname,thresh,oupr);
		visushrink2(inpg,rows,cols,J,wname,method,ext,thresh,level,oupg);
		//modwtshrink2(inpg,rows,cols,J,wname,thresh,oupg);
		visushrink2(inpb,rows,cols,J,wname,method,ext,thresh,level,oupb);
		//modwtshrink2(inpb,rows,cols,J,wname,thresh,oupb);
		

		p = img;
		pdiff = diffimg;
		//maxr = max_and_clean(oupr,N);
		//maxg = max_and_clean(oupg,N);
		//maxb = max_and_clean(oupb,N);

		rescale(oupr,N);
		printf("Rescale R \n");
		rescale(oupg,N);
		printf("Rescale G \n");
		rescale(oupb,N);
		printf("Rescale B \n");
		
		i = 0;
		for(pd = outimg; pd != outimg+img_size;pd+=channels,p+=channels,pdiff+=channels) {
			*pd = (uint8_t) oupr[i];
			*(pd+1) = (uint8_t) oupg[i];
			*(pd+2) = (uint8_t) oupb[i];
			*pdiff = (uint8_t) fabs(oupr[i] -inpr[i]);
			*(pdiff+1) = (uint8_t) fabs(oupg[i] -inpg[i]);
			*(pdiff+2) = (uint8_t) fabs(oupb[i] -inpb[i]);
			if (channels == 4) {
				*(pd+3) = *(p+3);
				*(pdiff+3) = *(p+3);
			}
			i++;
		}

		//stbi_write_bmp("building10_denoised.bmp",cols,rows,channels,outimg);
		stbi_write_png("building10_denoised.png", cols, rows, channels, outimg, cols*channels);
		stbi_write_png("building10_diff.png", cols, rows, channels, diffimg, cols*channels);

		printf("Image written i %d\n",i);

		free(inpr);
		free(inpg);
		free(inpb);
		printf("Inputs Freed \n");
		free(oupr);
		free(oupg);
		free(oupb);
		printf("Outputs freed \n");
	}

	stbi_image_free(img);
	stbi_image_free(outimg);
	stbi_image_free(diffimg);
}

void ebayesshrink2test() {
	int rows, cols, channels, N, K, J;
	int i;
	size_t img_size;
	double maxr,maxg,maxb;
	double *inp,*inpr,*inpg,*inpb, *oup, *oupr, *oupg, *oupb;
	unsigned char *p,*pd,*outimg,*diffimg,*pdiff;
	unsigned char *img = stbi_load("building1_10.png",&cols,&rows,&channels,0);
	char *wname = "sym4";
	char *method = "dwt";
	char *ext = "sym";
	char *thresrule = "median";
	char *prior = "laplace";
	char *level = "first";
	int isotone = 0;
	int stabadjustment = 1;
	int universalthresh = 1;
	int bayesfac = 0;
	double a = 0.5;
	double vscale = 0;

	
	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	J = 4;
    N = rows*cols;
	img_size = N * channels;
	outimg = malloc(img_size);
	diffimg = malloc(img_size);

	if (channels == 1 || channels == 2) {
		K = 1;
	} else if (channels == 3 || channels == 4) {
		K = 3;
	} else {
		printf("This example only supports images with 1,2,3 or 4 channels. \n");
		exit(-1);
	}

	if (K == 1) {
		inp = (double*)calloc(N, sizeof(double));
		oup = (double*)calloc(N, sizeof(double));
		p = img;

		for(i = 0; i < N;++i) {
			inp[i] = *p;
			p+=channels;
		}

		//sureshrink2(inp,rows,cols,J,wname,method,ext,thresh,level,oup);

		//modwtshrink2(inp,rows,cols,J,wname,thresh,oup);

		p = img;
		for(pd = outimg; pd != outimg +img_size;pd+=channels,p+=channels,oup++) {
			if (*oup > 255) {
				*oup = 255;
			}
			if (*oup < 0) {
				*oup = 0;
			}
			*pd = (uint8_t) *oup;
			if (channels == 2) {
				*(pd+1) = *(p+1);
			}
		}

		stbi_write_bmp("lena512visu.bmp",cols,rows,channels,outimg);

		free(inp);
		free(oup);
	} else {
		inpr = (double*)calloc(N, sizeof(double));
		inpg = (double*)calloc(N, sizeof(double));
		inpb = (double*)calloc(N, sizeof(double));
		oupr = (double*)calloc(N, sizeof(double));
		oupg = (double*)calloc(N, sizeof(double));
		oupb = (double*)calloc(N, sizeof(double));

		p = img;

		for (i = 0; i < N; ++i) {
			inpr[i] = (double) *p;
			inpg[i] = (double) *(p+1);
			inpb[i] = (double) *(p+2);
			p+=channels;
		}

		

		//visushrink2(inpr,rows,cols,J,wname,method,ext,thresh,level,oupr);
		ebayesthresh_wavelet2(inpr,rows,cols,J,wname,method,ext,level,vscale,prior,&a,bayesfac,thresrule,isotone,stabadjustment,universalthresh,oupr);
		ebayesthresh_wavelet2(inpg,rows,cols,J,wname,method,ext,level,vscale,prior,&a,bayesfac,thresrule,isotone,stabadjustment,universalthresh,oupg);
		ebayesthresh_wavelet2(inpb,rows,cols,J,wname,method,ext,level,vscale,prior,&a,bayesfac,thresrule,isotone,stabadjustment,universalthresh,oupb);
		//modwtshrink2(inpr,rows,cols,J,wname,thresh,oupr);
		//visushrink2(inpg,rows,cols,J,wname,method,ext,thresh,level,oupg);
		//modwtshrink2(inpg,rows,cols,J,wname,thresh,oupg);
		//visushrink2(inpb,rows,cols,J,wname,method,ext,thresh,level,oupb);
		//modwtshrink2(inpb,rows,cols,J,wname,thresh,oupb);
		

		p = img;
		pdiff = diffimg;
		//maxr = max_and_clean(oupr,N);
		//maxg = max_and_clean(oupg,N);
		//maxb = max_and_clean(oupb,N);

		rescale(oupr,N);
		printf("Rescale R \n");
		rescale(oupg,N);
		printf("Rescale G \n");
		rescale(oupb,N);
		printf("Rescale B \n");
		
		i = 0;
		for(pd = outimg; pd != outimg+img_size;pd+=channels,p+=channels,pdiff+=channels) {
			*pd = (uint8_t) oupr[i];
			*(pd+1) = (uint8_t) oupg[i];
			*(pd+2) = (uint8_t) oupb[i];
			*pdiff = (uint8_t) fabs(oupr[i] -inpr[i]);
			*(pdiff+1) = (uint8_t) fabs(oupg[i] -inpg[i]);
			*(pdiff+2) = (uint8_t) fabs(oupb[i] -inpb[i]);
			if (channels == 4) {
				*(pd+3) = *(p+3);
				*(pdiff+3) = *(p+3);
			}
			i++;
		}

		//stbi_write_bmp("building10_denoised.bmp",cols,rows,channels,outimg);
		stbi_write_png("building10_denoisedebayes.png", cols, rows, channels, outimg, cols*channels);
		stbi_write_png("building10_diffebayes.png", cols, rows, channels, diffimg, cols*channels);

		printf("Image written i %d\n",i);

		free(inpr);
		free(inpg);
		free(inpb);
		printf("Inputs Freed \n");
		free(oupr);
		free(oupg);
		free(oupb);
		printf("Outputs freed \n");
	}

	stbi_image_free(img);
	stbi_image_free(outimg);
	stbi_image_free(diffimg);
}

int main() {
    //visushrinktest();
	//fdrshrinktest();
	//sorttest();
	//qnormvectest();
	//shrinktest();
	//neightests();
	//neighcoeffshrinktest();
	//neighblockshrinktest();
	//wdenoisetest();
	//ebayesparamstest();
	
	//ebayesshrinktest();
	//denoisebayestests();
	//denoiseexponentialtest();
	//denoisecauchytest();
	//denoisepostmeantest();
	//denoiseexphardtest();
	//ebayesthresh_lctests();
	//stb_read();
	//convert_to_gray();
	//imagedwttest();
	//imagedwttest2();
	//imrecontest();
	//visushrink2test();
	//shrink2test();
	addNoisetest();
	//ebayesshrink2test();
    return 0;
}
