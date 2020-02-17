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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
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

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}

	visushrink(inp,N,J,wname,method,ext,thresh,level,oup);

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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
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

	N = 951;
	J = 4;

	sig = (double*)malloc(sizeof(double)* N);
	inp = (double*)malloc(sizeof(double)* N);
	oup = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
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

int main() {
	printf("VisuShrink \n");
	visushrinktest();
	printf("FDRShrink \n");
	fdrshrinktest();
	printf("NeighCoeffShrink \n");
	neighcoeffshrinktest();
	printf("NeighBlockShrink \n");
	neighblockshrinktest();
	return 0;
}
