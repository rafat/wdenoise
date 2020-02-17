#include "test2d.h"

#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#endif

void visushrink2test() {
	int rows, cols, channels, N, K, J;
	int i;
	size_t img_size;
	double maxr,maxg,maxb;
	double mse_noisy,psnr_noisy,mse_denoised,psnr_denoised;
	double *inp,*inpr,*inpg,*inpb, *oup, *oupr, *oupg, *oupb;
	unsigned char *p,*pd,*outimg;
	unsigned char *img = stbi_load("../data/Boats_noisy15.png",&cols,&rows,&channels,0);
	char *wname = "sym5";
	char *method = "dwt";
	char *ext = "sym";
	char *thresh = "soft";
	char *level = "all";
	char *oupimagegray = "denoised_visu_gray.png";
	char *oupimagergb = "denoised_visu_rgb.png";


	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	J = 4;
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

		rescale(oup,N);

		i = 0;
		for(pd = outimg; pd != outimg +img_size;pd+=channels,p+=channels) {
			*pd = (uint8_t) oup[i];
			if (channels == 2) {
				*(pd+1) = *(p+1);
			}
			i++;
		}

		stbi_write_png(oupimagegray,cols,rows,channels,outimg,cols*channels);

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

		stbi_write_png(oupimagergb,cols,rows,channels,outimg,cols*channels);

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

void ebayesshrink2test() {
	int rows, cols, channels, N, K, J;
	int i;
	size_t img_size;
	double maxr,maxg,maxb;
	double *inp,*inpr,*inpg,*inpb, *oup, *oupr, *oupg, *oupb;
	unsigned char *p,*pd,*outimg;
	unsigned char *img = stbi_load("../data/Boats_noisy15.png",&cols,&rows,&channels,0);
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
	char *oupimagegray = "denoised_ebayes_gray.png";
	char *oupimagergb = "denoised_ebayes_rgb.png";


	printf("Image Loaded. Width : %d Height : %d Channels : %d \n", cols, rows, channels);

	J = 4;
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

		//sureshrink2(inp,rows,cols,J,wname,method,ext,thresh,level,oup);

		ebayesthresh_wavelet2(inp,rows,cols,J,wname,method,ext,level,vscale,prior,&a,bayesfac,thresrule,isotone,stabadjustment,universalthresh,oup);

		//modwtshrink2(inp,rows,cols,J,wname,thresh,oup);

		p = img;
		rescale(oup,N);

		i = 0;
		for(pd = outimg; pd != outimg +img_size;pd+=channels,p+=channels) {
			*pd = (uint8_t) oup[i];
			if (channels == 2) {
				*(pd+1) = *(p+1);
			}
			i++;
		}

		stbi_write_png(oupimagegray,cols,rows,channels,outimg,cols*channels);

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

		//stbi_write_bmp("building10_denoised.bmp",cols,rows,channels,outimg);
		stbi_write_png(oupimagergb, cols, rows, channels, outimg, cols*channels);

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

int main() {
	visushrink2test();
	ebayesshrink2test();
	return 0;
}
