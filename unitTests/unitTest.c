/*  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License" */
/*  or "3-clause license") */
/*  Please read attached file License.txt */

#include "unitTest.h"

static double RMS_Error(double *data, double *rec, int N) {
	int i;
	double sum = 0;
	for (i = 0; i < N; ++i) {
		sum += (data[i] - rec[i])*(data[i] - rec[i]);
	}
	return N > 1 ? sqrt(sum / ((double)N - 1)) : sqrt(sum);
}

static double REL_Error(double *data, double *rec, int N) {
	int i;
	double sum1 = 0;
	double sum2 = 0;
	for (i = 0; i < N; ++i) {
		sum1 += (data[i] - rec[i])*(data[i] - rec[i]);
		sum2 += data[i] * data[i];
	}
	return sqrt(sum1) / sqrt(sum2);
}

static void beta_laplace_test() {
	int N;
	double a,s,err;
	N = 5;
	double x[5] = { -2, 1, 0, -4, 5 };
	double ey[5] = { 0.889852029651143, -0.380041716606011, -0.561817771773154, 285.459466672351, 15639.8849145429 };

	double y[5] = { 0, 0, 0, 0, 0 };
	a = 0.5;
	s = 1;

	beta_laplace(x, N, s, a, y);

	err = RMS_Error(y, ey, N);

	if (err > (double) UNIT_EPSILON) {
		printf("\n ERROR : Beta Laplace Test Failed. Exiting. \n");
		exit(-1);
	}
}

static void postmean_test() {
	int N;
	double a, s, w, err;
	const char* prior = "laplace";
	N = 5;
	double x[5] = { -2, 1, 0, -4, 5 };
	double ey[5] = { -1.01158962199946,0.270953305249239,0,-3.48800924041643,4.4997151290092 };

	double y[5] = { 0, 0, 0, 0, 0 };
	s = 1;
	a = 0.5;
	w = 0.5;
	postmean(x,N,s,prior,a,w,y);

	err = RMS_Error(y, ey, N);

	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : Postmean Test Failed. Exiting. \n");
		exit(-1);
	}
}

static void postmed_test() {
	int N;
	double a, s, w, err;
	const char* prior = "laplace";
	N = 5;
	double x[5] = { -2, 1, 0, -4, 5 };
	double ey[5] = { -0.829992882781227,0,0,-3.49568406354978,4.49992059554046 };

	double y[5] = { 0, 0, 0, 0, 0 };
	s = 1;
	a = 0.5;
	w = 0.5;
	postmed(x, N, s, prior, a, w, y);

	err = RMS_Error(y, ey, N);

	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : Postmed Test Failed. Exiting. \n");
		exit(-1);
	}
}

static void tfromw_test() {
	int N, bayesfac;
	double a, s, err;
	const char* prior = "laplace";

	N = 4;
	bayesfac = 0;
	double w[4] = { 0.2, 0.4, 0.6, 0.8 };
	double ey[4] = { 2.44873377028853, 1.92279064562172, 1.40956187155098, 0.767900790087879 };

	double y[4] = { 0, 0, 0, 0 };
	s = 1;
	a = 0.5;

	tfromw(w, N, s, prior, bayesfac, a, y);

	err = RMS_Error(y, ey, N);

	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : Tfromw Test Failed. Exiting. \n");
		exit(-1);
	}

}

static void tfromx_test() {
	int N, uthresh, bayesfac;
	double s, a, tt,err, y;
	const char* prior = "laplace";
	N = 100;
	s = 1.0;
	uthresh = 1;
	bayesfac = 0;
	a = 0.5;
	double x[100] = { -0.5604756, -0.2301775, 1.558708, 0.07050839, 0.1292877, 1.715065, 0.4609162, -1.265061, -0.6868529, -0.445662,
		1.224082, 0.3598138, 0.4007715, 0.1106827, -0.5558411, 1.786913, 0.4978505, -1.966617, 0.7013559, -0.4727914, -1.067824, -0.2179749,
		-1.026004, -0.7288912, -0.6250393, -1.686693, 0.837787, 0.1533731, -1.138137, 1.253815, 0.4264642, -0.2950715, 0.8951257, 0.8781335, 0.8215811,
		0.6886403, 0.5539177, -0.06191171, -0.3059627, -0.380471, -0.694707, -0.2079173, -1.265396, 2.168956, 1.207962, -1.123109, -0.4028848, -0.4666554,
		0.7799651, -0.08336907, 0.2533185, -0.02854676, -0.04287046, 1.368602, -0.225771, 1.516471, -1.548753, 0.5846137, 0.1238542, 0.2159416, 0.3796395,
		-0.5023235, -0.3332074, -1.018575, -1.071791, 0.3035286, 0.4482098, 0.05300423, 0.9222675, 2.050085, -0.4910312, -2.309169, 1.005739, -0.7092008,
		-0.6880086, 1.025571, -0.284773, -1.220718, 0.1813035, -0.1388914, 0.005764186, 0.3852804, -0.37066, 0.6443765, -0.2204866, 0.331782, 1.096839,
		0.4351815, -0.3259316, 1.148808, 0.9935039, 0.548397, 0.2387317, -0.6279061, 1.360652, -0.6002596, 2.187333, 1.532611, -0.2357004, -1.026421 };

	y = 3.03485425654799;

	tfromx(x, N, s, prior, bayesfac, a, uthresh, &tt);
	
	err = RMS_Error(&y, &tt, 1);
	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : Tfromx Test Failed. Exiting. \n");
		exit(-1);
	}
}

static void wfromt_test() {
	int N;
	N = 5;
	double err;
	double tt[5] = { 1,2,3,4,5 };
	double w[5] = { 0.734187187788918,0.368633767549335,0.0661925474440213,0.00348003803260551,6.39312743790131e-05 };
	double y[5] = { 0, 0, 0, 0, 0 };
	double a = 0.5;
	double s = 1;

	wfromt(tt, N, s, "laplace", a, y);
	
	err = RMS_Error(y, w, N);
	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : Wfromt Test Failed. Exiting. \n");
		exit(-1);
	}

}

static void wfromx_test() {
	int N, uthresh;
	double w,err,y,a,s;
	const char* prior = "laplace";
	uthresh = 1;
	N = 100;
	a = 0.5;
	s = 1;
	double x[100] = { -0.5604756, -0.2301775, 1.558708, 0.07050839, 0.1292877, 1.715065, 0.4609162, -1.265061, -0.6868529, -0.445662,
		1.224082, 0.3598138, 0.4007715, 0.1106827, -0.5558411, 1.786913, 0.4978505, -1.966617, 0.7013559, -0.4727914, -1.067824, -0.2179749,
		-1.026004, -0.7288912, -0.6250393, -1.686693, 0.837787, 0.1533731, -1.138137, 1.253815, 0.4264642, -0.2950715, 0.8951257, 0.8781335, 0.8215811,
		0.6886403, 0.5539177, -0.06191171, -0.3059627, -0.380471, -0.694707, -0.2079173, -1.265396, 2.168956, 1.207962, -1.123109, -0.4028848, -0.4666554,
		0.7799651, -0.08336907, 0.2533185, -0.02854676, -0.04287046, 1.368602, -0.225771, 1.516471, -1.548753, 0.5846137, 0.1238542, 0.2159416, 0.3796395,
		-0.5023235, -0.3332074, -1.018575, -1.071791, 0.3035286, 0.4482098, 0.05300423, 0.9222675, 2.050085, -0.4910312, -2.309169, 1.005739, -0.7092008,
		-0.6880086, 1.025571, -0.284773, -1.220718, 0.1813035, -0.1388914, 0.005764186, 0.3852804, -0.37066, 0.6443765, -0.2204866, 0.331782, 1.096839,
		0.4351815, -0.3259316, 1.148808, 0.9935039, 0.548397, 0.2387317, -0.6279061, 1.360652, -0.6002596, 2.187333, 1.532611, -0.2357004, -1.026421 };

	y = 0.0609124723599925;

	w = wfromx(x, N, s, prior, a, uthresh);

	err = RMS_Error(&w, &y, 1);

	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : Wfromx Test Failed. Exiting. \n");
		exit(-1);
	}

	
}

static void wandafromx_test() {
	int N, uthresh;
	double s, a, w,err;
	double y[2] = { 0, 0 };
	double ey[2] = { 0.371583145802847, 3 };
	N = 100;
	s = 1.0;
	uthresh = 1;
	double x[100] = { -0.5604756, -0.2301775, 1.558708, 0.07050839, 0.1292877, 1.715065, 0.4609162, -1.265061, -0.6868529, -0.445662,
		1.224082, 0.3598138, 0.4007715, 0.1106827, -0.5558411, 1.786913, 0.4978505, -1.966617, 0.7013559, -0.4727914, -1.067824, -0.2179749,
		-1.026004, -0.7288912, -0.6250393, -1.686693, 0.837787, 0.1533731, -1.138137, 1.253815, 0.4264642, -0.2950715, 0.8951257, 0.8781335, 0.8215811,
		0.6886403, 0.5539177, -0.06191171, -0.3059627, -0.380471, -0.694707, -0.2079173, -1.265396, 2.168956, 1.207962, -1.123109, -0.4028848, -0.4666554,
		0.7799651, -0.08336907, 0.2533185, -0.02854676, -0.04287046, 1.368602, -0.225771, 1.516471, -1.548753, 0.5846137, 0.1238542, 0.2159416, 0.3796395,
		-0.5023235, -0.3332074, -1.018575, -1.071791, 0.3035286, 0.4482098, 0.05300423, 0.9222675, 2.050085, -0.4910312, -2.309169, 1.005739, -0.7092008,
		-0.6880086, 1.025571, -0.284773, -1.220718, 0.1813035, -0.1388914, 0.005764186, 0.3852804, -0.37066, 0.6443765, -0.2204866, 0.331782, 1.096839,
		0.4351815, -0.3259316, 1.148808, 0.9935039, 0.548397, 0.2387317, -0.6279061, 1.360652, -0.6002596, 2.187333, 1.532611, -0.2357004, -1.026421 };
	wandafromx(x, N, s, uthresh, &a, &w);

	y[0] = w;
	y[1] = a;

	err = RMS_Error(y, ey, 2);

	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : Wandafromx Test Failed. Exiting. \n");
		exit(-1);
	}
}

static void ebayesthresh_test1() {
	int i, N;
	double *y;
	const char *prior = "laplace";
	double a, sdev,err;
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

	double x[100] = {0.253319, -0.028547, -0.042870, 1.368602, -0.225771, 1.516471, -1.548753,
	 0.584614, 0.123854, 0.215942, 0.379639, -0.502323, -0.333207, -1.018575, -1.071791, 0.303529, 
	 0.448210, 0.053004, 0.922267, 2.050085, -0.491031, -2.309169, 1.005739, -0.709201, -0.688009,
	  1.025571, -0.284773, -1.220718, 0.181303, -0.138891, 0.005764, 0.385280, -0.370660, 0.644377, 
	  -0.220487, 0.331782, 1.096839, 0.435181, -0.325932, 1.148808, 0.993504, 0.548397, 0.238732, 
	  -0.627906, 1.360652, -0.600260, 2.187333, 1.532611, -0.235700, -1.026421, -1.831358, -0.203471, 
	  2.870725, -0.206526, -0.693043, 3.385102, 0.136928, -4.198064, -1.753932, 0.027673, 1.872817, 
	  1.327592, -0.816340, 0.165803, -0.592275, 3.874980, 1.101377, -4.573940, 0.553007, -1.969712, 
	  -2.018001, -1.383424, -2.542566, -1.713875, 0.593783, -4.025337, 1.910961, 0.384707, -3.238131, 
	  2.436322, 2.297479, -0.138639, 1.831484, 1.333770, -0.410085, 2.508618, -0.352805, 0.616124, 
	  1.297178, -2.204835, -0.687630, -0.678032, -4.102937, 2.823244, 0.814388, -2.777124, -2.267525,
	   -0.245394, 3.660039, -1.453769};

	i = 0;
	N = 100;
	y = (double*)malloc(sizeof(double)* N);

	a = 0.5;
	bayesfac = 0;
	ls = 1;
	universalthresh = 1;
	stabadjustment = 1;
	isotone = 0;
	sdev = 1.0;

	ebayesthresh(x, N, prior, &a, bayesfac, &sdev, ls, thresrule, universalthresh, stabadjustment, isotone, y);

	err = RMS_Error(y, ey, N);

	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : EBayesThresh Test1 Failed. Exiting. \n");
		exit(-1);
	}

	free(y);
}

static void ebayesthresh_test2() {
	int i, N;
	double temp2[200];
	double *y;
	const char *prior = "laplace";
	double a, sdev, err;
	int bayesfac;
	int ls;
	const char* thresrule = "median";
	int universalthresh;
	int stabadjustment;
	int isotone;

	double ey[100] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.04121693067389,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.35978570570315, 0, 0, -0.240101857101125,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		3.83090910100726, 0.0712822438458154, 0, -3.59696604125997,
		3.86540675250451, 0.0897481774495907, 0, 0, 0, 0, -3.41370542189694,
		0, 0, 3.49859567543283, 0, 0, -1.81030820004674, 3.71037602480949,
		0, 0, 0, 2.81312247449243, 0, 0, 0, -2.00442016176308, -4.6661776361885,
		0, 0, 1.67416061651496, 0, 0, -1.53846644772252, 0, 5.8288454853458,
		0, 0, 0, 0, 0, 0, 3.00394257291387, 0, 0, -2.26047995381615, 0, 0 };

	double x[100] = {-1.331166, -0.358887, 0.492588, 0.230233, 1.073996, -0.658229, 0.887909,
	    0.947728, 0.265894, -0.523127, -0.730413, 0.903830, -2.255313, 1.582136, -0.217579, 1.741538,
	    -1.674214, -1.234555, -0.720547, 0.816514, -1.763461, -1.145425, -0.119363, -2.394993,
	   -1.064722, 0.770347, -1.999774, 1.080433, 1.441736, -0.682557, -1.141903, -0.117389, -0.579715,
	    -0.374162, -0.193691, 0.039320, 1.880382, -0.203125, -0.404977, 0.300170, 0.450253, 0.407507,
	    0.164689, -1.254269, 1.332580, 0.772674, -0.815840, 1.229292, -0.026717, 0.933490, 0.162605,
		-0.147363, -1.572679, 4.332919, 1.963001, -0.001016, -4.101710, 4.367169, 1.966740, -0.620258,
		0.482401, 1.260424, -0.677547, -3.922599, -0.729914, -1.177269, 4.005277, -1.053414, 1.761779,
		-2.632946, 4.213529, -0.074253, 0.110951, -1.026372, 3.364423, -0.207552, 0.384041, 1.220266,
		-2.751659, -5.166236, 1.507626, -0.129277, 2.555716, 0.895012, -0.014008, -2.483422, 
		-1.476252, 6.328846, 0.586972, -0.694025, 0.074610, 0.646651, -0.484505, 0.957027, 3.534982, 
		0.630696, -0.186338, -2.924747, 1.305072, -0.391821};

	i = 0;
	N = 100;
	

	y = (double*)malloc(sizeof(double)* N);

	a = 0.5;
	bayesfac = 0;
	ls = 1;
	universalthresh = 1;
	stabadjustment = 1;
	isotone = 0;
	sdev = 1.0;

	ebayesthresh(x, N, prior, &a, bayesfac, &sdev, ls, thresrule, universalthresh, stabadjustment, isotone, y);

	err = RMS_Error(y, ey, N);

	if (err > (double)UNIT_EPSILON) {
		printf("\n ERROR : EBayesThresh Test2 Failed. Exiting. \n");
		exit(-1);
	}

	free(y);
}

void runUnitTests() {
	printf("Running Unit Tests : \n \n");

	printf("Running Beta Laplace Test ... \n");
	beta_laplace_test();
	printf("PASSED : Beta Laplace Test \n");

	printf("Running Postmean Test ... \n");
	postmean_test();
	printf("PASSED : Postmean Test \n");

	printf("Running Postmed Test ... \n");
	postmed_test();
	printf("PASSED : Postmed Test \n");

	printf("Running Tfromw Test ... \n");
	tfromw_test();
	printf("PASSED : Tfromw Test \n");

	printf("Running Tfromx Test ... \n");
	tfromx_test();
	printf("PASSED : Tfromx Test \n");

	printf("Running Wfromt Test ... \n");
	wfromt_test();
	printf("PASSED : Wfromt Test \n");

	printf("Running Wfromx Test ... \n");
	wfromx_test();
	printf("PASSED : Wfromx Test \n");

	printf("Running Wandafromx Test ... \n");
	wfromx_test();
	printf("PASSED : Wandafromx Test \n");

	printf("Running EBayesThresh Test 1 ... \n");
	ebayesthresh_test1();
	printf("PASSED : EBayesThresh Test 1 \n");

	printf("Running EBayesThresh Test 2 ... \n");
	ebayesthresh_test2();
	printf("PASSED : EBayesThresh Test 2 \n");
}

int main() {
	runUnitTests();
	printf("\n All Unit Tests Passed \n");
	return 0;
}
