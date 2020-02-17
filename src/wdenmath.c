#include "wdenmath.h"

int compare_double(const void* a, const void* b)
{
	double arg1 = *(const double*)a;
	double arg2 = *(const double*)b;

	if (arg1 < arg2) return -1;
	if (arg1 > arg2) return 1;
	return 0;

}

int compare_double_desc(const void* a, const void* b)
{
	double arg1 = *(const double*)a;
	double arg2 = *(const double*)b;

	if (arg1 < arg2) return 1;
	if (arg1 > arg2) return -1;
	return 0;

}

void range(double *x, int N, double *max, double *min) {

	qsort(x, N, sizeof(double), compare_double);
	*max = x[N - 1];
	*min = x[0];
}

static int isign(int N) {
	int M;
	if (N >= 0) {
		M = 1;
	}
	else {
		M = -1;
	}

	return M;
}

static int iabs(int N) {
	if (N >= 0) {
		return N;
	}
	else {
		return -N;
	}
}


static double entropy_s(double *x,int N) {
  int i;
  double val,x2;

  val = 0.0;

  for(i = 0; i < N; ++i) {
    if (x[i] != 0) {
      x2 = x[i] * x[i];
      val -= x2 * log(x2);
    }
  }
  return val;
}

static double entropy_t(double *x,int N, double t) {
  int i;
  double val,x2;
  if (t < 0) {
    printf("Threshold value must be >= 0");
    exit(1);
  }
  val = 0.0;

  for(i = 0; i < N; ++i) {
    x2 = fabs(x[i]);
    if (x2 > t) {
      val += 1;
    }

  }

  return val;

}

static double entropy_n(double *x,int N,double p) {
  int i;
  double val,x2;
  if (p < 1) {
    printf("Norm power value must be >= 1");
    exit(1);
  }
  val = 0.0;
  for(i = 0; i < N; ++i) {
    x2 = fabs(x[i]);
    val += pow(x2,(double)p);

  }

  return val;
}

static double entropy_l(double *x,int N) {
  int i;
  double val,x2;

  val = 0.0;

  for(i = 0; i < N; ++i) {
    if (x[i] != 0) {
      x2 = x[i] * x[i];
      val += log(x2);
    }
  }
  return val;
}

static double exparg(int l)
{

	double lb = 0.69314718055995;
	int m = (l == 0) ? DBL_MAX_EXP : DBL_MIN_EXP - 1;

	return m * lb * 0.99999;
} 

double erf__(double x)
{
	/* -----------------------------------------------------------------------
	*             EVALUATION OF THE REAL ERROR FUNCTION
	* ----------------------------------------------------------------------- */

	/* Initialized data */

	double c = .564189583547756;
	double a[5] = { 7.7105849500132e-5, -.00133733772997339,
		.0323076579225834, .0479137145607681, .128379167095513 };
	double b[3] = { .00301048631703895, .0538971687740286,
		.375795757275549 };
	double p[8] = { -1.36864857382717e-7, .564195517478974,
		7.21175825088309, 43.1622272220567, 152.98928504694,
		339.320816734344, 451.918953711873, 300.459261020162 };
	double q[8] = { 1., 12.7827273196294, 77.0001529352295,
		277.585444743988, 638.980264465631, 931.35409485061,
		790.950925327898, 300.459260956983 };
	double r[5] = { 2.10144126479064, 26.2370141675169,
		21.3688200555087, 4.6580782871847, .282094791773523 };
	double s[4] = { 94.153775055546, 187.11481179959,
		99.0191814623914, 18.0124575948747 };

	/* System generated locals */
	double ret_val;

	/* Local variables */
	double t, x2, ax, bot, top;

	ax = fabs(x);
	if (ax <= 0.5) {
		t = x * x;
		top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
		bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;

		return x * (top / bot);
	}
	/* else: ax > 0.5 */

	if (ax <= 4.) { /*  ax in (0.5, 4] */

		top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
			+ p[5]) * ax + p[6]) * ax + p[7];
		bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
			+ q[5]) * ax + q[6]) * ax + q[7];
		ret_val = 0.5 - exp(-x * x) * top / bot + 0.5;
		if (x < 0.0) {
			ret_val = -ret_val;
		}
		return ret_val;
	}

	/* else: ax > 4 */

	if (ax >= 5.8) {
		return x > 0 ? 1 : -1;
	}
	x2 = x * x;
	t = 1.0 / x2;
	top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
	bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
	t = (c - top / (x2 * bot)) / ax;
	ret_val = 0.5 - exp(-x2) * t + 0.5;
	if (x < 0.0) {
		ret_val = -ret_val;
	}
	return ret_val;

} 

double erfc1__(int ind, double x)
{
	/* ----------------------------------------------------------------------- */
	/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

	/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
	/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
	/* ----------------------------------------------------------------------- */

	/* Initialized data */

	double c = .564189583547756;
	double a[5] = { 7.7105849500132e-5, -.00133733772997339,
		.0323076579225834, .0479137145607681, .128379167095513 };
	double b[3] = { .00301048631703895, .0538971687740286,
		.375795757275549 };
	double p[8] = { -1.36864857382717e-7, .564195517478974,
		7.21175825088309, 43.1622272220567, 152.98928504694,
		339.320816734344, 451.918953711873, 300.459261020162 };
	double q[8] = { 1., 12.7827273196294, 77.0001529352295,
		277.585444743988, 638.980264465631, 931.35409485061,
		790.950925327898, 300.459260956983 };
	double r[5] = { 2.10144126479064, 26.2370141675169,
		21.3688200555087, 4.6580782871847, .282094791773523 };
	double s[4] = { 94.153775055546, 187.11481179959,
		99.0191814623914, 18.0124575948747 };

	double ret_val;
	double e, t, w, bot, top;

	double ax = fabs(x);
	if (ax <= 0.5) {
		double t = x * x,
			top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0,
			bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
		ret_val = 0.5 - x * (top / bot) + 0.5;
		if (ind != 0) {
			ret_val = exp(t) * ret_val;
		}
		return ret_val;
	}
	if (ax <= 4.0) {
		top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
			+ p[5]) * ax + p[6]) * ax + p[7];
		bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
			+ q[5]) * ax + q[6]) * ax + q[7];
		ret_val = top / bot;

	}
	else { 
		if (x <= -5.6) {
			ret_val = 2.0;
			if (ind != 0) {
				ret_val = exp(x * x) * 2.0;
			}
			return ret_val;
		}
		if (ind == 0 && (x > 100.0 || x * x > -exparg(1))) {
			return 0.0;
		}
		t = 1. / (x * x);
		top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
		bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
		ret_val = (c - t * top / bot) / ax;
	}

	if (ind != 0) {
		if (x < 0.0)
			ret_val = exp(x * x) * 2.0 - ret_val;
	}
	else {
		w = x * x;
		t = w;
		e = w - t;
		ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
		if (x < 0.0)
			ret_val = 2.0 - ret_val;
	}
	return ret_val;

}

/*
double erf(double x) {
	return erf__(x);
}

double erfc(double x) {
	return erfc1__(0, x);
}
*/
double erfcx(double x) {
	return erfc1__(2, x);

}

double erfinv(double x) {
	double a, t, oup, num, den, xinf, sign, temp, pi;
	if (x > 0.0) {
		sign = 1.0;
	}
	else {
		sign = -1.0;
	}
	a = fabs(x);
	xinf = 1.79e308;
	pi = 3.1415926535897932384626434;

	if (a > 1.) {
		printf("The Input takes values between -1 and +1");
		exit(1);
	}

	if (a <= 0.75) {
		double p[5] = { 4.62680202125696, -16.6805947126248, 17.6230176190819, \
			- 5.4309342907266, 0.236997019142 };
		double q[5] = { 4.26606447606664, -17.5930726990431, 22.7331464544494, \
			- 9.9016863476727, 1.000000000000 };

		t = x*x - 0.75*0.75;
		num = p[0] + t * (p[1] + t * (p[2] + t * (p[3] + t * p[4])));
		den = q[0] + t * (q[1] + t * (q[2] + t * (q[3] + t * q[4])));
		oup = x * num / den;

	}
	else if (0.75 < a && a <= 0.9375) {
		double p[6] = { -0.041199817067782, 0.643729854003468, -3.28902674093993, \
			6.24518431579026, -3.65953596495117, 0.30260114943200 };
		double q[6] = { -0.029324540620124, 0.501148500527886, -2.90144687299145, \
			6.65393051963183, -5.40640580412825, 1.0000000000000 };

		t = x*x - 0.9375*0.9375;
		num = p[0] + t * (p[1] + t * (p[2] + t * (p[3] + t * (p[4] + t * p[5]))));
		den = q[0] + t * (q[1] + t * (q[2] + t * (q[3] + t * (q[4] + t * q[5]))));
		oup = x * num / den;
	}
	else if (0.9375 < a && a  < 1.) {
		double p[6] = { .1550470003116, 1.382719649631, .690969348887, \
			- 1.128081391617, .680544246825, -.16444156791 };
		double q[3] = { .155024849822, 1.385228141995, 1.000000000000 };

		t = 1.0 / sqrt(-log(1.0 - a));
		oup = sign*(p[0] / t + p[1] + t*(p[2] + t*(p[3] + t*(p[4] + t*p[5])))) /
			(q[0] + t*(q[1] + t*(q[2])));
	}
	else {
		oup = sign * xinf;
	}

	if (a < 1.) {
		temp = (erf(oup) - x) / (2 / sqrt(pi) * exp(-pow(oup, 2.0)));
		oup -= temp / (1 + temp*oup);
	}

	return oup;
}

double erfcinv(double x) {
	double oup, xinf, pi, t;

	xinf = 1.79e308;
	pi = 3.1415926535897932384626434;

	if (x > 2. || x < 0.) {
		printf("The Input takes values between 0.0 and 2.0");
		exit(1);
	}

	if (x >= 0.0625 && x <2.0) {
		oup = erfinv(1.0 - x);
	}
	else if (x < 0.0625 && x >= 1.0e-100) {
		double p[6] = { 0.1550470003116, 1.382719649631, 0.690969348887, \
			- 1.128081391617, 0.680544246825, -0.16444156791 };
		double q[3] = { 0.155024849822, 1.385228141995, 1.000000000000 };

		t = 1.0 / sqrt(-log(x));
		oup = (p[0] / t + p[1] + t*(p[2] + t*(p[3] + t*(p[4] + t*p[5])))) /
			(q[0] + t*(q[1] + t*(q[2])));

	}
	else if (x < 1.0e-100 && x > 1.0e-1000) {
		double p[4] = { 0.00980456202915, 0.363667889171, 0.97302949837, -0.5374947401 };
		double q[3] = { 0.00980451277802, 0.363699971544, 1.000000000000 };
		t = 1.0 / sqrt(-log(x));
		oup = (p[0] / t + p[1] + t*(p[2] + t*p[3])) / (q[0] + t*(q[1] + t*(q[2])));

	}
	else if (x <= 1.0e-1000) {
		oup = xinf;
	}
	else if (x == 2.) {
		oup = -xinf;
	}

	return oup;

}

double dnorm(double x, double mu, double sigma) {
	double oup, t;

	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}

	t = x - mu;
	t = t * t;

	oup = exp(-t / (2 * sigma *sigma)) / (sqrt(2 * PIVAL) * sigma);

	return oup;
}

double pnorm(double x, double mu, double sigma) {
	double oup, t;

	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}

	t = (x - mu) / (sqrt(2.0) * sigma);

	oup = 0.5 * (1. + erf(t));
	return oup;
}

double qnorm(double p, double mu, double sigma) {
	double oup, t;
	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}

	if (p < 0. || p > 1.0) {
		printf("Probablity Values can only take values between 0.0 and 1.0");
		exit(1);
	}

	t = sqrt(2.0) * erfinv(2.0 * p - 1.0);
	oup = t * sigma + mu;

	return oup;
}

void dnorm_vec(double *x,int N, double mu, double sigma,double *oup) {
	double t,den1,den2;
	int i;

	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}

	den1 = 2 * sigma *sigma;
	den2 = sqrt(2 * PIVAL) * sigma;

	for(i = 0;i < N; ++i) {
		t = x[i] - mu;
		t = t * t;

		oup[i] = exp(-t / den1) / den2;
	}

}

void pnorm_vec(double *x,int N, double mu, double sigma, double *oup) {
	double t,den1;
	int i;

	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}

	den1 = sqrt(2.0) * sigma;

	for(i = 0; i < N; ++i) {
		t = (x[i] - mu) / den1;

		oup[i] = 0.5 * (1. + erf(t));
	}	
	
}

void qnorm_vec(double *p,int N, double mu, double sigma,double *oup) {
	double t;
	int i;
	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}

	for (i = 0; i < N;++i) {
		if (p[i] < 0. || p[i] > 1.0) {
			printf("Probablity Values can only take values between 0.0 and 1.0");
			exit(1);
		}

		t = sqrt(2.0) * erfinv(2.0 * p[i] - 1.0);
		oup[i] = t * sigma + mu;
	}


}

void beta_cauchy(double *x, int N, double *beta) {
	int i;
	double d0;
	double *phix;

	phix = (double*)malloc(sizeof(double) *N);

	for (i = 0; i < N; ++i) {
		beta[i] = x[i];
		phix[i] = dnorm(x[i], 0, 1);
	}
	d0 = dnorm(0, 0, 1);
	for (i = 0; i < N; ++i) {
		if (beta[i] == 0.0) {
			beta[i] = -0.5;
		}
		else {
			beta[i] = (d0 / phix[i] - 1.0) / (x[i] * x[i]) - 1.0;
		}
	}

	free(phix);
}

void beta_laplace(double *x, int N, double s, double a, double *beta) {
	int i;
	double xpa, xma, temp1,temp2;


	for (i = 0; i < N; ++i) {
		xpa = fabs(x[i])/s + s*a;
		xma = fabs(x[i])/s - s*a;
		temp1 = 1 / xpa;
		if (xpa < 35) {
			temp1 = pnorm(-xpa, 0, 1)/dnorm(xpa,0,1) ;
		}
		temp2 = 1 / fabs(xma);
		if (xma > 35) {
			xma = 35;
		}
		if (xma > -35) {
			temp2 = pnorm(xma, 0, 1) / dnorm(xma, 0, 1);
		}

		beta[i] = (a * s / 2.0) * (temp1 + temp2) - 1.0;
	}

}

void postmean_cauchy(double *x, int N, double w, double *oup) {
	int i;
	double ex,z;

	for (i = 0; i < N; ++i) {
		oup[i] = x[i];
	}

	for (i = 0; i < N; ++i) {
		if (x[i] != 0.0) {
			ex = exp(-x[i] * x[i] / 2.0);
			z = w * (x[i] - (2.0 * (1 - ex)) / x[i]);
			z = z / (w * (1 - ex) + (1 - w) * ex * x[i] * x[i]);
			oup[i] = z;
		}
		
	}

}

void postmean_laplace(double *x, int N, double s, double ax, double w, double *oup) {
	double *wpost,*beta;
	double sx,cp1,cp2,ef,xpa,xma,temp,xf,pmc,a;
	int i;

	wpost = (double*)malloc(sizeof(double) *N);
	beta = (double*)malloc(sizeof(double) *N);
	a = ax;
	if (ax > 20) {
		a = 20;
	}

	beta_laplace(x,N,s,a,beta);

	for (i = 0; i < N; ++i) {
		wpost[i] = 1.0 - (1.0 - w)/(1 + w*beta[i]);
	}

	for (i = 0; i < N; ++i) {
		sx = x[i] > 0 ? 1 : -1;
		if (x[i] == 0) 
			sx = 0;
		xf = fabs(x[i]);
		xma = xf/s - s*a;
		xpa = xf/s + s*a;
		xpa = xpa > 35 ? 35 : xpa;
		xma = xma < -35 ? -35 : xma;
		cp1 = pnorm(xma, 0, 1);
		cp2 = pnorm(-xpa, 0, 1);
		temp = 2 * a * xf;
		ef = temp > 100 ? exp(100) : exp(temp);
		pmc = xf -  a*s*s*(2 * cp1 / (cp1 + ef * cp2) - 1);
		oup[i] = sx * pmc * wpost[i];
	}

	free(wpost);
	free(beta);
}

void postmean(double *x, int N, double s, const char* prior, double a, double w, double *oup) {

	if (!strcmp(prior, "laplace")) {
		postmean_laplace(x, N, s, a, w, oup);
	} else if (!strcmp(prior, "cauchy")) {
		postmean_cauchy(x, N, w, oup);
	}
	else {
		printf("Error : Only two prios are accepted - laplace and cauchy \n");
		exit(1);
	}
}

void cauchy_medzero(double *x, int N,  double *z, double *w, double *res,void *params) {
	int i;
	double hh,dhhh,yleft,yright2,z2;

	for (i = 0; i < N; ++i) {
		hh = z[i] - x[i];
		z2 = z[i] * z[i];
		dhhh = dnorm(hh, 0, 1);
		yleft = pnorm(hh,0,1) - z[i]*dhhh + ((z[i]*x[i] - 1)*dhhh*pnorm(-x[i],0,1)) / dnorm(x[i],0,1);
		yright2 = 1 + exp(-z2 / 2) * (z2 * (1 / w[i] - 1) - 1);
		res[i] = yright2 / 2 - yleft;
	}
}

void vecbinsolv(double *zf, int nz, custom_funcprior *func, double *tlo, double *thi, int nits,double *z,double *w,double *tsol) {
	int i,j;
	double *tmid,*fmid;

	tmid = (double*)malloc(sizeof(double) *nz);
	fmid = (double*)malloc(sizeof(double) *nz);
	for (j = 0; j < nits; ++j) {
		for (i = 0; i < nz; ++i) {
			tmid[i] = (thi[i] + tlo[i]) / 2;
		}
		FUNCPRIOR_EVAL(func, tmid, nz,z, w, fmid);
		for (i = 0; i < nz; ++i) {
			if (fmid[i] <= zf[i]) {
				tlo[i] = tmid[i];
			}
			else {
				thi[i] = tmid[i];
			}
		}
	}

	for (i = 0; i < nz; ++i) {
		tsol[i] = (tlo[i] + thi[i]) / 2;
	}

	free(tmid);
	free(fmid);
}

void postmed_cauchy(double *x, int nx, double wx, double *oup) {
	double *w,*z,*zf,*tlo,*thi,*tsol;
	double ax,amax,sgn;
	int *index;
	int i,sj,nits,j;
	custom_funcprior cauchymed = { cauchy_medzero, 0 };
	nits = 30;

	index = (int*)malloc(sizeof(int) * nx);

	sj = 0;
	amax = 0.0;
	for (i = 0; i < nx; ++i) {
		ax = fabs(x[i]);
		if (amax < ax) {
			amax = ax;
		}
		if (ax < 7) {
			index[i] = 1;
			sj += 1;
		}
		else {
			index[i] = 0;
			oup[i] = ax - 2 / ax;
		}
	}
	if (sj > 0) {
		zf = (double*)malloc(sizeof(double) *sj);
		w = (double*)malloc(sizeof(double) *sj);
		z = (double*)malloc(sizeof(double) *sj);
		tlo = (double*)malloc(sizeof(double) *sj);
		thi = (double*)malloc(sizeof(double) *sj);
		tsol = (double*)malloc(sizeof(double) *sj);


		for (i = 0; i < sj; ++i) {
			zf[i] = 0.0;
			tlo[i] = 0.0;
			thi[i] = amax;
			w[i] = wx;
		}
		j = 0;
		for (i = 0; i < nx; ++i) {
			if (index[i] == 1) {
				z[j] = fabs(x[i]);
				j += 1;
			}
		}

		vecbinsolv(zf, sj, &cauchymed, tlo, thi, nits, z, w, tsol);

		j = 0;
		for (i = 0; i < nx; ++i) {
			if (index[i] == 1) {
				oup[i] = tsol[j];
				j += 1;
			}

		}

		free(w);
		free(z);
		free(tlo);
		free(thi);
		free(zf);
		free(tsol);
	}

	for (i = 0; i < nx; ++i) {
		if (oup[i] < 1e-07) {
			oup[i] = 0;
		}

		sgn = x[i] > 0 ? 1 : -1;

		oup[i] = sgn * oup[i];
	}

	//vecbinsolv(zf,sj,&cauchymed,tlo,thi,nits,z,w,tsol);

	
	free(index);
}

void postmed_laplace(double *x, int N, double s, double ax, double w, double *oup) {
	double a,sx,xma,zz,temp,mucor;
	double *beta,*xabs;
	int i;

	beta = (double*)malloc(sizeof(double)*N);
	xabs = (double*)malloc(sizeof(double)*N);

	a = ax;
	if (ax > 20) {
		a = 20;
	}

	for (i = 0; i < N; ++i) {
		xabs[i] = x[i] >= 0 ? x[i] : -x[i];
	}

	beta_laplace(xabs, N, s, a, beta);

	for (i = 0; i < N; ++i) {
		sx = x[i] > 0 ? 1 : -1;
		if (x[i] == 0)
			sx = 0;
		xma = xabs[i]/s - s*a;
		zz = 1 / a * (1/s*dnorm(xma,0,1)) * (1 / w + beta[i]);
		if (xma > 25)
			zz = 0.5;
		temp = zz > 1 ? 1 : zz;
		mucor = qnorm(temp, 0, 1);
		temp = xma > mucor ? (xma - mucor) : 0;
		oup[i] = temp * sx * s;
	}

	free(beta);
	free(xabs);
}

void postmed(double *x, int N, double s, const char* prior, double a, double w, double *oup) {

	if (!strcmp(prior, "laplace")) {
		postmed_laplace(x, N, s, a, w, oup);
	}
	else if (!strcmp(prior, "cauchy")) {
		postmed_cauchy(x, N, w, oup);
	}
	else {
		printf("Error : Only two prios are accepted - laplace and cauchy \n");
		exit(1);
	}
}

void isotone(double *yn, int N, double *w, int increasing,double *z) {
	// Pool-adjacent-violators algorithm
	// increasing = 1 , the curve is set to increasing.
	// else, the curve is set to decreasing
	// Weight vector w is a vector of positive elements
	int i,j,k,l,N2;
	double w2;
	int *index;
	double *weight,*y;

	if (N == 1) {
		*z = *yn;
	}
	y = (double*)calloc(N, sizeof(double));
	if (increasing == 1) {
		for (i = 0; i < N; ++i) {
			y[i] = yn[i];
		}
	}
	else {
		for (i = 0; i < N; ++i) {
			y[i] = -yn[i];
		}
	}

	index = (int*)calloc(N, sizeof(int));
	weight = (double*)calloc(N, sizeof(double));

	for (i = 0; i < N; ++i) {
		z[i] = 0.0;
	}

	k = 0;
	j = 0;
	index[0] = 0;
	weight[0] = w[0];
	z[0] = y[0];

	for (j = 1; j < N; ++j) {
		k += 1;
		index[k] = j;
		weight[k] = w[j];
		z[k] = y[j];

		l = k - 1;

		if (l < 0) {
			l = 0;
		}

		while (k >= 1 && (z[l] >= z[k])) {
			w2 = weight[k] + weight[k - 1];
			z[k-1] += (weight[k] / w2) * (z[k] - z[k-1]);
			weight[k - 1] = w2;
			k -= 1;
			l = k - 1;

			if (l < 0) {
				l = 0;
			}
		}
	}

	N2 = N-1;

	while (N2 >= 0) {
		for (j = index[k]; j <= N2; ++j) {
			z[j] = z[k];
		}
		N2 = index[k] - 1;
		k -= 1;
	}

	if (increasing != 1) {
		for (i = 0; i < N; ++i) {
			z[i] = -z[i];
		}
	}
	
	free(index);
	free(weight);
	free(y);
}

void wfromt(double *tt, int N, double s, const char* prior, double a, double *w) {
	int i;
	double tma,dn,pn,pv2;
	double *wi, *beta;

	pv2 = sqrt((double)PIVAL / 2.0);
	wi = (double*)calloc(N, sizeof(double));
	beta = (double*)calloc(N, sizeof(double));
	if (!strcmp(prior, "laplace")) {
		for (i = 0; i < N; ++i) {
			tma = tt[i]/s - s*a;
			wi[i] = 1.0 / fabs(tma);
			if (tma > -35) {
				dn = dnorm(tma, 0, 1);
				pn = pnorm(tma, 0, 1);
				wi[i] = pn / dn;
			}
		}
		beta_laplace(tt, N, s, a, beta);
		for (i = 0; i < N; ++i) {
			wi[i] = a*s*wi[i] - beta[i];
			w[i] = 1 / wi[i];
		}

	}
	else if (!strcmp(prior, "cauchy")) {
		for (i = 0; i < N; ++i) {
			tma = tt[i];
			dn = dnorm(tma, 0, 1);
			pn = pnorm(tma, 0, 1);
			wi[i] = (pn - tma*dn - 0.5)/(pv2 * dn * tma * tma);
			if (wi[i] == wi[i]) {
				wi[i] += 1.0;
			}
			else {
				wi[i] = 1.0;
			}
			w[i] = 1.0 / wi[i];
		}
	}


	free(wi);
	free(beta);
}

double wfromx(double *x, int N, double s, const char* prior, double a, int uthresh) {
	int i,j;
	double tuniv,wlo,whi,wmid;
	double *tt,*beta;
	double w, shi,slo, smid;

	tt = (double*)calloc(1, sizeof(double));
	beta = (double*)calloc(N, sizeof(double));

	if (!strcmp(prior, "cauchy")) {
		s = 1.0;
	}

	if (uthresh) {
		tuniv = sqrt(2.0 * log((double)N)) * s;
		wfromt(&tuniv, 1, s, prior, a, tt);
		wlo = tt[0];
	}
	else {
		wlo = 0.0;
	}


	if (!strcmp(prior, "laplace")) {
		beta_laplace(x, N, s, a, beta);
	}
	else if (!strcmp(prior, "cauchy")) {
		beta_cauchy(x, N, beta);
	}

	whi = 1.0;
	shi = 0.0;
	slo = 0.0;
	for (i = 0; i < N; ++i) {
		if (beta[i] > 1.0e+20) {
			beta[i] = 1.0e+20;
		}
		shi += (beta[i] / (1.0 + beta[i]));
		slo += (beta[i] / (1.0 + wlo *beta[i]));
	}

	if (shi >= 0.0) {
		free(tt);
		free(beta);
		w = 1.0;
		return w;
	}

	if (slo <= 0.0) {
		free(tt);
		free(beta);
		return wlo;
	}

	for (j = 0; j < 30; ++j) {
		wmid = sqrt(wlo * whi);
		smid = 0.0;
		for (i = 0; i < N; ++i) {
			smid += (beta[i] / (1.0 + wmid *beta[i]));
		}
		if (smid == 0.0) {
			free(tt);
			free(beta);
			return wmid;
		}

		if (smid > 0.0) {
			wlo = wmid;
		}
		else {
			whi = wmid;
		}
	}

	w = sqrt(wlo * whi);

	free(tt);
	free(beta);

	return w;
}

int wmonfromx(double *x, int N, const char* prior, double a, double tol, int maxits,double *w) {
	/*
	Returns 1 if convergence is achieved 
	.
	Returns 0 if more iterations are needed
	*/
	int i, j, winit;
	double *wmin, *beta, *ps, *ww,*wnew,*rv;
	double aa,mx,mn,zinc,tuniv;

	wmin = (double*)calloc(1, sizeof(double));
	beta = (double*)calloc(N, sizeof(double));
	ps = (double*)calloc(N, sizeof(double));
	ww = (double*)calloc(N, sizeof(double));
	wnew = (double*)calloc(N, sizeof(double));
	rv = (double*)calloc(N, sizeof(double));

	tuniv = sqrt(2.0 * log((double)N));
	wfromt(&tuniv, 1, 1.0,prior, a, wmin);
	winit = 1;

	if (!strcmp(prior, "laplace")) {
		beta_laplace(x, N, 1.0, a, beta);
	}
	else if (!strcmp(prior, "cauchy")) {
		beta_cauchy(x, N, beta);
	}

	for (i = 0; i < N; ++i) {
		w[i] = winit * 1.0;
	}

	for (j = 0; j < maxits; ++j) {
		for (i = 0; i < N; ++i) {
			aa = w[i] + 1.0 / beta[i];
			ps[i] = w[i] + aa;
			ww[i] = 1.0 / (aa*aa);
		}
		isotone(ps, N, ww, 0, wnew);
		for (i = 0; i < N; ++i) {
			if (wnew[i] < *wmin) {
				wnew[i] = *wmin;
			}
		}
		for (i = 0; i < N; ++i) {
			if (wnew[i] > 1) {
				wnew[i] = 1;
			}
			rv[i] = wnew[i] - w[i];
		}

		range(rv, N, &mx, &mn);
		zinc = fabs(mx) > fabs(mn) ? fabs(mx) : fabs(mn);
		if (zinc < tol)
			return 1;
	}

	free(wmin);
	free(beta);
	free(ps);
	free(ww);
	free(wnew);
	free(rv);
	return 0;
}

void cauchy_threshzero(double *z, int N, double *w,double *y) {
	int i;
	double z2, zi;

	for (i = 0; i < N; ++i) {
		z2 = z[i] * z[i];
		zi = z[i];
		y[i] = pnorm(zi,0,1) - z[i] * dnorm(zi,0,1) - 0.5 - (z2 * exp( - z2/2) * (1.0/w[i] - 1.0 ))/2.0;
	}
}

void laplace_threshzero(double *x, int N, double s, double a, double *w, double *z) {
	int i;
	double xma;
	double *beta;

	beta = (double*)malloc(sizeof(double) * N);
	
	a = a < 20.0 ? a : 20.0;

	beta_laplace(x, N, s, a, beta);

	for (i = 0; i < N; ++i) {
		xma = x[i]/s - s*a;
		z[i] = pnorm(xma,0,1) - ((dnorm(xma,0,1)/s) * (1/w[i] + beta[i]))/a;
	}

	free(beta);
}

void beta_cauchy_func(double *x, int N, double *z, double *w, double *beta,void *params) {
	beta_cauchy(x, N, beta);
}

void beta_laplace_func(double *x, int N, double *z, double *w, double *beta, void *params) {
	double a,s;
	a = z[0];
	s = z[1];
	beta_laplace(x, N, s, a, beta);
}

void cauchy_threshzero_func(double *x, int N, double *z, double *w, double *beta, void *params) {
	cauchy_threshzero(x, N, w, beta);
}

void laplace_threshzero_func(double *x, int N, double *z, double *w, double *beta, void *params) {
	double a, s;
	a = z[0];
	s = z[1];
	laplace_threshzero(x, N, s, a,w, beta);
}

void tfromw(double *w, int N, double s, const char* prior,int bayesfac, double a, double *tt) {
	int i;
	double *z;
	double *thi, *tlo;
	double *xpar;
	int nits;
	custom_funcprior funcl1 = { beta_laplace_func, 0 };
	custom_funcprior funcc1 = { beta_cauchy_func, 0 };
	custom_funcprior funcl0 = { laplace_threshzero_func, 0 };
	custom_funcprior funcc0 = { cauchy_threshzero_func, 0 };

	z = (double*)malloc(sizeof(double)* N);
	tlo = (double*)malloc(sizeof(double)* N);
	thi = (double*)malloc(sizeof(double)* N);
	xpar = (double*)malloc(sizeof(double)* 2);

	nits = 30;

	if (bayesfac) {
		for (i = 0; i < N; ++i) {
			z[i] = (1.0 / w[i]) - 2.0;
			tlo[i] = 0.0;
			thi[i] = 10.0;
		}
		if (!strcmp(prior, "laplace")) {
			xpar[0] = a;
			xpar[1] = s;
			vecbinsolv(z,N, &funcl1, tlo, thi, nits, xpar, NULL, tt);
		}
		else if (!strcmp(prior, "cauchy")) {
			vecbinsolv(z, N, &funcc1, tlo, thi, nits, NULL, NULL, tt);
		}
	}
	else {
		for (i = 0; i < N; ++i) {
			z[i] = 0.0;
			tlo[i] = 0.0;
			thi[i] = 10.0;
		}
		if (!strcmp(prior, "laplace")) {
			for (i = 0; i < N; ++i) {
				thi[i] = s*(25.0 + s*a);
			}
			xpar[0] = a;
			xpar[1] = s;
			vecbinsolv(z, N, &funcl0, tlo, thi, nits, xpar, w, tt);
		}
		else if (!strcmp(prior, "cauchy")) {
			vecbinsolv(z, N, &funcc0, tlo, thi, nits, NULL, w, tt);
		}
	}

	free(xpar);
	free(z);
	free(tlo);
	free(thi);
}

void threshld(double *x, int N, double td, const char* thresh,double *z) {
	int i,sgn;

	if (!strcmp(thresh, "hard")) {
		for (i = 0; i < N; ++i) {
			z[i] = fabs(x[i]) >= td ? x[i] : 0;
		}
	}
	else if (!strcmp(thresh, "soft")) {
		for (i = 0; i < N; ++i) {
			if (fabs(x[i]) < td) {
				z[i] = 0.0;
			}
			else {
				sgn = x[i] > 0 ? 1 : -1;
				if (x[i] == 0) {
					sgn = 0;
				}
				z[i] = sgn * (fabs(x[i]) - td);
			}
		}
	}

}

void wpost_laplace(double *w, int N, double *x, double s, double a, double *z) {
	int i;
	double *beta;

	beta = (double*)malloc(sizeof(double) * N);

	beta_laplace(x, N, s, a, beta);

	for (i = 0; i < N; ++i) {
		z[i] = 1.0 - (1.0 - w[i]) / (1.0 + w[i] * beta[i]);
	}

	free(beta);
}

double negloglik_laplace(double *xx, int M, void *params) {
	int i,N;
	double a,sdev;
	double *beta;
	double wlomax, whimax;
	double loglik, temp;
	const char* prior = "laplace";
	double *xpar = (double*)params;

	a = xx[1];
	sdev = (double)xpar[0];
	N = (int)xpar[3];
	beta = (double*)malloc(sizeof(double) * N);
	
	//wlomax = xpar[3]; //tlo
	//whimax = xpar[4]; //thi
	wfromt(xpar + 1, 1, sdev, prior, a, &whimax);
	wfromt(xpar + 2, 1, sdev, prior, a, &wlomax);
	
	//printf("a %g wlo %g \n", a,wlomax);

	loglik = 0;
	beta_laplace(xpar+4, N,sdev, a, beta);

	for (i = 0; i < N; ++i) {
		temp = 1.0 + (xx[0] * (whimax - wlomax) + wlomax) * beta[i];
		//printf("whimax %g wlomax %g beta %g temp %g \n",whimax,wlomax,beta[i],temp);
		loglik += log(temp);
	}

	loglik = -loglik;
	//printf("loglik %g \n", loglik);

	free(beta);
	return loglik;
}

double negloglik_laplace2(double *xx, int M, void *params) {
	int i, N;
	double a, sdev;
	double *beta;
	double w;
	double loglik, temp, temp2;
	const char* prior = "laplace";
	double *xpar = (double*)params;

	a = xx[1];
	sdev = (double)xpar[0];
	N = (int)xpar[3];
	beta = (double*)malloc(sizeof(double) * N);

	//wlomax = xpar[3]; //tlo
	//whimax = xpar[4]; //thi
	wfromt(xx, 1, sdev, prior, a, &w);

	//printf("a %g wlo %g \n", a,wlomax);

	loglik = 0;
	beta_laplace(xpar + 4, N, sdev, a, beta);

	for (i = 0; i < N; ++i) {
		//temp = 1.0 + (xx[0] * (whimax - wlomax) + wlomax) * beta[i];
		//printf("whimax %g wlomax %g beta %g temp %g \n",whimax,wlomax,beta[i],temp);
		if (fabs(xpar[i+4]) < 35) {
			temp = dnorm(xpar[4+i], 0, 1);
			temp2 = w*beta[i] + 1;
		}
		loglik += log(temp) + log(temp2);
	}

	//printf("xx %g", xx[0]);

	loglik = -loglik;
	//printf("loglik %g \n", loglik);

	free(beta);
	return loglik;
}
/*
void wandafromx(double  *xx, int NX,double s,int universalthresh,double *a, double *w) {
	double thi,tlo,wlo,whi;
	int i,N,retval,ls;
	double lo[2] = { 0, 0.04 };
	double hi[2] = { 1, 3 };
	double x[2] = { 0.5, 0.5 };
	double *params;
	const char* prior = "laplace";
	N = 2;
	ls = 1;// Change ls when sdev is a vector
	if (universalthresh) {
		thi = s * sqrt(2.0 * log((double)NX));
	}
	else {
		thi = DBL_MAX;
	}
	printf("thi %g \n", thi);
	params = (double*)malloc(sizeof(double)*(2 * ls + 2 + NX));
	tlo = 0;
	params[0] = s;
	params[1] = tlo;
	params[2] = thi;
	params[3] = (double)NX;
	for (i = 0; i < NX; ++i) {
		params[i + 2 * ls + 2] = xx[i];
	}
	//// params = {s,tlo,thi,NX,xx}
	custom_function loglik_min = { &negloglik_laplace, params };

	//custom_function loglik_min = { &negloglik_laplace2, params };


	retval = lbfgsb_min(&loglik_min, NULL, x, N, 0, lo, hi);

	*a = x[1];
	wfromt(&tlo, 1, s, prior, *a, &whi);
	wfromt(&thi, 1, s, prior, *a, &wlo);
	*w = x[0] * (whi - wlo) + wlo;
	
	free(params);
}
*/

void wandafromx(double  *xx, int NX, double s, int universalthresh, double *a, double *w) {
	double thi, tlo, wlo, whi;
	int i, N, retval, ls;
	double lo[2] = { 0, 0.04 };
	double hi[2] = { 1, 3 };
	double x[2] = { 1.0, 0.5 };
	double *params;
	const char* prior = "laplace";
	N = 2;
	ls = 1;// Change ls when sdev is a vector
	if ( s != 1.0) {
		printf("This implementation of wandafromx only accepts scalar standard deviation value of 1.0. \n");
		printf("Please rescale the data to set the value to 1.0 \n");
		exit(-1);
	}
	if (universalthresh) {
		thi = s * sqrt(2.0 * log((double)NX));
	}
	else {
		thi = DBL_MAX;
	}
	hi[0] = thi;
	//printf("thi %g \n", thi);
	params = (double*)malloc(sizeof(double)*(2 * ls + 2 + NX));
	tlo = 0;
	params[0] = s;
	params[1] = tlo;
	params[2] = thi;
	params[3] = (double)NX;
	for (i = 0; i < NX; ++i) {
		params[i + 2 * ls + 2] = xx[i];
	}
	//// params = {s,tlo,thi,NX,xx}
	//custom_function loglik_min = { &negloglik_laplace, params };

	custom_function loglik_min = { &negloglik_laplace2, params };


	//retval = lbfgsb_min(&loglik_min, NULL, x, N, 0, lo, hi);

	minfun2d(&loglik_min, NULL, x, N, 0, lo, hi);

	*a = x[1];
	wfromt(x, 1, s, prior, *a, w);
	/*
	if (retval == 1) {
	*a = x[1];
	wfromt(&tlo, 1, s, prior, *a, &whi);
	wfromt(&thi, 1, s, prior, *a, &wlo);
	*w = x[0] * (whi - wlo) + wlo;
	}
	else {
	printf("Convergence Not Reached. \n");
	exit(1);
	}
	*/
	free(params);
}

void tfromx(double *x, int N, double s, const char* prior, int bayesfac,double a, int uthresh,double *tt) {
	double w;

	if (!strcmp(prior, "cauchy")) {
		s = 1.0;
	}
	if (!strcmp(prior, "laplace") && a < 0.0) {
		wandafromx(x, N, s, uthresh, &a, &w);
	}
	else {
		w = wfromx(x, N, s, prior, a,uthresh);
	}

	tfromw(&w, 1, s, prior, bayesfac, a, tt);
}

int unique_sort(double *x, int N) {
	int newl,i,iter;

	newl = N;

	qsort(x, N, sizeof(double), compare_double);
	iter = 0;

	for (i = 1; i < N; ++i) {
		if (x[i] == x[iter]) {
			x[i] = DBL_MAX;
			newl--;
		}
		else {
			iter = i;
		}
	}

	qsort(x, N, sizeof(double), compare_double);

	return newl;
}

static int array_max_index(double *array, int N) {
	int i, index;
	double m = array[0];
	index = 0;

	for (i = 1; i < N; ++i) {
		if (array[i] > m) {
			m = array[i];
			index = i;
		}
	}

	return index;
}

double zetafromx(double *xd, int nx, double *cs,const char* prior, double a, double *pilo ) {
	int i,j,nx2,mz,nlmin,nit,sgn,nlmax, index;
	double zmax, tt,ze,ld,rd,zlo,zhi,zstep,likd,tmp,ptemp;
	double *beta,*zs1,*zs2,*zj,*cb,*zlmax,*lmin,*cbil,*cbir,*cbi,*zlmin,*zm,*pz;

	beta = (double*)malloc(sizeof(double) * nx);
	zs1 = (double*)malloc(sizeof(double) * nx);
	zs2 = (double*)malloc(sizeof(double) * nx);
	zj = (double*)malloc(sizeof(double) * 2 * nx);
	cb = (double*)malloc(sizeof(double) * nx);
	zlmax = (double*)malloc(sizeof(double) * 2 * nx);
	zlmin = (double*)malloc(sizeof(double) * 2 * nx);
	pz = (double*)malloc(sizeof(double) * nx);
	

	//  Find the beta values and the minimum weight if necessary

	if (pilo == NULL) {
		tt = sqrt(2.0 * log((double)nx));
		wfromt(&tt, 1, 1.0,prior, a, &ptemp);
		pilo = (double*)calloc(1, sizeof(double));
		*pilo = ptemp;
	}
	
	if (!strcmp(prior, "laplace")) {
		beta_laplace(xd, nx, 1.0,a, beta);
	} else if (!strcmp(prior, "cauchy")) {
		beta_cauchy(xd, nx, beta);
	}

	//  Find jump points zj in derivative of log likelihood as function
	//    of z, and other preliminary calculations

	for (i = 0; i < nx; ++i) {
		zs1[i] = *pilo / cs[i];
		zs2[i] = 1.0 / cs[i];
		zj[i] = zs1[i];
		zj[nx + i] = zs2[i];
		cb[i] = cs[i] * beta[i];
	}

	nx2 = 2 * nx;

	mz = unique_sort(zj, nx2);
	//mdisplay(zj, 1, mz);
	//  Find left and right derivatives at each zj
	//  and check which are local minima
	//  Check internal zj first

	lmin = (double*)calloc(mz, sizeof(double));

	for (j = 1; j < mz - 1; ++j) {
		ze = zj[j];
		cbil = (double*)calloc(nx, sizeof(double));

		for (i = 0; i < nx; ++i) {
			if ((ze > zs1[i]) && (ze <= zs2[i])) {
				cbil[i] = cb[i];
			}
		}
		ld = 0.0;

		for (i = 0; i < nx; ++i) {
			ld += (cbil[i]/(1.0 + ze*cbil[i]));
		}

		if (ld <= 0) {
			cbir = (double*)calloc(nx, sizeof(double));

			for (i = 0; i < nx; ++i) {
				if ((ze >= zs1[i]) && (ze < zs2[i])) {
					cbir[i] = cb[i];
				}

				rd = 0.0;

				for (i = 0; i < nx; ++i) {
					rd += (cbir[i] / (1.0 + ze*cbir[i]));
				}

				if (rd >= 0) {
					lmin[j] = rd;
				}
			}

			free(cbir);
		}

		free(cbil);
	}
	//mdisplay(lmin, 1, mz);

	//  Deal with the two end points in turn, finding right deriv
	//   at lower end and left deriv at upper

	//  In each case the corresponding end point is either a local min
	//   or a local max depending on the sign of the relevant deriv

	cbir = (double*)calloc(nx, sizeof(double));
	cbil = (double*)calloc(nx, sizeof(double));
	nlmax = 0;
	for (i = 0; i < nx; ++i) {
		if (zj[0] == zs1[i]) {
			cbir[i] = cb[i];
		}
	}

	rd = 0.0;

	for (i = 0; i < nx; ++i) {
		rd += (cbir[i] / (1.0 + zj[0]*cbir[i]));
	}

	if (rd > 0) {
		lmin[0] = 1;
	}
	else {
		*zlmax = zj[0];
		nlmax = 1;
	}

	for (i = 0; i < nx; ++i) {
		if (zj[mz-1] == zs2[i]) {
			cbil[i] = cb[i];
		}
	}

	ld = 0.0;

	for (i = 0; i < nx; ++i) {
		ld += (cbil[i] / (1.0 + zj[mz-1]*cbil[i]));
	}

	if (ld < 0) {
		lmin[mz - 1] = 1;
	}
	else {
		*zlmax = zj[mz-1];
		nlmax = 1;
	}

	//mdisplay(zlmax, 1, mz);

	// Flag all local minima and do a binary search between them to find the local maxima
	nlmin = 0;


	for (i = 0; i < mz; ++i) {
		if (lmin[i] == 1) {
			zlmin[nlmin] = zj[i];
			nlmin += 1;
		}
	}
	//mdisplay(zlmin, 1, mz);

	if (nlmin >= 2) {
		for (j = 1; j < nlmin;++j) {
			zlo = zlmin[j - 1];
			zhi = zlmin[j];
			ze = (zlo + zhi) / 2;
			//printf("ze %g \n", ze);
			zstep = (zhi - zlo) / 2;
			for (nit = 0; nit < 10; ++nit) {
				cbi = (double*)calloc(nx, sizeof(double));
				for (i = 0; i < nx; ++i) {
					if ((ze >= zs1[i]) && (ze <= zs2[i])) {
						cbi[i] = cb[i];
					}
				}
				likd = 0.0;
				for (i = 0; i < nx; ++i) {
					likd += (cbi[i] / (1.0 + ze * cbi[i]));
				}
				zstep /= 2.0;
				sgn = likd > 0.0 ? 1 : -1;
				if (likd == 0.0) {
					sgn = 0;
				}
				ze += (sgn*zstep);
				//printf("ze %g \n", ze);
			}
			zlmax[nlmax] = ze;
			nlmax += 1;
		}
	}
	//mdisplay(zlmax, 1, mz);
	// Evaluate all local maxima and find global max

	zm = (double*)malloc(sizeof(double) * nlmax);

	for (j = 0; j < nlmax; ++j) {
		for (i = 0; i < nx; ++i) {
			pz[i] = zlmax[j] < zs2[i] ? zlmax[j] : zs2[i];
		}
		for (i = 0; i < nx; ++i) {
			pz[i] = pz[i] > zs1[i] ? pz[i] : zs1[i];
		}
		tmp = 0.0;

		for (i = 0; i < nx; ++i) {
			tmp += log(1.0 + cb[i] * pz[i]);
		}
		zm[j] = tmp;
	}
	
	index = array_max_index(zm,nlmax);
	zmax = zlmax[index];

	free(beta);
	free(zs1);
	free(zs2);
	free(zj);
	free(cb);
	free(lmin);
	free(cbir);
	free(cbil);
	free(zlmax);
	free(zlmin);
	free(zm);
	free(pz);
	return zmax;
}

double macheps() {
	double macheps;
	macheps = 1.0;

	while ((macheps + 1.0) > 1.0) {
		macheps = macheps / 2.0;
	}

	macheps = macheps * 2;

	return macheps;
}

void mdisplay(double *A, int row, int col) {
	int i, j;
	printf("\n MATRIX Order : %d X %d \n \n", row, col);

	for (i = 0; i < row; i++) {
		printf("R%d: ", i);
		for (j = 0; j < col; j++) {
			printf("%g ", A[i*col + j]);
		}
		printf(":R%d \n", i);
	}
}

int mvalue(int N) {
	int mval;

	if (N <= 10) {
		mval = N;
	}
	else if (N > 10 && N <= 20) {
		mval = 10;
	}
	else if (N > 20 && N <= 200) {
		mval = 15;
	}
	else if (N > 200) {
		mval = 20;
	}

	return mval;
}

int grad_cd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
	double eps3, double *f) {
	int retval;
	retval = 0;
	if (funcgrad == NULL) {
		//printf("FD Gradient \n");
		retval = grad_calc2(funcpt, x, N, dx, eps3, f);
	}
	else {
		//printf("Analytic gradient \n");
		FUNCGRAD_EVAL(funcgrad, x, N, f);
	}
	return retval;

}

double signx(double x) {
	double sgn;
	if (x > 0.) {
		sgn = 1.0;
	}
	else if (x < 0.) {
		sgn = -1.0;
	}
	else {
		sgn = 0.0;
	}

	return sgn;
}

int grad_calc2(custom_function *funcpt, double *x, int N, double *dx, double eps3, double *f) {
	int j, retval;
	double stepsize, stepmax, temp;
	double fp, fm;

	retval = 0;

	for (j = 0; j < N; ++j) {
		if (fabs(x[j]) >= 1.0 / fabs(dx[j])) {
			stepmax = x[j];
		}
		else {
			stepmax = signx(x[j]) * 1.0 / fabs(dx[j]);
		}

		stepsize = stepmax * eps3;
		temp = x[j];
		x[j] += stepsize;
		stepsize = x[j] - temp;
		fp = FUNCPT_EVAL(funcpt, x, N);
		if (fp >= DBL_MAX || fp <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			return 15;
		}
		if (fp != fp) {
			printf("Program Exiting as the function returns NaN");
			return 15;
		}
		x[j] = temp - stepsize;
		fm = FUNCPT_EVAL(funcpt, x, N);
		if (fm >= DBL_MAX || fm <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			return 15;
		}
		if (fm != fm) {
			printf("Program Exiting as the function returns NaN");
			return 15;
		}
		f[j] = (fp - fm) / (2 * stepsize);
		x[j] = temp;
	}

	return retval;
}


int grad_calc(custom_function *funcpt, double *x, int N, double *dx, double eps2, double *f) {
	int i, j, retval;
	double step, fd, stepmax;
	double *xi;

	fd = eps2; // square root of macheps
	retval = 0;
	xi = (double*)malloc(sizeof(double)*N);

	for (i = 0; i < N; ++i) {
		if (fabs(x[i]) >= 1.0 / fabs(dx[i])) {
			stepmax = x[i];
		}
		else {
			stepmax = signx(x[i]) * 1.0 / fabs(dx[i]);
		}
		step = fd * stepmax;
		for (j = 0; j < N; ++j) {
			xi[j] = x[j];
		}
		xi[i] += step;
		f[i] = (FUNCPT_EVAL(funcpt, xi, N) - FUNCPT_EVAL(funcpt, x, N)) / step;
		if (f[i] >= DBL_MAX || f[i] <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(xi);
			return 15;
		}
		if (f[i] != f[i]) {
			printf("Program Exiting as the function returns NaN");
			free(xi);
			return 15;
		}
		//xi[i] -= step;
	}

	free(xi);
	return retval;

}

int grad_fd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
	double eps2, double *f) {
	int retval;
	retval = 0;
	if (funcgrad == NULL) {
		//printf("FD Gradient \n");
		retval = grad_calc(funcpt, x, N, dx, eps2, f);
	}
	else {
		//printf("Analytic gradient \n");
		FUNCGRAD_EVAL(funcgrad, x, N, f);
	}

	return retval;

}

int grad_rextp(custom_function *funcpt, double *x, int N, int r, int v, double *a, double *h, double *xp, double *xm, double meps, double *f) {
	int i, k, j, m, retval, m1, i1;
	double eps, ztol, d, fx, fn, pm1;
	retval = 0;

	eps = 1.0e-04;
	d = 0.0001;
	ztol = sqrt(meps / (7.0e-07));

	for (i = 0; i < N; ++i) {
		fx = fabs(x[i]);
		h[i] = fx*d;

		if (fx < ztol) {
			h[i] += eps;
		}
	}

	//mdisplay(h, 1, N);

	for (k = 0; k < r; ++k) {
		if (N == 1) {
			for (i = 0; i < N; ++i) {
				xp[i] = x[i] + h[i];
				xm[i] = x[i] - h[i];
			}
			fn = FUNCPT_EVAL(funcpt, xp, N) - FUNCPT_EVAL(funcpt, xm, N);
			for (i = 0; i < N; ++i) {
				a[k*N + i] = fn / (2 * h[i]);
				if (a[k*N + i] != a[k*N + i]) {
					printf("Program Exiting as the function returns NaN");
					return 15;
				}
			}
		}
		else {
			if ((k != 0) && (fabs(a[(k - 1)*N + i]) < 1e-20)) {
				a[k*N + i] = 0;
			}
			else {
				for (i = 0; i < N; ++i) {
					for (j = 0; j < N; ++j) {
						xp[j] = xm[j] = x[j];
						if (j == i) {
							xp[j] += h[j];
							xm[j] -= h[j];
						}
					}
					fn = FUNCPT_EVAL(funcpt, xp, N) - FUNCPT_EVAL(funcpt, xm, N);
					a[k*N + i] = fn / (2 * h[i]);
					if (a[k*N + i] != a[k*N + i]) {
						printf("Program Exiting as the function returns NaN");
						return 15;
					}
				}
			}
		}
		for (i = 0; i < N; ++i) {
			h[i] /= (double)v;
		}
	}

	//mdisplay(a, r, N);

	for (m = 0; m < r - 1; ++m) {
		m1 = m + 1;
		pm1 = pow(4.0, (double)m1);
		for (i = 1; i < r - m; ++i) {
			i1 = i - 1;
			for (k = 0; k < N; ++k) {
				a[i1*N + k] = (a[i*N + k] * pm1 - a[i1*N + k]) / (pm1 - 1.0);
			}
		}

	}

	for (i = 0; i < N; ++i) {
		f[i] = a[i];
	}


	return retval;
}

/* Using minfun2d instead

int lbfgsb_min(custom_function *funcpt, custom_gradient *funcgrad, double *x, int NN, int MM, double *l, double *u) {
	int i, retval;
	double eps, f;
	double gtol;
	integer *nbd, *iwa;
	//     static char task[60]; 
	static integer taskValue;
	static integer *task = &taskValue; // must initialize !! 
	//      http://stackoverflow.com/a/11278093/269192 
	static double factr;
	//     static char csave[60]; 
	static integer csaveValue;
	static integer *csave = &csaveValue;
	static double dsave[29];
	static integer isave[44];
	static logical lsave[4];
	static double pgtol;
	static integer iprint;
	double *g, *dx, *wa;
	int maxiter, nmax, mmax;
	int r, v;
	double *a, *h, *xp, *xm;
	integer N, M;

	N = (integer)NN;
	M = (integer)MM;
	r = 4;
	v = 2;

	retval = 0;

	nmax = 1024;
	mmax = 20;

	if (N > nmax) {
		nmax = N;
	}

	if (M <= 0) {
		M = mvalue(N);
	}
	else if ((M > mmax) | (M > N)) {
		printf("M should be less than or equal to 20 or N whichever is smaller. \n");
	}


	//(2mmax + 5)nmax + 12mmax^2 + 12mmax. 

	eps = macheps();
	gtol = 1.0e-5;// pow(eps, 1.0 / 3.0); //pgtol
	factr = 1.0e7;
	pgtol = 1.0e-5;

	// Intialize bounds parameter nbd

	nbd = (integer*)malloc(sizeof(integer) * N);
	g = (double*)malloc(sizeof(double)*N);
	dx = (double*)malloc(sizeof(double)*N);
	h = (double*)malloc(sizeof(double)*N);
	xp = (double*)malloc(sizeof(double)*N);
	xm = (double*)malloc(sizeof(double)*N);
	a = (double*)malloc(sizeof(double)*r*N);

	wa = (double*)calloc((2 * mmax + 5)*nmax + 12 * mmax*mmax + 12 * mmax, sizeof(double));
	iwa = (integer*)calloc(3 * nmax, sizeof(integer));


	for (i = 0; i < N; ++i) {
		dx[i] = 1.0;
	}

	if ((l == NULL) && (u == NULL)) {
		nbd[0] = 0;
	}
	else if (u == NULL) {
		nbd[0] = 1;
	}
	else if (l == NULL) {
		nbd[0] = 3;
	}
	else {
		nbd[0] = 2;
	}

	for (i = 1; i < N; ++i) {
		nbd[i] = nbd[0];
	}

	maxiter = 1000;

	iprint = -1;// Set iprint to 1 for verbose output 

	//     We start the iteration by initializing task. 

	*task = (integer)START;
	//     s_copy(task, "START", (ftnlen)60, (ftnlen)5); 
	//        ------- the beginning of the loop ---------- 
L111:
	//     This is the call to the L-BFGS-B code. 
	setulb(&N, &M, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
		iprint, csave, lsave, isave, dsave);
	//     if (s_cmp(task, "FG", (ftnlen)2, (ftnlen)2) == 0) { 
	if (IS_FG(*task)) {
		f = FUNCPT_EVAL(funcpt, x, NN);
		grad_rextp(funcpt, x, NN, r, v, a, h, xp, xm, eps, g);
		//         go back to the minimization routine. 
		goto L111;
	}

	//     if (s_cmp(task, "NEW_X", (ftnlen)5, (ftnlen)5) == 0) { 
	if (*task == NEW_X) {
		goto L111;
	}

	free(nbd);
	free(g);
	free(dx);
	free(wa);
	free(iwa);
	free(xp);
	free(xm);
	free(h);
	free(a);

	return retval;
}

*/

static void minfun2d1(custom_function *funcpt, double *x0, int N, int I, double *xl, double *xu,double *a, double *b,double *v, double *u, double *xmin, double *fmin) {
	/*
	This is a C Translation of the code by Anestis Antoniadis and Maarten Jansen. The matlab version is available
	at http://www-ljk.imag.fr/membres/Anestis.Antoniadis/EBayesThresh/down.html
	*************************************************************************************************************
	COPYRIGHT NOTICE
	% Copyright Anestis Antoniadis and Maarten Jansen
	% antonia@imag.fr and mjansen@win.tue.nl
	% Part of this procedure is a modified copy of the function mingcv.m from
	% PiefLab (see www.cs.kuleuven.ac.be/~maarten/software/pieflab.html)
	*/
	int i, k, F0, F1, F2, Niter;
	double eps,f12,f02,ba,fv,fu;
	F0 = F1 = Niter = 1;
	F2 = 2;

	eps = 1e-07;

	while ((xu[I] - xl[I]) / (double) F2 > eps) {
		F0 = F1;
		F1 = F2;
		F2 = F0 + F1;
		Niter += 1;
	}

	for (i = 0; i < N; ++i) {
		if (i == I) {
			a[i] = xl[i];
			b[i] = xu[i];
		}
		else {
			a[i] = b[i] = x0[i];
		}
	}

	f12 = (double)F1 / (double)F2;
	f02 = (double)F0 / (double)F2;

	for (i = 0; i < N; ++i) {
		ba = b[i] - a[i];
		v[i] = a[i] + f12 * ba;
		u[i] = a[i] + f02 * ba;
	}

	fv = FUNCPT_EVAL(funcpt, v, N);

	fu = FUNCPT_EVAL(funcpt, u, N);

	for (k = 0; k < Niter-1; ++k) {
		if (fu > fv) {
			for (i = 0; i < N; ++i) {
				a[i] = u[i];
				u[i] = v[i];
			}
			fu = fv;
			for (i = 0; i < N; ++i) {
				v[i] = b[i] - u[i] + a[i];
			}
			fv = FUNCPT_EVAL(funcpt, v, N);
		}
		else {
			for (i = 0; i < N; ++i) {
				b[i] = v[i];
				v[i] = u[i];
			}
			fv = fu;
			for (i = 0; i < N; ++i) {
				u[i] = a[i] + b[i] - v[i];
			}
			fu = FUNCPT_EVAL(funcpt, u, N);
		}
	}

	for (i = 0; i < N; ++i) {
		xmin[i] = (a[i] + b[i]) / 2;
	}

	*fmin = FUNCPT_EVAL(funcpt, xmin, N);

}

double l2norm(double *vec, int N) {
	double l2, sum;
	int i;
	sum = 0.;
	for (i = 0; i < N; ++i) {
		sum += vec[i] * vec[i];
	}
	l2 = sqrt(sum);
	return l2;
}

void minfun2d(custom_function *funcpt, custom_gradient *funcgrad, double *x, int NN, int MM, double *xl, double *xu) {
	/*
	This is a C Translation of the code by Anestis Antoniadis and Maarten Jansen. The matlab version is available
	at http://www-ljk.imag.fr/membres/Anestis.Antoniadis/EBayesThresh/down.html
	*************************************************************************************************************
	COPYRIGHT NOTICE
	% Copyright Anestis Antoniadis and Maarten Jansen
	% antonia@imag.fr and mjansen@win.tue.nl
	% Part of this procedure is a modified copy of the function mingcv.m from
	% PiefLab (see www.cs.kuleuven.ac.be/~maarten/software/pieflab.html)
	*/
	int i, I;
	double eps,relerr,fmin;
	double *a, *b, *u, *v,*xmin,*x2;

	a = (double*)malloc(sizeof(double)*NN);
	b = (double*)malloc(sizeof(double)*NN);
	u = (double*)malloc(sizeof(double)*NN);
	v = (double*)malloc(sizeof(double)*NN);
	xmin = (double*)malloc(sizeof(double)*NN);
	x2 = (double*)malloc(sizeof(double)*NN);

	I = 0;
	eps = 1e-07;
	relerr = 1.0;
	fmin = 0;

	while (relerr > eps) {
		minfun2d1(funcpt, x, NN, I, xl, xu, a, b, v, u, xmin, &fmin);
		for (i = 0; i < NN; ++i) {
			x2[i] = xmin[i] - x[i];
		}
		relerr = l2norm(x2,NN) / (l2norm(xmin, NN) - l2norm(x, NN));
		I = (I + 1) % NN;
		for (i = 0; i < NN; ++i) {
			x[i] = xmin[i];
		}
	}

	free(a);
	free(b);
	free(u);
	free(v);
	free(xmin);
	free(x2);
}
