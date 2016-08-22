/**
 * Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
 * Signal Analysis and Machine Perception Laboratory,
 * Department of Electrical, Computer, and Systems Engineering,
 * Rensselaer Polytechnic Institute, Troy, NY 12180, USA
 */

/**
 * This is the C/MEX code of dynamic time warping of two signals
 *
 * compile:
 *     mex dtw_c.c
 *
 * usage:
 *     d=dtw_c(s,t)  or  d=dtw_c(s,t,w)
 *     where s is signal 1, t is signal 2, w is window parameter
 */

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double cubicRoot(double d)
{
	  if(d<0.0)
		        return -cubicRoot(-d);
	    else
			      return pow(d,1.0/3.0);
}

/* This function solves the following problem:
 min_{x>=0} x^3+ax+b */
double root_c(double a, double b)
{
    double x=0, y=0;
    double a3=4*pow(a,3), b2=27*pow(b,2);
    double delta = a3+b2;
	int k;
    if(delta<=0) /* 3 distinct real roots or 1 real multiple solution */
    {
	    double r3  = 2*sqrt(-a/3);
        double th3 = atan2(sqrt(-delta/108),-b/2)/3;
        double ymax=0, xopt=0;
        for(k=0;k<=4;k=k+2)
        {
            x = r3*cos(th3+((k*3.14159265)/3));
            y=pow(x,4)/4+a*pow(x,2)/2+b*x;
	 	    if(y<ymax)
	               {ymax=y; xopt=x;}
        }
        return xopt;
    }
    else /* 1 real root and two complex */
    {
         double z = sqrt(delta/27);
         x = cubicRoot(0.5*(-b+z))+cubicRoot(0.5*(-b-z));
         y = pow(x,4)/4+a*pow(x,2)/2+b*x;
         return x;
    }
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
	double * dp;
    /*  create a pointer to the input matrix s */
    double a = mxGetScalar(prhs[0]);

    /*  create a pointer to the input matrix t */
    double b = mxGetScalar(prhs[1]);

    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);

    /*  create a C pointer to a copy of the output matrix */
    dp = mxGetPr(plhs[0]);

    /*  call the C subroutine */
    dp[0]=root_c(a,b);

    return;

}
