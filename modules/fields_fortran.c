#include "../include/fields_fortran.h"
#include <math.h>

double B_h = 0.005; //holding field strength
double B = 1.4; //remnant magnet field strength
double d = 0.0254; //thickness of layer of PM array
double L = 0.02; //characteristic spacing of magnets
//double mu = -9.662364e-27/1.674927351e-27; //mu in units where m=1
double mu = -9.6623647e-27;
double m_n = 1.674927471e-27;
double grav = 9.80665e0;
double N = 3.0; //how far out to go in field ripple expansion

void force_(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t) //analytical form of halbach field force, mu*del(mod(B))
{
	double A = sqrt(8.0)*B/M_PI; //parameter related to B -- shows up in expansion

	double n,m;

	double x = *x_in;
	double y = *y_in;
	double z = *z_in;

	double gx=0.0, gy=0.0, gz=0.0, R, r, B_tot;

	if (x > 0.0)
	{
		R = 1.0;
		r = 0.5;
	}
	else
	{
		R = 0.5;
		r = 1.0;
	}

	double rho = sqrt(y*y+z*z);
	double r_zeta = sqrt((rho-R)*(rho-R)+x*x);

	if (z < -1.0 && r_zeta < r)
	{

		double rho = sqrt(y*y+z*z);
		double r_zeta = sqrt((rho-R)*(rho-R)+x*x);
		double zeta = r-r_zeta;
		double eta = r*atan(x/(rho-R));
//		double Bsum = 0.0, B_zeta = 0.0 , B_eta = 0.0;
		double B_zeta = 0.0 , B_eta = 0.0;
		
		double k_n;
		
		for(n = 1.0; n <= N; n+=1.0) {
			k_n = 2*M_PI*(4.0*n-3.0)/L;
			B_zeta += pow(-1.0,n)/(4.0*n-3.0)*(1-exp(-k_n*d))*exp(-(k_n)*zeta)*cos((k_n)*eta);
			B_eta += pow(-1.0,n)/(4.0*n-3.0)*(1-exp(-k_n*d))*exp(-(k_n)*zeta)*sin((k_n)*eta);
		}
		
		double Hold = B_h*B_h*(r+R)*(r+R)/rho/rho/rho/rho;
		B_tot = sqrt(B_h*B_h*(r+R)*(r+R)/rho/rho + A*A*B_zeta*B_zeta + A*A*B_eta*B_eta);

		gx = 0.5/B_tot*mu*(A*A*(B_zeta*B_zeta)*x/r_zeta - A*A*(B_eta*B_eta)*r*(rho-R)/r_zeta/r_zeta);
		gy = 0.5/B_tot*mu*(A*A*(B_zeta*B_zeta)*y*(1.0-R/rho)/r_zeta + A*A*(B_eta*B_eta)*x*y*r/rho/r_zeta/r_zeta - 2.0*Hold*y);
		gz = 0.5/B_tot*mu*(A*A*(B_zeta*B_zeta)*z*(1.0-R/rho)/r_zeta + A*A*(B_eta*B_eta)*r*x*z/rho/r_zeta/r_zeta - 2.0*Hold*z);

//		double k_m,k_n;
//
//		for (m = 1.0;m<=N;m+=1.0)
//		{
//			k_m = 2*M_PI*(4.0*m-3.0)/L;
//
//			for (n = 1.0;n<=N;n+=1.0)
//			{
//				k_n = 2*M_PI*(4.0*n-3.0)/L;
//
//				B_zeta += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*d))*(1-exp(-k_n*d))*(k_n+k_m)*exp(-(k_n+k_m)*zeta)*cos((k_n-k_m)*eta);
//				B_eta += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*d))*(1-exp(-k_n*d))*(k_n-k_m)*exp(-(k_n+k_m)*zeta)*sin((k_n-k_m)*eta);
//
//				Bsum += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*d))*(1-exp(-k_n*d))*exp(-(k_n+k_m)*zeta)*cos((k_n-k_m)*eta);
//			}
//		}
//
//		double Hold = B_h*B_h*(r+R)*(r+R)/rho/rho/rho/rho;
//		B_tot = sqrt(B_h*B_h*(r+R)*(r+R)/rho/rho + A*A*Bsum);
//
//		gx = 0.5/B_tot*mu*(A*A*B_zeta*x/r_zeta - A*A*B_eta*r*(rho-R)/r_zeta/r_zeta);
//		gy = 0.5/B_tot*mu*(A*A*B_zeta*y*(1.0-R/rho)/r_zeta + A*A*B_eta*x*y*r/rho/r_zeta/r_zeta - 2.0*Hold*y);
//		gz = 0.5/B_tot*mu*(A*A*B_zeta*z*(1.0-R/rho)/r_zeta + B_eta*r*x*z/rho/r_zeta/r_zeta - 2.0*Hold*z);

	}

	gz -= grav*m_n;

	*totalU = -mu*B_tot+grav*m_n*z;

	*fx = gx;
	*fy = gy;
	*fz = gz;
}

void fieldstrength_(double *x_in, double *y_in, double *z_in, double *totalB, double* t) //-mu*mod(B)
{
	double A = sqrt(8.0)*B/M_PI; //parameter related to B -- shows up in expansion

	double x = *x_in;
	double y = *y_in;
	double z = *z_in;

	double R,r;

	if (x > 0.0)
	{
		R = 1.0;
		r = 0.5;
	}
	else
	{
		R = 0.5;
		r = 1.0;
	}

	double rho = sqrt(y*y+z*z);
	double r_zeta = sqrt((rho-R)*(rho-R)+x*x);
	double B_tot = 0.0;

	if (z < -1.0 && r_zeta < r)
	{
		double zeta = r-r_zeta;
		double eta = r*atan(x/(rho-R));
		double Bsum = 0.0;
		double k_m,k_n,m,n;

		for (m = 1.0;m<=N;m+=1.0)
		{
			k_m = 2*M_PI*(4.0*m-3.0)/L;
			for (n = 1.0;n<=N;n+=1.0)
			{
				k_n = 2*M_PI*(4.0*n-3.0)/L;

				Bsum += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*d))*(1-exp(-k_n*d))*exp(-(k_n+k_m)*zeta)*cos((k_n-k_m)*eta);
			}
		}

		B_tot = sqrt(B_h*B_h*(r+R)*(r+R)/(y*y+z*z)+A*A*Bsum);

	}

	*totalB = B_tot;
}

void potential_(double *x_in, double *y_in, double *z_in, double *totalU, double* t) //-mu*mod(B) + g*z. remember that mu is already negative.
{
	double A = sqrt(8.0)*B/M_PI; //parameter related to B -- shows up in expansion (see Walsrom, et al).

	double x = *x_in;
	double y = *y_in;
	double z = *z_in;

	double R,r;

	if (x > 0.0)
	{
		R = 1.0;
		r = 0.5;
	}
	else
	{
		R = 0.5;
		r = 1.0;
	}

	double rho = sqrt(y*y+z*z);
	double r_zeta = sqrt((rho-R)*(rho-R)+x*x);
	double B_tot = 0.0;

	if (z < -1.0 && r_zeta < r)
	{
		double zeta = r-r_zeta;
		double eta = r*atan(x/(rho-R));
		double Bsum = 0.0;
		double k_m,k_n,m,n;

		for (m = 1.0;m<=N;m+=1.0)
		{
			k_m = 2*M_PI*(4.0*m-3.0)/L;
			for (n = 1.0;n<=N;n+=1.0)
			{
				k_n = 2*M_PI*(4.0*n-3.0)/L;

				Bsum += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*d))*(1-exp(-k_n*d))*exp(-(k_n+k_m)*zeta)*cos((k_n-k_m)*eta);
			}
		}

		B_tot = sqrt(B_h*B_h*(r+R)*(r+R)/(y*y+z*z)+A*A*Bsum);

	}

	*totalU = -mu*B_tot+grav*m_n*z;
}
