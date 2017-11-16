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
int N_terms = 3; //how far out to go in field ripple expansion

void force_(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t) //analytical form of halbach field force, mu*del(mod(B))
{
//	double A = sqrt(8.0)*B/M_PI; //parameter related to B -- shows up in expansion
	double A = 4*B/(M_PI*sqrt(2));

//	double n;

	double x = *x_in;
	double y = *y_in;
	double z = *z_in;

	double gx=0.0, gy=0.0, gz=0.0, R, r;

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
//		double eta = ;
//		double zeta = ;
		double eta = r*atan(x/(sqrt(y*y + z*z) - R));
		double zeta = r - sqrt(x*x + (sqrt(y*y + z*z) - R)*(sqrt(y*y + z*z) - R));
		double sum_cos=0.0, sum_sin=0.0, sum_k_cos=0.0, sum_k_sin=0.0;
		double cos_term=0.0, sin_term=0.0;
		
		double k_n;

		for (int n = 1.0; n <= N_terms; n += 1.0)
		{
			k_n = 2*M_PI*(4.0*n-3.0)/L;
			
			cos_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*d))*exp(-k_n*zeta)*cos(k_n*eta);
			sin_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*d))*exp(-k_n*zeta)*sin(k_n*eta);
			
			sum_cos += cos_term;
			sum_k_cos += k_n*cos_term;
			sum_sin += sin_term;
			sum_k_sin += k_n*sin_term;
		}
		
		double b_zeta = A*sum_cos;
		double b_eta = A*sum_sin;
		double b_hold = B_h*(r+R)
			                /
			         (sqrt(y*y + z*z));
		
		double b_tot = sqrt(b_zeta*b_zeta + b_eta*b_eta + b_hold*b_hold);
	
		//zeta, eta
		double d_BZeta[2] = {-1*A*sum_k_cos, -1*A*sum_k_sin};
		double d_BEta[2] = {-1*A*sum_k_sin, A*sum_k_cos};

		//x,y,z
		double d_Bh[3] = {0.0,
						  -B_h*(r+R)*y/(pow(y*y + z*z, 3.0/2.0)),
						  -B_h*(r+R)*z/(pow(y*y + z*z, 3.0/2.0))};
		double d_Zeta[3] = {-(x
							  /
							  sqrt(
								  x*x + ((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z)))
							  )
							 ),
							-y*(sqrt(y*y + z*z) - R)
								/(
								  sqrt(y*y + z*z)
								  *sqrt(x*x + ((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
							     ),
							-z*(sqrt(y*y + z*z) - R)
								/(
								  sqrt(y*y + z*z)
								  *sqrt(x*x + ((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
							     )};
		double d_Eta[3] = {r
			               /(
			                 (sqrt(y*y + z*z) - R)
			                 *(1.0 + x*x/((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
		                    ),
						   -r*x*y
							/(
							   sqrt(y*y + z*z)
							   *((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z)))
							   *(1.0 + x*x/((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
						     ),
						   -r*x*z
							/(
							   sqrt(y*y + z*z)
							   *((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z)))
							   *(1.0 + x*x/((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
						     )};
		
		gx = mu*(1.0/b_tot)*(
			b_zeta*(d_BZeta[0]*d_Zeta[0] + d_BZeta[1]*d_Eta[0])
			+ b_eta*(d_BEta[0]*d_Zeta[0] + d_BEta[1]*d_Eta[0]));
		
		gy = mu*(1.0/b_tot)*(
			b_zeta*(d_BZeta[0]*d_Zeta[1] + d_BZeta[1]*d_Eta[1])
			+ b_eta*(d_BEta[0]*d_Zeta[1] + d_BEta[1]*d_Eta[1])
			+ b_hold*d_Bh[1]);
		
		gz = mu*(1.0/b_tot)*(
			b_zeta*(d_BZeta[0]*d_Zeta[2] + d_BZeta[1]*d_Eta[2])
			+ b_eta*(d_BEta[0]*d_Zeta[2] + d_BEta[1]*d_Eta[2])
			+ b_hold*d_Bh[2]);
		
		*totalU = -mu*b_tot + grav*m_n*z;
	}

	gz -= grav*m_n;

//	*totalU = -mu*b_tot + grav*m_n*z;

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

		for (m = 1.0;m<=3.0;m+=1.0)
		{
			k_m = 2*M_PI*(4.0*m-3.0)/L;
			for (n = 1.0;n<=3.0;n+=1.0)
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

		for (m = 1.0;m<=3.0;m+=1.0)
		{
			k_m = 2*M_PI*(4.0*m-3.0)/L;
			for (n = 1.0;n<=3.0;n+=1.0)
			{
				k_n = 2*M_PI*(4.0*n-3.0)/L;

				Bsum += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*d))*(1-exp(-k_n*d))*exp(-(k_n+k_m)*zeta)*cos((k_n-k_m)*eta);
			}
		}

		B_tot = sqrt(B_h*B_h*(r+R)*(r+R)/(y*y+z*z)+A*A*Bsum);

	}

	*totalU = -mu*B_tot+grav*m_n*z;
}
