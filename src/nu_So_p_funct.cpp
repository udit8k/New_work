#include"main.h"


// POP scattering calculation
void nu_So_p_funct()
{	double k[],k_pp,k_mm;
 	
	double C_mm[i]= (k[i]*k[i]+k_mm*k_mm)/(2*k_mm*k[i]);
 	double C_pp[i]= (k[i]*k[i]+k_pp*k_pp)/(2*k_pp*k[i]);
	double B_mm[i]= ((1+3*C_mm[i]*C_mm[i])/2)*log(abs((1+C_mm[i])/(1-C_mm[i]))) - 3*C_mm[i];
	double B_pp[i]= ((1+3*C_pp[i]*C_pp[i])/2)*log(abs((1+C_pp[i])/(1-C_pp[i]))) - 3*C_pp[i];
	double pre_mat[i]= B_mm[i]*(N_poph(double omega, double T)+1)*(1-f0(double E1, double e_f, double T)) + B_mm[i]*(N_poph(double omega, double T)+0)*(f0(double E1, double e_f, double T))
		+B_pp[i]*(N_poph(double omega, double T)+1)*(1-f0(double E1, double e_f, double T)) + B_pp[i]*(N_poph(double omega, double T)+1)*(f0(double E1, double e_f, double T));
 	
}
