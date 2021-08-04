#include"main.h"


// POP scattering calculation
void nu_So_p_funct(double T, double efef, int ii)
{	double k[ii],counter;
 	double k_pp;
 	double k_mm;
 	
	double C_mm[ii];
 	double C_pp[ii];
	double B_mm[ii];
 	double B_pp[ii];
	double pre_mat[ii];
 		
 for(i=0;i<){
 	k_pp= kplus(int counter, double omega, int points, double energy[]);
 	k_mm= kminus(int counter, double omega, int points, double energy[]);
 	
	C_mm[i]= (k[i]*k[i]+k_mm*k_mm)/(2*k_mm*k[i]);
 	C_pp[i]= (k[i]*k[i]+k_pp*k_pp)/(2*k_pp*k[i]);
	B_mm[i]= ((1+3*C_mm[i]*C_mm[i])/2)*log(abs((1+C_mm[i])/(1-C_mm[i]))) - 3*C_mm[i];
	B_pp[i]= ((1+3*C_pp[i]*C_pp[i])/2)*log(abs((1+C_pp[i])/(1-C_pp[i]))) - 3*C_pp[i];
	pre_mat[i]= B_mm[i]*(N_poph(double omega, double T)+1)*(1-f0(double E1, double e_f, double T)) + B_mm[i]*(N_poph(double omega, double T)+0)*(f0(double E1, double e_f, double T))
		+B_pp[i]*(N_poph(double omega, double T)+1)*(1-f0(double E1, double e_f, double T)) + B_pp[i]*(N_poph(double omega, double T)+1)*(f0(double E1, double e_f, double T));
 }	
}
for (int aa=0;aa<iv_number;aa++)
{
	N_e[aa] = N_poph(we[aa],T);
}
