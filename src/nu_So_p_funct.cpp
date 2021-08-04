#include"main.h"


// POP scattering calculation
void nu_So_p_funct(double T, double efef, int points)
{	double k[points],counter;
 	double k_pp;
 	double k_mm;
 	
	double C_mm[points];
 	double C_pp[points];
	double B_mm[points];
 	double B_pp[points];
	double pre_mat[points];
	
 	 double arr[points];
 
 	 double k_minus = kminus(int counter, double omega, int points, double energy[]);

	    for (int i=0;i<points;i++)
		arr[i] = abs(k_grid[i] - k_minus);
	    int minus_index =FindMinInd(arr,points);


	double k_plus = kplus(int counter, double omega, int points, double energy[]);
	    //cout<<"k_plus = "<<k_plus<<endl;

	    for (int i=0;i<points;i++)
		arr[i] = abs(k_grid[i] - k_plus);
	int plus_index =FindMinInd(arr,points);
 
	double k_pp= kplus(int counter, double omega, int points, double energy[]);
 	double k_mm= kminus(int counter, double omega, int points, double energy[]);
 	
 for (int i=0;i<points;i++)
	{
		S_o_grid[i] = 1;
		S_o_grid_total[i]=0;
	}	

 for(i=0;i<points;i++){
 		
	C_mm[i]= (k[i]*k[i]+k_mm*k_mm)/(2*k_mm*k[i]);
 	C_pp[i]= (k[i]*k[i]+k_pp*k_pp)/(2*k_pp*k[i]);
	B_mm[i]= ((1+3*C_mm[i]*C_mm[i])/2)*log(abs((1+C_mm[i])/(1-C_mm[i]))) - 3*C_mm[i];
	B_pp[i]= ((1+3*C_pp[i]*C_pp[i])/2)*log(abs((1+C_pp[i])/(1-C_pp[i]))) - 3*C_pp[i];
	pre_mat[i]= B_mm[i]*(N_poph(double omega, double T)+1)*(1-f0(energy_n[minus_index],efefp,T)) + B_mm[i]*(N_poph(double omega, double T)+0)*(f0(energy_n[minus_index],efefp,T))
			+B_pp[i]*(N_poph(double omega, double T)+1)*(1-f0f0(energy_n[plus_index],efefp,T)) + B_pp[i]*(N_poph(double omega, double T)+1)*(f0f0(energy_n[plus_index],efefp,T));
		

	if (scattering_mechanisms[1] == 1)
{

	if ((((1+3*C_mm[i]*C_mm[i])/2)*log(abs((1+C_mm[i])/(1-C_mm[i]))) - 3*C_mm[i]==0)&&(((1+3*C_pp[i]*C_pp[i])/2)*log(abs((1+C_pp[i])/(1-C_pp[i]))) - 3*C_pp[i]) ==0)
	S_o_grid[counter1] = average_dummy;
			    // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
	else
	S_o_grid[counter1] = (1/(epsilon_inf*epsilon_0)-1/(epsilon_s*epsilon_0))*pre_mat[i]*(e*e*omega)/(16*pi*h_bar*v_p[i]);

}
 }	
	    
	 
	for (int counter1=0;counter1<points;counter1++)
        	S_o_grid_total[counter1] = S_o_grid[counter1];

	

}
