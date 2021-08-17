#include"main.h"

// POP scattering calculation
void nu_So_p_funct( double T) {
cout<<"working"<<endl;
double k[points];
double k_pp;
double k_mm;
int i,counter;
double C_mm[points];
double C_pp[points];
double B_mm[points];
double B_pp[points];
double pre_mat[points],p_kminus_grid_pop[points],nu_So_p[points][2][2];
//double omega  =omega_LO;
double arr[points],p_S_o_grid[points],p_S_o_grid_total[points],p_k_grid[points],p_k[points],p_kplus_grid_pop[points];
int minus_index;
int plus_index;
double k_minus = kminus(counter,  omega_LO,  points,  energy_p);
double k_plus = kplus(counter,  omega_LO,  points,  energy_p);
double sum=0;
double average_dummy = sum/points;
double k_dum; 
	
	k_pp= k_plus;     //redudant term  
	k_mm= k_minus;
 	
 		for(i=0;i<points;i++){
			p_S_o_grid[i] = 1;  			 
			p_S_o_grid_total[i]=0;			
					}	

		for(counter = 0;counter < points;counter++)
	    		sum = sum +  p_S_o_grid[counter];		

	

		for(counter = 0;counter < points;counter++){	
			k_dum =p_k_grid[counter];			//3
		
		//grid_
				p_kplus_grid_pop[counter] = kplus(counter,omega_LO,points,energy_p);   	//4
				p_kminus_grid_pop[counter] = kminus(counter,omega_LO,points,energy_p);		//5

		
			for ( i=0;i<points;i++)
				arr[i] = abs(p_k_grid[i] - p_kminus_grid_pop[counter]);			//4
		
			minus_index =FindMinInd(arr,points);

			for ( i=0;i<points;i++)
				arr[i] = abs(p_k_grid[i] - p_kplus_grid_pop[counter]);			//5
		
		 	plus_index =FindMinInd(arr,points);
//for(i=0;i<points;i++){		
			
			C_mm[counter]= (p_k[counter]*p_k[counter]+k_mm*k_mm)/(2*k_mm*p_k[counter]);		//6
 			C_pp[counter]= (p_k[counter]*p_k[counter]+k_pp*k_pp)/(2*k_pp*p_k[counter]);		//6
			
			B_mm[counter]= ((1+3*C_mm[counter]*C_mm[counter])/2)*log(abs((1+C_mm[counter])/(1-C_mm[counter]))) - 															3*C_mm[counter];
			B_pp[counter]= ((1+3*C_pp[counter]*C_pp[counter])/2)*log(abs((1+C_pp[counter])/(1-C_pp[counter]))) - 															3*C_pp[counter];
			
			pre_mat[counter]= B_mm[counter]*(N_poph(omega_LO,T)+1)*(1-f0(energy_p[minus_index],efef_p,T)) + 								B_mm[counter]*(N_poph(omega_LO,T)+0)*(f0(energy_p[minus_index],efef_p,T))
							+B_pp[counter]*(N_poph(omega_LO,T)+1)*(1-f0(energy_p[plus_index],efef_p,T)) + 									B_pp[counter]*(N_poph(omega_LO,T)+1)*(f0(energy_p[plus_index],efef_p,T));
		

			if (scattering_mechanisms[1] == 1){
				if ((((1+3*C_mm[counter]*C_mm[counter])/2)*log(abs((1+C_mm[counter])/(1-C_mm[counter]))) - 3*C_mm[counter]==0)&&(((1+3*C_pp[counter]*C_pp[counter])/2)*log(abs((1+C_pp[counter])/(1-C_pp[counter]))) - 3*C_pp[counter]) ==0)
					
					p_S_o_grid[counter] = average_dummy;		
			    // Just to avoid instability since S_o is the only denominator in g_LO and cannot be zero
				else
				p_S_o_grid[counter] = (1/(epsilon_inf[counter]*epsilon_0)-1/            	(epsilon_s[counter]*epsilon_0))*pre_mat[counter]*(e*e*omega_LO)/(16*pi*h_bar*v_p[counter]);

}
// }	
 }
	    
	 
for ( counter=0;counter<points;counter++)
	nu_So_p[counter][2][2] = p_S_o_grid[counter];				
		
		FILE *fid1;
		    fid1 = fopen("nu_So_pop.txt","w");
		    for (int i = 0; i < points; i++)
		        fprintf(fid1,"%d    %f  \n", i+1, nu_So_p[i][2][2]);
		fclose(fid1);
}
