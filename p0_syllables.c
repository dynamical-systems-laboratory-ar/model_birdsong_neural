//-------------------------------------------
/*To compile this file:
gcc p0_syllables.c -lm -o p0_syllables
(rk4.c file should be in the same folder)
*/
//-------------------------------------------
#include <stdio.h> 
#include <math.h>
#include "rk4.c"
//-------------------------------------------
#define PI 3.14159265359
//-------------------------------------------
  FILE *ptr1,*ptr2,*ptr3;
//-------------------------------------------
struct Par_uno {double Vk1,Vk2;} a;
struct Par_dos {double Vk1,Vk2;} b;
//-------------------------------------------

//-------------------------------------------
void rewind(FILE *f); 

void takens_p2(int n, double V[], double dV[], double t){
	double x,y,z,w;
                                                                                     
	x=V[0]; y=V[1]; z=V[2]; w=V[3];

    dV[0]  = 249.5*(-x + 1/(1+exp(-(-7.5+9*x-1*y+a.Vk1+9*z ))));  	//ERe
    dV[1]  = 249.5*(-y + 1/(1+exp(-(-11.5+10*x-(-2)*y+a.Vk2+0*z ))));	//ERi
    dV[2]  = 20*(-z + 1/(1+exp(-(-3+6*z-3*w+b.Vk1 ))));   		//RA exc
    dV[3]  = 20*(-w + 1/(1+exp(-(-6+6*z+6*w+b.Vk2 ))));  		//RA inh


return;
}


//-------------------------------------------
void kick(double t,double t0,double dur_k, double amp1,double amp2,double *para1,double *para2){
if( (t>=t0) && (t<=(t0+dur_k)) ){  *para1 = amp1;    *para2 = amp2; }
}

void clock(double t,int *k,int N0,double T[]){
	if ( (t>=T[*k+1]) && (*k<N0) ){
		*k=*k+1;
	}
}

void fast_kick(double t,double dt,double t0,double dur_k, double amp1,double amp2,double *para1,double *para2,double *t_fire){
 if( (t>=t0) && (t<=(t0+dur_k)) ){
	*para1 = amp1;    *para2 = amp2;
	if( *t_fire<0 ){ *t_fire=t0; }
 }
 if( (t<=(t0+dur_k)) && (t<(t0+dur_k+dt)) ){ *t_fire=t0; }
}


void slow_kick(double t,double dt,double delay,double dur_k, double amp1,double amp2,double *para1,double *para2,double t_fire,double *Tslow,double *check){
if( (t>=t_fire+delay) && (*check<0.0) && (t<=(t_fire+delay+dur_k)) ){ *check=t_fire; }
if( (*check>0.0) && (*Tslow<=dur_k) ){	 *para1 = amp1;    *para2 = amp2;	*Tslow=*Tslow+dt;	}
if( *Tslow>dur_k ){ *Tslow=0.0; *check=-100.0; }
}
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
main(int argc, char *argv[]){
  double V[4],t,dt,t_max,t_scan,dur;
  double tau1,tau2,t_ini,Para1[20],Para2[20],T[20],w_f[20],a_f[20];
  int i,k,N0;
  double duracion[20],t_wait,Tick[20],forma[20],delay[20],estira[20],periodo[20];
  double t_fire[20],Tslow[20],check[20],Tack[20];
  double TTick,TTack,di[50],dh[50],delta;
  double Xscan[4],Tscan[4],listo_abs,x_old1,t_old1,x_old2,t_old2,t_old3,x_old3,t_old4,x_old4;
  int listo_max,listo_min,listo_reset,listo_print,listo_end,Gest;
  double bondad_x,bondad_y,bondad_dur,prod,RdM[50],RdP[50],RD[50];
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
	//Parameters:
	delta=0.005;
	T[1]=1.15261;	di[1]=0.011;	dh[1]=0.007;
	T[2]=1.1825;	di[2]=0.014;	dh[2]=0.01;	N0=2;	T[3]=1.21534;
//--------------------------------------------------------------------------------------------
	ptr3=fopen("gesture_p2t3.dat","w");

	for(i=0;i<5;i++){ V[i]=0; }

	t = 0;		dt = 0.0001;   	k=1;
	t_ini = 0.0;	t_max = 2.97;
	listo_max=0;	listo_min=0;	listo_reset=0;	listo_print=0;	listo_end=0;	

	//------------------
	while(t < t_max){

		a.Vk1=0;	a.Vk2=0;	b.Vk1=0;	b.Vk2=0;
		clock(t,&k,N0,T);
        //----------
        // P0:
		kick(t,T[k],di[k],20,0,&a.Vk1,&a.Vk2);		      //IA activity
		kick(t,T[k]+delta,dh[k],20,0,&b.Vk1,&b.Vk2);      //HVC activity
		//------------------
		kick(t,T[1],0.011,14,0,&a.Vk1,&a.Vk2);		     //First pulse (IA)
		kick(t,T[1]+delta,0.007,14,10,&b.Vk1,&b.Vk2);	 //First pulse (HVC)

		//------------------------------------------------------------
		rk4(takens_p2, V, 4, t, dt);			
		t = t+dt;					// time scaled by 1/2.

fprintf(ptr3,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",V[0],V[2],t,0.1*a.Vk1/(a.Vk1+0.001),a.Vk2/(a.Vk2+0.001),0.1*b.Vk1/(b.Vk1+0.001),b.Vk2/(b.Vk2+0.001));
//Reference: $1 Pressure | $2 RAe | $3 time | $4 IAe | $5 IAi = 0| $6 HVCe | &7 HVCi	
		
	}
close(ptr3);
//-------------------------------------------

  return(0); 
}



