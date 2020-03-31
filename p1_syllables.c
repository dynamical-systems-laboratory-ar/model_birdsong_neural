//-------------------------------------------
/*To compile this file:
gcc p1_syllables.c -lm -o p1_syllables
(rk4.c file should be in the same folder)
*/
//-------------------------------------------
//-------------------------------------------
#include <stdio.h> 
#include <math.h>
#include "rk4.c"
//-------------------------------------------
#define PI 3.14159265359
//-------------------------------------------
  FILE *ptr1,*ptr2,*ptr3,*ptr4,*ptr5;
//-------------------------------------------
struct Par_uno {double vk1,vk2,Vk1,Vk2,uk1,uk2,Uk1,Uk2;} a;
struct Par_dos {double vk1,vk2,Vk1,Vk2,uk1,uk2,Uk1,Uk2;} b;
//-------------------------------------------

//-------------------------------------------
void takens_p1(int n, double u[], double du[], double t){
	double x,y,z,w;
                                                                                     
	x=u[0]; y=u[1]; z=u[2]; w=u[3];

  du[0]  = 249.5*(-x + 1/(1+exp(-(-6+10*x-10*y+a.uk1+12*z ))));  	//RAm exc
  du[1]  = 249.5*(-y + 1/(1+exp(-(-8+10*x-(-2)*y+a.uk2+4*z ))));	//RAm inh
  du[2]  = 20*(-z + 1/(1+exp(-(-5.25+10*z-10*w+b.uk1 ))));   		//RA exc
  du[3]  = 20*(-w + 1/(1+exp(-(-5+10*z-(-2)*w+b.uk2 ))));  		    //RA inh
return;
}

//-------------------------------------------
void kick(double t,double t0,double dur_k, double amp1,double amp2,double *para1,double *para2){
if( (t>=t0) && (t<=(t0+dur_k)) ){  *para1 = amp1;    *para2 = amp2; }
}

void fast_kick(double t,double dt,double t0,double dur_k, double amp1,double amp2,double *para1,double *para2,double *t_fire){
 if( (t>=t0) && (t<=(t0+dur_k)) ){
	*para1 = amp1;    *para2 = amp2;
	if( *t_fire<0 ){ *t_fire=t0; }
 }
 if( (t<=(t0+dur_k)) && (t<(t0+dur_k+dt)) ){ *t_fire=t0; }
}
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
main(int argc, char *argv[]){
  double v[4],u[4],V[4],U[4],t,dt,t_max,t_scan;
  double tau1,tau2,t_ini,Para1[20],Para2[20],T[20],w_f[20],a_f[20];
  int i,k,N0;
  double dur[20],t_wait,Tick[20],Tuck[20],forma[20],delay[20],estira[20],period[20];
  double t_fire[20],Tslow[20],check[20],Tack[20];
  double TTick,TTack,dh;
  double cambio[20],cambio_t[20],RoX,RoY,RoZ,RoW; int llm=1;
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
//-------------------------------------------------------------------------------------------

	k=0;	t_ini=0;
	T[0]=0;	Para1[0]=0;	Para2[0]=0;
	//---------------------------------------
	//Parameters:
	period[1]=0.0193355;	T[1]=0.5;	dur[1]=0.74;	
    
	//---------------------------------------
	ptr4=fopen("gesture_p1.dat","w");

	for(i=1;i<5;i++){ v[i]=0; u[i]=0; V[i]=0; U[i]=0;}

	t = 0;		dt = 0.0001;   
	t_scan = 0.0;	t_max = 4.5;

	dh=0.0075;
	Tick[1]=T[1]+0.02+period[1]-dh;

for(i=0;i<20;i++){ 	t_fire[i]=0; Tslow[i]=0; check[i]=-100; }
//------------------

	while(t < t_max){
		a.uk1=0;	a.uk2=0;	b.uk1=0;	b.uk2=0;

		clock(t,&k,N0,T);

		//----------
		//P1:
			kick(t,T[1]-0.014,0.02,0.5,0,&a.uk1,&a.uk2);	      //First pulse (IA)
			kick(t,T[1],0.02,5.5,0,&b.uk1,&b.uk2);		          //First pulse (HVC)
			if( (t>T[1]+0.02+period[1]-dh) && (t<dur[1]) ){
				b.uk2=0;
				fast_kick(t,dt,Tick[1],dh,6,0,&b.uk1,&b.uk2,&t_fire[2]);	//HVC activity
			}
			if ( (t>=Tick[1]+period[1]) && (t<t_max) ){ Tick[1]=Tick[1]+period[1]; }

		//----------
		rk4(takens_p1, u, 4, t, dt);
		t = t+dt;			// time scaled by 1/2.

		if(t > t_scan){
fprintf(ptr4,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",u[0],u[2],t,0.1*a.uk1/(a.uk1+0.001),a.uk2/(a.uk2+0.001),0.1*b.uk1/(b.uk1+0.001),b.uk2/(b.uk2+0.001));
//Reference: $1 Pressure | $2 RAe | $3 time | $4 IAe | $5 IAi = 0| $6 HVCe | &7 HVCi	
		}
		
	}
close(ptr4);
//-------------------------------------------

  return(0); 
}
