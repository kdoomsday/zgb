
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


# define L 	120
# define N	L*L
# define Cal	50000
# define pasos	100000
# define yx	0.1

//int vecinoaleatorio(int I,int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N]);
void reaccionCO( int Cell[N], int I,int FE[1],int O[1],int CO2[1],int u[N],int d[N],int l[N],int r[N]);
void reaccionO2(int Cell[N], int I, int k,int FE[1],int CO[1],int CO2[1],int u[N],int d[N],int l[N],int r[N]);
void warm(double y,int Cell[N], int FE[1],int CO[1],int O[1], int CO2[1],int X[1],int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N]);
int rdtsc()
{
asm volatile ("rdtsc");
}

main(){
int Cell[N],i,j,k,k1,k2,p,I,J,x,v,sem,u[N],d[N],l[N],r[N],ur[N],ul[N],dr[N],dl[N],bd,br,bl,vec,CO2[1],PMC,FE[1],X[1],Otot,COtot,CO[1],O[1];
double CO2p,Oc,COc,Xc,suma,produc,q,y,yoc,yo,kco,kx,sO,s2O,sCO,s2CO,sCO2,s2CO2,s2X,sigma[4],resta,invN,invNP,inv;
char sal[99];

srand(rdtsc());
FILE *fp;

sprintf(sal,"ZGBHa1%d %d %d %lf.dat",L,Cal,pasos,yx);
fp=fopen(sal,"w");



i=0; j=0; k=0;  
invN=(double)1/(N);
invNP=(double)(invN/(pasos));
inv=((pasos*N)-1);
inv=(double)(1/(inv));

//___________________creamos la matriz de vecinos_____________

	
for(i=0;i<(L-1);++i){
	
	for (j=((i+1)*L+1);j<((i+2)*L-1);++j){
	
//___________________Vecinos primarios________________________
	
	u[j]=j-L;
	d[j]=j+L;
	l[j]=j-1;
	r[j]=j+1;
//___________________Vecinos secundarios___________________________________	
	
	ur[j]=u[j]+1;
	ul[j]=u[j]-1;
	dr[j]=d[j]+1;
	dl[j]=d[j]-1;


}

}
i=1;
for(i=1;i<L;++i){

//___________________Borde superior_______________________________________

	u[i]=N-L+i;
	d[i]=i+L;
	l[i]=i-1;
	r[i]=i+1;

	ur[i]=u[i]+1;
	ul[i]=u[i]-1;
	dr[i]=d[i]+1;
	dl[i]=d[i]-1;

//____________________Borde Derecho_______________________________________	
	
	br=i*L-1;
	
	u[br]=br-L;
	d[br]=br+L;
	l[br]=br-1;
	r[br]=br-L+1;

	ur[br]=r[br]-L;
	ul[br]=u[br]-1;
	dr[br]=r[br]+L;
	dl[br]=d[br]-1;

//___________________Borde Inferior_______________________________________

	bd=N-L+i-1;
	
	u[bd]=bd-L;
	d[bd]=i-1;
	l[bd]=bd-1;
	r[bd]=bd+1;

	ur[bd]=u[bd]+1;
	ul[bd]=u[bd]-1;
	dr[bd]=d[bd]+1;
	dl[bd]=d[bd]-1;

	
//___________________Borde izquierdo_______________________________________

	bl=i*L;
	
	u[bl]=bl-L;
	d[bl]=bl+L;
	l[bl]=(i+1)*L-1;
	r[bl]=bl+1;

	ur[bl]=u[bl]+1;
	ul[bl]=l[bl]-L;
	dr[bl]=d[bl]+1;
	dl[bl]=l[bl]+L;

}

//____________________Esquinas______________________________________

        u[0]=N-L;
	d[0]=L;
	l[0]=L-1;
	r[0]=1;

	ur[0]=u[0]+1;
	ul[0]=N-1;
	dr[0]=d[0]+1;
	dl[0]=2*L-1;
//______________________esquina derecha sup
	u[L-1]=N-1;
	d[L-1]=2*L-1;
	l[L-1]=L-2;
	r[L-1]=0;

	ur[L-1]=N-L;
	ul[L-1]=u[L-1]-1;
	dr[L-1]=L;
	dl[L-1]=d[L-1]-1;	
//______________________esquina inferior izq

	u[N-L]=N-2*L;
	d[N-L]=0;
	l[N-L]=N-1;
	r[N-L]=N-L+1;

	ur[N-L]=u[N-L]+1;
	ul[N-L]=N-L-1;
	dr[N-L]=d[N-L]+1;
	dl[N-L]=L-1;

	u[N-1]=N-L-1;
	d[N-1]=L-1;
	l[N-1]=N-2;
	r[N-1]=N-L;

	ur[N-1]=N-2*L;
	ul[N-1]=u[N-1]-1;
	dr[N-1]=0;
	dl[N-1]=d[N-1]-1;
	
	

y=yx;
 //____________________________Probabilidad Yco_________________________

 	for (y=yx;y<1;y+=0.01){

 	yo=1-y-yx;
 	yoc=yo+y;

 	p=0; CO[0]=0; CO2[0]=0; O[0]=0; FE[0]=N; X[0]=0; s2O=0; s2CO=0; s2X=0; sigma[0]=0;sigma[1]=0;sigma[2]=0; sigma[3]=0;
 	Oc=0; COc=0; produc=0; s2CO2=0; Xc=0; 


 //____________________________Iniciamos la matriz y calentamos_______________________

 	warm(y,Cell,FE,CO,O,CO2,X,u,d,l,r,ur,ul,dr,dl);


 //____________________________Llenamos la matriz________________________

 	p=0;

		for(p=0;p<(pasos);++p){
		J=0; 
			for(J=0;J<N;++J){
				
		
			I=N*((double)rand()/RAND_MAX);
 		if (Cell[I]==0) {
 //____________________________Desorpción________________________________

	 /*			if (Cell[I]==1){ 
				q=((double)rand()/RAND_MAX);
	
						if (q<kco){ 
		
						CO[0]=CO[0]-1;
						Cell[I]=0;
			}
		}
	
				if (Cell[I]==2){ 
				q=((double)rand()/RAND_MAX);
	
						if (q<kx){ 
		
						CO[0]=CO[0]-1;
						Cell[I]=0;
			}
		}
*/	
 //______________________________Absorcion___________________________________

			q=((double)rand()/RAND_MAX);
				
 //__________________________________ CO=1
				
			if (q<=y) {
										
				Cell[I]=1; 
				reaccionCO(Cell,I,FE,O,CO2,u,d,l,r);
				if (Cell[I]==1) {++CO[0]; --FE[0];}
				else  if (Cell[I]==0) ++CO2[0];
				
			}
 //__________________________________  O=-1			
			else if (q<yoc) {
				
				v=2*((double)rand()/RAND_MAX);
				if (v==0) {
					k2=l[I]; k1=r[I];
					if ((Cell[k1]==0)&&(Cell[k2]==0)){
						Cell[l[I]]=-1; Cell[r[I]]=-1;
						
						k=k1;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						else if (Cell[I]==0) ++CO2[0];
							
						k=k2;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						else if (Cell[I]==0) ++CO2[0];
											
					}
				}
									
				if (v==1) {
										
					k2=u[I]; k1=d[I];
					if ((Cell[u[I]]==0)&&(Cell[d[I]]==0)){
						Cell[u[I]]=-1; Cell[d[I]]=-1;
						
						k=k1;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						else if (Cell[I]==0) ++CO2[0];
									
						k=k2;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						else if (Cell[I]==0) ++CO2[0];
										
					}
				}
												
										
			}
			
			else { Cell[I]=2; ++X[0]; --FE[0];
			}
		}
					
	}	
			//}
			

	

		Oc+=O[0];			
		COc+=CO[0];
		Xc+=X[0];
		//produc=CO2[0];
		s2O+=O[0]*O[0];
			
		s2CO+=CO[0]*CO[0];
		s2CO2+=CO2[0]*CO2[0];;
		suma=FE[0]+O[0]+CO[0];
			
			
			
			}
	 produc=CO2[0];
 //________________________________Promedios___________________
	 suma=((double)suma*invNP); 
	 Oc=((double)Oc*invNP); 
	 COc=((double)COc*invNP);
	 produc=((double)produc*(invNP));
	 Xc=((double)Xc*(invNP));
	
	
 //__________________________Promedios del cuadrado
	 s2O=((double)s2O*invNP);
	 s2CO=((double)s2CO*invNP);
	 s2CO2=((double)s2CO2*invNP);
	 
	
 //________________________________Error del O_________________
	
	 sigma[0]=s2O-(Oc*Oc); 
	 sigma[0]=(double)(sigma[0])*inv;
	 sigma[0]=sqrt(sigma[0]);
	
 //_______________________________Error del CO
	 sigma[1]=s2CO-(COc*COc);
	 sigma[1]=(double)(sigma[1])*inv;
	 sigma[1]=sqrt(sigma[1]);
	
 //_______________________________Error del CO2
	 sigma[2]=s2CO2-(produc*produc);
	 sigma[2]=(double)(sigma[2])*inv;
	 sigma[2]=sqrt(sigma[2]);
	 
   //_______________________________Error de X
	 sigma[3]=s2X-(Xc*Xc);
	 sigma[3]=(double)(sigma[3])*inv;
	 sigma[3]=sqrt(sigma[3]);
	
	fprintf(fp,"%.4lf  %.4lf  %.4lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",y,Oc,COc,produc,Xc,suma,sigma[0],sigma[1],sigma[2],sigma[3]);

	}
}
//---------------------------FUNCIONES


/*int vecinoaleatorio(int I,int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N]){
int z,v;

	v=4*((double)rand()/RAND_MAX);
	switch(v){
	/*case 0:
	J=u[I];
	break;
	case 1:
	J=d[I];
	break;
	case 2:
	J=l[I];
	break;
	case 3:
	J=r[I];
	break;
	case 0:
	z=ur[I];
	break;
	case 1:
	z=ul[I];
	break;
	case 2:
	z=dr[I];
	break;
	case 3:
	z=dl[I];
	break;
	return (z); }
}*/

//______________REACCION DEL CO

void reaccionCO( int Cell[N], int I,int FE[1],int O[1],int CO2[1],int u[N],int d[N],int l[N],int r[N]){
	int j;


j=0; 

	j=4*((double)rand()/RAND_MAX);
	switch(j){
	case 0:
	
	if (Cell[u[I]]==-1){
	
	
		Cell[I]=0; 
		Cell[u[I]]=0;
		--O[0];
		++FE[0]; 
		
		}
	
	else { 
	
	if (Cell[d[I]]==-1){
		
		Cell[I]=0;
		Cell[d[I]]=0;
		--O[0];
		++FE[0];
		
		}
								
	else{ 
				
	if (Cell[l[I]]==-1){	
	
		Cell[I]=0;
		Cell[l[I]]=0;
		--O[0];
		++FE[0];
			
		}
	
	else { 
	
	if (Cell[r[I]]==-1){
	
		Cell[I]=0;
		Cell[r[I]]=0;
		--O[0];
		++FE[0];
			
		}
		}
		}
		}
	break;
	case 1:
	if (Cell[d[I]]==-1){

		Cell[I]=0;
		Cell[d[I]]=0;
		--O[0];
		++FE[0];
			
		}

	else { 
	
	if (Cell[l[I]]==-1){

		Cell[I]=0;
		Cell[l[I]]=0;
		--O[0];
		++FE[0];
			
		}
							
	else{ 
				
	if (Cell[r[I]]==-1){
	
		Cell[I]=0;
		Cell[r[I]]=0;
		--O[0];
		++FE[0];

		}
	
	else { 
	
	if (Cell[u[I]]==-1){
	
		Cell[I]=0;
		Cell[u[I]]=0;
		--O[0];
		++FE[0];

				
		}
		}
		}
		}
	
	break;
	
	case 2:
	
	if (Cell[l[I]]==-1){
	
		Cell[I]=0;
		Cell[l[I]]=0;
		--O[0];
		++FE[0];

		}
	
	else { 
	
	if (Cell[r[I]]==-1){
	
		Cell[I]=0;
		Cell[r[I]]=0;
		--O[0];
		++FE[0];

		}
							
	else{ 
				
	if (Cell[u[I]]==-1){
	
		Cell[I]=0;
		Cell[u[I]]=0;
		--O[0];
		++FE[0];

		}
	
	else { 
	
	if (Cell[d[I]]==-1){
	
		Cell[I]=0;
		Cell[d[I]]=0;
		--O[0];
		++FE[0];
			
		}
		}
		}
		}

	break;

	case 3:
	if (Cell[r[I]]==-1){
				Cell[I]=0;
				Cell[r[I]]=0;
				--O[0];
				++FE[0];

				
				
			}
			else { 
	
	if (Cell[u[I]]==-1){
				Cell[I]=0;
				Cell[u[I]]=0;
				--O[0];
				++FE[0];

				}
								
				else{ 
				
	if (Cell[d[I]]==-1){
				Cell[I]=0;
				Cell[d[I]]=0;
				--O[0];
				++FE[0];

				
				}
				else { 
	
	if (Cell[l[I]]==-1){
				Cell[I]=0;
				Cell[l[I]]=0;
				--O[0];
				++FE[0];

				
			}
			}
			}
			}
	break;
	
	}
}

//_____________________REACCION DEL O2

void reaccionO2(int Cell[N],int I, int k,int FE[1],int CO[1],int CO2[1],int u[N],int d[N],int l[N],int r[N]){
int j;


j=0; 
	j=4*((double)rand()/RAND_MAX);
			
switch(j){
	case 0:
	if (Cell[u[k]]==1){
				Cell[k]=0;
				Cell[u[k]]=0;
				--CO[0];
				++FE[0];

				
				}
				else { 
	
	if (Cell[d[k]]==1){
				Cell[k]=0;
				Cell[d[k]]=0;
				--CO[0];
				++FE[0];

				
				}
								
				else{ 
				
	if (Cell[l[k]]==1){
				Cell[k]=0;
				Cell[l[k]]=0;
				--CO[0];
				++FE[0];

				
				}
				else { 
	
	if (Cell[r[k]]==1){
				Cell[k]=0;
				Cell[r[k]]=0;
				--CO[0];
				++FE[0];

			
			}
			}
			}
			}
	break;
	case 1:
	if (Cell[d[k]]==1){
				Cell[k]=0;
				Cell[d[k]]=0;
				--CO[0];
				++FE[0];

				}
				else { 
	
	if (Cell[l[k]]==1){
				Cell[k]=0;
				Cell[l[k]]=0;
				--CO[0];
				++FE[0];
				
				}
							
				else{ 
				
	if (Cell[r[k]]==1){
				Cell[k]=0;
				Cell[r[k]]=0;
				--CO[0];
				++FE[0];

				
				}
				else { 
	
	if (Cell[u[k]]==1){
				Cell[k]=0;
				Cell[u[k]]=0;
				--CO[0];
				++FE[0];

			
			}
			}
			}
			}
	break;
	case 2:
	if (Cell[l[k]]==1){
				Cell[k]=0;
				Cell[l[k]]=0;
				--CO[0];
				++FE[0];

				
				}
				else {
	
	if (Cell[r[k]]==1){
				Cell[k]=0;
				Cell[r[k]]=0;
				--CO[0];
				++FE[0];

				
				}
							
				else{ 
				
	if (Cell[u[k]]==1){
				Cell[k]=0;
				Cell[u[k]]=0;
				--CO[0];
				++FE[0];

				
				}
				else { 
	
	if (Cell[d[k]]==1){
				Cell[k]=0;
				Cell[d[k]]=0;
				--CO[0];
				++FE[0];

			
			}
			}
			}
			}
	break;
	case 3:
	if (Cell[r[k]]==1){
	
		Cell[k]=0;
				Cell[r[k]]=0;
				--CO[0];
				++FE[0];

				
				
			}
			else {
	
	if (Cell[u[k]]==1){
				Cell[k]=0;
				Cell[u[k]]=0;
				--CO[0];
				++FE[0];

				}
								
				else{ 
				
	if (Cell[d[k]]==1){
				Cell[k]=0;
				Cell[d[k]]=0;
				--CO[0];
				++FE[0];

				
				}
				else { 
	
	if (Cell[l[k]]==1){
	
	
				Cell[k]=0;
				Cell[l[k]]=0;
				--CO[0];
				++FE[0];

				
			}
			}
			}
			}
	break;
	
	}

}

//_______________CALENTAMIENTO

void warm(double y, int Cell[N], int FE[1],int CO[1],int O[1],int CO2[1],int X[1],int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N]){

int i,p,J,I,k,k1,k2,x,v;
double q,yoc,yo,invNp,invn,co2,oc,coc;
char sal[99];
/*FILE *fw;

sprintf(sal, "Calentamiento MCS %lf .dat",y);
fw=fopen(sal,"w");*/

yo=1-y-yx;
yoc=yo+y;
/*invn=((double)1/(N));

co2=0;
oc=0;
coc=0;*/
//_______________________________Iniciamos la matriz

i=0;
for(i=0;i<N;++i) Cell[i]=0;

//_________________________________Calentamos
p=1;
for(p=1;p<(Cal+1);++p){
	
	//invNp=((double) (invn/(p)));	
	J=0; 
	for(J=0;J<N;++J){		
	I=N*((double)rand()/RAND_MAX);

if (Cell[I]==0) {
 //____________________________Desorpción________________________________

	 /*			if (Cell[I]==1){ 
				q=((double)rand()/RAND_MAX);
	
						if (q<kco){ 
		
						CO[0]=CO[0]-1;
						Cell[I]=0;
			}
		}
	
				if (Cell[I]==2){ 
				q=((double)rand()/RAND_MAX);
	
						if (q<kx){ 
		
						CO[0]=CO[0]-1;
						Cell[I]=0;
			}
		}
*/	
 //______________________________Absorcion___________________________________

			q=((double)rand()/RAND_MAX);
				
 //__________________________________ CO=1
				
			if (q<=y) {
										
				Cell[I]=1; 
				reaccionCO(Cell,I,FE,O,CO2,u,d,l,r);
				if (Cell[I]==1) {++CO[0]; --FE[0];}
				else if (Cell[I]==0) ++co2;
				
			}
 //__________________________________  O=-1			
			else if (q<yoc) {
				
				v=2*((double)rand()/RAND_MAX);
				if (v==0) {
					k2=l[I]; k1=r[I];
					if ((Cell[k1]==0)&&(Cell[k2]==0)){
						Cell[l[I]]=-1; Cell[r[I]]=-1;
						
						k=k1;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						//else if (Cell[k]==0) ++co2;
							
						k=k2;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						//else if (Cell[k]==0) ++co2;
											
					}
				}
									
				if (v==1) {
										
					k2=u[I]; k1=d[I];
					if ((Cell[u[I]]==0)&&(Cell[d[I]]==0)){
						Cell[u[I]]=-1; Cell[d[I]]=-1;
						
						k=k1;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						//else if (Cell[k]==0) ++co2;
									
						k=k2;
						reaccionO2(Cell,I,k,FE,CO,CO2,u,d,l,r);
						if (Cell[k]==-1) {++O[0]; --FE[0];}
						//else if (Cell[k]==0) ++co2;
										
					}
				}
												
										
			}
			
			else { Cell[I]=2; ++X[0]; --FE[0];
			}
		}
	}
/*	oc+=O[0];
	coc+=CO[0];
	fprintf (fw,"%d  %.4lf  %.4lf  %lf  \n",p,co2*invNp,oc*invNp,coc*invNp);*/
}
/*		v=0;
for(v=0;v<N;++v){ 
if (v%L==0) printf("\n");

printf(" %d ",Cell[v]);

}
printf("\n");*/

		}

	
