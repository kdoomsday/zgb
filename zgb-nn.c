#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

# define L 	40
# define N	L*L
# define Cal	500
# define pasos	800
# define yx	0

int vecinoaleatorio(int I,int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N]);

void reaccionCO( int Cell[N], int I,int FE[1],int O[1],int CO2[1],int u[N],int d[N],int l[N],int r[N]);

void reaccionO2(int Cell[N], int I, int vec, int k,int FE[1],int CO[1],int CO2[1],int u[N],int d[N],int l[N],int r[N]);

void warm(double y,int Cell[N], int FE[1],int CO[1],int O[1], int CO2[1], int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N],double X[1]);

// Esto obtiene una medida de tiempo relacionada con el procesador. Es utilizada como semilla del random
int rdtsc()
{
	asm volatile ("rdtsc");
}

int main() {
	int Cell[N],i,j,k,p,I,J,x,v,sem,u[N],d[N],l[N],r[N],ur[N],ul[N],dr[N],dl[N],bd,br,bl,vec,CO2[1],PMC,FE[1],Otot,COtot,CO[1],O[1];
	double CO2p,Oc,COc,Xc,suma,produc,a,y,yoc,yo,kco,kx,sO,s2O,sCO,s2X,s2CO,sCO2,s2CO2,sigma[3],resta,invN,invNP,inv,X[1];
	char sal[99];

	/*printf("introduzca la semilla=");
	scanf("%d",&sem); printf("\n");
	srand (sem); */ //time(NULL));

	// Inicializar un generador de aleatorios con una semilla que varia
	srand(rdtsc());
	FILE *fp;

	sprintf(sal,"ZGBnn%d %d %d %d.dat",L,Cal,pasos,yx);

	fp=fopen(sal,"w");

	i=0; j=0; k=0; y=0;
	invN=(double)1/(N);
	invNP=(double)(invN/(pasos));
	inv=((pasos*N)-1);
	inv=(double)(1/(inv));
	yo=1-y-yx;
	yoc=yo+y;


	//___________________creamos la matriz de vecinos_____________


	for(i=0;i<(L-1);++i) {
		for (j=((i+1)*L+1);j<((i+2)*L-1);++j) {

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




	//____________________________Probabilidad Yco_________________________

	for (y=0;y<1;y+=0.01){
	yo=1-y-yx;
	yoc=yo+y;

	p=0; CO[0]=0;  O[0]=0; FE[0]=N; s2O=0; s2CO=0; sigma[0]=0; sigma[1]=0;sigma[2]=0; sigma[3]=0;
	Oc=0; COc=0; produc=0; s2CO2=0; X[0]=0;


	//____________________________Iniciamos la matriz y calentamos_______________________

	warm(y,Cell,FE,CO,O,CO2,u,d,l,r,ur,ul,dr,dl,X);
	CO2[0]=0;

	//____________________________Llenamos la matriz________________________


	for(p=0;p<(pasos);++p){
		J=0;
		for(J=0;J<N;++J){
			I=N*((double)rand()/RAND_MAX);

	//____________________________Desorpción________________________________

	/*			if (Cell[I]==1){
				a=((double)rand()/RAND_MAX);

						if (a<kco){

						CO[0]=CO[0]-1;
						Cell[I]=0;
			}
		}

				if (Cell[I]==2){
				a=((double)rand()/RAND_MAX);

						if (a<kx){

						CO[0]=CO[0]-1;
						Cell[I]=0;
			}
		}
	*/

			if (Cell[I]==0) {


	numalazar:	  a=((double)rand()/RAND_MAX);
			if (a==1) { goto numalazar;
			}
	//___________________ CO=1

				else if (a<y) {
				Cell[I]=1;
				reaccionCO(Cell,I,FE,O,CO2,u,d,l,r);
				if (Cell[I]==1) {++CO[0]; --FE[0];
				}
				else if (Cell[I]==0)  ++CO2[0];
				}
	//__________________  O=-1
			//printf("CO=%d\n",CO[0]);
			else if (a<yoc){

				vec=vecinoaleatorio(I,u,d,l,r,ur,ul,dr,dl);
				if (Cell[vec]==0) {

					Cell[I]=-1; Cell[vec]=-1;
					i=2*((int)rand()/RAND_MAX);
					if (i==0) k=I;
					else k=vec;
					reaccionO2(Cell,I,vec,k,FE,CO,CO2,u,d,l,r);
					if (Cell[k]==-1) {++O[0]; --FE[0];}
					else if (Cell[k]==0) ++CO2[0];

					if (k==I) k=vec;
					else k=I;
					reaccionO2(Cell,I,vec,k,FE,CO,CO2,u,d,l,r);
					if (Cell[k]==-1) {++O[0]; --FE[0];}
						else if (Cell[k]==0) ++CO2[0];
				}
			}
			else { Cell[I]=2; ++X[0]; --FE[0];
			}
			// printf("CO=%d\n",CO[0]);
			}
		}

		Oc+=O[0];
		COc+=CO[0];
		Xc+=X[0];
		//produc=CO2[0];
		s2O+=O[0]*O[0];

		s2CO+=CO[0]*CO[0];
		s2X+=X[0]*X[0];
		s2CO2+=CO2[0]*CO2[0];
		suma=FE[0]+O[0]+CO[0]+X[0];

	}
	produc=CO2[0];
	//________________________________Promedios___________________
	suma=((double)suma*invN);
	Oc=((double)Oc*invNP);
	COc=((double)COc*invNP);
	produc=((double)produc*(invNP));
	Xc=((double)Xc*(invNP));

	//__________________________Promedios del cuadrado
	s2O=((double)s2O*invNP);
	s2CO=((double)s2CO*invNP);
	s2CO2=((double)s2CO2*invNP);
	s2X=((double)s2X*invNP);

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


int vecinoaleatorio(int I,int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N]){
int z,v;

	v=4*((double)rand()/RAND_MAX);
	switch(v){
	/*case 0:
	z=u[I];
	break;
	case 1:
	z=d[I];
	break;
	case 2:
	z=l[I];
	break;
	case 3:
	z=r[I];
	break;*/
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
}

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

void reaccionO2(int Cell[N],int I, int vec, int k,int FE[1],int CO[1],int CO2[1],int u[N],int d[N],int l[N],int r[N]){
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

void warm(double y, int Cell[N], int FE[1],int CO[1],int O[1], int CO2[1],int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N],double X[1]){

int i,p,J,I,k,vec,x,v,co2,xc,oc,coc;
double a,yo,yoc,invNp;
char name[99];

FILE *fw;

yo=1-y-yx;
yoc=yo+y;


//_______________________________Iniciamos la matriz

i=0;
for(i=0;i<N;++i) Cell[i]=0;

//_________________________________Calentamos
p=0;
for(p=0;p<(Cal);++p){
	invNp=((double) 1/p);
	invNp= ((double)invNp/N);
	J=0; 
	for(J=0;J<N;++J){		
	I=N*((double)rand()/RAND_MAX);

//____________________________Desorpción________________________________

/*			if (Cell[I]==1){ 
			a=((double)rand()/RAND_MAX);

					if (a<kco){ 
	
					CO[0]=CO[0]-1;
					Cell[I]=0;
		}
	}

			if (Cell[I]==2){ 
			a=((double)rand()/RAND_MAX);

					if (a<kx){ 
	
					CO[0]=CO[0]-1;
					Cell[I]=0;
		}
	}
*/

//__________________________Absorcion________________________

	if (Cell[I]==0) {
		  
		  a=((double)rand()/RAND_MAX);
//___________________ CO=1
		  if (a<y) {
			
			Cell[I]=1; 
			
			reaccionCO(Cell,I,FE,O,CO2,u,d,l,r);
			
			if (Cell[I]==1) {++CO[0]; --FE[0];
			}
			else if (Cell[I]==0) ++co2;
			
		  }
//  O=-1		  
		else if (a<yoc){										
			
			vec=vecinoaleatorio(I,u,d,l,r,ur,ul,dr,dl);
			
			if (Cell[vec]==0) { 
				
				Cell[I]=-1; Cell[vec]=-1;
				i=2*((int)rand()/RAND_MAX);
				if (i==0) k=I;
				else k=vec;
				reaccionO2(Cell,I,vec,k,FE,CO,CO2,u,d,l,r);
				if (Cell[k]==-1) {++O[0]; --FE[0];}
				else if (Cell[I]==0) ++co2;
				
									
				if (k==I) k=vec;
				else k=I;
				reaccionO2(Cell,I,vec,k,FE,CO,CO2,u,d,l,r);
				if (Cell[k]==-1) {++O[0]; --FE[0];}
				else if (Cell[I]==0) ++co2;
			   	
			 }	
			}
		
			else { Cell[I]=2; ++X[0]; --FE[0];
		
		}
		}
	}


oc+=O[0];
coc+=CO[0];

sprintf(name,"ZGBnn%d %d %d %lf Calentamiento y=%lf.dat",L,Cal,pasos,yx,y);
fw=fopen(name,"w");
fprintf(fw,"%d  %.4lf	%.4lf	%.4lf  \n",p,co2*invNp,oc*invNp,coc*invNp);
}
/*		v=0;
for(v=0;v<N;++v){ 
if (v%L==0) printf("\n");

printf(" %d ",Cell[v]);

}
printf("\n");*/
	}
