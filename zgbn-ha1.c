
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "helpers.h"

//int vecinoaleatorio(int I,int u[N],int d[N],int l[N],int r[N],int ur[N],int ul[N],int dr[N],int dl[N]);
void reaccionCO(MATRIZ m, int I, int *FE,int *O);

// Hubo reaccion en la posicion pos y la posicion posReaccion
void reaccionoCO(MATRIZ m, int pos, int posReaccion, int *O, int *FE);
void revisarReaccionCO(MATRIZ m, int pos, int posRevision[], int len, int *FE, int *O);

void reaccionO2(MATRIZ m, int I, int k,int *FE,int *CO,int *CO2);
void revisarReaccionO2(MATRIZ m, int pos, int orden[], int len, int *CO, int *CO2);

void warm(double y, MATRIZ m, int *FE,int *CO,int *O, int *CO2,int *X);

MATRIZ crearMatriz();


int rdtsc() {
	asm volatile ("rdtsc");
}


MATRIZ crearMatriz() {
  MATRIZ m;
  int i, j, ind;
  int il, ir, ju, jd; // left of i, right of i, up of j, down of j
  
  for (i=0; i < L; i++) {
    for (j=0; j < L; j++) {
      ind = POS(i, j);
      il = DECR(i);
      ir = INCR(i);
      ju = DECR(j);
      jd = INCR(j);
      
      m.u[ind] = POS(i, ju);
      m.d[ind] = POS(i, jd);
      m.l[ind] = POS(il, j);
      m.r[ind] = POS(ir, j);
      m.ul[ind] = POS(il, ju);
      m.ur[ind] = POS(ir, ju);
      m.dl[ind] = POS(il, jd);
      m.dr[ind] = POS(ir, jd);
    }
  }
  
  return m;
}

int main(int argc, char *argv[]) {
  MATRIZ m;
  int i,j,k,k1,k2,p,I,J,x,v,sem,bd,br,bl,vec,
    CO2,	// Numero de produccion de CO2
    PMC,
    FE,
    X,	// Numero de impurezas
    Otot,
    COtot,
    CO,
    O;
  double CO2p,Oc,COc,Xc,suma,produc,q,y,yoc,yo,kco,kx,sO,s2O,sCO,s2CO,sCO2,s2CO2,s2X,sigma[4],resta,invN,invNP,inv;
  char sal[99];

  srand(rdtsc());
  FILE *fp;

  sprintf(sal,"ZGBHa1%d %d %d %lf.dat",L,Cal,pasos,yx);
  fp = fopen(sal,"w");



  i=0; j=0; k=0;  
  invN=(double)1/(N);
  invNP=(double)(invN/(pasos));
  inv=((pasos*N)-1);
  inv=(double)(1/(inv));

//___________________creamos la matriz de vecinos_____________
	m = crearMatriz();
 //____________________________Probabilidad Yco_________________________
	y=yx;
 	for (y=yx;y<1;y+=0.01) {
    yo  = 1-y-yx;
    yoc = yo+y;

    p=0; CO=0; CO2=0; O=0; FE=N; X=0; s2O=0; s2CO=0; s2X=0;
    sigma[0]=0;sigma[1]=0;sigma[2]=0; sigma[3]=0;
    
    Oc=0; COc=0; produc=0; s2CO2=0; Xc=0; 
 //____________________________Iniciamos la matriz y calentamos_______________________
    warm(y, m, &FE, &CO, &O, &CO2, &X);
 //____________________________Llenamos la matriz________________________
    p=0;

		for(p=0;p<(pasos);++p) {
      J=0; 
			for(J=0; J < N; ++J) {
        I = N * ALEATORIO;
        if (m.cell[I] == 0) {
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
          q = ALEATORIO;
 //__________________________________ CO=1
				
          if (q <= y) { // Probabilidad de absorcion de CO
            m.cell[I] = MONOXIDO; 
            reaccionCO(m, I, &FE, &O);
            if (m.cell[I] == MONOXIDO) { ++CO; --FE; }
            else  if (m.cell[I] == LIBRE) ++CO2;
          }
 //__________________________________  O=-1			
          else if (q < yoc) { // Probabilidad de absorcion de O2
            v = 2 * ALEATORIO;
            if (v == 0) {
              k2 = m.l[I];
              k1 = m.r[I];
              if (m.cell[k1]==LIBRE && m.cell[k2]==LIBRE){
                m.cell[m.l[I]] = OXIGENO;
                m.cell[m.r[I]] = OXIGENO;
                
                k=k1;
                reaccionO2(m, I, k, &FE, &CO, &CO2);
                if (m.cell[k] == OXIGENO) {++O; --FE;}
                else if (m.cell[I] == LIBRE) ++CO2;
                  
                k=k2;
                reaccionO2(m, I, k, &FE, &CO, &CO2);
                if (m.cell[k] == OXIGENO) {++O; --FE;}
                else if (m.cell[I] == LIBRE) ++CO2;
                          
              }
            }
            else {
              k2 = m.u[I];
              k1 = m.d[I];
              if (m.cell[m.u[I]]==0 && m.cell[m.d[I]]==0) {
                m.cell[m.u[I]] = OXIGENO;
                m.cell[m.d[I]] = OXIGENO;
                
                k = k1;
                reaccionO2(m,I,k,&FE, &CO, &CO2);
                if (m.cell[k] == OXIGENO) {++O; --FE;}
                else if (m.cell[I] == LIBRE) ++CO2;
                      
                k = k2;
                reaccionO2(m, I, k, &FE, &CO, &CO2);
                if (m.cell[k] == OXIGENO) { ++O; --FE; }
                else if (m.cell[I] == LIBRE) ++CO2;
              }
            }
          }
          else {
          m.cell[I] = IMPUREZA;
          ++X; --FE;
        }
        }
      }
			//}

      Oc += O;			
      COc += CO;
      Xc += X;
      //produc=CO2[0];
      s2O += O*O;
        
      s2CO += CO*CO;
      s2CO2 += CO2*CO2;
      suma = FE+O+CO;
    }
    
    produc = CO2;
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

void revisarReaccionCO(MATRIZ m, int pos, int posRevision[], int len, int *FE, int *O) {
  int i;
  for (i=0; i < len; i++) {
    if (m.cell[posRevision[i]] == OXIGENO) {
      m.cell[i] = LIBRE;
      m.cell[posRevision[i]] = LIBRE;
      --*O;
      ++*FE;
      return;
    }
  }
}

/* Si hay oxigeno cerca de la posicion, reacciona con una al azar. Si */
void reaccionCO(MATRIZ m, int I, int *FE, int *O) {
  int j;
  int len = 4;
  int posiciones[len];
  
  j = 4 * ALEATORIO;
  
  switch(j) {
    case 0:
      posiciones[0] = m.u[I]; posiciones[1] = m.d[I]; posiciones[2] = m.l[I]; posiciones[3] = m.r[I];
      break;
      
    case 1:
      posiciones[1] = m.d[I]; posiciones[2] = m.l[I]; posiciones[3] = m.r[I]; posiciones[0] = m.u[I];
      break;
      
    case 2:
      posiciones[2] = m.l[I]; posiciones[3] = m.r[I]; posiciones[0] = m.u[I]; posiciones[1] = m.d[I];
      break;
      
    case 3:
      posiciones[3] = m.r[I]; posiciones[0] = m.u[I]; posiciones[1] = m.d[I]; posiciones[2] = m.l[I];
      break;
  }
  
  revisarReaccionCO(m, I, posiciones, len, FE, O);
}

//_____________________REACCION DEL O2

void revisarReaccionO2(MATRIZ m, int pos, int orden[], int len, int *CO, int *FE) {
  int i;
  for (i=0; i < len; i++) {
    if (m.cell[orden[i]] == MONOXIDO) {
      m.cell[pos] = LIBRE;
      m.cell[orden[i]] = LIBRE;
      --*CO;
      ++*FE;
      return;
    }
  }
}

void reaccionO2(MATRIZ m, int I, int k, int *FE, int *CO, int *CO2) {
	int j = 4 * ALEATORIO;
  int len = 4;
  int orden[len];
  
  switch(j) {
    case 0:
      orden[1] = m.d[k]; orden[2] = m.l[k]; orden[3] = m.r[k]; orden[0] = m.u[k];
      break;
      
    case 1:
      orden[2] = m.l[k]; orden[3] = m.r[k]; orden[0] = m.u[k]; orden[1] = m.d[k];
      break;
      
    case 2:
      orden[3] = m.r[k]; orden[0] = m.u[k]; orden[1] = m.d[k]; orden[2] = m.l[k];
      break;
      
    case 3:
      orden[0] = m.u[k]; orden[1] = m.d[k]; orden[2] = m.l[k]; orden[3] = m.r[k];
      break;
  }
  
  revisarReaccionO2(m, k, orden, len, CO, FE);
}

//_______________CALENTAMIENTO

// @param y Probabilidad de aceptar un CO
void warm(double y, MATRIZ m, int *FE,int *CO,int *O,int *CO2,int *X) {

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
	for(i=0; i < N; ++i) m.cell[i]=0;

	//_________________________________Calentamos
	p=1;
	for(p=1; p < (Cal+1); ++p) {

		//invNp=((double) (invn/(p)));
		J=0;
		for(J=0; J < N; ++J) {
			I = N * ALEATORIO;	// Posicion aleatoria en la matriz

			if (m.cell[I] == LIBRE) {
	//____________________________Desorpción________________________________

		/*			if (Cell[I]==1){
					q=((double)rand()/RAND_MAX);

							if (q<kco){

							*CO=*CO-1;
							Cell[I]=0;
				}
			}

					if (Cell[I]==2){
					q=((double)rand()/RAND_MAX);

							if (q<kx){

							*CO=*CO-1;
							Cell[I]=0;
				}
			}
	*/
	//______________________________Absorcion___________________________________

        q = ALEATORIO;

	//__________________________________ CO=1

				if (q <= y) {	// Absorbe CO, quizas reacciona
					m.cell[I] = MONOXIDO; // 1 = CO
					reaccionCO(m,I,FE,O);
					if (m.cell[I] == MONOXIDO) {
						++*CO;
						--*FE;
					}
					else if (m.cell[I] == LIBRE) ++co2;
				}
	//__________________________________  O=-1
				else if (q < yoc) {	// Posiblemente absorbe oxigeno
					v = 2 * ALEATORIO;	// Elegir horizontal o vertical

					if (v == 0) {	// Horizontal
						k2 = m.l[I];
            k1 = m.r[I];
						if (m.cell[k1]==LIBRE && m.cell[k2]==LIBRE) {
							m.cell[m.l[I]] = OXIGENO;
              m.cell[m.r[I]] = OXIGENO;

							k = k1;
							reaccionO2(m,I,k,FE,CO,CO2);
							if (m.cell[k] == OXIGENO) { ++*O; --*FE; }
							//else if (Cell[k]==0) ++co2;

							k = k2;
							reaccionO2(m,I,k,FE,CO,CO2);
              
							if (m.cell[k] == OXIGENO) { ++*O; --*FE; }
							//else if (Cell[k]==0) ++co2;
						}
					}
					else {	// Vertical
						k2 = m.u[I];
            k1 = m.d[I];
						if (m.cell[m.u[I]]==LIBRE && m.cell[m.d[I]]==LIBRE) {
							m.cell[m.u[I]] = OXIGENO;
              m.cell[m.d[I]] = OXIGENO;

							k=k1;
							reaccionO2(m,I,k,FE,CO,CO2);
							if (m.cell[k] == OXIGENO) { ++*O; --*FE; }
							//else if (Cell[k]==0) ++co2;

							k=k2;
							reaccionO2(m, I, k, FE, CO, CO2);
							if (m.cell[k] == OXIGENO) { ++*O; --*FE; }
							//else if (Cell[k]==0) ++co2;
						}
					}
				}
				else {	// Impureza
					m.cell[I] = IMPUREZA;
					++*X;
					--*FE;
				}
			}
		}
/*	oc+=*O;
	coc+=*CO;
	fprintf (fw,"%d  %.4lf  %.4lf  %lf  \n",p,co2*invNp,oc*invNp,coc*invNp);*/
	}
/*		v=0;
for(v=0;v<N;++v){ 
if (v%L==0) printf("\n");

printf(" %d ",Cell[v]);

}
printf("\n");*/

}

	
