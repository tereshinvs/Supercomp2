// This is the explicit conjugate gradient method for descrete Puasson problem 
// on uniform mesh. 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Domain size.
const double A = 120.0;
const double B =  90.0;

int NX, NY;							        // the number of internal points on axes (ox) and (oy).
double hx, hy;                              // mesh steps on (0x) and (0y) axes

#define Test
#define Print
#define Step 10

#define Max(A,B) ((A)>(B)?(A):(B))
#define R2(x,y) ((x)*(x)+(y)*(y))
#define Cube(x) ((x)*(x)*(x))
#define x(i) ((i)*hx)
#define y(j) ((j)*hy)

#define LeftPart(P,i,j)\
((-(P[NX*(j)+i+1]-P[NX*(j)+i])/hx+(P[NX*(j)+i]-P[NX*(j)+i-1])/hx)/hx+\
 (-(P[NX*(j+1)+i]-P[NX*(j)+i])/hy+(P[NX*(j)+i]-P[NX*(j-1)+i])/hy)/hy)

#ifdef Test
    double Solution(double x,double y)
    {    
        double kappa2 = (16.0/R2(A,B));
        
        return sqrt(R2(A,B)) / (1.0 + kappa2*R2(x-A/2,y-B/2));
    }
#endif

double BoundaryValue(double x, double y)
{
    #ifdef Test
        return Solution(x,y);
    #else
        return exp(1.0 - x*x * y*y);
    #endif
}

int RightPart(double * rhs)
{
    int i, j;
    double kappa2 = (16.0/R2(A,B));

    #ifdef Test
        for(j=0; j<NY; j++)
            for(i=0; i<NX; i++)
                rhs[j*NX+i] = 4.0*sqrt(R2(A,B))*kappa2*(1.0-kappa2*R2(x(i)-A/2,y(j)-B/2))\
                            /Cube(1.0 + kappa2*R2(x(i)-A/2,y(j)-B/2));
        return 0;
    #else
        memset(rhs,0,NX*NY*sizeof(double));
    
    	for (j = 0; j < NY; ++j)
    		for (i = 0; i < NX; ++i) {
    			double sqrX = x(i) * x(i);
    			double sqrY = y(i) * y(i);
    			rhs[j*NX + i] = 2.0 * (sqrX + sqrY) * (1.0 - 2.0 * sqrX * sqrY) * exp(1.0 - sqrX * sqrY);
    		}
    
        return 0;
    #endif
}

int main(int argc, char * argv[])
{
	double * SolVect;					    // the solution array
    double * ResVect;					    // the residual array
	double * BasisVect;					    // the vector of A-orthogonal system in CGM
	double * RHS_Vect;					    // the right hand side of Puasson equation.
	double sp, alpha, tau, NewValue, err;	// auxiliary values
	int SDINum, CGMNum;					    // the number of steep descent and CGM iterations.
	int counter;						    // the current iteration number.

	int i,j;
	char str[127];
	FILE * fp;

    // command line analizer
	switch (argc)
	{
	case 4:{
				SDINum = 1;
				CGMNum = atoi(argv[3]);
				break;
		   }
	case 5:{
				SDINum = Max(atoi(argv[3]),1);		// SDINum >= 1
				CGMNum = atoi(argv[4]);
				break;
		   }
	default:{
				printf("Wrong number of parameters in command line.\nUsage: <ProgName> "
                       "<Nodes number on (0x) axis> <Nodes number on (0y) axis> "
                       "[the number of steep descent iterations] "
                       "<the number of conjugate gragient iterations>\nFinishing...\n");
				return(-1);
			}
	}
    
    NX = atoi(argv[1]); NY = atoi(argv[2]);
    hx = A / (NX-1);    hy = B / (NY-1);

	sprintf(str,"PuassonSerial_ECGM_%dx%d.log", NX, NY);
	fp = fopen(str,"w");
	fprintf(fp,"The Domain: [0,%f]x[0,%f], number of points: N[0,A] = %d, N[0,B] = %d;\n"
			   "The steep descent iterations number: %d\n"
			   "The conjugate gradient iterations number: %d\n",
			    A,B, NX,NY, SDINum,CGMNum);
	
	SolVect   = (double *)malloc(NX*NY*sizeof(double));
	ResVect   = (double *)malloc(NX*NY*sizeof(double));
	RHS_Vect  = (double *)malloc(NX*NY*sizeof(double));

// Initialization of Arrays
	memset(ResVect,0,NX*NY*sizeof(double));
    RightPart(RHS_Vect);
    
	for(i=0; i<NX; i++)
	{
		SolVect[i] = BoundaryValue(x(i),0.0);
		SolVect[NX*(NY-1)+i] = BoundaryValue(x(i),B);
	}
	for(j=0; j<NY; j++)
	{
		SolVect[NX*j] = BoundaryValue(0.0,y(j));
		SolVect[NX*j+(NX-1)] = BoundaryValue(A,y(j));
	}

// Iterations ...
	#ifdef Test
		err = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				err = Max(err, fabs(Solution(x(i),y(j))-SolVect[NX*j+i]));
		fprintf(fp,"\nNo iterations have been performed. The residual error is %.12f\n", err);
	#endif

// Steep descent iterations begin ...
	#ifdef Print
		printf("\nSteep descent iterations begin ...\n");
	#endif

	for(counter=1; counter<=SDINum; counter++)
	{
// The residual vector r(k) = Ax(k)-f is calculating ...
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				ResVect[NX*j+i] = LeftPart(SolVect,i,j)-RHS_Vect[NX*j+i];

// The value of product (r(k),r(k)) is calculating ...
		sp = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				sp += ResVect[NX*j+i]*ResVect[NX*j+i]*hx*hy;
		tau = sp;

// The value of product sp = (Ar(k),r(k)) is calculating ...
		sp = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				sp += LeftPart(ResVect,i,j)*ResVect[NX*j+i]*hx*hy;
		tau = tau/sp;

// The x(k+1) is calculating ...
		err = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
			{
				NewValue = SolVect[NX*j+i]-tau*ResVect[NX*j+i];
				err = Max(err, fabs(NewValue-SolVect[NX*j+i]));
				SolVect[NX*j+i] = NewValue;
			}

        if(counter%Step == 0)
        {
            printf("The Steep Descent iteration %d has been performed.\n",counter);
            
    #ifdef Print
            fprintf(fp,"\nThe Steep Descent iteration k = %d has been performed.\n"
					"Step \\tau(k) = %f.\nThe difference value is estimated by %.12f.\n",\
					counter, tau, err);
    #endif

    #ifdef Test
			err = 0.0;
			for(j=1; j < NY-1; j++)
				for(i=1; i < NX-1; i++)
					err = Max(err, fabs(Solution(x(i),y(j))-SolVect[NX*j+i]));
			fprintf(fp,"The Steep Descent iteration %d have been performed. "
					   "The residual error is %.12f\n", counter, err);
    #endif
        }
    }
// the end of steep descent iteration.

	BasisVect = ResVect;    // g(0) = r(k-1).
    ResVect = (double *)malloc(NX*NY*sizeof(double));
    memset(ResVect,0,NX*NY*sizeof(double));

// CGM iterations begin ...
// sp == (Ar(k-1),r(k-1)) == (Ag(0),g(0)), k=1.
	#ifdef Print
		printf("\nCGM iterations begin ...\n");
	#endif

	for(counter=1; counter<=CGMNum; counter++)
	{
	// The residual vector r(k) is calculating ...
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				ResVect[NX*j+i] = LeftPart(SolVect,i,j)-RHS_Vect[NX*j+i];

	// The value of product (Ar(k),g(k-1)) is calculating ...
		alpha = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				alpha += LeftPart(ResVect,i,j)*BasisVect[NX*j+i]*hx*hy;
		alpha = alpha/sp;

	// The new basis vector g(k) is being calculated ...
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				BasisVect[NX*j+i] = ResVect[NX*j+i]-alpha*BasisVect[NX*j+i];

	// The value of product (r(k),g(k)) is being calculated ...
		tau = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				tau += ResVect[NX*j+i]*BasisVect[NX*j+i]*hx*hy;
		
	// The value of product sp = (Ag(k),g(k)) is being calculated ...
		sp = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				sp += LeftPart(BasisVect,i,j)*BasisVect[NX*j+i]*hx*hy;
		tau = tau/sp;

	// The x(k+1) is being calculated ...
		err = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
			{
				NewValue = SolVect[NX*j+i]-tau*BasisVect[NX*j+i];
				err = Max(err, fabs(NewValue-SolVect[NX*j+i]));
				SolVect[NX*j+i] = NewValue;
			}

		if(counter%Step == 0)
		{
			printf("The %d iteration of CGM method has been carried out.\n", counter);
            
#ifdef Print            
            fprintf(fp,"\nThe iteration %d of conjugate gradient method has been finished.\n"
                       "The value of \\alpha(k) = %f, \\tau(k) = %f. The difference value is %f.\n",\
                        counter, alpha, tau, err);
#endif

#ifdef Test
			err = 0.0;
			for(j=1; j < NY-1; j++)
				for(i=1; i < NX-1; i++)
					err = Max(err, fabs(Solution(x(i),y(j))-SolVect[NX*j+i]));
			fprintf(fp,"The %d iteration of CGM have been performed. The residual error is %.12f\n",\
                        counter, err);
#endif
		}
	}
// the end of CGM iterations.

// printing some results ...
	fprintf(fp,"\nThe %d iterations are carried out. The error of iterations is estimated by %.12f.\n",
                SDINum+CGMNum, err);
	fclose(fp);

	sprintf(str,"PuassonSerial_ECGM_%dx%d.dat", NX, NY);
	fp = fopen(str,"w");
		fprintf(fp,"# This is the conjugate gradient method for descrete Puasson equation.\n"
				"# A = %f, B = %f, N[0,A] = %d, N[0,B] = %d, SDINum = %d, CGMNum = %d.\n"
				"# One can draw it by gnuplot by the command: splot 'MyPath\\FileName.dat' with lines\n",\
				A, B, NX, NY, SDINum, CGMNum);
		for (j=0; j < NY; j++)
		{
			for (i=0; i < NX; i++)
				fprintf(fp,"\n%f %f %f", x(i), y(j), SolVect[NX*j+i]);
			fprintf(fp,"\n");
		}
	fclose(fp);

	free(SolVect); free(ResVect); free(BasisVect); free(RHS_Vect);
	return(0);
}
