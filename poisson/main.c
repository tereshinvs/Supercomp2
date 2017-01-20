#include <mpi.h>

#include <omp.h>

// This is the explicit conjugate gradient method for descrete Puasson problem 
// on uniform mesh. 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Domain size.
const long double A = 2.0;
const long double B = 2.0;

int NX, NY;							        // the number of internal points on axes (ox) and (oy).
long double hx, hy;                              // mesh steps on (0x) and (0y) axes

#define Max(A,B) ((A)>(B)?(A):(B))
#define R2(x,y) ((x)*(x)+(y)*(y))

#define x(i) (-A + (i)*hx)
#define y(j) (-B + (j)*hy)

#define LeftPart(P,i,j)\
((-(P[NX*(j)+i+1]-P[NX*(j)+i])/hx+(P[NX*(j)+i]-P[NX*(j)+i-1])/hx)/hx+\
 (-(P[NX*(j+1)+i]-P[NX*(j)+i])/hy+(P[NX*(j)+i]-P[NX*(j-1)+i])/hy)/hy)

#define ResVect(i,j)\
(((i)<=0||(j)<=0||(i)>=NX-1||(j)>=NY-1) ? 0.0 : \
	(LeftPart(SolVect,i,j)-RHS_Vect[NX*(j)+i]))

#define LeftPartResVect(i,j)\
((-(ResVect(i+1,j)-ResVect(i,j))/hx+(ResVect(i,j)-ResVect(i-1,j))/hx)/hx+\
 (-(ResVect(i,j+1)-ResVect(i,j))/hy+(ResVect(i,j)-ResVect(i,j-1))/hy)/hy)


#define ParallelCycle(expression)\
for (int q = displs[myid], i = displs[myid]%NX, j = displs[myid]/NX;\
	q < displs[myid] + size[myid];\
	++q, j = i == NX-1 ? j+1 : j, i = i == NX-1 ? 0 : i+1) {\
	if (i == 0 || i == NX-1 || j == 0 || j == NY-1)\
		continue;\
	{expression;}\
}

long double BoundaryValue(long double x, long double y)
{
    return 1 + sin(x * y);
}

void RightPart(long double* rhs)
{
    int i, j;
    memset(rhs, 0, NX*NY * sizeof(long double));

	for (j = 0; j < NY; ++j)
		for (i = 0; i < NX; ++i)
			rhs[j*NX + i] = (x(i)*x(i) + y(j)*y(j)) * sin(x(i) * y(j));
}

void clearArray(long double* array, int size)
{
	for (int i = 0; i < size; ++i)
		array[i] = 0;
}

void fillBound(long double *array)
{
	for(int i=0; i<NX; i++)
	{
		array[i] = BoundaryValue(x(i),-B);
		array[NX*(NY-1)+i] = BoundaryValue(x(i),B);
	}
	for(int j=0; j<NY; j++)
	{
		array[NX*j] = BoundaryValue(-A,y(j));
		array[NX*j+(NX-1)] = BoundaryValue(A,y(j));
	}	
}

int main(int argc, char* argv[])
{
    NX = 100; NY = 100;
    hx = 2 * A / (NX-1);
    hy = 2 * B / (NY-1);

	long double* SolVect = (long double *)malloc(NX*NY*sizeof(long double));
	long double* RHS_Vect = (long double *)malloc(NX*NY*sizeof(long double));			    // the right hand side of Puasson equation.
	long double* TmpSolVect = (long double *)malloc(NX*NY*sizeof(long double));
	long double* BasisVect = (long double *)malloc(NX*NY*sizeof(long double));
	long double* TmpBasisVect = (long double *)malloc(NX*NY*sizeof(long double));
	long double sp, alpha, tau, NewValue, err=1e6;	// auxiliary values

// Initialization of Arrays
    RightPart(RHS_Vect);

    clearArray(SolVect, NX*NY);
    fillBound(SolVect);
    fillBound(TmpSolVect);
    clearArray(BasisVect, NX*NY);
    clearArray(TmpBasisVect, NX*NY);

	int rc;
	if (rc = MPI_Init(&argc, &argv)) {
		printf("Error\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		return 0;
	}

//	MPI_Comm comm;
//	MPI_Cart_comm_create(&comm);

	int numprocs, myid;
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	long double startWTime = 0.0;
	if (myid == 0)
		startWTime = MPI_Wtime();

	int* size = (int*)malloc(numprocs * sizeof(int));
	int* displs = (int*)malloc(numprocs * sizeof(int));
	for (int i = 0; i < numprocs; ++i) {
		size[i] = i < (numprocs - 1)
			? (NX*NY) / numprocs
			: (NX*NY - (((NX*NY) / numprocs) * (numprocs - 1)));
		displs[i] = i == 0 ? 0 : (displs[i - 1] + size[i - 1]);
	}

	printf("Start: %d %d\n", myid, numprocs);
	printf("Start: %d %d %d %lf %lf\n", myid, NX, NY, (double)hx, (double)hy);
	fflush(stdout);

	long double tmp;

// Working
// The value of product (r(k),r(k)) is calculating ...
	sp = 0.0;
	tmp = 0.0;
	ParallelCycle({tmp += ResVect(i,j)*ResVect(i,j)*hx*hy;});

	MPI_Reduce(&tmp, &sp, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	tau = sp;

// The value of product sp = (Ar(k),r(k)) is calculating ...
	sp = 0.0;
	tmp = 0.0;
	ParallelCycle({tmp += LeftPartResVect(i,j)*ResVect(i,j)*hx*hy;});

	MPI_Reduce(&tmp, &sp, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	tau = tau/sp;

// The x(k+1) is calculating ...
	err = 0.0;
	MPI_Bcast(&tau, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	tmp = 0.0;
	ParallelCycle({
		NewValue = SolVect[NX*j+i]-tau*ResVect(i,j);
		tmp += (NewValue-SolVect[NX*j+i])*(NewValue-SolVect[NX*j+i]);
		TmpSolVect[NX*j+i] = NewValue;
	});

	MPI_Reduce(&tmp, &err, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	err = sqrt(err);
	MPI_Bcast(&err, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Allgatherv(TmpSolVect + displs[myid], size[myid], MPI_LONG_DOUBLE,
		SolVect, size, displs, MPI_LONG_DOUBLE, MPI_COMM_WORLD);

	for (int i = 1; i < NX-1; ++i)
		for (int j = 1; j < NY-1; ++j)
		BasisVect[NX*j+i] = ResVect(i,j);

// the end of steep descent iteration.
	while (err > 1e-4)
	{
	// The value of product (Ar(k),g(k-1)) is calculating ...
		alpha = 0.0;
		tmp = 0.0;
		ParallelCycle({tmp += LeftPartResVect(i,j)*BasisVect[NX*j+i]*hx*hy;});

		MPI_Reduce(&tmp, &alpha, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		alpha = alpha/sp;
		MPI_Bcast(&alpha, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

	// The new basis vector g(k) is being calculated ...
		ParallelCycle({TmpBasisVect[NX*j+i] = ResVect(i,j)-alpha*BasisVect[NX*j+i];});
		MPI_Allgatherv(TmpBasisVect + displs[myid], size[myid], MPI_LONG_DOUBLE,
			BasisVect, size, displs, MPI_LONG_DOUBLE, MPI_COMM_WORLD);

	// The value of product (r(k),g(k)) is being calculated ...
		tau = 0.0;
		tmp = 0.0;
		ParallelCycle({tmp += ResVect(i,j)*BasisVect[NX*j+i]*hx*hy;});
		MPI_Reduce(&tmp, &tau, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
	// The value of product sp = (Ag(k),g(k)) is being calculated ...
		sp = 0.0;
		ParallelCycle({tmp += LeftPart(BasisVect,i,j)*BasisVect[NX*j+i]*hx*hy;});
		MPI_Reduce(&tmp, &sp, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		tau = tau/sp;
		MPI_Bcast(&tau, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);

	// The x(k+1) is being calculated ...
		err = 0.0;
		tmp = 0.0;
		ParallelCycle({
			NewValue = SolVect[NX*j+i]-tau*BasisVect[NX*j+i];
			tmp += (NewValue-SolVect[NX*j+i])*(NewValue-SolVect[NX*j+i]);
			TmpSolVect[NX*j+i] = NewValue;
		});
		MPI_Reduce(&tmp, &err, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		err = sqrt(err);
		MPI_Bcast(&err, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Allgatherv(TmpSolVect + displs[myid], size[myid], MPI_LONG_DOUBLE,
			SolVect, size, displs, MPI_LONG_DOUBLE, MPI_COMM_WORLD);
	}
// the end of CGM iterations.

    if (myid == 0) {
    	printf("%lf\n", (double)err);
	    for (int i = 1; i < NX-1; ++i) {
	    	for (int j = 1; j < NY-1; ++j)
	    		printf("%d", fabs(RHS_Vect[NX*j+i] - LeftPart(SolVect, i, j)) < 1e-3);
		    printf("\n");
		}
	    printf("\n");
/*
	    for (int i = 0; i < NX; ++i) {
	    	for (int j = 0; j < NY; ++j)
	    		printf("%d", fabs(SolVect[NX*j+i] - BoundaryValue(x(i), y(j))) < 1e-4);
		    printf("\n");
	    }*/

	    printf("Time = %lf\n", (double)(MPI_Wtime() - startWTime));
	}

	free(SolVect); free(RHS_Vect);

	MPI_Finalize();
	return 0;
}
