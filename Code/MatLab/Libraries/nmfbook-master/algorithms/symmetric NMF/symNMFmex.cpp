#include "math.h"
#include "mex.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

struct mat {
	double *values;
	int *row_index;
	int *col_begin;
	int n;
	int m;
};

// This function computes d^(1/3)
// for either positive or negative value d.
double cubicRoot(double d)
{
  if(d<0.0)
      return -cubicRoot(-d);
  else
      return pow(d,1.0/3.0);
}

// This function solves the following problem:
// min_{x>=0} x^3+ax+b
double BestPolynomialRoot(double a, double b)
{
    double x=0, y=0;
    double a3=4*pow(a,3), b2=27*pow(b,2);
    double delta = a3+b2;
    
    if(delta<=0) // 3 distinct real roots or 1 real multiple solution
    {
        double r3  = 2*sqrt(-a/3);
        double th3 = atan2(sqrt(-delta/108),-b/2)/3;
        double ymax=0, xopt=0;
        for(int k=0;k<=4;k=k+2)
        {
            x = r3*cos(th3+((k*3.14159265)/3));
            if(x>=0)
            {
                y=pow(x,4)/4+a*pow(x,2)/2+b*x;
                if(y<ymax)
                    {ymax=y; xopt=x;}
            }
        }
        return xopt;
    }
    else // 1 real root and two complex
    {
        double z = sqrt(delta/27);
        x = cubicRoot(0.5*(-b+z))+cubicRoot(0.5*(-b-z));
        y = pow(x,4)/4+a*pow(x,2)/2+b*x;
        if(y<0 && x>=0)
            return x;
        else
            return 0;
    }
}

// Precomputations before the iterations over the variables of H:
// nL[i] is norm(H(i,:),'fro')^2
// nC[j] is norm(H(:,j),'fro')^2
// HH is H^TH
void precomputations(int n, int r, double *H, double *nL, double *nC, double *HH)
{
    for(int i=0; i<n; i++)
    {
        nL[i]=0;
        for(int j=0; j<r; j++)
            nL[i] += H[j*n+i]*H[j*n+i];
    }
    
    for(int j=0; j<r; j++)
    {
        nC[j]=0;
        for(int i=0; i<n; i++)
            nC[j] += H[j*n+i]*H[j*n+i];
    }
    
    for(int j1=0; j1<r; j1++)
        for(int j2=j1; j2<r; j2++)
        {
            HH[j2*r+j1]=0;
            for(int i=0; i<n; i++)
                HH[j2*r+j1] += H[j1*n+i]*H[j2*n+i];
        }
}

// In the algorithm, we need the values A(i,i)
// This function stores these values for sparse matrix A
void diagA(int n, struct mat *A, double *dA)
{
    bool flag=false;
    double valdiag=-1;
    for(int i=0; i<n; i++)
    {
        int begin=A->col_begin[i], end=A->col_begin[i+1];
        flag=false;
        for(int q=begin;q<end;q++)
        {
            if(i==A->row_index[q])
                {flag=true; valdiag=A->values[q];}
            if(A->row_index[q]>=i)
                break;
        }
        dA[i]=flag?valdiag:0;           
    }
}

// This function checks if the stopping criteria (maxiter or timelimit) are satisfied
bool stoppingcondition(double *e, double *t, int iter, double timelimit, int maxiter)
{
    if(t[iter]>=timelimit)
    {
        for(int i=iter+1;i<maxiter;i++)
        {
            t[i]=-1;
            e[i]=-1;
        }
        return true;
    }
    return false;
}

// Main algorithm for dense matrices
void iterations(double *A, double *H, int n, int r, int maxiter, double timelimit, double *e, double *t, int shuffle_columns)
{
    clock_t t0  = clock();
       
    double *nL  = (double *)malloc(sizeof(double)*(n));
    double *nC  = (double *)malloc(sizeof(double)*(r));
    double *HH  = (double *)malloc(sizeof(double)*(r*r));
    double *HHH = (double *)malloc(sizeof(double)*(n));
    double *AH  = (double *)malloc(sizeof(double)*(n));
    int *icol   = (int *)malloc(sizeof(int)*(r));

    precomputations(n,r,H,nL,nC,HH);

    double *Hkn,*HHkr,*Ain,fdecrease=0,hold=0,hnew=0,s1=0,s2=0,a=0,b=0;
    int k=0,indt=0,indal=0;
    
    for(int iter=0;iter<maxiter;iter++)
    {
        fdecrease=0;
        for(int i=0;i<r;i++)
            icol[i]=i;
        //Shuffle columns
        if(shuffle_columns==1)
            for(int i=r-1;i>=1;i--)
            {
                indal       = rand()%(i+1);
                indt        = icol[indal];
                icol[indal] = icol[i];
                icol[i]     = indt;
            }
        // Iteration over the columns
        for(int kk=0; kk<r; kk++)
        {
            k=icol[kk]; //columns are shuffled
            
            Hkn        = &(H[k*n]);
            HHkr       = &(HH[k*r]);
            
            // Precomputations of AH[i] values
            for(int i=0; i<n; i++)
            { 
                Ain   = &(A[i*n]);
                AH[i] = 0;
                for(int j=0;j<n;j++)
                    AH[i] += Hkn[j]*Ain[j];
            }
            
            // Iteration over the n variables of column k
            for(int i=0; i<n; i++)
            {
                HHH[i]=0;
                for(int j=0;j<r;j++)
                    HHH[i]+=j<=k?(HHkr[j]*H[j*n+i]):(HH[j*r+k]*H[j*n+i]);
                
                hold = Hkn[i];
                a    = nC[k]+nL[i]-A[i*n+i]-2*pow(hold,2);
                b    = HHH[i]-AH[i]-Hkn[i]*a-pow(hold,3);
                hnew = BestPolynomialRoot(a,b);
                s1   = hold-hnew;
                
                // Update
                if(s1!=0)
                {
                    Hkn[i] = hnew;
                    s2     = pow(hnew,2) - pow(hold,2);
                    nC[k] += s2;
                    nL[i] += s2;
                    Ain = &(A[i*n]);
                    for(int j=i;j<n;j++)
                        AH[j] -= Ain[j]*s1;
                    for(int j=0;j<k;j++)
                        HHkr[j] -= s1*H[j*n+i];
                    HH[k*r+k] += s2;
                    for(int j=k+1;j<r;j++) 
                        HH[j*r+k] -= s1*H[j*n+i];
                    fdecrease += 4*b*s1-2*a*s2 + pow(hold,4)-pow(hnew,4);
                }
            }
        }
        e[iter]=fdecrease;
        t[iter]=(double) (clock()-t0)/CLOCKS_PER_SEC;
        if(stoppingcondition(e,t,iter,timelimit,maxiter))
            break;
    }

    free(nL);
    free(nC);
    free(HH);
    free(HHH);
    free(AH);
    free(icol);
}

// Main algorithm for sparse matrices
void iterations_sp(struct mat *A, double *H, int n, int r, int maxiter, double timelimit, double *e, double *t, int shuffle_columns)
{
	clock_t t0  = clock();
       
    double *nL  = (double *)malloc(sizeof(double)*(n));
    double *nC  = (double *)malloc(sizeof(double)*(r));
    double *HH  = (double *)malloc(sizeof(double)*(r*r));
    double *HHH = (double *)malloc(sizeof(double)*(n));
    double *AH  = (double *)malloc(sizeof(double)*(n));
    double *dA  = (double *)malloc(sizeof(double)*(n));
    int *icol   = (int *)malloc(sizeof(int)*(r));
    
    precomputations(n,r,H,nL,nC,HH);
    diagA(n,A,dA);
    
    double *Hkn, *HHkr,fdecrease=0,hold=0,hnew=0,s1=0,s2=0,a=0,b=0;
    int k=0,indt=0,indal=0;
    
    for(int iter=0;iter<maxiter;iter++)
    {
        fdecrease=0;
        for(int i=0;i<r;i++)
            icol[i]=i;
        //Shuffle columns
        if(shuffle_columns==1)
            for(int i=r-1;i>=1;i--)
            {
                indal       = rand()%(i+1);
                indt        = icol[indal];
                icol[indal] = icol[i];
                icol[i]     = indt;
            }
        // Iteration over the columns
        for(int kk=0; kk<r; kk++)
        {
            k    = icol[kk];
            Hkn  = &(H[k*n]);
            HHkr = &(HH[k*r]);
            
            // Iteration over the n variables of column k
            for(int i=0; i<n; i++)
            {
                double HHHi=0;
                for(int j=0;j<r;j++)
                    HHHi+=j<=k?(HHkr[j]*H[j*n+i]):(HH[j*r+k]*H[j*n+i]);
                double AHi=0;
                for(int q=A->col_begin[i];q<A->col_begin[i+1];q++)
                    AHi+=Hkn[A->row_index[q]]*(A->values[q]);
                
                hold = Hkn[i];
                a    = nC[k]+nL[i]-dA[i]-2*pow(hold,2);
                b    = HHHi-AHi-Hkn[i]*a-pow(hold,3);
                hnew = BestPolynomialRoot(a,b);
                s1   = hold-hnew;
                
                // Update
                if(s1!=0)
                {
                    Hkn[i] = hnew;
                    s2     = pow(hnew,2) - pow(hold,2);
                    nC[k] += s2;
                    nL[i] += s2;
                    for(int j=0;j<k;j++)
                        HHkr[j] -= s1*H[j*n+i];
                    HH[k*r+k] += s2;
                    for(int j=k+1;j<r;j++) 
                        HH[j*r+k] -= s1*H[j*n+i];
                    fdecrease += 4*b*s1-2*a*s2 + pow(hold,4)-pow(hnew,4);
                }
            }
        }
        e[iter]=fdecrease;
        t[iter]=(double) (clock()-t0)/CLOCKS_PER_SEC;
        if(stoppingcondition(e,t,iter,timelimit,maxiter))
            break;
    }
    free(nL);
    free(nC);
    free(HH);
    free(HHH);
    free(AH);
    free(dA);
    free(icol);
}

// Reading a sparse matrix
struct mat *read_matrix(const mxArray *matrix)
{
	int nnz = mxGetNzmax(matrix);
	
    double *elements;
	elements = mxGetPr(matrix);
    
	struct mat *V = (struct mat *)malloc(sizeof(struct mat));
	V->n = mxGetM(matrix);
	V->m = mxGetN(matrix);
    
	V->values = mxGetPr(matrix);
	V->row_index = (int *)malloc(sizeof(int)*nnz);
	V->col_begin = (int *)malloc(sizeof(int)*(V->m+1));
    
    mwIndex *ir,*jc;
	ir = mxGetIr(matrix);
	jc = mxGetJc(matrix);
	for (int i=0; i<nnz; i++)
        V->row_index[i] = (int)ir[i];
	for (int i=0; i<=V->m; i++)
        V->col_begin[i] = (int)jc[i];
    return V;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *values;
    double timelimit;
    int maxiter, n, r, shuffle_columns;
    
    srand(time(NULL));
    
    // Input argument #1
    double *H0 = mxGetPr(prhs[1]);
    n = mxGetM(prhs[1]);
    r = mxGetN(prhs[1]);
    
    // Input argument #2
    values  = mxGetPr(prhs[2]);
	maxiter = values[0];
    
    // Input argument #3
	values    = mxGetPr(prhs[3]);
	timelimit = values[0];
                
    // Input argument #4
    values          = mxGetPr(prhs[4]);
    shuffle_columns = values[0];

    // Output argument #0
    plhs[0]   = mxCreateDoubleMatrix(n,r,mxREAL);
	double *H = mxGetPr(plhs[0]);
    for(int i=0;i<n*r;i++)
        H[i] = H0[i];

    // Output argument #1
    plhs[1]   = mxCreateDoubleMatrix(maxiter,1,mxREAL);
    double *e = mxGetPr(plhs[1]);
    
    // Output argument #2
    plhs[2]   = mxCreateDoubleMatrix(maxiter,1,mxREAL);
    double *t = mxGetPr(plhs[2]);

    // Algorithm and Input argument #0    
    if (!mxIsSparse(prhs[0]))
	{
		// The matrix A is dense
		double *A = mxGetPr(prhs[0]);
        iterations(A,H,n,r,maxiter,timelimit,e,t,shuffle_columns);
	}
	else
    {
		// The matrix A is sparse
        struct mat *Asp = read_matrix(prhs[0]);
        iterations_sp(Asp,H,n,r,maxiter,timelimit,e,t,shuffle_columns);
        free(Asp->row_index);
        free(Asp->col_begin);
        free(Asp);
    }
}