#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"./../RosenBrock/test_f.h"

//int gsl_blas_dgemv(CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix *A, const gsl_vector *x, double beta, gsl_vector *y)

double	get_rosenbrock( const double a, const double b, gsl_vector* v )
{
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	return RosenBrock(a,b,x,y);
}
void	get_drosenbrock( const double a, const double b, gsl_vector* v, gsl_vector* dv )
{
	double x = gsl_vector_get(v,0);
	double y = gsl_vector_get(v,1);
	double w[2];
	D_RosenBrock(a,b,x,y,&w[0]);
	
	gsl_vector_set(dv,0,w[0]);
	gsl_vector_set(dv,1,w[1]);
	return;
}

double	bfgs_ls_subproblem( const double a, const double b, const double c, const gsl_vector* v, const gsl_vector* p )
{
	double ret = 0;
	gsl_vector* vt = gsl_vector_calloc(2);

	double t = 1.; double beta = 0.84;
	double tr_std;
	double bt_std;

	for(int i=0;i<10000;i++)
	{
	//int gsl_blas_daxpy(double alpha, const gsl_vector *x, gsl_vector *y)
		gsl_vector_set_zero(vt);
		gsl_blas_daxpy(1.,v,vt);	// vt += v
		gsl_blas_daxpy( t,p,vt);	// vt += t*p;
		
		bt_std = c - t*0.5*pow(gsl_blas_dnrm2(p),2.);
		tr_std = get_rosenbrock(a,b,vt);
	
		if( tr_std < bt_std )
		{
			gsl_vector_free(vt);
			return t;
		}
		else
		{
			t = t*beta;
		}
	}

	gsl_vector_free(vt);
	return -1;
}

gsl_matrix* vwT( const gsl_vector* v, const gsl_vector* w );
void get_bfgs_H( const gsl_vector* s, const gsl_vector* y, gsl_matrix* H );

int main( int argc, char* argv[] )
{
	const double a = 1.;
	const double b = 100.;
	const double x = -3;
	const double y = -4;

	int iter;
	double norm_sq;

	double cost;
	double bfgs_a;	// step-size within LS-subproblem
	gsl_vector* bfgs_v = gsl_vector_calloc(2);
	gsl_vector* bfgs_tv = gsl_vector_calloc(2);
	gsl_vector* bfgs_gv = gsl_vector_calloc(2);
	gsl_vector* bfgs_tgv = gsl_vector_calloc(2);
	gsl_vector* bfgs_y = gsl_vector_calloc(2);
	gsl_vector* bfgs_p = gsl_vector_calloc(2);
	gsl_vector* bfgs_s = gsl_vector_calloc(2);
	gsl_matrix* bfgs_inv_b = gsl_matrix_calloc(2,2);


	// set
	gsl_vector_set(bfgs_v,0,x);
	gsl_vector_set(bfgs_v,1,y);
	gsl_matrix_set(bfgs_inv_b,0,0,1);
	gsl_matrix_set(bfgs_inv_b,1,1,1);



	iter = 0;
	for(;;)
	{

		cost = get_rosenbrock(a,b,bfgs_v);
		if ( iter == 0 ) {	get_drosenbrock(a,b,bfgs_v,bfgs_gv);	}
		
		norm_sq = gsl_blas_dnrm2(bfgs_gv)/(bfgs_v->size);
		if( norm_sq < 10.E-8 )
		{
			printf("Termination Condition Achieved!\n");
			printf("%.4s\t%8.5d\t%12.6lf\t%12.6lf%12.4e\n","CYC",iter,gsl_vector_get(bfgs_v,0),gsl_vector_get(bfgs_v,1),norm_sq);
			break;
		}
		iter++;
		printf("%.4s\t%8.5d\t%12.6lf\t%12.6lf%12.6lf\n","CYC",iter,gsl_vector_get(bfgs_v,0),gsl_vector_get(bfgs_v,1),norm_sq);


		// p = - B_inv*gv
		gsl_vector_set_zero(bfgs_p);
		gsl_blas_dgemv(CblasNoTrans,-1.,bfgs_inv_b,bfgs_gv,1.,bfgs_p);
		
		// LS-subproblem
		bfgs_a = bfgs_ls_subproblem(a,b,cost,bfgs_v,bfgs_p);
		//printf("%lf\n",bfgs_a);
		//break;

		// set 's_(k)'
		gsl_vector_set_zero(bfgs_s);
		gsl_blas_daxpy(bfgs_a,bfgs_p,bfgs_s);
		
		// set 'v_(k+1)'
		gsl_vector_set_zero(bfgs_tv);
		gsl_blas_daxpy(1.,bfgs_v,bfgs_tv);	// tv += v
		gsl_blas_daxpy(1.,bfgs_s,bfgs_tv);	// tv += s

		// get grad f( v_(k+1) )
		get_drosenbrock(a,b,bfgs_tv,bfgs_tgv);
		// set 'y_(k)'
		gsl_vector_set_zero(bfgs_y);
		gsl_blas_daxpy(-1.,bfgs_gv,bfgs_y);	// y = -bfgs_gv	
		gsl_blas_daxpy(1.,bfgs_tgv,bfgs_y);	// y = bfgs_tgv - bfgs-gv

		// update v_(k+1) <- v_(k)
		// update grad v_(k+1) <- grad v_(k)
		gsl_vector_memcpy(bfgs_v,bfgs_tv);
		gsl_vector_memcpy(bfgs_gv,bfgs_tgv);

		// GET INVERSE B
		get_bfgs_H(bfgs_s,bfgs_y,bfgs_inv_b);


//printf("%.4s\t%8.5d\t%12.6lf\t%12.6lf%12.6lf\n","CYC",iter,x,y,sqrt(norm_sq));


		// END INVERSE B

	}
	gsl_vector_free(bfgs_v);
	gsl_vector_free(bfgs_tv);
	gsl_vector_free(bfgs_gv);
	gsl_vector_free(bfgs_tgv);
	gsl_vector_free(bfgs_y);
	gsl_vector_free(bfgs_p);
	gsl_vector_free(bfgs_s);
	gsl_matrix_free(bfgs_inv_b);
	return 0;
}

gsl_matrix* vwT( const gsl_vector* v, const gsl_vector* w )
{
	gsl_matrix* ret = NULL;
	if( v->size == w->size )
	{
		size_t len = v->size;
		ret = gsl_matrix_calloc(len,len);

		for(int i=0;i<len;i++)
		{	for(int j=0;j<len;j++)
			{	gsl_matrix_set(ret,i,j, gsl_vector_get(v,i)*gsl_vector_get(w,j));
			}
		}
		return ret;
	}
	else {
		return ret;
	}
}

void get_bfgs_H( const gsl_vector* s, const gsl_vector* y, gsl_matrix* H )
{
	double c_sy;	gsl_blas_ddot(s,y,&c_sy);

	double c_yHy;
	gsl_vector* Hy = gsl_vector_calloc(H->size1);
	gsl_blas_dgemv(CblasNoTrans,1.,H,y,0.,Hy);
	gsl_blas_ddot(y,Hy,&c_yHy);

	gsl_matrix* ss = vwT( s , s );
	gsl_matrix* Hys = vwT( Hy , s );
	gsl_matrix* sy = vwT( s , y ); 

	gsl_matrix* syH = gsl_matrix_calloc(H->size1,H->size2);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,sy,H,0.,syH);	// C = a * op(A) * op(B) + b * C


	// H
	// gsl_matrix* H

	// const * ss
	gsl_matrix_scale(ss, ((c_sy+c_yHy)/c_sy/c_sy) );
	// const * ( Hys + syH )
	gsl_matrix_add(Hys,syH);
	gsl_matrix_scale(Hys,(-1./c_sy));

	// final calculation
	gsl_matrix_add(H,ss);
	gsl_matrix_add(H,Hys);

	gsl_vector_free(Hy);
	gsl_matrix_free(ss);
	gsl_matrix_free(Hys);
	gsl_matrix_free(sy);
	gsl_matrix_free(syH);

	return;
}


/* Possible Output

CYC	   00069	    1.000257	    1.000534    0.004288
CYC	   00070	    1.000181	    1.000375    0.003011
CYC	   00071	    1.000127	    1.000264    0.002121
CYC	   00072	    1.000089	    1.000186    0.001495
CYC	   00073	    1.000058	    1.000121    0.000969
CYC	   00074	    1.000038	    1.000078    0.000628
CYC	   00075	    1.000024	    1.000051    0.000408
CYC	   00076	    1.000016	    1.000033    0.000264
CYC	   00077	    1.000010	    1.000021    0.000171
CYC	   00078	    1.000007	    1.000014    0.000111
CYC	   00079	    1.000004	    1.000009    0.000072
CYC	   00080	    1.000003	    1.000006    0.000047
CYC	   00081	    1.000002	    1.000004    0.000030
CYC	   00082	    1.000001	    1.000002    0.000020
CYC	   00083	    1.000001	    1.000002    0.000013
CYC	   00084	    1.000000	    1.000001    0.000008
CYC	   00085	    1.000000	    1.000001    0.000005
CYC	   00086	    1.000000	    1.000000    0.000003
CYC	   00087	    1.000000	    1.000000    0.000002
CYC	   00088	    1.000000	    1.000000    0.000001
CYC	   00089	    1.000000	    1.000000    0.000001
CYC	   00090	    1.000000	    1.000000    0.000001
CYC	   00091	    1.000000	    1.000000    0.000000
CYC	   00092	    1.000000	    1.000000    0.000000
CYC	   00093	    1.000000	    1.000000    0.000000
CYC	   00094	    1.000000	    1.000000    0.000000
Termination Condition Achieved!
CYC	   00094	    1.000000	    1.000000  7.0947e-08

*/	
