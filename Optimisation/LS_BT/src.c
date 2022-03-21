#include<stdio.h>
#include<math.h>
#include"./../RosenBrock/test_f.h"

int main( int argc, char* argv[] )
{
    const double a = 1.;
    const double b = 100.;
    double x = -3;
    double y = -4;

    double tx,ty;

    int iter = 0;

    double cost;
    double dcost[2];
    double trial_dcost[2];
    double LC,RC;
   
    double norm_sq;
    double t;
    const double beta = 0.84;
    double trial_cost;
    double backtracking_std;

    for(;;)
    {
        cost = RosenBrock(a,b,x,y);
        D_RosenBrock(a,b,x,y,&dcost[0]); 
        norm_sq  = dcost[0]*dcost[0] + dcost[1]*dcost[1];

        if( sqrt(norm_sq)/2. < 10.E-8 )
	{ printf("Termination Condition Achieved!\n");
          printf("%.4s\t%8.5d\t%12.6lf\t%12.6lf%12.4e\n","CYC",iter,x,y,sqrt(norm_sq)/2.);
	  break;
	}
        
        t = 1.;
        // sub solver
        for(int i=0;i<1000;i++)
        {
            tx = x - t*dcost[0];
            ty = y - t*dcost[1];
            trial_cost = RosenBrock(a,b,tx,ty);
            backtracking_std = cost - 0.5*t*norm_sq;

            if( trial_cost < backtracking_std )
            {   
	/*
		// trial_dcost
		D_RosenBrock(a,b,tx,ty,&trial_dcost[0]);
		LC = (dcost[0]*trial_dcost[0] + dcost[1]*trial_dcost[1]);
		RC = 0.9*norm_sq; // dcost[0]*dcost[0] + dcost[1]*dcost[1];
 		if ( LC <= RC )
		{	
			x = tx;
			y = ty;
			break;
		}
	*/
		x = tx;
		y = ty;
		break;
            }
            else
            { t = beta*t;
            }
        }
        iter++;

        printf("%.4s\t%8.5d\t%12.6lf\t%12.6lf%12.6lf\n","CYC",iter,x,y,sqrt(norm_sq));
    }

    return 0;
}


/* Possible Output

CYC	   03798	    1.000000	    1.000000    0.000001
CYC	   03799	    1.000000	    1.000000    0.000000
CYC	   03800	    1.000000	    1.000000    0.000001
CYC	   03801	    1.000000	    1.000000    0.000000
CYC	   03802	    1.000000	    1.000000    0.000001
CYC	   03803	    1.000000	    1.000000    0.000000
CYC	   03804	    1.000000	    1.000000    0.000001
CYC	   03805	    1.000000	    1.000000    0.000000
CYC	   03806	    1.000000	    1.000000    0.000001
CYC	   03807	    1.000000	    1.000000    0.000000
CYC	   03808	    1.000000	    1.000000    0.000001
CYC	   03809	    1.000000	    1.000000    0.000000
CYC	   03810	    1.000000	    1.000000    0.000001
CYC	   03811	    1.000000	    1.000000    0.000000
CYC	   03812	    1.000000	    1.000000    0.000001
CYC	   03813	    1.000000	    1.000000    0.000000
CYC	   03814	    1.000000	    1.000000    0.000001
CYC	   03815	    1.000000	    1.000000    0.000000
CYC	   03816	    1.000000	    1.000000    0.000001
CYC	   03817	    1.000000	    1.000000    0.000000
CYC	   03818	    1.000000	    1.000000    0.000001
CYC	   03819	    1.000000	    1.000000    0.000000
CYC	   03820	    1.000000	    1.000000    0.000001
CYC	   03821	    1.000000	    1.000000    0.000000
CYC	   03822	    1.000000	    1.000000    0.000001
CYC	   03823	    1.000000	    1.000000    0.000000
CYC	   03824	    1.000000	    1.000000    0.000001
CYC	   03825	    1.000000	    1.000000    0.000000
CYC	   03826	    1.000000	    1.000000    0.000001
CYC	   03827	    1.000000	    1.000000    0.000000
CYC	   03828	    1.000000	    1.000000    0.000001
Termination Condition Achieved!
CYC	   03828	    1.000000	    1.000000  9.7729e-08

*/
