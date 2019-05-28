#include "Quantum.h"

//This performs a wigner transformation from a wavefunction into phase space distribution in atomic unit, dp needs to be specified
void WignerTransform_au(gsl_matrix_complex * dest, gsl_vector_complex * source, int grades, double dx, double dp)
{
    gsl_complex P;
    gsl_complex temp;
    gsl_complex power;

    int i;
    int j;
    int k;
    
    //Start converting. In this code, i stands for the x grid points, j for the p grid points, and k for the sum step cursor representing integration \int dy \psi^\dagger(x-y) \psi(x+y) exp(2py/h).
    for(i=0;i<2*grades+1;i++)
    {   
        for(j=0;j<2*grades+1;j++)
        {
            GSL_SET_COMPLEX(&P,0,0);
            for(k=-grades;k<=grades;k++)
            {
                if(i-k>=0 && i+k<2*grades+1)
                {
                    GSL_SET_COMPLEX(&power,0,2.0*(j-grades)*dp*dx*k);
                    //Here is the integrated part, namely dy \psi^\dagger(x-y) \psi(x+y) exp(2py/h). The operation steps are:
                    //dx *( ( (\psi(x-y))^dagger * \psi(x+y) ) * exp{2py} )
                    temp = gsl_complex_mul_real(gsl_complex_mul(gsl_complex_mul(gsl_complex_conjugate(gsl_vector_complex_get(source,i-k)),gsl_vector_complex_get(source,i+k)),gsl_complex_exp(power)),dx);
                    P = gsl_complex_add(P,gsl_complex_div_real(temp,M_PI));
                }
            }

            gsl_matrix_complex_set(dest,j,i,P);
        }
    }
}


//This performs a wigner transformation from a wavefunction into phase space distribution in SI unit, dp needs to be specified
void WignerTransform_SI(gsl_matrix_complex * dest, gsl_vector_complex * source, int grades, double dx, double dp)
{
    gsl_complex P;
    gsl_complex temp;
    gsl_complex power;

    int i;
    int j;
    int k;
    
    double h = 6.626070040818181818e-34;
    //Start converting. In this code, i stands for the x grid points, j for the p grid points, and k for the sum step cursor representing integration \int dy \psi^\dagger(x-y) \psi(x+y) exp(2py/h).
    for(i=0;i<2*grades+1;i++)
    {   
        for(j=0;j<2*grades+1;j++)
        {
            GSL_SET_COMPLEX(&P,0,0);
            for(k=-grades;k<=grades;k++)
            {
                if(i-k>=0 && i+k<2*grades+1)
                {
                    GSL_SET_COMPLEX(&power,0,2.0*(j-grades)*dp*dx*k/h);
                    //Here is the integrated part, namely dy \psi^\dagger(x-y) \psi(x+y) exp(2py/h). The operation steps are:
                    //dx *( ( (\psi(x-y))^dagger * \psi(x+y) ) * exp{2py/h} )
                    temp = gsl_complex_mul_real(gsl_complex_mul(gsl_complex_mul(gsl_complex_conjugate(gsl_vector_complex_get(source,i-k)),gsl_vector_complex_get(source,i+k)),gsl_complex_exp(power)),dx);
                    P = gsl_complex_add(P,gsl_complex_div_real(temp,M_PI*h));
                }
            }

            gsl_matrix_complex_set(dest,j,i,P);
        }
    }
}
