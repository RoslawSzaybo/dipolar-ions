#include <stdio.h>
#include "input.h"

int test_input(int n1, int n3, int n5, int j1, int j2, 
                float omega_rho, float omega_z)
{
    const int jmax = 2; // maximal implemented j
    int is_OK = 1;
    if ( n1<1 || n3<1 || n5<1 )
    {
            printf(" Warning.\n");
            printf(" There has to be at least one mode of each allowed.\n\n");
            is_OK = 0;
    }
    
    if (j1 < 0 || j1 > jmax || j2 < 0 || j2 > jmax)
    {
            printf(" Warning.\n");
            printf(" Jmax has to be no smaller than 0"
            " and no larger than %d.\n\n", jmax);
            is_OK =0;
    }

    // this two will usually differ by an order of magnitude
    if ( omega_z > omega_rho)
    {
            printf(" Omega_z > omega_rho.\n");
            is_OK = 0;
    }

    return is_OK;
}