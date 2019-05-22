#include "../versor.h"
#include "../input.h"
#include "../state.h"
#include "../hamiltonian.h"
#include <stdio.h>
#include <math.h>

void choose_test_versor(basis b, parameters p);
void present_braH(versor psi0, parameters pars, basis b);

int main(int argc, char *argv[])
{
    int n1=20, 
    n3 = 4,
    n5 = 4, 
    j1 = 2,
    j2 = 2;

    if(argc > 1)
    {
        if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'b')
        {
            printf("Define basis truncation:\n");
            scanf("%d %d %d %d %d", &n1, &n3, &n5, &j1, &j2);
            printf("\n");
        }
    }

    basis b = {n1, n3, n5, j1, j2};
    check_basis_truncation(b);

    printf(" Basis truncation: |%d,%d,%d;%d,%d>\n", 
        b.n1, b.n3, b.n5, b.j1, b.j2);

    // SrYb^+
    // $ \omega_1 = \sqrt{3} 0.16MHz \approx 0.28 MHz $
    // $ \omega_3 \approx 1.39 MHz $
    // $ \omega_3 / \omega_1 \approx 5$
    float omega_rho = 1.4f;
    float omega_z = 0.16f;
    float omega_1 = sqrt(3.0)*omega_z;
    float omega_3 = sqrt(omega_rho*omega_rho - omega_z*omega_z);

    hamiltonian activate = {0, 0, 0, 0, 0};
    activate.normal_modes = 1;
    
    const parameters pars = 
    {261.f, 1.f, 4.75f, 503.7f, omega_1, omega_3, activate};

    print_active_terms_of_Hamiltonian(pars);

    while(1)
        choose_test_versor(b, pars);

    printf("\n");
    return 0;
}

void choose_test_versor(basis b, parameters p)
{
    int versor_idx;
    int basis_size = get_basis_size(b);

    printf("Give an index of the test versor: ");
    scanf("%d", &versor_idx);
    while( !(versor_idx < basis_size && versor_idx >= 0) )
    {
        printf("Incorrect index.\n"
        "It must be less than %d and no less than 0.\n", basis_size);
        printf("Give an index of the test versor: ");
        scanf("%d", &versor_idx);
    }
    versor psi0 = get_versor_from_index(versor_idx, b);
    
    printf("The input versor is: ");
    show_versor(psi0);
    printf("\t\t");
    printf("index = %d (sanity check)", get_index_from_versor(psi0, b));
    printf("\n");

    // acts from left with Hamiltonian on the bra `psi0`
    // prints the result
    present_braH(psi0, p, b);

    printf("\n");
}


void present_braH(versor psi0, parameters pars, basis b)
{
    // initialisation of the `state` structure 
    // intially it contains only a singe versor, the test versor.
    state state0;
    state_init(&state0);
    state_add(&state0, psi0, (fcomplex){1.0f, 0.0f});

    // apply the Hamiltonian and save the output bra state
    state sps = bra_H(&state0, pars);

    // present the resut state
    printf("\n <psi0|H = \n");
    versor loop_versor;
    fcomplex loop_amplitude;
    for(int l=0; l < sps.length; l++)
    {
            loop_versor = state_get_versor(&sps, l);
            loop_amplitude = state_get_amplitude(&sps, l);
            printf("\t(%8.2f,%8.2f)", loop_amplitude.re, loop_amplitude.im);
            show_bra_versor(loop_versor);
            printf("+\t\tidx = %7d\n", get_index_from_versor(loop_versor, b));
    }

    // present only the important part
    printf("\n But the physical versors (within truncation) give\n");
    printf(" <psi0|H = \n");
    for(int l=0; l < sps.length; l++)
    {
            loop_versor = state_get_versor(&sps, l);
            loop_amplitude = state_get_amplitude(&sps, l);
            if(!valid_versor(loop_versor, b))
                    continue;

            printf("\t(%8.2f,%8.2f)", loop_amplitude.re, loop_amplitude.im);
            show_bra_versor(loop_versor);
            printf("+\t\tidx = %7d\n", get_index_from_versor(loop_versor, b));
    }

    // free memory
    state_free(&sps);
    state_free(&state0);
}