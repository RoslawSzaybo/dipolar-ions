#include "../versor.h"
#include "../input.h"
#include "../state.h"
#include "../hamiltonian.h"
#include <stdio.h>
#include <math.h>

void choose_test_versor_idx(basis b, parameters p);
void choose_test_versor_quantum_numbers(basis b, parameters p);
void present_braH(versor psi0, parameters pars, basis b);
void present_versor_contribution(dcomplex amp, versor bra, basis b);

int main(int argc, char *argv[])
{
    int n1=20, 
    n3 = 20,
    n5 = 20, 
    j1 = 2,
    j2 = 2;

    if (argc > 1)
    {
        if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'b')
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
    double omega_rho = 1.4;
    double omega_z = 0.16;
    double omega_1 = sqrt(3.0)*omega_z;
    double omega_3 = sqrt(omega_rho*omega_rho - omega_z*omega_z);

    hamiltonian activate;
    activate.normal_modes = 0;
    activate.T_rot = 1;
    activate.Vqd_zeroth = 0;
    activate.Vqd_first = 0;
    activate.Vdd_zeroth = 0;
    
    /* MgH^+
    const parameters pars = 
    {25.0, 1.0, 3.0, 1.9e5, omega_1, omega_3, activate};
    */

    const parameters pars = 
    // SrYb^+
    {261.0, 1.0, 4.745, 503.7, omega_1, omega_3, activate};

    print_system_parameters(pars);
    print_active_terms_of_Hamiltonian(pars);

    while(1)
        choose_test_versor_quantum_numbers(b, pars);

    printf("\n");
    return 0;
}

void choose_test_versor_idx(basis b, parameters p)
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


/*
Select the initial state (one versor) by exclicitly listing 
the quantum numbers of the state.
*/
void choose_test_versor_quantum_numbers(basis b, parameters p)
{
    versor psi0;

    printf("Please give the numbers of the test state: ");
    scanf("%d %d %d %d %d %d %d", &psi0.n1, &psi0.n3, &psi0.n5,
    &psi0.j1, &psi0.m1, &psi0.j2, &psi0.m2);

    while( !valid_versor(psi0, b) )
    {
        printf("Incorrect versor.\n");
        printf("Please give the numbers of the test state: ");
        scanf("%d %d %d %d %d %d %d", &psi0.n1, &psi0.n3, &psi0.n5,
        &psi0.j1, &psi0.m1, &psi0.j2, &psi0.m2);
    }
    
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
    state_add(&state0, psi0, (dcomplex){1.0, 0.0});

    // apply the Hamiltonian and save the output bra state
    state sps = bra_H(&state0, pars);

    // present the resut state
    printf("\n <psi0|H = \n");
    versor loop_versor;
    dcomplex loop_amplitude;
    for(int l=0; l < sps.length; l++)
    {
        loop_versor = state_get_versor(&sps, l);
        loop_amplitude = state_get_amplitude(&sps, l);
        present_versor_contribution(loop_amplitude, loop_versor, b);
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

        present_versor_contribution(loop_amplitude, loop_versor, b);
    }

    // free memory
    state_free(&sps);
    state_free(&state0);
}

void present_versor_contribution(dcomplex amp, versor bra, basis b)
{
    printf("\t(%14.11f,%14.11f)", amp.re, amp.im);
    show_bra_versor(bra);
    printf("+\t\t");
    printf("idx = %7d", get_index_from_versor(bra, b));
    printf("\n");
}
