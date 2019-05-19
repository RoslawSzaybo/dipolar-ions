#include "../versor.h"
#include "../input.h"
#include "../state.h"
#include "../hamiltonian.h"
#include <stdio.h>

void choose_test_versor(basis b, parameters p);
void present_braH(versor psi0, parameters pars, basis b);

int main(int argc, char *argv[])
{
    int n1=20, n35=4, j=2;

    if(argc > 1)
    {
        if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'b')
        {
            printf("Define basis truncation:\n");
            scanf("%d %d %d", &n1, &n35, &j);
            printf("\n");
        }
    }

    basis b = {n1, n35, n35, j, j};
    check_basis_truncation(b);

    printf(" Basis truncation: |%d,%d,%d;%d,%d>\n", 
        b.n1, b.n3, b.n5, b.j1, b.j2);

    // SrYb^+
    const parameters pars = {261.f, 1.f, 4.75f, 503.7f, 1.4f, 0.16f};

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

    state sps = bra_H(&state0, pars);
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

    state_free(&sps);
    state_free(&state0);
}