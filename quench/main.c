#include "../spectrum/versor.h"
#include "../spectrum/input.h"
#include "../spectrum/state.h"
#include "../spectrum/hamiltonian.h" // valid_versor()
#include <stdio.h>
#include <math.h>

versor choose_test_versor_quantum_numbers(basis b, parameters p);
float get_dt();
int get_no_of_steps();
void state_print(state *psi);
void print_limited_state(state *psi);
void print_the_big_four(state *psi);

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
    float omega_rho = 1.4f;
    float omega_z = 0.16f;
    float omega_1 = sqrt(3.0)*omega_z;
    float omega_3 = sqrt(omega_rho*omega_rho - omega_z*omega_z);

    hamiltonian activate;
    activate.normal_modes = 1;
    activate.T_rot = 1;
    activate.Vqd_zeroth = 1;
    activate.Vqd_first = 1;
    activate.Vdd_zeroth = 1;
    
    const parameters pars = 
    // SrYb$^+$
    {261.f, 1.f, 4.745f, 503.7f, omega_1, omega_3, activate};

    print_active_terms_of_Hamiltonian(pars);

    /* 
    versor psi0 = choose_test_versor_quantum_numbers(b, pars);
    float time_step = get_dt();
    */
    versor psi0 = (versor){8, 3, 0, 0, 0, 0, 0};
    float time_step = 100.f;
    int steps = get_no_of_steps();
    float ns = 1e-9f;
    float dt = time_step*ns;
    float N2 = 0.f;
    float threshold = 1e-20f;

    state bra;
    state braH;
    state_init(&bra);
    state_init(&braH);
    state_add(&bra, psi0, (fcomplex){1.f, 0.f});

    printf("# Attention please!\n");
    printf("#  All states below are if fact <psi|,"
    " so if you want |ket> complex conjugate amplitudes.\n");

    printf("t = %f ns\n", 0.f);
    state_sort(&bra);
    N2 = state_normalisation(&bra);
    printf("N^2 = %f\n", N2);
    state_times_float(&bra, 1./sqrt(N2));
    state_print(&bra);

    int i=1;
    for (; i<steps; i++)
    {
        braH = bra_H(&bra, pars);
        // i\hbar \ket{psi} = H \ket{psi}
        // bra{psi} = \frac{i}{\hbar}\bra{psi}H
        // H is in unit \hbar \omega_1 so you have to multiply it 
        // \bra{psi} = i \omega_1 bra_H()
        state_times_i_float(&braH, dt*pars.omega_1); 
        state_add_state(&bra, &braH);

        printf("t = %f ns\t", (float)(i)*dt/ns);
        state_sort(&bra);
        state_keep_only_first_max_versors(&bra, 200);
        printf("bra->length = %d\t", bra.length);
        N2 = state_normalisation(&bra);
        printf("N^2 = %f\n", N2);
        // normalisation
        state_times_float(&bra, 1./sqrt(N2));
        print_the_big_four(&bra);
    }

    state_free(&bra);
    state_free(&braH);
    return 0;
}

void state_print(state *psi)
{
    printf("|psi> = ");
    int i = 0;
    fcomplex amp;
    versor ket;
    for(; i < psi->length; i++)
    {
        if(i!=0)
            printf("        ");
        amp = psi->amplitudes[i];
        printf( "(%10.7f,%10.7f)", amp.re, amp.im );
        ket = psi->kets[i];
        show_versor( ket );
        printf( " +\n" );
    }
    printf( "         ...\n" );
}

void print_limited_state(state *psi)
{
    printf("|psi> = ");
    int i = 0;
    fcomplex amp;
    versor ket;
    float total_amp2 = 0.0;
    int cnt = 0;

    int print_more = 1;
    
    while (print_more)
    {
        if (i!=0)
            printf("        ");
        amp = psi->amplitudes[i];
        printf( "(%10.7f,%10.7f)", amp.re, amp.im );
        ket = psi->kets[i];
        show_versor( ket );
        printf( " +\n" );

        cnt++;
        total_amp2 += fcomplex_amplitude_sqr(&amp);
        if ( (1.0f - total_amp2 < 1e-6f) || cnt > 20 )
            print_more = 0;
    }
    printf( "         ...\n" );
}

void print_the_big_four(state *psi)
{
    printf("|psi> = ");
    int i = 0;
    fcomplex amp;
    versor ket;

    for (; i< 4; i++)
    {
        if (i!=0)
            printf("        ");
        amp = psi->amplitudes[i];
        printf( "(%10.7f,%10.7f)", amp.re, amp.im );
        ket = psi->kets[i];
        show_versor( ket );
        printf( " +\n" );
    }
    printf( "         ...\n" );
}

/*
* Select the initial state (one versor) by listing all the quantum numbers of
* the state explicitly.
*/
versor choose_test_versor_quantum_numbers(basis b, parameters p)
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
    printf("index = %d", get_index_from_versor(psi0, b));
    printf("\n");

    return psi0;
}

/* 
* Read a time step from the input
*/
float get_dt()
{
    float dt=0.f;
    printf(" Specify time step dt: ");
    scanf("%f", &dt);
    return dt;
}

/* 
* Ask the user for a number of steps in the propagation
*/
int get_no_of_steps()
{
    int N =0;
    printf(" How many steps do you want to have in the propagation: ");
    scanf("%d", &N);
    if(N < 0)
    {
        printf("Error!\n");
        printf("Negative number of the propagation steps N = %d\n", N);
        exit(0);
    }
    return N;
}