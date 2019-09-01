#include "../spectrum/versor.h"
#include "../spectrum/input.h"
#include "../spectrum/state.h"
#include "../spectrum/hamiltonian.h" // valid_versor()
#include <stdio.h>
#include <math.h>

versor choose_test_versor_quantum_numbers(basis b, parameters p);
double get_dt();
int get_no_of_steps();
void state_print(state *psi);
void print_limited_state(state *psi, double threshold);
void print_the_big_four(state *psi);
parameters SrYb_parameters();
parameters excited_SrYb_parameters();

int main(int argc, char *argv[])
{
    //const parameters pars = SrYb_parameters();
    const parameters pars = excited_SrYb_parameters();

    print_system_parameters(pars);
    print_active_terms_of_Hamiltonian(pars);
    /* 
    versor psi0 = choose_test_versor_quantum_numbers(b, pars);
    double time_step_ns = get_dt();
    */
    versor psi0 = (versor){0, 0, 0, 1, 0, 0, 0};
    double dt = 10.0e-9;
    double print_every = 1.0e-6;
    double propagation_time = 5.0e-3;
    double print_threshold = 1.0e-4;
    double N2 = 0.0;


    printf("# Attention please!\n");
    printf("#  All printed states below are in fact duals i.e. <psi|,"
    " so if you want vectors i.e |psi> you have to complex conjugate amplitudes.\n");

    state bra;
    state braH;
    state_init(&bra);
    state_init(&braH);
    state_add(&bra, psi0, (dcomplex){1.0, 0.0});

    printf("t = %f us\n", 0.0);
    state_sort(&bra);
    N2 = get_state_normalisation(&bra);
    state_times_double(&bra, 1.0/sqrt(N2));
    state_print(&bra);

    int i=1;
    int steps = propagation_time/dt;
    int print_divisor = print_every/dt;
    for (; i<steps; i++)
    {
        braH = bra_H(&bra, pars);
        // i\hbar \ket{psi} = H \ket{psi}
        // bra{psi} = \frac{i}{\hbar}\bra{psi}H
        // H is in unit \hbar \omega_1 so you have to multiply it 
        // \bra{psi} = i \omega_1 bra_H()
        state_times_i_double(&braH, dt*pars.omega_1); 
        state_add_state(&bra, &braH);
        state_sort(&bra);
        state_keep_only_first_max_versors(&bra, 200);
        N2 = get_state_normalisation(&bra);
        state_times_double(&bra, 1.0/sqrt(N2));

        if (!(i % print_divisor))
        {
            printf("t = %f us\t", (double)(i)*dt*1.0e6);
            printf("N^2 = %f\n", N2);
            //print_limited_state(&bra, print_threshold);
            print_the_big_four(&bra);
        }
    }

    state_free(&bra);
    state_free(&braH);
    return 0;
}

parameters SrYb_parameters()
{
    // SrYb^+
    // $ \omega_1 = \sqrt{3} \omega_z 
    // $ \omega_3 = \sqrt{omega_rho^2 - omega_z^2}$
    double omega_rho = 1.4;
    double omega_z = 0.16;
    double omega_1 = sqrt(3.0)*omega_z;
    double omega_3 = sqrt(omega_rho*omega_rho - omega_z*omega_z);

    hamiltonian activate;
    activate.normal_modes = 1;
    activate.T_rot = 1;
    activate.Vqd_zeroth = 1;
    activate.Vqd_first = 1;
    activate.Vdd_zeroth = 1;

    double mass = 261.0;
    double charge = 1.0;
    double dipole = 4.745;
    double B = 503.7;
    
    return (parameters){mass, charge, dipole, B, omega_1, omega_3, activate};
}

parameters excited_SrYb_parameters()
{
    // SrYb^+
    // $ \omega_1 = \sqrt{3} \omega_z 
    // $ \omega_3 = \sqrt{omega_rho^2 - omega_z^2}$
    double omega_rho = 1.4;
    double omega_z = 0.16;
    double omega_1 = sqrt(3.0)*omega_z;
    double omega_3 = sqrt(omega_rho*omega_rho - omega_z*omega_z);

    hamiltonian activate;
    activate.normal_modes = 1;
    activate.T_rot = 1;
    activate.Vqd_zeroth = 1;
    activate.Vqd_first = 1;
    activate.Vdd_zeroth = 1;

    double mass = 261.0;
    double charge = 1.0;
    double dipole = 150.0;
    double B = 12.0;
    
    return (parameters){mass, charge, dipole, B, omega_1, omega_3, activate};
}

void state_print(state *psi)
{
    printf("|psi> = ");
    int i = 0;
    dcomplex amp;
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

/*  
 * Print all the states of an amplitude larger than a threshold
 */
void print_limited_state(state *psi, double threshold)
{
    printf("|psi> = ");
    dcomplex amp;
    versor ket;
    double amp_modulus_2;

    int i=0;
    for (;; i++)
    {
        // we can not print more versors than 
        // the number of versors available in the state
        if ( i>= psi->length )
            break;

        // if the versors are too insignificant, we skip them
        amp = psi->amplitudes[i];
        amp_modulus_2 = dcomplex_amplitude_sqr(&amp);

        if (amp_modulus_2 < threshold)
            break;

        // spectial formatng for the first line
        if (i!=0)
            printf("        ");

        // printing of a versors
        printf( "(%10.7f,%10.7f)", amp.re, amp.im );
        ket = psi->kets[i];
        show_versor( ket );
        printf( " +\n" );
    }
    printf( "         ...\n" );
}

void print_the_big_four(state *psi)
{
    printf("|psi> = ");
    int i = 0;
    dcomplex amp;
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
double get_dt()
{
    double dt=0.0;
    printf(" Specify time step dt: ");
    scanf("%lf", &dt);
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
