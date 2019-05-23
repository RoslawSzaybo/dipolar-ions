#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "input.h"

/* the code takes a specific number of parameters. */
void check_argc(int argc, char *argv[])
{
        const int parameters_no = 12;
        if (argc != parameters_no)
        {
                printf(" Error!\n");
                printf("  argc = %d\n", argc);
                printf("  argc should be equal to %d\n", parameters_no);
                printf("This is how to execute the program:\n" );
                printf("\t$ %s n1 n3 n5 j1 j2 "
                "mass[u] charge[e] dipole[D] B[MHz]"
                " omega_rho[MHz] omega_z[MHz]\n", argv[0]);
                exit( 0 );
        }
}

/* omega_rho has to be significantly larger than omega_z */
void check_physical_parameters(char *argv[])
{
    float omega_rho = atof(argv[10]);
    float omega_z = atof(argv[11]);

    int is_OK = 1;
    // this two will usually differ by an order of magnitude
    // this is mostly to alert when the input is typed in a 
    // wrong order
    if (omega_z > omega_rho)
    {
            printf(" Error!\n");
            printf(" Omega_z > omega_rho.\n");
            exit(0);
    }
}

/* not all j are implemented. Not all ns are allowed. 
Incorrect basis truncation prints error message and
stops the program execution. */
void check_basis_truncation(basis b)
{
    int is_OK = 1;

    const int jmax = 2; // maximal implemented j

    int j1 = b.j1;
    int j2 = b.j2;
    if (j1 < 0 || j1 > jmax || j2 < 0 || j2 > jmax)
    {
            printf(" Warning.\n");
            printf(" Jmax has to be no smaller than 0"
            " and no larger than %d.\n\n", jmax);
            is_OK = 0;
    }

    int n1 = b.n1;
    int n3 = b.n3;
    int n5 = b.n5;
    if ( n1<1 || n3<1 || n5<1 )
    {
            printf(" Warning.\n");
            printf(" There has to be at least one mode (the zeroth mode) "
            "allowed in each normal mode oscillations.\n\n");
            is_OK = 0;
    }

    if(!is_OK)
    {
        exit(0);
    }
}

/* Returns a number of versors, which can be described by
 the set of quantum numbers {n1,n2,n3,j1,m1,j2,m2}, 
 and which satisfy the truncation conditions stored in `b`, i.e. 
 `n1 < b.n1, ..., -b.j1 < j1 < b,j1, ... ` */
int get_basis_size(basis b)
{
        int n1 = b.n1, n3 = b.n3, n5 = b.n5;
        int n = n1*n3*n5;

        int j1 = b.j1;
        n *= j1*j1+2*j1+1;

        int j2 = b.j2;
        n *= j2*j2+2*j2+1;

        return n;
}

/* Converts the command line argumets into numbers and 
saves it in a convenient structure. */
basis get_basis_truncation(char *argv[])
{
        // truncation parametrs
        int n1 = atoi(argv[1]),
        n3 = atoi(argv[2]),
        n5 = atoi(argv[3]),
        j1 = atoi(argv[4]),
        j2 = atoi(argv[5]);

        return (basis){n1,n3,n5,j1,j2};
}

parameters get_system_parameters(char *argv[])
{
        // physical parameters
        float mass = atof(argv[6]),
        charge = atof(argv[7]),
        dipole = atof(argv[8]),
        B = atof(argv[9]),
        omega_rho = atof(argv[10]),
        omega_z = atof(argv[11]);

        // frequencies
        // omega_1 - the asymmetric modes along the axis
        // omega_3 - the asymmetric modes in the radial direction
        float omega_1 = sqrt(3.f)*omega_z;
        float omega_3 = sqrt(omega_rho*omega_rho - omega_z*omega_z);

        hamiltonian active_terms;
        active_terms.normal_modes = 1;
        active_terms.T_rot = 1;
        active_terms.Vqd_first = 1;
        active_terms.Vqd_zeroth = 1;
        active_terms.Vdd_zeroth = 1;

        return (parameters)
        {mass, charge, dipole, B, omega_1, omega_3, active_terms};
}

void print_input(basis b, parameters p, char *argv[])
{
        float omega_rho = atof(argv[10]);
        float omega_z = atof(argv[11]);

        // pure input
        printf("# Parameters of the model\n");
        printf("#  mass:\t\t%10.2f u\n", p.mass);
        printf("#  charge:\t%10.2f e\n", p.charge);
        printf("#  dipole:\t%10.2f D\n", p.dipole);
        printf("#  B:\t\t%10.2f MHz\n", p.B);
        printf("#  omega_rho:\t%10.2f MHz\n", omega_rho);
        printf("#  omega_z:\t%10.2f MHz\n\n", omega_z);
        printf("# Basis truncation:\t|%d,%d,%d;%d,%d>\n", 
        b.n1, b.n3, b.n5, b.j1, b.j2);
        // values which are functions of the input
        printf("# omega_1:\t%10.2f MHz\n", p.omega_1);
        printf("# omega_3/5:\t%10.2f MHz\n", p.omega_3);
        int basis_size = get_basis_size(b);
        printf("# Basis size:\t\t%d\n\n", basis_size);
}

void print_active_terms_of_Hamiltonian(parameters p)
{
        // Hamiltonian contains many terms
        // not always all of them are active
        printf("# Active terms of the Hamiltionian:\n");
        printf("# \tNomal modes:\t%d\n", p.active_terms.normal_modes);
        printf("# \tT_rot:      \t%d\n", p.active_terms.T_rot);
        printf("# \tV_qd^{(0)}: \t%d\n", p.active_terms.Vqd_zeroth);
        printf("# \tV_qd^{(1)}: \t%d\n", p.active_terms.Vqd_first);
        printf("# \tV_dd^{(0)}: \t%d\n", p.active_terms.Vdd_zeroth);
        printf("\n");
}

void *my_malloc(size_t size, const char *name)
{
        void *array = malloc(size);
        if (array == NULL)
        {
                printf("Not enough memory to allocate %s.\n", name);
                exit(0);
        }
        else
        {
                return array;
        }
}