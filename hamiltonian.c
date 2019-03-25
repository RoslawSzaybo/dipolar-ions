#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian.h"
#include "state.h"

int test_input(int n1, int n3, int n5, int j1, int j2)
{
        const int jmax = 2; // maximal implemented j
        int is_OK = 1;
        if( !(n1*n3*n5) )
        {
                printf("There has to be at least one mode of each allowed.\n");
                is_OK =0;
        }
        
        if(j1<0 || j1 > jmax || j2 < 0 || j2 > jmax)
        {
                printf("At least j=0 has to be in the basis.\n");
                is_OK =0;
        }

        return is_OK;
}

int jm_jump(int j, int m)
{
        /* returns the number of rotational states behind $\ket{j,m}$ */
        int jump=0;
        for(int i=0; i<j; i++)
                jump += 2*i+1;

        jump += m+j;

        return jump;
}

int get_index_from_versor(versor psi, const basis b)
{
        /*
        product basis 
        $$
        \ket{\psi} = \ket{n1}\ket{n3}\ket{n5}\ket{j1,m1}\ket{j2,m2}
        $$
        order of my basis is as follows:
        - vibrations (nx) (b.n1, b.n3, b.n5)
        - rotations (j1, m1), (j2, m2)
        $$
        \sum_{j=0}^{j=b.j1}(2j+1) = (b.j1)^2 + b.j1 + 1
        $$
        */
        int no_vib35 = b.n3 * b.n5;
        int no_vib5 = b.n5;
        int no_j1_rot = b.j1*b.j1 + 2*b.j1 + 1;
        int no_j2_rot = b.j2*b.j2 + 2*b.j2 + 1;
        int no_rot = no_j1_rot * no_j2_rot; 

        int idx=0;
        idx += psi.n1 * no_vib35 * no_rot;
        idx += psi.n3 * no_vib5 * no_rot;
        idx += psi.n5 * no_rot;
        idx += jm_jump(psi.j1, psi.m1) * no_j2_rot;
        idx += jm_jump(psi.j2, psi.m2);

        return idx;
}

int getj(int idx)
{
        int j=0;
        while(idx/(2*j+1))
        {
                idx -= 2*j+1;
                j++;
        }

        return j;
}

int getm(int idx)
{
        int j=0;
        while(idx/(2*j+1))
        {
                idx -= 2*j+1;
                j++;
        }

        return idx-j;
}

versor get_versor_from_index(int idx, const basis b)
{
        int no_vib35 = b.n3 * b.n5;
        int no_vib5 = b.n5;
        int no_j1_rot = b.j1*b.j1 + 2*b.j1 + 1;
        int no_j2_rot = b.j2*b.j2 + 2*b.j2 + 1;
        int no_rot = no_j1_rot * no_j2_rot; 

        versor psi = {0, 0, 0, 0, 0, 0, 0};

        psi.n1 = idx/ no_vib35 / no_rot;
        idx %= no_vib35 * no_rot;
        psi.n3 = idx/no_vib5/no_rot;
        idx %= no_vib5 * no_rot;
        psi.n5 = idx/ no_rot;
        idx %= no_rot;

        psi.j1 = getj(idx/no_j2_rot);
        psi.m1 = getm(idx/no_j2_rot);

        idx %= no_j2_rot;

        psi.j2 = getj(idx);
        psi.m2 = getm(idx);

        return psi;
}

void print_versor(const versor psi, const basis b)
{
        printf("|%4d,%2d,%2d;%d,%3d;%d,%3d>\t\t"
        "index = %d\n",
        psi.n1, psi.n3, psi.n5, psi.j1, psi.m1, psi.j2, psi.m2,
        get_index_from_versor(psi, b));
}



/* test */
void test_idx_to_versor_translation() {
        const basis b = {1000, 2, 3, 0, 0};

        for(int idx=0; idx<20; idx++)
        {
                versor psi1 = get_versor_from_index(idx, b);
                print_versor(psi1, b);
        }
}


void apply_harmonic_oscillator(state * input, state * output)
{
        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                // Harminic oscillator terms
                const float homega1 = 1.1f;
                const float homega3 = 100.0f;
                const float homega5 = 100.0f;

                versor psi_n1 = psi;
                fcomplex A_psi_n1 = (fcomplex){
                        (float)psi.n1 + 0.5f,
                        0.0f
                };
                A_psi_n1.re *= homega1;
                state_add(output, psi_n1, fcomplex_multiply(&A_psi_n1, &A));

                versor psi_n3 = psi;
                fcomplex A_psi_n3 = (fcomplex){
                        (float)psi.n3 + 0.5f,
                        0.0f
                };
                A_psi_n3.re *= homega3;
                state_add(output, psi_n3, fcomplex_multiply(&A_psi_n3, &A));

                versor psi_n5 = psi;
                fcomplex A_psi_n5 = (fcomplex){
                        (float)psi.n5 + 0.5f,
                        0.0f
                };
                A_psi_n5.re *= homega5;
                state_add(output, psi_n5, fcomplex_multiply(&A_psi_n5, &A));
        }
}


/*
$$
\bra{output} = \bra{input} B(\hat{L}_1^2 + hat{L}_2^2 )
$$
*/
void apply_rotational_kinetic_energy(state * input, state * output)
{
        const float B1 = 1.0f;
        const float B2 = 1.0f;

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                versor psi_L2_1 = psi;
                float j1 = (float)(psi.j1);
                float BJJ1 =  B1*j1*(j1+1.f);
                state_add(output, psi_L2_1, fcomplex_times_float(&A, BJJ1));

                versor psi_L2_2 = psi;
                float j2 = (float)(psi.j2);
                float BJJ2 =  B2*j2*(j2+1.f);
                state_add(output, psi_L2_2, fcomplex_times_float(&A, BJJ2));
        }
}

/*
$$
\bra{output} = \bra{input} \omega_1^{-1/2}(\hat{a}_1 + hat{a}_1^\dagger )
$$
*/
void apply_a1_plus_a1dagger(state * input, state * output)
{
        float omega_1 = 0.1;
        float sqrt_omega_1 = sqrt(omega_1);
        float sqrt_omega_1_inv = 1.0f/sqrt_omega_1;

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                versor psi_a1 = psi;
                psi_a1.n1++;
                float a1_factor = sqrt(psi.n1 + 1.0);
                state_add(output, psi_a1, fcomplex_times_float(&A, a1_factor));

                versor psi_a1dagger = psi;
                psi_a1dagger.n1--;
                float a1_dagger_factor;
                if( psi.n1 >= 0)
                {
                        a1_dagger_factor = sqrt(psi.n1);
                }
                else
                {
                        a1_dagger_factor = 0.f;
                }
                
                state_add(output, psi_a1dagger, fcomplex_times_float(&A, a1_dagger_factor));
        }
        state_times_float(output, sqrt_omega_1_inv);
}


void apply_a3_plus_a3dagger(state * input, state * output)
{
        float omega_3 = 0.1;
        float sqrt_omega_3 = sqrt(omega_3);
        float sqrt_omega_3_inv = 1.0f/sqrt_omega_3;

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                versor psi_a3 = psi;
                psi_a3.n3++;
                float a3_factor = sqrt(psi.n3 + 1.0);
                state_add(output, psi_a3, fcomplex_times_float(&A, a3_factor));

                versor psi_a3dagger = psi;
                psi_a3dagger.n3--;
                float a3_dagger_factor;
                if( psi.n3 >= 0)
                {
                        a3_dagger_factor = sqrt(psi.n3);
                }
                else
                {
                        a3_dagger_factor = 0.f;
                }
                state_add(output, psi_a3dagger, fcomplex_times_float(&A, a3_dagger_factor));
        }
        state_times_float(output, sqrt_omega_3_inv);
}


void apply_a5_plus_a5dagger(state * input, state * output)
{
        float omega_5 = 0.1;
        float sqrt_omega_5 = sqrt(omega_5);
        float sqrt_omega_5_inv = 1.0f/sqrt_omega_5;

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                versor psi_a5 = psi;
                psi_a5.n5++;
                float a5_factor = sqrt(psi.n5 + 1.0);
                state_add(output, psi_a5, fcomplex_times_float(&A, a5_factor));

                versor psi_a5dagger = psi;
                psi_a5dagger.n5--;
                float a5_dagger_factor;
                if( psi.n5 >= 0)
                {
                        a5_dagger_factor = sqrt(psi.n5);
                }
                else
                {
                        a5_dagger_factor = 0.f;
                }
                state_add(output, psi_a5dagger, fcomplex_times_float(&A, a5_dagger_factor));
        }
        state_times_float(output, sqrt_omega_5_inv);
}


/*
$$
\bra{output} = \bra{input} hat{d}_1^z 
$$
*/
void apply_dz1(state * input, state * output)
{
        const float d = 1.0f;

        // list of all the constant factors
        // \bra{ja,m}\hat{d}^z\ket{jb,m}.
        // Each of them is a real number which don't depend on the sign of m.
        // So all of those products can be referred to by the  
        // smaller j and by the absolute value of the corresponding m.
        float dz[6] = {
                // j=0, |m| = 0 
                sqrt(1.0/3.0),                        
                // j=1, |m| = 0,1 
                sqrt(4.0/15.0), sqrt(1.0/5.0),        
                // j=2, |m| = 0,1,2 
                sqrt(9.0/35.0), sqrt(8.0/35.0), sqrt(1.0/7.0)
        };

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                int j = psi.j1;
                int m = abs(psi.m1);

                // If the input is incorrect then skip it.
                // We limit the computations to max(j) = 2;
                // m value is tested below
                if(j > 2 || j < 0 )
                        continue;
                
                versor psi_dz1 = psi;
                if( m == j )
                {
                        int idx_j = j+1;
                        int idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j,|j|} in the dz list
                        psi_dz1.j1 = j+1;
                        state_add(output, psi_dz1, fcomplex_times_float(&A, dz[idx]));
                }
                else if(m < j)
                {
                        // go down with j
                        psi_dz1.j1 = j-1;
                        int idx_j = j;
                        int idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j-1,|j-1|} in the dz list
                        idx = idx - (j-1) + m; // a trick I like which returns the \ket{j-1,|m|}
                        state_add(output, psi_dz1, fcomplex_times_float(&A, dz[idx]));

                        // go up with j
                        psi_dz1.j1 = j+1;
                        idx_j = j+1;
                        idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j,|j|}
                        idx = idx - j + m ; // \ket{j, |m|}
                        state_add(output, psi_dz1, fcomplex_times_float(&A, dz[idx]));
                }
                else
                        continue; // incorrect value of m;
        }
}

/*
$$
\bra{output} = \bra{input} hat{d}_2^z 
$$
*/
void apply_dz2(state * input, state * output)
{
        const float d = 1.0f;

        // list of all the constant factors
        // \bra{ja,m}\hat{d}^z\ket{jb,m}.
        // Each of them is a real number which don't depend on the sign of m.
        // So all of those products can be referred to by the  
        // smaller j and by the absolute value of the corresponding m.
        float dz[6] = {
                // j=0, |m| = 0 
                sqrt(1.0/3.0),                        
                // j=1, |m| = 0,1 
                sqrt(4.0/15.0), sqrt(1.0/5.0),        
                // j=2, |m| = 0,1,2 
                sqrt(9.0/35.0), sqrt(8.0/35.0), sqrt(1.0/7.0)
        };

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                int j = psi.j2;
                int m = abs(psi.m2);

                // If the input is incorrect then skip it.
                // We limit the computations to max(j) = 2;
                // m value is tested below
                if(j > 2 || j < 0 )
                        continue;
                
                versor psi_dz2 = psi;
                if( m == j )
                {
                        int idx_j = j+1;
                        int idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j,|j|} in the dz list
                        psi_dz2.j2 = j+1;
                        state_add(output, psi_dz2, fcomplex_times_float(&A, dz[idx]));
                }
                else if(m < j)
                {
                        // if(j!= 0) - it's useless because of what is above
                        // for m<j, j has to be larger than 0;
                        // go down with j
                        psi_dz2.j2 = j-1;
                        int idx_j = j;
                        // the index of the element ket{j-1,|j-1|} in the dz list
                        int idx = (idx_j*idx_j+idx_j)/2-1; 
                        // a trick I like. It returns the index of \ket{j-1,|m|}
                        idx = idx - (j-1) + m; 
                        state_add(output, psi_dz2, fcomplex_times_float(&A, dz[idx]));

                        // go up with j
                        psi_dz2.j2 = j+1;
                        idx_j = j+1;
                        idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j,|j|}
                        idx = idx - j + m ; // \ket{j, |m|}
                        state_add(output, psi_dz2, fcomplex_times_float(&A, dz[idx]));
                }
                else
                        continue; // incorrect value of m; double check xd
        }
}


void apply_dy1(state * input, state * output)
{
        const float d = 1.0f;
        // WARNING, all of them have to be multiplied by 
        // the imaginary unit i
        // list of all the constant factors
        float dy[6] = {
                // 0
                // \ket{0,0} down to \ket{1,-1}, \ket{1,1}
                sqrt(1.0/6.0),sqrt(1.0/6.0),
                // 1
                // \ket{1,-1} down to \ket{2,-2}, \ket{2,0}
                sqrt(1.0/5.0), sqrt(1.0/30.0),
                // \ket{1,0} down to \ket{2,-1}, \ket{2,1}
                sqrt(1.0/10.0), sqrt(1.0/10.0),
                // \ket{1,1} down to \ket{2,0}, \ket{2,2}
                sqrt(1.0/30.0), sqrt(1.0/5.0),
                // 2
                // \ket{2,-2} down to \ket{3,-3} \ket{3,-1}
                sqrt(3.0/14.0), sqrt(1.0/70.0),
                // \ket{2,-1} down to \ket{3,-2} \ket{3,0}
                sqrt(1.0/7.0), sqrt(3.0/70.0),
                // \ket{2,0} down to \ket{3,-1} \ket{3,1}
                sqrt(3.0/35.0), sqrt(3.0/35.0),
                // \ket{2,1} down to \ket{3,0} \ket{3,2}
                sqrt(3.0/70.0), sqrt(1.0/7.0), 
                // \ket{2,2} down to \ket{3,1} \ket{3,3}
                sqrt(1.0/70.0), sqrt(3.0/14.0)
        };

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                // minus comes from the complex conjugation
                // going up add another minus from the 
                // state_add(output, psi_dy1, fcomplex_times_i_float(&A, -dy[idx]));
        }
}



void apply_dx1(state * input, state * output)
{
        const float d = 1.0f;
        // WARNING, all of them have to be multiplied by 
        // the imaginary unit i
        // list of all the constant factors
        float dx[18] = {
                // 0
                // \ket{0,0} down to \ket{1,-1}, \ket{1,1}
                sqrt(1.0/6.0), -sqrt(1.0/6.0),
                // 1
                // \ket{1,-1} down to \ket{2,-2}, \ket{2,0}
                sqrt(1.0/5.0), -sqrt(1.0/30.0),
                // \ket{1,0} down to \ket{2,-1}, \ket{2,1}
                sqrt(1.0/10.0), -sqrt(1.0/10.0),
                // \ket{1,1} down to \ket{2,0}, \ket{2,2}
                sqrt(1.0/30.0), -sqrt(1.0/5.0),
                // 2
                // \ket{2,-2} down to \ket{3,-3} \ket{3,-1}
                sqrt(1.0/14.0), -sqrt(1.0/70.0),                 // HERE \ket{2,-2} to \ket{3,-3} should be checked
                // \ket{2,-1} down to \ket{3,-2} \ket{3,0}
                sqrt(1.0/7.0), -sqrt(3.0/70.0),
                // \ket{2,0} down to \ket{3,-1} \ket{3,1}
                sqrt(3.0/35.0), -sqrt(3.0/35.0),
                // \ket{2,1} down to \ket{3,0} \ket{3,2}
                sqrt(3.0/70.0), -sqrt(1.0/7.0), 
                // \ket{2,2} down to \ket{3,1} \ket{3,3}
                sqrt(1.0/70.0), -sqrt(1.0/14.0)                 // here is the same point
        };

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                // state_add(output, psi_dy1, fcomplex_times_float(&A, dx[idx]));
        }
}

/*
$$
\bra{output} 
= 
\bra{input} \frac{q}{4\pi\epsilon_0} \frac{hat{d}_1^z -\hat{d}_2^z}{D_0^2}
$$
*/
void apply_charge_dipole_zero(state *input, state *output)
{
        const float q = 1.f;
        const float pie = 1.f; // \frac{1}{4\pi\epsilon_0}
        const float D0_2_inv = 1.f/2.f/2.f;

        state work_state;
        state_init(&work_state);
        
        // dz_1
        apply_dz1(input, &work_state);
        state_times_float(&work_state, q*pie*D0_2_inv);
        state_add_state(output, &work_state);
        state_free(&work_state);

        // -dz_2
        apply_dz2(input, &work_state);
        state_times_float(&work_state, -q*pie*D0_2_inv);
        state_add_state(output, &work_state);
        state_free(&work_state);
}

/*
This function returns a BRA state (it's not a ket)
which is a result of acting with bra called psi
on the Hamiltonian.
$$
\bra{return} = \bra{psi}\hat{H}
$$
*/
state bra_H(state* psi)
{
        state output_bra;
        state_init(&output_bra);

        apply_harmonic_oscillator(psi, &output_bra);
        apply_rotational_kinetic_energy(psi, &output_bra);
        apply_charge_dipole_zero(psi, &output_bra);

        //apply_dz1_m_dz2(psi, &output_bra);
        /*state aadagger;
        state_init(&aadagger);
        apply_a1_plus_a1dagger(psi, &aadagger);
        state_add_state(&output_bra, &aadagger);
        state_free(&aadagger);
        */

       /*
        state dz;
        state_init(&dz);
        apply_dz1(psi, &dz);
        state_add_state(&output_bra, &dz);
        state_free(&dz);
        */

        return output_bra;
}


int valid_versor(versor psi, basis b)
{
        if( 0 > psi.n1 || psi.n1 >= b.n1)
                return 0;
        
        if( 0 > psi.n3 || psi.n3 >= b.n3)
                return 0;

        if( 0 > psi.n5 || psi.n5 >= b.n5)
                return 0;

        if( psi.j1 < 0 || psi.j1 > b.j1)
                return 0;

        if( psi.m1 < -b.j1 || psi.m1 > b.j1 )
                return 0;

        if( psi.j2 < 0 || psi.j2 > b.j2)
                return 0;

        if( psi.m2 < -b.j2 || psi.m2 > b.j2 )
                return 0;

        return 1;
}


void construct_Hamiltonian(fcomplex* a, basis b)
{
        // fill-in with zeroes
        int basis_size = b.n1*b.n3*b.n5*(b.j1*b.j1+2*b.j1+1)*(b.j2*b.j2+2*b.j2+1);
        fcomplex zero = {0.0f, 0.0f};

        for(int i = 0; i<basis_size*basis_size; i++)
                a[i] = zero;

        // scan through the matrix rows and replace all the non-zero elements
        for(int i = 0; i<basis_size; i++)
        {
                // scan through rows
                versor psi0 = get_versor_from_index(i, b);
                state state0;
                state_init(&state0);
                state_add(&state0, psi0, (fcomplex){1.0f, 0.0f});

                // act with H from the left
                state psiH = bra_H(&state0);

                versor loop_versor;
                fcomplex loop_amplitude;
                for(int l=0; l<psiH.length; l++)
                {
                        loop_versor = state_get_versor(&psiH, l);
                        // not sure but I have already put reutrn NULL earlier
                        /* do something like
                        if(!loop_versor)
                                continue;
                        // or check it another way
                        */

                        // check if the versor is a valid versor
                        if(!valid_versor(loop_versor, b))
                                continue;

                        // check if it is above the diagonal
                        int loc_idx = get_index_from_versor(loop_versor, b);
                        if( loc_idx < i )
                                continue;

                        // apply amplitude;
                        fcomplex amp = state_get_amplitude(&psiH, l);
                        /* again
                        if(!amp)
                                continue;
                        */

                        // add amplitudes
                        a[i*basis_size+loc_idx].re += amp.re;
                        a[i*basis_size+loc_idx].im += amp.im;
                }
               state_free(&psiH);
        }
}


void print_only_versor(const versor psi)
{
        printf("|%5d,%5d,%5d;%d,%3d;%d,%3d>",
        psi.n1, psi.n3, psi.n5, psi.j1, psi.m1, psi.j2, psi.m2);
}

void test_bra_H() 
{
        const basis b = {1, 1, 1, 0, 2};

        versor psi1 = get_versor_from_index(0, b);

        state state0;
        state_init(&state0);

        state_add(&state0, psi1, (fcomplex){1.0f, 0.0f});
        printf("Initial state is only a single versor.\n");
        printf("Amplitude = %4.2f+i%4.2f\t",1.0f,0.0f);
        print_only_versor(psi1);
        printf("\t\tidx = %7d\n", get_index_from_versor(psi1, b));

        state sps = bra_H(&state0);
        versor loop_versor;
        fcomplex loop_amplitude;
        for(int l=0; l < sps.length; l++)
        {
                loop_versor = state_get_versor(&sps, l);
                loop_amplitude = state_get_amplitude(&sps, l);
                printf("Amplitude = %4.2f+i%4.2f\t", 
                loop_amplitude.re, loop_amplitude.im);
                print_only_versor(loop_versor);
                printf("\t\tidx = %7d\n", get_index_from_versor(loop_versor, b));
        }
}