#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hamiltonian.h"

void apply_harmonic_oscillator(state * input, state * output, 
                                const parameters pars)
{
        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                // Harmonic oscillator energy scale factors
                // everything is in the units of 
                // $\hbar \omega_1$
                const float homega1 = 1.0;
                const float homega3 = pars.omega_3/pars.omega_1;
                const float homega5 = pars.omega_3/pars.omega_1;

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
\bra{output} = \bra{input} \frac{2\pi}{\hbar} B(\hat{L}_1^2 + hat{L}_2^2 )
$$
*/
void apply_rotational_kinetic_energy(state * input, state * output, 
                                        const parameters pars)
{
    // The kinetic rotational energy of a state of angular quantum 
    // number j is equal to 2\pi B \hbar j(j+1).
    // 
    // Everything is in the units of $\hbar \omega_1$
    // after dividing the above by $\hbar \omega_1$ 
    // one ends up with 2.pi*B/\omega_1
    const float B1 = 2*M_PI*pars.B/pars.omega_1;
    const float B2 = 2*M_PI*pars.B/pars.omega_1;

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
\bra{output} = \bra{input} (\hat{a}_1 + hat{a}_1^\dagger )
$$
*/
void apply_a1_plus_a1dagger(state * input, state * output,
                                    const parameters pars)
{
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
}


void apply_a3_plus_a3dagger(state * input, state * output,
                                    const parameters pars)
{
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
}


void apply_a5_plus_a5dagger(state * input, state * output,
                                    const parameters pars)
{
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
}


/*
$$
\bra{output} = \bra{input} hat{d}_1^z 
$$
*/
void apply_dz1(state * input, state * output,
                const parameters pars)
{
        const float d = pars.dipole;

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
                if( m == j ) // maximal/minimal m so only jump up in j
                {
                        int idx_j = j+1;
                        int idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j,|j|} in the dz list
                        psi_dz1.j1 = j+1;
                        state_add(output, psi_dz1, fcomplex_times_float(&A, d*dz[idx]));
                }
                else if(m < j)
                {
                        // go down with j
                        psi_dz1.j1 = j-1;
                        int idx_j = j;
                        int idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j-1,|j-1|} in the dz list
                        idx = idx - (j-1) + m; // a trick I like which returns the \ket{j-1,|m|}
                        state_add(output, psi_dz1, fcomplex_times_float(&A, d*dz[idx]));

                        // go up with j
                        psi_dz1.j1 = j+1;
                        idx_j = j+1;
                        idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j,|j|}
                        idx = idx - j + m ; // \ket{j, |m|}
                        state_add(output, psi_dz1, fcomplex_times_float(&A, d*dz[idx]));
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
void apply_dz2(state * input, state * output,
                const parameters pars)
{
        const float d = pars.dipole;

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
                        state_add(output, psi_dz2, fcomplex_times_float(&A, d*dz[idx]));
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
                        state_add(output, psi_dz2, fcomplex_times_float(&A, d*dz[idx]));

                        // go up with j
                        psi_dz2.j2 = j+1;
                        idx_j = j+1;
                        idx = (idx_j*idx_j+idx_j)/2-1; // the index of the element ket{j,|j|}
                        idx = idx - j + m ; // \ket{j, |m|}
                        state_add(output, psi_dz2, fcomplex_times_float(&A, d*dz[idx]));
                }
                else
                        continue; // incorrect value of m;
        }
}


int get_dy_idx(int j, int m)
{
        // function which determines the index
        // of the factor $\bra{j+1,m-1}\hat{d}^y\ket{j,m}$
        // in the list dy

        // \sum_j=0^J (2j+1) = (J+1) + 2 \sum_j=0^J j 
        // = (J+1) + (J+1)J = (J+1)(J+1);
        // number of $|ket{j,j}$ in the states list
        int idx = (j+1)*(j+1); 
        // number of $\ket{j,m}$ in the rotations states list
        // number of states in given J 
        // minus 
        // the positon of the $\ket{j,m}$ state in a list of J=const states
        idx -= 2*j+1 - (m+j+1);
        // there are twice as many entires in dy
        idx *= 2;
        // but the one which takes to j+1, m-1 is first, and we count from 0
        idx -= 2;
        return idx;
}

// highest implemented input j
const int jmax = 2;

void apply_dy1(state * input, state * output,
                        const parameters pars)
{
        const float d = pars.dipole;
        // WARNING, all of them have to be multiplied by 
        // the imaginary unit i

        // list of all the constant factors
        // 2*(jmax+1)^2 = 18
        float dy[18] = {
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
                int j = psi.j1;
                int m = psi.m1;

                // check if the input state is correct
                if(j<0 || j > jmax || m > j || m < -j)
                        continue;

                // always add j+1, m-1 
                versor psi_work = psi;
                psi_work.j1=j+1;
                psi_work.m1=m-1;

                int idx = get_dy_idx(j, m);
                // I work with bra so there is one conjugation -> -i
                state_add(output, psi_work, fcomplex_times_i_float(&A, -d*dy[idx]));

                // always add j+1, m+1;
                // almost like in the previous case but I go to m+1
                psi_work.m1=m+1;
                idx += 1;
                state_add(output, psi_work, fcomplex_times_i_float(&A, -d*dy[idx]));

                // sometimes j-1, m-1 is also included
                if(m > -(j-1))
                {
                        psi_work.j1=j-1;
                        psi_work.m1=m-1;

                        idx = get_dy_idx(j-1,m-1)+1; // going to larger m so +1
                        // here there is a second conjugation as the 
                        // bracket goes in the other direction
                        state_add(output, psi_work, fcomplex_times_i_float(&A, d*dy[idx]));
                }

                // similarily, j-1, m+1 is also sometimes included
                if(m < j-1)
                {
                        psi_work.j1=j-1;
                        psi_work.m1=m+1;

                        idx = get_dy_idx(j-1,m+1);
                        state_add(output, psi_work, fcomplex_times_i_float(&A, d*dy[idx]));
                }
        }
}

void apply_dy2(state * input, state * output,
                        const parameters pars)
{
        const float d = pars.dipole;
        // WARNING, all of them have to be multiplied by 
        // the imaginary unit i

        // list of all the constant factors
        // 2*(jmax+1)^2 = 18
        float dy[18] = {
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

                int j = psi.j2;
                int m = psi.m2;
                versor psi_work = psi;

                // check if the input state is correct
                if(j<0 || j > jmax || m > j || m < -j)
                        continue;

                // always add j+1, m-1 
                psi_work.j2=j+1;
                psi_work.m2=m-1;
                int idx = get_dy_idx(j, m);
                // I work with bra so there is one conjugation -> -i
                state_add(output, psi_work, fcomplex_times_i_float(&A, -d*dy[idx]));

                // always add j+1, m+1;
                // almost like in the previous case but I go to m+1
                psi_work.m2=m+1;
                idx += 1;
                state_add(output, psi_work, fcomplex_times_i_float(&A, -d*dy[idx]));

                // sometimes j-1, m-1 is also included
                if(m > -(j-1))
                {
                        psi_work.j2=j-1;
                        psi_work.m2=m-1;

                        idx = get_dy_idx(j-1,m-1)+1;// going to larger m so +1
                        // here there is a second conjugation as the 
                        // bracket goes in the other direction
                        state_add(output, psi_work, fcomplex_times_i_float(&A, d*dy[idx]));
                }

                // similarily, j-1, m+1 is also sometimes included
                if(m < j-1)
                {
                        psi_work.j2=j-1;
                        psi_work.m2=m+1;

                        idx = get_dy_idx(j-1,m+1); 
                        state_add(output, psi_work, fcomplex_times_i_float(&A, d*dy[idx]));
                }
        }
}



void apply_dx1(state * input, state * output,
                        const parameters pars)
{
        const float d = pars.dipole;
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
                sqrt(3.0/14.0), -sqrt(1.0/70.0),
                // \ket{2,-1} down to \ket{3,-2} \ket{3,0}
                sqrt(1.0/7.0), -sqrt(3.0/70.0),
                // \ket{2,0} down to \ket{3,-1} \ket{3,1}
                sqrt(3.0/35.0), -sqrt(3.0/35.0),
                // \ket{2,1} down to \ket{3,0} \ket{3,2}
                sqrt(3.0/70.0), -sqrt(1.0/7.0), 
                // \ket{2,2} down to \ket{3,1} \ket{3,3}
                sqrt(1.0/70.0), -sqrt(3.0/14.0) 
        };

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                int j = psi.j1;
                int m = psi.m1;
                versor psi_work = psi;

                // check if the input state is correct
                if(j<0 || j > jmax || m > j || m < -j)
                        continue;

                // always add j+1, m-1 
                psi_work.j1=j+1;
                psi_work.m1=m-1;
                int idx = get_dy_idx(j, m);
                state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));

                // always add j+1, m+1 
                psi_work.j1=j+1;
                psi_work.m1=m+1;
                idx += 1;
                state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));

                // j-1,m-1
                if(m > -(j-1))
                {
                        idx = get_dy_idx(j-1,m-1)+1;
                        psi_work.j1=j-1;
                        psi_work.m1=m-1;
                        state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));
                }
                // j-1,m+1
                if(m < j-1)
                {
                        idx = get_dy_idx(j-1,m+1);
                        psi_work.j1=j-1;
                        psi_work.m1=m+1;
                        state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));
                }
        }
}


void apply_dx2(state * input, state * output,
                        const parameters pars)
{
        const float d = pars.dipole;
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
                sqrt(3.0/14.0), -sqrt(1.0/70.0),
                // \ket{2,-1} down to \ket{3,-2} \ket{3,0}
                sqrt(1.0/7.0), -sqrt(3.0/70.0),
                // \ket{2,0} down to \ket{3,-1} \ket{3,1}
                sqrt(3.0/35.0), -sqrt(3.0/35.0),
                // \ket{2,1} down to \ket{3,0} \ket{3,2}
                sqrt(3.0/70.0), -sqrt(1.0/7.0), 
                // \ket{2,2} down to \ket{3,1} \ket{3,3}
                sqrt(1.0/70.0), -sqrt(3.0/14.0)
        };

        for (int p=0; p < input->length; p++)
        {
                versor psi = input->kets[p];
                fcomplex A = input->amplitudes[p];

                int j = psi.j2;
                int m = psi.m2;
                versor psi_work = psi;

                // check if the input state is correct
                if(j<0 || j > jmax || m > j || m < -j)
                        continue;

                // always add j+1, m-1 
                psi_work.j2=j+1;
                psi_work.m2=m-1;
                int idx = get_dy_idx(j, m);
                state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));

                // always add j+1, m+1 
                psi_work.j2=j+1;
                psi_work.m2=m+1;
                idx += 1;
                state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));

                // j-1,m-1
                if(m > -(j-1))
                {
                        idx = get_dy_idx(j-1,m-1)+1;
                        psi_work.j2=j-1;
                        psi_work.m2=m-1;
                        state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));
                }
                // j-1,m+1
                if(m < j-1)
                {
                        idx = get_dy_idx(j-1,m+1);
                        psi_work.j2=j-1;
                        psi_work.m2=m+1;
                        state_add(output, psi_work, fcomplex_times_float(&A, d*dx[idx]));
                }
        }
}


/*
$$
\bra{output} = \bra{input}(\hat{d}^y_1 - \hat{d}^y_2)
$$
*/
void apply_dx1_m_dx2(state *input, state *output,
                        const parameters pars)
{
        state work_state;
        state_init(&work_state);
        
        apply_dx1(input, &work_state, pars);
        state_add_state(output, &work_state);
        state_free(&work_state);

        apply_dx2(input, &work_state, pars);
        state_times_float(input, -1.f);
        state_add_state(output, &work_state);
        state_free(&work_state);
}

/*
$$
\bra{output} = \bra{input}(\hat{d}^y_1 - \hat{d}^y_2)
$$
*/
void apply_dy1_m_dy2(state *input, state *output,
                        const parameters pars)
{
        state work_state;
        state_init(&work_state);
        
        apply_dy1(input, &work_state, pars);
        state_add_state(output, &work_state);
        state_free(&work_state);

        apply_dy2(input, &work_state, pars);
        state_times_float(input, -1.f);
        state_add_state(output, &work_state);
        state_free(&work_state);
}

/*
$$
\bra{output} = \bra{input}(\hat{d}^y_1 - \hat{d}^y_2)
$$
*/
void apply_dz1_m_dz2(state *input, state *output,
                        const parameters pars)
{
        state work_state;
        state_init(&work_state);
        
        apply_dz1(input, &work_state, pars);
        state_add_state(output, &work_state);
        state_free(&work_state);

        apply_dz2(input, &work_state, pars);
        state_times_float(input, -1.f);
        state_add_state(output, &work_state);
        state_free(&work_state);
}

/*
$$
\bra{output} 
= 
\bra{input} \frac{q}{4\pi\epsilon_0} \frac{hat{d}_1^z -\hat{d}_2^z}{|2z_0^2}
$$

the factor divided by \hbar\omega_1 "reduces" to
$
(\omega_1/q)^{1/3}
(\frac{m}{6})^{2/3}
(d/\hbar)
(frac{1}{4\pi\epsilon_0})^{1/3}
$

pie := \frac{1}{4\pi\epsilon_0}
= by def = \frac{\mu_0 c^2}{4\pi}
= 10^{-7} c^2 N/A^2
= 8987551787.368176 (kg*m^3/s^4/A^2)
= 0.138931 (m^2 u / s^2 /e^2)
= (0.517925 (m^2 u / s^2 /e^2)^{1/3})^{3}

1kg = 6.022 * 10^{26} u

\hbar \approx 1.0545718 J\cdot s

\frac{1}{\hbar} \approx 
1.574555188 \cdot 10^7 s/(m^2 u)

1/\hbar (pie)^{1/3}
\approx
8.155 * 10^6 [s^{1/3} m/ (u^{2/3} e^{2/3} )]

[electric dipole moment unit]
1D 
\approx 
0.20819 *10^{-10} e m

so the final formula for the prefactor
which allows to use the units described in the main funcion
gives a numeric value 
$$
100(1/(se)^(1/3))*
8.155 * 10^6 s^{1/3} m/ (u^{2/3} e^{2/3} )
0.20819*10^{-10}em*
(1/6)^(2/3)
u^{2/3}
=
2.06*10^{-2}
$$

and the formula for the prefactor becomes
(|omega_1|/|q|)^{1/3}
|mass|^{2/3}
|dipole|
*/
void apply_charge_dipole_zero(state *input, state *output,
                        const parameters pars)
{
        const float units_and_constants = 0.0206f;

        float factor = pow(pars.omega_1/pars.charge, 1.0/3.0);
        factor *= pow(pars.mass, 2.0/3.0);

        state work_state;
        state_init(&work_state);
        
        // dz_1
        apply_dz1(input, &work_state, pars);
        state_times_float(&work_state, factor*units_and_constants);
        state_add_state(output, &work_state);
        state_free(&work_state);

        // -dz_2
        apply_dz2(input, &work_state, pars);
        state_times_float(&work_state, factor*units_and_constants);
        state_add_state(output, &work_state);
        state_free(&work_state);
}

/*
$$
\bra{output} 
= 
\bra{input} 
\frac{q}{4\pi\epsilon_0}
\frac{1}{|2z_0|^3}
(\frac{\hbar}{m})^(1/2)
\begin{pmatrix}
\omega_3^{-1/2}(a3 + a3^dagger)\\
\omega_5^{-1/2}(a5 + a5^dagger)\\
\omega_1^{-1/2}(a1 + a1^dagger)
\end{pmatrix}
\cdot
(\vec{hat{d}}_1 -\vec{\hat{d}_2})
$$

The units and constants collected together give
$$
(1 + \delta_{iz})
\frac{
        |omega_1||mass|^{1/2}|dipole|
}{
        |omega_i|^{1/2}|charge|
}
1.37689e-5;
$$
*/
void apply_charge_dipole_first(state *input, state *output,
                        const parameters pars)
{
        float factor = 1.37689e-5f; 
        factor *= pars.omega_1/pars.charge;
        factor *= pow(pars.mass, 1.0/2.0);

        state sum_container;
        state_init(&sum_container);

        state work_state;

        /* z-part */
        state_init(&work_state);
        apply_dz1_m_dz2(input, &work_state, pars);
        state_times_float(&work_state, factor);
        float loc_omega = 2.0/pow(pars.omega_1, 1.0/2.0);
        state_times_float(&work_state, loc_omega);
        apply_a1_plus_a1dagger(&work_state, &sum_container, pars);
        state_free(&work_state);

        /* y-part */
        state_init(&work_state);
        apply_dy1_m_dy2(input, &work_state, pars);
        state_times_float(&work_state, factor);
        loc_omega = 1.0/pow(pars.omega_3, 1.0/2.0);
        state_times_float(&work_state, loc_omega);
        apply_a5_plus_a5dagger(&work_state, &sum_container, pars);
        state_free(&work_state);

        /* x-part */
        state_init(&work_state);
        apply_dx1_m_dx2(input, &work_state, pars);
        state_times_float(&work_state, factor);
        loc_omega = 1.0/pow(pars.omega_3, 1.0/2.0);
        state_times_float(&work_state, loc_omega);
        apply_a3_plus_a3dagger(&work_state, &sum_container, pars);
        state_free(&work_state);

        state_add_state(output, &sum_container);
        state_free(&sum_container);
}



/*
$$
\bra{output} 
= 
\bra{input} 
\frac{1}{4\pi\epsilon_0}
\frac{1}{|2z_0|^3}
\vec{\hat{d}}_1
\cdot
\vec{\hat{d}}_2
-
3
\hat{d}^z_1\hat{d}^z_2
=
\bra{input} 
\frac{1}{4\pi\epsilon_0}
\frac{1}{|2z_0|^3}
(
\hat{d}^x_1
\hat{d}^x_2
+
\hat{d}^y_1
\hat{d}^y_2
-2
\hat{d}^z_1
\hat{d}^z_2
)

the and units give and constants give the following prefactor
|d|^2/|q|^2 |m| |omega_1| 
*1.136*10^{-9}
$$
*/
void apply_dipole_dipole_zero(state *input, state *output,
                        const parameters pars)
{
        float factor = 1.136e-09f;
        factor *= pars.mass * pars.omega_1 / pars.charge / pars.charge;
        state sum_container;
        state_init(&sum_container);

        state work_state;
        state second_work_state;

        /* x-part */
        state_init(&work_state);
        apply_dx1(input, &work_state, pars);
        state_init(&second_work_state);
        apply_dx2(&work_state, &second_work_state, pars);
        state_free(&work_state);
        state_add_state(&sum_container, &second_work_state);
        state_free(&second_work_state);

        /* y-part */
        state_init(&work_state);
        apply_dy1(input, &work_state, pars);
        state_init(&second_work_state);
        apply_dy2(&work_state, &second_work_state, pars);
        state_free(&work_state);
        state_add_state(&sum_container, &second_work_state);
        state_free(&second_work_state);

        /* z-part */
        state_init(&work_state);
        apply_dz1(input, &work_state, pars);
        state_init(&second_work_state);
        apply_dz2(&work_state, &second_work_state, pars);
        state_free(&work_state);
        state_times_float(&second_work_state, -2.0f);
        state_add_state(&sum_container, &second_work_state);
        state_free(&second_work_state);

        state_times_float(&sum_container, factor);

        state_add_state(output, &sum_container);
        state_free(&sum_container);
}

/*
This function returns a BRA state (it's not a ket)
which is a result of acting with bra called psi
on the Hamiltonian.
$$
\bra{return} = \bra{psi}\hat{H}
$$
*/
state bra_H(state* psi, const parameters pars)
{
    state output_bra;
    state_init(&output_bra);

    // quantised $T_tr + V_qq + V_trap$
    //apply_harmonic_oscillator(psi, &output_bra, pars);

    // quantised $T_rot$
    //apply_rotational_kinetic_energy(psi, &output_bra, pars);

    // quantised $V_{qd}$ in the zeroth order
    apply_charge_dipole_zero(psi, &output_bra, pars);
    //apply_charge_dipole_first(psi, &output_bra, pars);

    //apply_dipole_dipole_zero(psi, &output_bra, pars);

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


/*
The Hamiltonian Matix is expressed in the product basis 
$\{\ket{n1,n3,n5,j1,m1,j2,m2}\}$. Its element H_{bra,ket} 
are scalar products $\bra{...} \hat{H} \ket{...'}$, where 
dots indicate set of quantum numbers. 

The basis is orthonormal, and the Hamiltonian do not mix too many states
(a rule of thumb), thus the matrix H has very a few non-zero elements. 

The non-zero elements are found by decomposing the covector 
$\bra{...}\hat{H}$ in the basis. 

In a loop over all $\bra{...}$ the decomposition of $\bra{...}\hat{H}$ 
is repetedly computed, and non-zero contributions are added to the 
matrix H.
*/
void construct_Hamiltonian(fcomplex* a, const basis b, const parameters pars)
{
        // fill-in with zeroes
        long basis_size = b.n1*b.n3*b.n5*(b.j1*b.j1+2*b.j1+1)*(b.j2*b.j2+2*b.j2+1);
        fcomplex zero = {0.0f, 0.0f};

        for(long i = 0; i<basis_size*basis_size; i++)
                a[i] = zero;

        // scan through the matrix rows and insert all the non-zero elements
        for(int i = 0; i<basis_size; i++)
        {
                versor psi0 = get_versor_from_index(i, b);
                state state0;
                state_init(&state0);
                state_add(&state0, psi0, (fcomplex){1.0f, 0.0f});

                // act with H from the left on the bra of the state state0
                state psiH = bra_H(&state0, pars);

                versor loop_versor;
                fcomplex loop_amplitude;
                for(int l=0; l<psiH.length; l++)
                {
                        loop_versor = state_get_versor(&psiH, l);
                        // it could be that state_get_versor returns something
                        // strange when &psiH[l] is strange.
                        //
                        // Here in an attempt to avoid versors which are 
                        // out of the truncation or don't make sense like
                        // |-23,1,1;1,0,1,0> (-23 doesn't make any sense)
                        // check if the versor is a valid versor
                        if(!valid_versor(loop_versor, b))
                                continue;

                        int loc_idx = get_index_from_versor(loop_versor, b);
                        // check if it is above the diagonal
                        // the input for cheev is a hermitian matrix 
                        // so onnly the part below diagonal matters
                        /*
                        if( loc_idx < i )
                                continue;
                        */

                        // apply amplitude;
                        fcomplex amp = state_get_amplitude(&psiH, l);

                        // add amplitudes
                        a[i*basis_size+loc_idx].re += amp.re;
                        a[i*basis_size+loc_idx].im += amp.im;
                }
               state_free(&psiH);
        }
}