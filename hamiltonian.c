/*
this is a file which used to serve as an example on how to use the
cheev function. 

*/

#include <stdlib.h>
#include <stdio.h>
#include "hamiltonian.h"

int jm_jump(int j, int m)
{
        /* returns the number of rotational states behind $\ket{j,m}$ */
        int jump=0;
        for(int i=0; i<j; i++)
                jump += 2*i+1;

        jump += m+j;

        return jump;
}

int get_index_from_state(state psi, const basis b)
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

state get_state_from_index(int idx, const basis b)
{
        int no_vib35 = b.n3 * b.n5;
        int no_vib5 = b.n5;
        int no_j1_rot = b.j1*b.j1 + 2*b.j1 + 1;
        int no_j2_rot = b.j2*b.j2 + 2*b.j2 + 1;
        int no_rot = no_j1_rot * no_j2_rot; 

        state psi = {0, 0, 0, 0, 0, 0, 0};

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

void print_state(const state psi, const basis b)
{
        printf("State n1=%d, n3=%d, n5=%d, j1=%d, m1 =%2d, j2=%d, m2=%2d\t\t"
        "index = %d\n",
        psi.n1, psi.n3, psi.n5, psi.j1, psi.m1, psi.j2, psi.m2,
        get_index_from_state(psi, b));
}


/* test */
void test_idx_to_state_translation() {
        const basis b = {1000, 2, 3, 0, 0};

        for(int idx=0; idx<20; idx++)
        {
                state psi1 = get_state_from_index(idx, b);
                print_state(psi1, b);
        }
}

