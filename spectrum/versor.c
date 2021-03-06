#include <stdio.h>
#include "versor.h"

void show_versor(const versor psi)
{
        printf("|%4d,%2d,%2d;%d,%3d;%d,%3d>",
        psi.n1, psi.n3, psi.n5, psi.j1, psi.m1, psi.j2, psi.m2);
}

void show_bra_versor(const versor psi)
{
        printf("<%4d,%2d,%2d;%d,%3d;%d,%3d|",
        psi.n1, psi.n3, psi.n5, psi.j1, psi.m1, psi.j2, psi.m2);
}

// returns 1 if all the components of two versors are the same
int versor_equals_versor(versor ket0, versor ket1)
{
        int the_same = 1;
        if ( ket0.n1 != ket1.n1 )
                return 0;
        if ( ket0.n3 != ket1.n3 )
                return 0;
        if ( ket0.n5 != ket1.n5 )
                return 0;

        if ( ket0.j1 != ket1.j1 )
                return 0;
        if ( ket0.m1 != ket1.m1 )
                return 0;

        if ( ket0.j2 != ket1.j2 )
                return 0;
        if ( ket0.m2 != ket1.m2 )
                return 0;

        return the_same;
}

/* For tests.
Print "|n1,n3,n5,j1,m1,j2,m2>\t index = idx"
where `idx` is the position of the versor in the 
list of the basis elements. */
void print_versor(const versor psi, const basis b)
{
        show_versor(psi);
        printf("\t\t" "index = %d\n", get_index_from_versor(psi, b));
}


/* This function tells the index of the $\ket{j,m}$ state on the 
list of all rotational states $\ket{J,M}.*/
int jm_jump(int j, int m)
{
        /* returns the number of rotational states behind $\ket{j,m}$ */
        int jump=0;
        for(int i=0; i<j; i++)
                jump += 2*i+1;

        jump += m+j;

        return jump;
}

/* 
We work in a product basis
$$
\ket{\psi} = \ket{n1}\ket{n3}\ket{n5}\ket{j1,m1}\ket{j2,m2}.
$$
This basis is sorted, so that to each versor there corresponds 
an index $idx$. The ordering is possible only thanks to the basis 
is truncated i.e. $n1 < N1, \ldots, j1 <= J1, \ldots, m2 \le M2$.
This function returns the index of a given set of quantum number 
and a truncation.

The product basis is undersood in the following way
$$
\ket{\psi} = \ket{n1}\ket{n3}\ket{n5}\ket{j1,m1}\ket{j2,m2}
$$
order of this basis is as follows:
- vibrations (nx) (b.n1, b.n3, b.n5)
- rotations (j1, m1), (j2, m2)
$$
\sum_{j=0}^{j=b.j1}(2j+1) = (b.j1)^2 + b.j1 + 1
$$ */
int get_index_from_versor(versor psi, const basis b)
{
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

/*
* Returns $j$ or $m$ quantum number of a versor $\ket{j,m}$
* which is in the position `idx` on a list of all
* $\ket{j,m}$. The list is sorted so that $j$ 
* is non-decreasing, and within the same $j$,
* $m$ increases: $\ket{0,0}, \ket{1,-1},
* \ket{1,0}, \ket{1,1}, \ket{2,-2}, \ldots$.
*
The two functions defined below serve as a supporting 
functions for the get_versor_from_index function. 
They output respectively `j` and `m` quantum numbers from 
the versor of an index `idx`. All versors are ordered 
(indexed):
- their `j` number don't decrease, 
- within the same `j`, `m` numbers increase.
It is a unique mapping, and here the maps 
from `idx` to `j,m` are implemented. */
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

/* Returns true if:
 * - any of the normal modes is in a negative excitation
 * - |mx| > jx, for x in 1,2
 * - j1 or j2 is negative
 * 
 * Returns false if:
 * - j1 or j2 is larger than jmax
 */
int versor_is_unphysical(const versor psi)
{
    if ( psi.n1<0 || psi.n3<0 || psi.n5<0 )
        return 1;

    if ( psi.j1 < 0 || psi.j2 < 0 )
        return 1;

    if ( psi.m1 > psi.j1 || psi.m1 < -psi.j1 )
        return 1;

    if ( psi.m2 > psi.j2 || psi.m2 < -psi.j2 )
        return 1;

    return 0;
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