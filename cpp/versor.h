#ifndef versor_h
#define versor_h

/* Baisc elements which allow translation between the
* quantum mechanial terms and linear algebra. */

/* Structure designed to store the basis truncation parameters. */
typedef struct {
    int n1, n3, n5, j1, j2;
} basis;

/* Structure allows to call different states which can be 
represented in the product basis */
typedef struct {
    int n1, n3, n5, j1, m1, j2, m2;
} versor;

/* 
We work in a product basis
$$
\ket{\psi} = \ket{n1}\ket{n3}\ket{n5}\ket{j1,m1}\ket{j2,m2}.
$$
This basis is sorted, so that to each versor there corresponds 
an index $idx$. The ordering is possible only because the basis 
is truncated i.e. $n1 < N1, \ldots, j1 <= J1, \ldots, m2 \le M2$.
This function returns the index of a given set of quantum number 
and a truncation.
*/
int get_index_from_versor(versor psi, const basis b);
int jm_jump(int j, int m);
versor get_versor_from_index(int idx, const basis b);
/*
* Returns $j$ or $m$ quantum number of a versor $\ket{j,m}$
* which is in the positin `idx` on a list of all
* $\ket{j,m}$. The list is sorted so that $j$ 
* is non-decreasing, and within the same $j$,
* $m$ is increasing: $\ket{0,0}, \ket{1,-1},
* \ket{1,0}, \ket{1,1}, \ket{2,-2}, \ldots$.
*/
int getj(int idx);
int getm(int idx);

void test_idx_to_versor_translation();

/* prints |n1,n3,n5,j1,m1,j2,m2> */
void show_versor(const versor psi);
void show_bra_versor(const versor psi);
int versor_equals_versor(versor ket0, versor ket1);

#endif // versor_h