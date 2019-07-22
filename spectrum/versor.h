#ifndef versor_h
#define versor_h

/* Baisc elements which allow translation between the
* quantum mechanial terms and linear algebra. */

/* Structure designed to store the basis truncation parameters. 
 * The numbers nx, x=1,3,5 denote number of modes allowed i.e. basis.n1=1 
 * means that only one mode can be excited, the mode withe versor.n1=0. 
 */
typedef struct {
    int n1, n3, n5, j1, j2;
} basis;

/* Structure allows to call different states which can be 
represented in the product basis */
typedef struct {
    int n1, n3, n5, j1, m1, j2, m2;
} versor;

int get_index_from_versor(versor psi, const basis b);
int jm_jump(int j, int m);
versor get_versor_from_index(int idx, const basis b);
int getj(int idx);
int getm(int idx);

void test_idx_to_versor_translation();

/* prints |n1,n3,n5,j1,m1,j2,m2> */
void show_versor(const versor psi);
void show_bra_versor(const versor psi);
int versor_equals_versor(versor ket0, versor ket1);
int versor_is_unphysical(const versor psi);

#endif // versor_h