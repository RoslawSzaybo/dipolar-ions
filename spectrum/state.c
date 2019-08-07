#include "state.h"
#include <stdlib.h>
#include <stdio.h>

void state_init(state* psi)
{
    psi->kets = NULL;
    psi->amplitudes = NULL;
    psi->length = 0;
    psi->cap = 0;
}

void state_free(state* psi)
{
    free(psi->kets);
    free(psi->amplitudes);
    psi->length = 0;
    psi->cap = 0;
}

int state_cnt(state *psi)
{
    return psi->length;
}

/*
* computes a sum of all amplitudes square. This is known as the
* normalisation factor. The normalisation factor should be equal
* to unity in for all states that we consider here.
*
* If the normalisation facor is different from unity, then the state
* sould be divided by the sqrt(get_state_normalisation).
*/
double get_state_normalisation(state *psi)
{
    double N2 = 0.;

    // we add first the small amplitudes
    // we assume that the state is sorted i.e. 
    // state_sort(&psi);
    // but if it is not then nothing happens
    // sort would give a better numerical stability
    int i=psi->length-1;
    for (;i>=0; i--)
        N2 += dcomplex_amplitude_sqr(&psi->amplitudes[i]);

    return N2;
}

// returns 1 when the state `psi` has as its member versor 
// a versor exactly the same as `ket`
int state_contains_versor(state *psi, const versor *ket)
{
    int contains = 0;
    int i;
    for (i=0; i < psi->length; i++)
    {
        if ( versor_equals_versor(*ket, psi->kets[i]) )
            return 1;
    }
    return 0;
} 

/* 
* Rearranges all of the versors and amplitudes stored in the local arrays 
* of the state so that their amplitude don't decrease i.e.
*  for each 0 <= i <= j < psi->length
*  |psi->amplitudes[i]|^2 >= |psi->amplitudes[j]|^2
*/
void state_sort(state *psi)
{
    versor tmp_ket;
    dcomplex tmp_amp;

    // assume that the state is not sored
    int is_sorted = 0;

    int i;
    double amp0, amp1;
    while (!is_sorted)
    {
        // suppose that after some changes the state is sorted
        is_sorted = 1;

        for (i=0; i<psi->length-1; i++)
        {
            amp0 = dcomplex_amplitude_sqr(&psi->amplitudes[i]);
            amp1 = dcomplex_amplitude_sqr(&psi->amplitudes[i+1]);

            // check if for all pairs amplitudes do not increase
            if (amp1 > amp0) 
            {
                // swap positions so that the state is ordered 
                is_sorted = 0;
                // replace versors
                tmp_ket = psi->kets[i];
                psi->kets[i] = psi->kets[i+1];
                psi->kets[i+1] = tmp_ket;

                // replace amplitudes
                tmp_amp = psi->amplitudes[i];
                psi->amplitudes[i] = psi->amplitudes[i+1];
                psi->amplitudes[i+1] = tmp_amp;
            }
        }
    }
}

/*
* Removes from the state an element which is stored 
* in the local arrays of the state under the index idx.
*/
void state_rm(state *psi, const int idx)
{
    // indices of the local arrays of the states are
    // non-negative integers and are no larger than
    // the number of stored elements
    if ( idx < 0 || idx >= psi->length )
    {
        printf("Error in state_rm(state *psi, %d).\n", idx);
        printf("Requested index do not point to any element.\n");
        exit(0);
    }

    // we shift all the vectors with index larger than idx 
    // to a slot in the local array of in index smaller by one

    // versors
    int i = idx;
    for(; i < psi->length-1; i++)
    {
        psi->kets[i] = psi->kets[i+1];
    }

    // amplitudes
    i = idx;
    for(; i < psi->length-1; i++)
    {
        psi->amplitudes[i] = psi->amplitudes[i+1];
    }

    // we forget about the last element
    //psi->kets[psi->length-1] = NULL;
    //psi->amplitudes[psi->length-1] = NULL;

    // the number of stored elements decreased
    psi->length--;
}

// returns the location of a versor `ket` on the list
// of all member versors of the state `psi`
// WARNING:
// the code assumes that the versor `ket` is a member
// versor of the state `psi`. 
int state_versor_location(state *psi, const versor *ket)
{
    int i;
    for (i=0; i < psi->length; i++)
    {
        if ( versor_equals_versor(*ket, psi->kets[i]) )
            return i;
    }

    printf("ERROR!");
    printf("A versor location was not found.");
    printf("state_versor_location() failed.");
    exit(0);
}

void state_add(state* psi, const versor ket, const dcomplex amplitude)
{
    // the versor `ket` is already present in the state;
    // it is enough to sum amplitudes
    if ( state_contains_versor(psi, &ket) ) 
    {
        int loc = state_versor_location(psi, &ket);
        psi->amplitudes[loc].re += amplitude.re;
        psi->amplitudes[loc].im += amplitude.im;
    }
    // otherwise a new versor has to be added to the state
    else
    {
        // `ket` is the first versos in the state
        if (psi->length == 0)
        {
            psi->cap = 25;
            psi->kets = (versor*)malloc(sizeof(versor)*psi->cap);
            psi->amplitudes = (dcomplex*)malloc(sizeof(dcomplex)*psi->cap);
        }
        // check if the extra versor can fit into the state
        if (psi->length == psi->cap)
        {
            psi->cap *= 2;
            psi->kets = (versor*)realloc(psi->kets, sizeof(versor)*psi->cap);
            if(psi->kets == NULL)
            {
                printf("The state cannot fit into the memory.\n");
                printf("The requested number of versors: %d\n", psi->cap);
                exit(0);
            }
            psi->amplitudes = (dcomplex*)realloc(psi->amplitudes, sizeof(dcomplex)*psi->cap);
            if(psi->kets == NULL)
            {
                printf("The state cannot fit into the memory.\n");
                printf("The requested number of versors: %d\n", psi->cap);
                exit(0);
            }
        }
        psi->kets[psi->length] = ket;
        psi->amplitudes[psi->length].re = amplitude.re;
        psi->amplitudes[psi->length].im = amplitude.im;
        psi->length++;
    }
}

versor state_get_versor(state* psi, int idx)
{
    if (idx >= psi->length || idx <0)
    {
        //return NULL;
        printf("Error!");
        printf("state_get_amplitude is attmeptig to"
        " reach index larger than the state length or"
        " a negative index.\n");
        exit(0);
    }
    
    return psi->kets[idx];
}

dcomplex state_get_amplitude(state* psi, int idx)
{
    if (idx >= psi->length || idx <0)
    {
        //return NULL;
        printf("Error!");
        printf("state_get_amplitude is attmeptig to"
        " reach index larger than the state length or"
        " a negative index.\n");
        exit(0);
    }

    return psi->amplitudes[idx];
}

void state_times_double(state *input, const double factor)
{
    for (int l=0; l < input->length; l++)
    {
        input->amplitudes[l].re *= factor;
        input->amplitudes[l].im *= factor;
    }
}

void state_times_i_double(state *input, const double factor)
{
    // (re + i im) * i factor = - im*factor + i re*factor
    double re;
    double im;

    for(int l=0; l < input->length; l++)
    {
        re = input->amplitudes[l].re;
        im = input->amplitudes[l].im;

        input->amplitudes[l].re = - im*factor;
        input->amplitudes[l].im =   re*factor;
    }
}


void state_add_state(state *sum, const state *component)
{
    for (int l=0; l<component->length; l++)
        state_add(sum, component->kets[l], component->amplitudes[l]);
}

/*
* This function removes the from the state all the versors which possess:
* - negative excitations of modes
* - |mx| larger than jx
* - negative j
*
* Warrning!
* This function do not remove versors with angular number j larger than 
* the maximal implemented jmax.
*/
void state_clean_unphysical_versors(state *psi)
{
    /* it is faster to rm a versor from the end of the state */
    int i = psi->length-1;

    versor subject;
    for (; i>= 0; i--)
    {
        subject = psi->kets[i];

        if (versor_is_unphysical(subject))
            state_rm(psi, i);
    }
}

/* 
 * This function removes from the state all the versors which have 
 * the square of its amplitude smaller than the threshold value.
 */
void state_cut(state *psi, double threshold)
{
    /* Rm works faster when applied from the end. */
    int i = psi->length-1;
    double amp = 0.;

    for (; i>=0; i--)
    {
        // square of an amplitude of the versor in a state
        amp = dcomplex_amplitude_sqr(&psi->amplitudes[i]);
        if ( amp < threshold )
            state_rm(psi, i);
    }
}

/* 
 * This function keeps only the first `max` versors of the sorted state
 */
void state_keep_only_first_max_versors(state *psi, int max)
{
    // the state is short enough
    if ( psi->length <= max )
        return;

    // otherwise we cut
    /* Rm works faster when applied from the end. */
    int i = psi->length-1;
    for (; i>=max; i--)
        state_rm(psi, i);
}