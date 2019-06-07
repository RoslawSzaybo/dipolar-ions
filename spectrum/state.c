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
    free (psi->kets);
    free (psi->amplitudes);
    psi->length = 0;
    psi->cap = 0;
}

int state_cnt(state *psi)
{
    return psi->length;
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
// the code assumes that the versos `ket` is a member
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
    exit(1);
}

void state_add(state* psi, const versor ket, const fcomplex amplitude)
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
            psi->amplitudes = (fcomplex*)malloc(sizeof(fcomplex)*psi->cap);
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
            psi->amplitudes = (fcomplex*)realloc(psi->amplitudes, sizeof(fcomplex)*psi->cap);
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
    if(idx >= psi->length)
    {
        //return NULL;
        printf("Error!");
        printf("state_get_amplitude is attmeptig to"
        " reach index larger than the state length.");
    }
    
    return psi->kets[idx];
}

fcomplex state_get_amplitude(state* psi, int idx)
{
    if(idx >= psi->length)
    {
        //return NULL;
        printf("Error!");
        printf("state_get_amplitude is attmeptig to"
        " reach index larger than the state length.");
    }

    return psi->amplitudes[idx];
}

void state_times_float(state * input, const float factor)
{
    for(int l=0; l < input->length; l++)
    {
        input->amplitudes[l].re *= factor;
        input->amplitudes[l].im *= factor;
    }
}

void state_times_i_float(state * input, const float factor)
{
    // (re + i im) * i factor = - im*factor + i re*factor
    float re;
    float im;

    for(int l=0; l < input->length; l++)
    {
        re = input->amplitudes[l].re;
        im = input->amplitudes[l].im;

        input->amplitudes[l].re = - im*factor;
        input->amplitudes[l].im =   re*factor;
    }
}


void state_add_state(state * sum, const state * component)
{
    for(int l=0; l<component->length; l++)
        state_add(sum, component->kets[l], component->amplitudes[l]);
}

void state_clean(state* psi, basis b)
{
    state *clean;
    state_init(clean);
}
