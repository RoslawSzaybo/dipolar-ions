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
            psi->amplitudes = (fcomplex*)realloc(psi->amplitudes, sizeof(fcomplex)*psi->cap);
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
