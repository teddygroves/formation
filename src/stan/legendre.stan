/*
  Functions for adjusting thermodynamic measurements given conditions
*/

real dh(real temperature, real ionic_strength){
    /* get the Debeye-Hukel number for a temperature and ionic strength */
    real a1 = 9.20483e-3;  // source for these numbers: https://gitlab.com/equilibrator/equilibrator-cache/-/blob/develop/src/equilibrator_cache/thermodynamic_constants.py
    real a2 = 1.284668e-5;
    real a3 = 4.95199e-8;
    real alpha = a1 * temperature + a2 * temperature^2 + a3 * temperature^3;
    real sqrtI = sqrt(ionic_strength);
    return alpha * sqrtI / (1 + 1.6 * sqrtI);
}
real lt(real dgf_std,   // standard condition formation energy of compound
        vector pkas,   // acid dissociation constants
        vector pkmgs,  // magnesium dissociation constants
        int nH,
        real pH,
        real pMg,
        real I,
        real charge,
        real t,
        real R){
    /* 
    Get the condition-specific formation energy of a microspecies using a
    Legendre transform.
    */
    real RT = R * t;
    real dgfmg = -455.3;  
    real dhfmg = -467.0;  // source for these numbers: https://gitlab.com/equilibrator/equilibrator-cache/-/blob/develop/src/equilibrator_cache/thermodynamic_constants.py
    real dgfmg_prime = (t / 298.15) * dgfmg + (1.0 - t / 298.15) * dhfmg;
    int nMg = rows(pkmgs) + 1;
    real ddg_over_rt =
      nH * RT * log(10) * pH
      + nMg * (RT * log(10) * pMg - dgfmg_prime)
      + dh(t, I) * (nH + 4 * nMg - charge^2);
    if (rows(pkas) > 0)
      ddg_over_rt -= sum(log(10) * pkas);
    if (rows(pkmgs) > 0)
      ddg_over_rt += rows(pkmgs) * dgfmg/RT - sum(log(10) * pkmgs);
    return dgf_std + ddg_over_rt * RT;
}

vector ms_pks(int ms, vector pk, int[] fst_pk, int[] n_pk){
  vector[n_pk[ms]] out;
  if (n_pk[ms] > 0) out = pk[fst_pk[ms]:fst_pk[ms]+n_pk[ms]-1];
  return out;
}

real get_ms_dgf(int ms,
                real dgf_cpd,
                int[] fst_pka,
                int[] n_pka,
                int[] fst_pkmg,
                int[] n_pkmg,
                int  nH,
                real charge,
                real I,
                real t,
                real pH,
                real pMg,
                real R,
                vector pka,
                vector pkmg){
  vector[n_pka[ms]] pkas = ms_pks(ms, pka, fst_pka, n_pka);
  vector[n_pkmg[ms]] pkmgs = ms_pks(ms, pkmg, fst_pkmg, n_pkmg);
  return lt(dgf_cpd, pkas, pkmgs, nH, pH, pMg, I, charge, t, R);
}

real get_dgr_prime(vector stoic_coefs,  // stoic of each reactant in the reaction
                    int [] rxn_cpds,     // compound of each reactant in the reaction 
                    int [] n_ms,         // number of microspecies for every compound
                    int [] fst_ms,       // first microspecies for every compound
                    int [] n_pka,        // number of pkas for every microspecies
                    int [] fst_pka,      // first value in pka for every microspecies 
                    int [] n_pkmg,       // number of pkmgs for every microspecies
                    int [] fst_pkmg,     // first value in pkmg for every microspecies
                    int [] nH,           // number of H+ ions for every microspecies
                    vector z,            // charge for each microspecies
                    real I,              // ionic strength
                    real t,              // temperature
                    real pH,
                    real pMg,
                    real R,
                    vector pka,
                    vector pkmg,
                    vector dgf){
/* Get the Gibbs energy change of a reaction in given conditions. */
  real out = 0;
  real RT = R * t;
  for (i in 1:size(rxn_cpds)){
    int cpd = rxn_cpds[i];
    real dgf_prime;
    vector[n_ms[cpd]] dgf_primes;
    for (n in 1:n_ms[cpd]){
      int ms = fst_ms[cpd] + n - 1;
      vector[n_pka[ms]] pkas = ms_pks(ms, pka, fst_pka, n_pka);
      vector[n_pkmg[ms]] pkmgs = ms_pks(ms, pkmg, fst_pkmg, n_pkmg);
      dgf_primes[n] = lt(dgf[cpd], pkas, pkmgs, nH[ms], pH, pMg, I, z[ms], t, R);
    }
    dgf_prime = -RT * log_sum_exp(-1 * dgf_primes / RT);
    out += stoic_coefs[i] * dgf_prime;
  }
  return out;
}
