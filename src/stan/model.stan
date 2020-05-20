/*
  A joing model of K prime measurements, PKa measurements and PKMg
  measurements, ignoring any group contribution.
*/

functions {
#include legendre.stan
#include ordered_ragged_array.stan
}
data {
  int<lower=1> N_measurement_kpr;
  int<lower=1> N_measurement_pka;
  int<lower=1> N_measurement_pkmg;
  int<lower=1> N_microspecies;
  int<lower=1> N_pka;
  int<lower=1> N_pkmg;
  int<lower=1> N_compound;
  int<lower=1> N_stoic;
  int<lower=1> N_reaction;
  int<lower=1> N_Hable_compound;
  int<lower=1> N_Mgable_compound;
  int<lower=1> n_cpd[N_reaction];
  vector[N_microspecies] charge;
  int<lower=0> nH[N_microspecies];
  int<lower=0> nMg[N_microspecies];
  int<lower=1,upper=N_compound> stoic_cpd[N_stoic];
  vector[N_stoic] stoic_coef;
  int<lower=1,upper=N_stoic> fst_stoic[N_reaction];
  int<lower=0,upper=N_microspecies> n_pka[N_microspecies];
  int<lower=1> n_pka_cpd[N_Hable_compound];
  int<lower=0,upper=N_pka> fst_pka[N_microspecies];
  int<lower=0,upper=N_microspecies> n_pkmg[N_microspecies];
  int<lower=1> n_pkmg_cpd[N_Mgable_compound];
  int<lower=0,upper=N_pkmg> fst_pkmg[N_microspecies];
  int<lower=1,upper=N_microspecies> n_ms[N_compound];
  int<lower=1,upper=N_microspecies> fst_ms[N_compound];
  int<lower=1,upper=N_microspecies> pka_obs_ix[N_measurement_pka];
  vector<lower=0>[N_measurement_pka] pka_obs;
  int<lower=1,upper=N_microspecies> pkmg_obs_ix[N_measurement_pkmg];
  vector<lower=0>[N_measurement_pkmg] pkmg_obs;
  vector[N_measurement_kpr] kpr_obs;
  int rxn[N_measurement_kpr];
  vector[N_measurement_kpr] temperature;
  vector[N_measurement_kpr] I;
  vector[N_measurement_kpr] pH;
  vector[N_measurement_kpr] pMg;
  vector[N_compound] prior_loc_dgf;
  int<lower=1,upper=2> prior_regime_dgf[N_compound];
}
transformed data {
  real R = 8.314e-3;
  vector[N_measurement_kpr] RT = R * temperature;
  vector[N_measurement_kpr] dgr_obs = -RT .* log(kpr_obs);
  real<lower=0> sigma_dgr = 10;
  real<lower=0> sigma_pka = 0.2;
  real<lower=0> sigma_pkmg = 0.2;
  vector<lower=0>[2] sigma_dgf = [30, 500]';
  int<lower=1> N_pka_diff = N_pka - N_Hable_compound;
  int<lower=1> N_pkmg_diff = N_pkmg - N_Mgable_compound;
}
parameters {
  vector[N_compound] dgf_z;
  vector[N_Hable_compound] first_pka;
  vector[N_Mgable_compound] first_pkmg;
  vector<lower=0>[N_pka_diff] pka_diffs;
  vector<lower=0>[N_pkmg_diff] pkmg_diffs;
}
transformed parameters {
  vector[N_pka] pka = ordered_ragged_array(first_pka, n_pka_cpd, pka_diffs);
  vector[N_pkmg] pkmg = ordered_ragged_array(first_pkmg, n_pkmg_cpd, pkmg_diffs);
  vector[N_compound] dgf = prior_loc_dgf + dgf_z .* sigma_dgf[prior_regime_dgf];
  vector[N_measurement_kpr] dgr_hat;
  for (m in 1:N_measurement_kpr) {
    int first = fst_stoic[rxn[m]];
    int last = first + n_cpd[rxn[m]] - 1;
    dgr_hat[m] = get_dgr_prime(stoic_coef[first:last],
                               stoic_cpd[first:last],
                               n_ms, fst_ms,
                               n_pka, fst_pka,
                               n_pkmg, fst_pkmg,
                               nH,
                               charge,
                               I[m], temperature[m], pH[m], pMg[m], R,
                               pka, pkmg, dgf);
  }
}
model {
  // priors
  target += std_normal_lpdf(dgf_z|);
  target += normal_lpdf(first_pka | 7, 3);
  target += normal_lpdf(first_pkmg | 5, 2);
  target += normal_lpdf(pka_diffs | 3, 1);
  target += normal_lpdf(pkmg_diffs | 3, 1);
  // likelihood
  target += normal_lpdf(dgr_obs | dgr_hat, sigma_dgr);
  target += normal_lpdf(pka_obs | pka[pka_obs_ix], sigma_pka);
  target += normal_lpdf(pkmg_obs | pkmg[pkmg_obs_ix], sigma_pkmg);
}
