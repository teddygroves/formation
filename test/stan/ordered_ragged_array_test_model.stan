functions {
#include ordered_ragged_array.stan
}
data {
  int N;
  int O;
  vector[N] first_vals;
  int n_elems[N];
  vector[O-N] diffs;
}
generated quantities {
  vector[O] actual = ordered_ragged_array(first_vals, n_elems, diffs);
}
