vector ordered_ragged_array(vector first_vals,  // first value of each subarray
                            int[] n_elems,      // number of elements of each subarray
                            vector diffs){      // differences between subarray elements
  /* A vector containing the values of a ragged array each of whose subarrays
   * are in descending order. 
   *
   * The reason for doing this is that the ordering constraint can be enforced
   * by ensuring that the elements of `diffs` are all positive.
   *
   */
  if (rows(first_vals) != size(n_elems)) {
    reject("Different number of first vals and element counts");
  }
  int N = rows(first_vals);  // length of the array
  int O = sum(n_elems);      // length of the output
  vector[O] out;
  int i_out = 1;
  int i_diff = 1;
  for (n in 1:N){
    real first_val = first_vals[n];
    int n_elem = n_elems[n];
    out[i_out] = first_val;
    if (n_elem > 1){
      vector[n_elem-1] row_diffs = diffs[i_diff:i_diff+n_elem-2];
      out[i_out+1:i_out+n_elem-1] = first_val - cumulative_sum(row_diffs);
      i_diff += n_elem - 1;
    }
    i_out += n_elem;
  }
  return out;
}
