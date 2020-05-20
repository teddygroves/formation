from cmdstanpy import CmdStanModel
import os


def test_ordered_ragged_array():
    """Test that the ordered_ragged_array function behaves as expected.
    
    The test ragged array: [[2], [3, 2], [5, 3, 2]]
    """
    here = os.path.dirname(os.path.realpath(__file__))
    stan_file = os.path.join(here, "stan/ordered_ragged_array_test_model.stan")
    include_path = os.path.join(here, "../src/stan")
    data = {
        "N": 3,
        "O": 6,
        "first_vals": [2, 3, 5],
        "n_elems": [1, 2, 3],
        "diffs": [1, 2, 1],
    }
    expected = [2., 3., 2., 5., 3., 2.]
    model = CmdStanModel(
        stan_file=stan_file, stanc_options={"include_paths": [include_path]}
    )
    samples = model.sample(data=data, fixed_param=True, iter_sampling=1)
    actual = samples.get_drawset(params=["actual"]).loc[0].tolist()
    assert actual == expected
