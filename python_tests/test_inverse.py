import pytest
from py_laplace_inversion import laplace_inversion


@pytest.mark.parametrize("rate", [.5, 1., 2., 3., 4., 20.])
def test_exponential_density(rate, exponential_density, exponential_transform):
    tol = pow(10, -13)
    delta = .05 
    m_exp = 11
    m = pow(2, m_exp)
    
    result = laplace_inversion(exponential_transform(rate), delta, m_exp, 48)[:m]
    inverse = exponential_density(rate)
    expected_values  = [inverse(x * delta) for x in range(m)]
    
    assert all([abs(result[i] - expected_values[i]) < tol for i in range(1, m)])
    assert abs(result[0] - .5 * expected_values[0]) < tol
