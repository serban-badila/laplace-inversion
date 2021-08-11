import pytest
from python_bindings.demo import laplace_inversion


@pytest.mark.parametrize("rate", [.5, 1., 2., 3., 4., 20.])
def test_exponential_distribution(rate, exponential_distribution, exponential_transform):
    tol = pow(10, -13)
    delta = .05 
    m_exp = 11
    m = pow(2, m_exp)
    result = laplace_inversion(exponential_transform(rate), delta, m_exp, 48)[1:m]
    inverse = exponential_distribution(rate)
    expected  = [inverse(x * delta) for x in range(1, m)]
    assert all([abs(result[i] - expected[i]) < tol for i in range(m-1)])
