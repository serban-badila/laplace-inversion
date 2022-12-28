import pytest 

from math import exp

@pytest.fixture()
def exponential_transform():
    def _transform(rate):
        def _inner(s):
            return rate / (s + rate)
        return _inner
    return _transform

@pytest.fixture()
def exponential_density():
    def _density(rate):
        def _inner(x): 
            return rate * exp( -rate * x) if x >= 0 else 0
        return _inner
    return _density
