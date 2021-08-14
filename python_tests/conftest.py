import pytest 

from math import exp

@pytest.fixture()
def exponential_transform():
    def _transform(rate):
        def __inner(s):
            return rate / (s + rate)
        return __inner
    return _transform

@pytest.fixture()
def exponential_density():
    def _density(rate):
        def __innner(x): 
            return rate * exp( -rate * x) if x >= 0 else 0
        return __innner
    return _density
