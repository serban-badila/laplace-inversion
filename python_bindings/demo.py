from build.python_bindings.py_laplace import laplace_inversion

from math import e, exp

def exponential_transform(s):
    return 1. / (s + .5)

def exponential(x):
    return pow(e, -.5*x)

if __name__ == "__main__":
    result = laplace_inversion(exponential_transform, .1, 9, 48)

    print(result[0:100])

    expected = [exponential((i) * .1) for i in range(101)] 
    # print(expected) 

    assert all([abs(expected[i+1] - result[i+1]) < pow(10,-12) for i in range(1, 100)]), f"{expected[1]} {result[1]}" 
