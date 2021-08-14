from py_laplace_inversion import laplace_inversion

from math import e, exp


def exponential_transform(s):
    return .5 / (s + .5)


def exponential(x):
    return .5 * pow(e, -.5*x)


if __name__ == "__main__":
    result = laplace_inversion(exponential_transform, .1, 9, 48)

    print(f"Inverse of an exponential rate={.5} transform:\n", result[0:10])

    expected = [exponential((i) * .1) for i in range(101)] 

    assert all([abs(expected[i+1] - result[i+1]) < pow(10,-12) for i in range(1, 100)]), f"{expected[1]} {result[1]}" 
