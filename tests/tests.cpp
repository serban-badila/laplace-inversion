#include <complex>
#include <vector>

#include <gtest/gtest.h>

#include "laplaceInversion.h"


/* Exponential distributions */

std::complex<double> exp_transform(std::complex<double> s) {
    /*Laplace transform of an exponential distribution.*/
    return .5 / (s + .5);
}

double exponential_df(double x) {
    /* Density function (DF) of an exponential distribution. */
    if (x >= 0.)
        return .5 * std::exp(-.5*x);
    else return 0.;    
}


/* Normal distributions */ 

std::complex<double> normal_transform(std::complex<double> s){
    /* Transform of a standard normal CDF.*/
    return std::exp(.5 * s * s);
}

double centered_normal(double x) {
    /* DF of a standard normal distribution.  */
    return 1./(std::sqrt(2 * PI)) * std::exp(- .5 * x * x);    
}


/* scaled exponential function */

std::complex<double> scaled_exp_transform(std::complex<double> s){
    return std::pow(s + 1., -2);
}

double scaled_exp(double x){
    return x * std::exp(-x);
}


/*  sine function */

std::complex<double> sine_transform(std::complex<double> s){
    return 1. / (s * s + 1.);
}


/* scaled cosine */

std::complex<double> scaled_cosine_transform(std::complex<double> s){
    return std::pow(s * s + 1., -2) * (s * s - 1.);
}


double scaled_cosine(double x){
    return x * std::cos(x);
}


/************** Tests ****************/

class oneDimensionalTest : public ::testing::Test {
    /* Set up the fixtures. */
    protected:
        int n {48};
        double delta {.1};
        unsigned int mexp {15};
        size_t m;
        double err;

        void SetUp() override {
            m = pow(2, mexp);
            err = pow(10, -14);
        }
};


TEST_F(oneDimensionalTest, InvertExponentialTransform){ 

    auto inverse = laplaceInversion::oneDimensionalInverse(exp_transform, delta, mexp, n);
 
    for (size_t j = 1; j < m; j++){
        // set the evaluation points away from the origin because the inverse has a discontinuity there
        EXPECT_NEAR (inverse[j], exponential_df(double(j)*delta), err);
    }
    EXPECT_EQ (inverse[0], .5 * (exponential_df(0) + exponential_df(-err))); // the inverse is "normalized" at discontinuities 
    // in the same way as the inverse Fourier transform.
}



TEST_F(oneDimensionalTest, InvertNormalTransform){

    auto inverse = laplaceInversion::oneDimensionalInverse(normal_transform, delta, mexp, n);

    for (size_t j = 0; j < m; j++){
        EXPECT_NEAR (inverse[j], centered_normal(double(j)*delta), err);
    }
}


TEST_F(oneDimensionalTest, InvertScaledExponentialTransform){

    auto inverse = laplaceInversion::oneDimensionalInverse(scaled_exp_transform, delta, mexp, n);

    for (size_t j = 0; j < m; j++){
        EXPECT_NEAR (inverse[j], scaled_exp(double(j)*delta), err);
    }
}

TEST_F(oneDimensionalTest, InvertSineTransform){

    auto inverse = laplaceInversion::oneDimensionalInverse(sine_transform, delta, mexp, n);

    for (size_t j = 0; j < m; j++){
        EXPECT_NEAR (inverse[j], std::sin(double(j)*delta), err * 100); // losing some precision here; could be the sine implementation ?
    }
}

TEST_F(oneDimensionalTest, InvertScaledCosineTransform){

    auto inverse = laplaceInversion::oneDimensionalInverse(scaled_cosine_transform, delta, mexp, n);

    for (size_t j = 0; j < m; j++){
        EXPECT_NEAR (inverse[j], scaled_cosine(double(j)*delta), err * pow(10, 6)); // also losing precision here; 
        // this is the worst case tolerance among these iterations.
    }
}
