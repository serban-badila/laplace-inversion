#include <vector>
#include <complex>

#include <fftw3.h>

#define PI 3.141592653589793238462643383279502884
#define OVRSMPL 8

namespace laplaceInversion {

void loadGaussQuadRule (const int , std::vector<double>& , std::vector<double>& );


std::vector<double> oneDimensionalInverse(
    std::complex<double> (*f_hat)(std::complex<double>), 
    double delta, 
    unsigned int mexp,
    int n,
    bool use_generic
);


std::vector<double> __oneDimensionalInverse(
    std::complex<double> (*f_hat)(std::complex<double>), 
    double delta, 
    int mexp, 
    int n, 
    std::vector<double> &lambda, 
    std::vector<double> &beta
)
    /*
    Private function. 

    Args:
        delta: mesh size;
    
    */
{       
    // approximate the initial transform with its quadrature
    std::complex<double> I {0.0 , 1.0}; // complex literals require >=c++14
    unsigned int m {static_cast<unsigned int>(pow(2, mexp))};
    unsigned int m2 {OVRSMPL * m};
    std::complex<double> quad_approx[m2+1];
    double a { 44.0 / static_cast<double> (m2) };

    std::complex<double> temp;
    for (int k=0; k <= m2; k++){
        for (int j = 1; j <= n/2; j++){
            temp = f_hat((a + lambda[j] * I + 2.0 * PI * k/m2 * I)/delta);
            quad_approx[k] += beta[j] * temp; 
        }
        quad_approx[k] *= 2.0 / delta;
    }
    quad_approx[0] = .5 * (quad_approx[0] + quad_approx[m2]);

    // now run the inverse Fast Fourier transform on the quadrature
    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) (quad_approx); // the last index will not be used in the fft
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m2);

    p = fftw_plan_dft_1d(m2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);    
    
    std::complex<double>* inverse = (std::complex<double>*) *out; // safe casting to complex<double> is guaranteed

    std::vector<double> v_inv_d(m2);
    for (int k=0; k < m; k++){
        v_inv_d[k] = exp(a*k)/static_cast<double>(m2) * inverse[k].real(); // the ifft doesn't normalize the inverse with <m2> so do it here
    }

    fftw_destroy_plan(p);
    fftw_free(out);

    return v_inv_d;
}


void loadGaussQuadRule(const int n, std::vector<double>& beta, std::vector<double>& lambda)
{   
    switch (n)
    {
    case 16:
        lambda.resize(8); beta.resize(8);
        lambda[0] = 4.44089209850063e-016; beta[0] = 1.00000000000000;
        lambda[1] = 6.28318530717958;      beta[1] = 1.00000000000004;
        lambda[2] = 12.5663706962589;      beta[2] = 1.00000015116847;
        lambda[3] = 18.8502914166954;      beta[3] = 1.00081841700481;
        lambda[4] = 25.2872172156717;      beta[4] = 1.09580332705189;
        lambda[5] = 34.296971663526;       beta[5] = 2.00687652338724;
        lambda[6] = 56.1725527716607;      beta[6] = 5.94277512934943;
        lambda[7] = 170.533131190126;      beta[7] = 54.9537264520382;
        break;
    case 32:
        lambda.resize(16); beta.resize(16);
        lambda[0]  = 0;                beta[0]  = 1.00000000000000;
        lambda[1]  = 6.28318530717958; beta[1]  = 1.00000000000000;
        lambda[2]  = 12.5663706143592; beta[2]  = 1.00000000000000;
        lambda[3]  = 18.8495559215388; beta[3]  = 1.00000000000000;	
        lambda[4]  = 25.1327412287184; beta[4]  = 1.00000000000000;
        lambda[5]  = 31.4159265359035; beta[5]  = 1.00000000000895;
        lambda[6]  = 37.6991118820067; beta[6]  = 1.00000004815464;
        lambda[7]  = 43.9823334683971; beta[7]  = 1.00003440685547;
        lambda[8]  = 50.2716029125234; beta[8]  = 1.00420404867308;
        lambda[9]  = 56.7584358919044; beta[9]  = 1.09319461846681;
        lambda[10] = 64.7269529917882; beta[10] = 1.51528642466058;
        lambda[11] = 76.7783110023797; beta[11] = 2.4132076646714;
        lambda[12] = 96.7780294888711; beta[12] = 4.16688127092229;
        lambda[13] = 133.997553190014; beta[13] = 8.3777001312961;
        lambda[14] = 222.527562038705; beta[14] = 23.6054680083019;
        lambda[15] = 669.650134867713; beta[15] = 213.824023377988;
        break;
    case 48:
        lambda.resize(24); beta.resize(24);
        lambda[1]  = 0;                beta[1]  = 1.00000000000000;
        lambda[2]  = 6.28318530717957; beta[2]  = 1.00000000000000;
        lambda[3]  = 12.5663706143592; beta[3]  = 1.00000000000000;  
        lambda[4]  = 18.8495559215388; beta[4]  = 1.00000000000000;  	
        lambda[5]  = 25.1327412287183; beta[5]  = 1.00000000000000;
        lambda[6]  = 31.4159265358979; beta[6]  = 1.00000000000000;
        lambda[7]  = 37.6991118430775; beta[7]  = 1.00000000000000;
        lambda[8]  = 43.9822971502571; beta[8]  = 1.00000000000000;
        lambda[9]  = 50.2654824574367; beta[9]  = 1.00000000000000;
        lambda[10] = 56.5486677646182; beta[10] = 1.00000000000234;
        lambda[11] = 62.8318530747628; beta[11] = 1.00000000319553;
        lambda[12] = 69.1150398188909; beta[12] = 1.00000128757818;
        lambda[13] = 75.3984537709689; beta[13] = 1.00016604436873;
        lambda[14] = 81.6938697567735; beta[14] = 1.00682731991922;
        lambda[15] = 88.1889420301504; beta[15] = 1.08409730759702;
        lambda[16] = 95.7546784637379; beta[16] = 1.3631917322868;
        lambda[17] = 105.767553649199; beta[17] = 1.85773538601497;
        lambda[18] = 119.58751936774;  beta[18] = 2.59022367414073;
        lambda[19] = 139.158762677521; beta[19] = 3.73141804564276;
        lambda[20] = 168.156165377339; beta[20] = 5.69232680539143;
        lambda[21] = 214.521886792255; beta[21] = 9.54600616545647; 
        lambda[22] = 298.972429369901; beta[22] = 18.8912132110256;
        lambda[23] = 497.542914576338; beta[23] = 52.7884611477405; 
        lambda[24] = 1494.71066227687; beta[24] = 476.448331869636; 
        break;
    default:
        throw std::invalid_argument("Uspecified quadrature size. Choose from 16, 32 and 48");
    }
}


std::vector<double> oneDimensionalInverse(
    std::complex<double> (*f_hat)(std::complex<double>), 
    double delta, 
    unsigned int mexp,
    int n
)
{
    /*
    Compute the inverse of a one-dimensional Laplace-Stieltjes Transform using den Iseger's algorithm.

    Args:
        f_hat: function pointer to the Laplace tansform to be inverted.
        mexp: the exponent defining the number of evaluation points (2^mexp) for the inverse where precision is guaranteed; 
            The transform evaluation uses an oversampling factor, <OVRSMPL> * 2^(mexp).
        delta: mesh size; 
        mexp: The number of evaluation/sampling points is 2^(mexp); Applies to the inverse; 
        n: The size of the Gauss quadrature; Choose between 16, 32 and 48; Defaults to 16. 
    
    Returns:
            A std::vector<double> containing the inverse function evaluated at {0, delta, 2*delta, ...}
    */

    int quad_size { int(n/2) };
    std::vector<double> beta;
    std::vector<double> lambda;
    loadGaussQuadRule(n, beta, lambda);
    std::vector<double> inverse_array;


    inverse_array = __oneDimensionalInverse(f_hat, delta, mexp, n, lambda, beta);
    return inverse_array;
}

} // namespace laplaceInversion