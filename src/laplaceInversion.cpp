#include "laplaceInversion.h"

namespace laplaceInversion {

std::vector<double> __ifft(std::vector<std::complex<double>>& input, size_t len, double a)
    /* Wrap the Fast Fourier Transform Cxx routine. */
{   
    size_t m = len / OVRSMPL;
    pocketfft::shape_t shape{len};
    pocketfft::stride_t strided(shape.size());
    size_t tmpd=sizeof(std::complex<double>);
    
    for (int i = shape.size() - 1; i >= 0; i--)
      {
      strided[i] = tmpd;
      tmpd *= shape[i];
      }
    size_t ndata {1};
    for (size_t i=0; i<shape.size(); i++)
      ndata*=shape[i];
    
    pocketfft::shape_t axes;
    for (size_t i = 0; i < shape.size(); i++)
      axes.push_back(i);
    
    std::vector<std::complex<double>> output(ndata);
    pocketfft::c2c(shape, strided, strided, axes, pocketfft::BACKWARD, input.data(), output.data(), 1.);

    std::vector<double> out (m);
    for (size_t k = 0; k < m; k++){
        out[k] = exp(a*k) / len * output[k].real(); // the ifft doesn't normalize the inverse with <len> so do it here
    }
    return out;
}

std::vector<double> __oneDimensionalInverse(
    std::function<std::complex<double>(std::complex<double>)> f_hat,
    double delta, 
    int mexp, 
    unsigned int n, 
    std::vector<double> &lambda, 
    std::vector<double> &beta
)
    /*
    Private function. 

    Args:
        mexp: the exponent defining the number of evaluation points (2^mexp).
        delta: mesh width.
        n: The size of the Gauss quadrature;
    
    */
{       
    // approximate the initial transform with its quadrature
    std::complex<double> I {0.0 , 1.0}; // complex literals require >=c++14
    size_t m {static_cast<size_t> (pow(2, mexp))};
    size_t m2 {OVRSMPL * m};
    std::vector<std::complex<double>> quad_approx (m2+1);
    double a { 44.0 / static_cast<double> (m2) };

    std::complex<double> temp;
    for (size_t k = 0; k <= m2; k++){
        for (size_t j = 0; j < n/2; j++){
            temp = f_hat((a + lambda[j] * I + 2.0 * PI * k/m2 * I)/delta);
            quad_approx[k] += beta[j] * temp; 
        }
        quad_approx[k] *= 2.0 / delta;
    }
    quad_approx[0] = .5 * (quad_approx[0] + quad_approx[m2]);

    // now run the inverse Fast Fourier transform on the quadrature
    auto output = __ifft(quad_approx, m2, a); // the last index of <quad_approx> will not be used in the ifft
    
    return output;
}


void __loadGaussQuadRule(const int n, std::vector<double>& beta, std::vector<double>& lambda)
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
        lambda[0]  = 0;                beta[0]  = 1.00000000000000;
        lambda[1]  = 6.28318530717957; beta[1]  = 1.00000000000000;
        lambda[2]  = 12.5663706143592; beta[2]  = 1.00000000000000;  
        lambda[3]  = 18.8495559215388; beta[3]  = 1.00000000000000;  	
        lambda[4]  = 25.1327412287183; beta[4]  = 1.00000000000000;
        lambda[5]  = 31.4159265358979; beta[5]  = 1.00000000000000;
        lambda[6]  = 37.6991118430775; beta[6]  = 1.00000000000000;
        lambda[7]  = 43.9822971502571; beta[7]  = 1.00000000000000;
        lambda[8]  = 50.2654824574367; beta[8]  = 1.00000000000000;
        lambda[9]  = 56.5486677646182; beta[9]  = 1.00000000000234;
        lambda[10] = 62.8318530747628; beta[10] = 1.00000000319553;
        lambda[11] = 69.1150398188909; beta[11] = 1.00000128757818;
        lambda[12] = 75.3984537709689; beta[12] = 1.00016604436873;
        lambda[13] = 81.6938697567735; beta[13] = 1.00682731991922;
        lambda[14] = 88.1889420301504; beta[14] = 1.08409730759702;
        lambda[15] = 95.7546784637379; beta[15] = 1.3631917322868;
        lambda[16] = 105.767553649199; beta[16] = 1.85773538601497;
        lambda[17] = 119.58751936774;  beta[17] = 2.59022367414073;
        lambda[18] = 139.158762677521; beta[18] = 3.73141804564276;
        lambda[19] = 168.156165377339; beta[19] = 5.69232680539143;
        lambda[20] = 214.521886792255; beta[20] = 9.54600616545647; 
        lambda[21] = 298.972429369901; beta[21] = 18.8912132110256;
        lambda[22] = 497.542914576338; beta[22] = 52.7884611477405; 
        lambda[23] = 1494.71066227687; beta[23] = 476.448331869636; 
        break;
    default:
        throw std::invalid_argument("Uspecified quadrature size. Choose from 16, 32 and 48");
    }
}


std::vector<double> oneDimensionalInverse(
    std::function<std::complex<double>(std::complex<double>)> f_hat,
    double delta, 
    unsigned int mexp,
    int n = 48
)
{
    /*
    Compute the inverse of a one-dimensional Laplace-Stieltjes Transform using den Iseger's algorithm.

    Args:
        f_hat: the Laplace tansform to be inverted.
        mexp: the exponent defining the number of evaluation points (2^mexp) for the inverse where precision is guaranteed; 
            The transform evaluation uses an oversampling factor, <OVRSMPL> * 2^(mexp).
        delta: mesh width.
        n: The size of the Gauss quadrature; Choose between 16, 32 and 48; Defaults to 48 (recommended for higher precision). 
    
    Returns:
            A std::vector containing the inverse function evaluated at {0, delta, 2*delta, ...}.
    */
    std::vector<double> lambda, beta;
    __loadGaussQuadRule(n, beta, lambda);

    auto inverse_array = __oneDimensionalInverse(f_hat, delta, mexp, n, lambda, beta);
    return inverse_array;
}

} // namespace
