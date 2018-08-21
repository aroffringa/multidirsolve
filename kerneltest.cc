#include <iostream>
#include "KernelSmoother.h"
#include "PieceWisePhaseFitter.h"

#include <complex>
#include <random>

int main(int argc, char* argv[])
{
  std::normal_distribution<double> gaus(0.0, 1.0);
  std::uniform_real_distribution<double> uniform(0.2, 1.0);
  std::mt19937 rnd;
  
  size_t n = 1000;
  std::vector<std::complex<double>> function(n), values(n);
  std::vector<double> frequencies(n), weights(n);
  for(size_t i=0; i!=n; ++i)
  {
    frequencies[i] = 150e6 + double(i)*50e6/n; // 150 - 200 MHz
    double x = double(i)/n;
    double x3 = 3.3*(x-0.4);
    double x2 = x-0.2;
    double y = x3*x3*x3-12.0*x2*x2+2.0;
    double stddev = uniform(rnd);
    
    function[i] = y;
    weights[i] = 1.0/(stddev*stddev);
    values[i] = y + gaus(rnd) * stddev;
  }

  std::vector<std::complex<double>> rectangular(values), triangular(values), gaussian(values), quadratic(values);
  
  typedef KernelSmoother<std::complex<double>, double> Smoother;
  
  Smoother rectangularSmoother(frequencies.data(), n, Smoother::RectangularKernel, 5.0e6);
  rectangularSmoother.Smooth(rectangular.data(), weights.data());
  
  Smoother triangularSmoother(frequencies.data(), n, Smoother::TriangularKernel, 5.0e6);
  triangularSmoother.Smooth(triangular.data(), weights.data());
  
  Smoother gaussianSmoother(frequencies.data(), n, Smoother::GaussianKernel, 5.0e6);
  gaussianSmoother.Smooth(gaussian.data(), weights.data());
  
  Smoother quadraticSmoother(frequencies.data(), n, Smoother::TriangularKernel, 5.0e6);
  quadraticSmoother.Smooth(quadratic.data(), weights.data());
  
  for(size_t i=0; i!=n; ++i)
  {
    std::cout << frequencies[i]*1e-6 << '\t' << function[i].real() << '\t' << values[i].real() << '\t'
      << rectangular[i].real() << '\t'
      << triangular[i].real() << '\t'
      << gaussian[i].real() << '\t'
      << quadratic[i].real() << '\n';
  }
  
  for(size_t i=0; i!=10000; ++i)
    quadraticSmoother.Smooth(quadratic.data(), weights.data());
  
  return 0;
}

