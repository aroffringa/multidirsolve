#include "Constraint.h"
#include "KernelSmoother.h"

#ifndef SMOOTHNESS_CONSTRAINT_H
#define SMOOTHNESS_CONSTRAINT_H

class SmoothnessConstraint : public Constraint
{
public:
  typedef std::complex<double> dcomplex;
  
  SmoothnessConstraint(const double* frequencies, size_t n, double bandwidthHz);
  
  std::vector<Constraint::Result> Apply(
    std::vector<std::vector<dcomplex> >& solutions, double, std::ostream* statStream) final override;
  
  void SetWeights(std::vector<double> &weights) final override {
    _weights = weights;
  }
  
  virtual void InitializeDimensions(size_t nAntennas,
                                    size_t nDirections,
                                    size_t nChannelBlocks) final override;
                                    
  struct FitData
  {
    FitData(const double* frequencies, size_t n, double kernelBandwidth)
      : smoother(frequencies, n, kernelBandwidth),
      data(n), weight(n)
    { }
    
    KernelSmoother<dcomplex, double> smoother;
    std::vector<dcomplex> data;
    std::vector<double> weight;
  };
  std::vector<FitData> _fitData;
  std::vector<double> _frequencies, _weights;
  double _bandwidth;
};

#endif

