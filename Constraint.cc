#ifdef AOPROJECT
#include "Constraint.h"
#include <omp.h> // for tec constraints
#else
#include <DPPP/Constraint.h>
#include <Common/OpenMP.h>
#endif


Constraint::Result PhaseConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions)
{
  for (uint ch=0; ch<solutions.size(); ++ch) {
    for (uint solIndex=0; solIndex<solutions[ch].size(); ++solIndex) {
      solutions[ch][solIndex] /= std::abs(solutions[ch][solIndex]);
    }
  }

  return Constraint::Result();
}

TECConstraint::TECConstraint(Mode mode) :
  _mode(mode),
  _nAntennas(0)
  ,_nDirections(0),
  _nChannelBlocks(0),
  _phaseFitters()
{
}

TECConstraint::TECConstraint(Mode mode, size_t nAntennas, size_t nDirections, 
                             size_t nChannelBlocks, const double* frequencies) :
  _mode(mode)
{
  init(nAntennas, nDirections, nChannelBlocks, frequencies);
}

void TECConstraint::init(size_t nAntennas, size_t nDirections, 
                         size_t nChannelBlocks, const double* frequencies) {
  _nAntennas = nAntennas;
  _nDirections = nDirections;
  _nChannelBlocks = nChannelBlocks;
  _phaseFitters.resize(
#ifdef AOPROJECT
      omp_get_max_threads()
#else
      LOFAR::OpenMP::maxThreads()
#endif
   );

  for(size_t i=0; i!=_phaseFitters.size(); ++i)
  {
    _phaseFitters[i].SetChannelCount(_nChannelBlocks);
    std::memcpy(_phaseFitters[i].FrequencyData(), frequencies, sizeof(double) * _nChannelBlocks);
  }
}

Constraint::Result TECConstraint::Apply(std::vector<std::vector<dcomplex> >& solutions)
{
  TECConstraint::TECResult res;
  res.tecvals.resize(_nAntennas*_nDirections);
#pragma omp parallel for
  for(size_t solutionIndex = 0; solutionIndex<_nAntennas*_nDirections; ++solutionIndex)
  {
    size_t thread =
#ifdef AOPROJECT
        omp_get_thread_num();
#else
        LOFAR::OpenMP::threadNum();
#endif

    for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
      _phaseFitters[thread].PhaseData()[ch] = std::arg(solutions[ch][solutionIndex]);
    }
    
    double alpha, beta=0.0;
    if(_mode == TECOnlyMode) {
      _phaseFitters[thread].FitDataToTEC1Model(alpha);
    } else {
      _phaseFitters[thread].FitDataToTEC2Model(alpha, beta);
    }

    res.tecvals[solutionIndex].first = alpha;
    res.tecvals[solutionIndex].second = beta;
    
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
     solutions[ch][solutionIndex] = std::polar<double>(1.0, _phaseFitters[thread].PhaseData()[ch]);
    }
  }

  return res;
}
