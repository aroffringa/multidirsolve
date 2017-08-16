#ifdef AOPROJECT
#include "TECConstraint.h"
#include <omp.h> // for tec constraints
#else
#include <DPPP_DDECal/TECConstraint.h>
#include <Common/OpenMP.h>
#endif

TECConstraintBase::TECConstraintBase(Mode mode) :
  _mode(mode),
  _nAntennas(0),
  _nDirections(0),
  _nChannelBlocks(0),
  _phaseFitters()
{
}

void TECConstraintBase::initialize(size_t nAntennas, size_t nDirections, size_t nChannelBlocks, const double* frequencies) {
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
  
  initializeChild();
}

void ApproximateTECConstraint::initializeChild()
{
  _pwFitters.resize(
#ifdef AOPROJECT
      omp_get_max_threads()
#else
      LOFAR::OpenMP::maxThreads()
#endif
   );
  _threadData.resize(_pwFitters.size());
  _threadFittedData.resize(_pwFitters.size());
  
  if(_fittingChunkSize == 0)
  {
    size_t
      n = _phaseFitters.front().Size();
    const double
      startFreq = _phaseFitters.front().FrequencyData()[0],
      endFreq = _phaseFitters.front().FrequencyData()[n-1];
      _fittingChunkSize = PieceWisePhaseFitter::CalculateChunkSize(startFreq, endFreq, n);
  }
}

std::vector<Constraint::Result> TECConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double)
{
  std::vector<Constraint::Result> res(2);

  res[0].vals.resize(_nAntennas*_nDirections);
  res[0].axes="ant,dir,freq";
  res[0].name="tec";
  res[0].dims.resize(3);
  res[0].dims[0]=_nAntennas;
  res[0].dims[1]=_nDirections;
  res[0].dims[2]=1;
  res[1]=res[0];
  res[1].name="scalarphase";

  // TODO chose this more cleverly?
  size_t refAntenna = 0;

  // Divide out the reference antenna
  for(size_t ch=0; ch!=_nChannelBlocks; ++ch)
  {
    for(size_t antennaIndex=0; antennaIndex!=_nAntennas; ++antennaIndex)
    {
      for(size_t d=0; d!=_nDirections; ++d)
      {
        size_t solutionIndex = antennaIndex*_nDirections + d;
        size_t refAntennaIndex = d + refAntenna*_nDirections;
        if(antennaIndex != refAntenna)
        {
          solutions[ch][solutionIndex] = solutions[ch][solutionIndex] / solutions[ch][refAntennaIndex];
        }
      }
    }
    for(size_t d=0; d!=_nDirections; ++d)
      solutions[ch][refAntenna*_nDirections + d] = 1.0;
  }
  
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
      if(std::isfinite(solutions[ch][solutionIndex].real()) &&
        std::isfinite(solutions[ch][solutionIndex].imag()))
      {
        _phaseFitters[thread].PhaseData()[ch] = std::arg(solutions[ch][solutionIndex]);
        _phaseFitters[thread].WeightData()[ch] = 1.0;
      }
      else {
        _phaseFitters[thread].PhaseData()[ch] = 0.0;
        _phaseFitters[thread].WeightData()[ch] = 0.0;
      }
    }
    
    double alpha, beta=0.0;
    if(_mode == TECOnlyMode) {
      _phaseFitters[thread].FitDataToTEC1Model(alpha);
    } else {
      _phaseFitters[thread].FitDataToTEC2Model(alpha, beta);
    }

    res[0].vals[solutionIndex] = alpha / -8.44797245e9;
    res[1].vals[solutionIndex] = beta;
    
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch) 
    {
      solutions[ch][solutionIndex] = std::polar<double>(1.0, _phaseFitters[thread].PhaseData()[ch]);
    }
  }

  return res;
}

std::vector<Constraint::Result> ApproximateTECConstraint::Apply(
    std::vector<std::vector<dcomplex> >& solutions, double time)
{
  if(_finishedApproximateStage)
    return TECConstraint::Apply(solutions, time);
  else {
    // TODO chose this more cleverly?
    size_t refAntenna = 0;

    // Divide out the reference antenna
    for(size_t ch=0; ch!=_nChannelBlocks; ++ch)
    {
      for(size_t antennaIndex=0; antennaIndex!=_nAntennas; ++antennaIndex)
      {
        for(size_t d=0; d!=_nDirections; ++d)
        {
          size_t solutionIndex = antennaIndex*_nDirections + d;
          size_t refAntennaIndex = d + refAntenna*_nDirections;
          if(antennaIndex != refAntenna)
          {
            solutions[ch][solutionIndex] = solutions[ch][solutionIndex] / solutions[ch][refAntennaIndex];
          }
        }
      }
      for(size_t d=0; d!=_nDirections; ++d)
        solutions[ch][refAntenna*_nDirections + d] = 1.0;
    }
    
  #pragma omp parallel for
    for(size_t solutionIndex = 0; solutionIndex<_nAntennas*_nDirections; ++solutionIndex)
    {
  #ifdef AOPROJECT
      size_t thread = omp_get_thread_num();
  #else
      size_t thread = LOFAR::OpenMP::threadNum();
  #endif
      std::vector<double>& data = _threadData[thread];
      std::vector<double>& fittedData = _threadFittedData[thread];
      
      data.resize(_nChannelBlocks);

      for(size_t ch=0; ch!=_nChannelBlocks; ++ch) {
        if(std::isfinite(solutions[ch][solutionIndex].real()) &&
          std::isfinite(solutions[ch][solutionIndex].imag()))
        {
          data[ch] = std::arg(solutions[ch][solutionIndex]);
          // TODO weight
        }
        else {
          data[ch] = 0.0;
          // TODO weight
        }
      }
      
      _pwFitters[thread].SlidingFit(_phaseFitters[thread].FrequencyData(), data, fittedData);

      for(size_t ch=0; ch!=_nChannelBlocks; ++ch) 
      {
        solutions[ch][solutionIndex] = std::polar<double>(1.0, _phaseFitters[thread].PhaseData()[ch]);
      }
    }

    return std::vector<Constraint::Result>();
  }
}