
#ifdef AOPROJECT
#include "MultiDirSolver.h"
#include "Matrix2x2.h"
#else
#include <DPPP_DDECal/MultiDirSolver.h>
#include <DPPP_DDECal/Matrix2x2.h>
#endif

using namespace arma;

MultiDirSolver::MultiDirSolver(size_t maxIterations, double accuracy, double stepSize) :
  _nAntennas(0),
  _nDirections(0),
  _nChannels(0),
  _nChannelBlocks(0),
  _maxIterations(maxIterations),
  _accuracy(accuracy),
  _stepSize(stepSize),
  _phaseOnly(false)
{
}

void MultiDirSolver::init(size_t nAntennas,
                          size_t nDirections,
                          size_t nChannels,
                          const std::vector<int>& ant1,
                          const std::vector<int>& ant2)
{
  _nAntennas = nAntennas;
  _nDirections = nDirections;
  _nChannels = nChannels;
  _nChannelBlocks = nChannels;
  _ant1 = ant1;
  _ant2 = ant2;
}

MultiDirSolver::SolveResult MultiDirSolver::processScalar(std::vector<Complex *>& data,
  std::vector<std::vector<Complex *> >& modelData,
  std::vector<std::vector<DComplex> >& solutions, double time) const
{
  const size_t nTimes = data.size();
  SolveResult result;
  
  std::vector<std::vector<DComplex> > nextSolutions(_nChannelBlocks);

#ifndef NDEBUG
  if (solutions.size() != _nChannelBlocks) {
    cout << "Error: 'solutions' parameter does not have the right shape" << endl;
    result.iterations = 0;
    return result;
  }
#endif

  result._results.resize(_constraints.size());
  
  // Model matrix ant x [N x D] and visibility vector ant x [N x 1],
  // for each channelblock
  // The following loop allocates all structures
  std::vector<std::vector<cx_mat> > gTimesCs(_nChannelBlocks);
  std::vector<std::vector<cx_vec> > vs(_nChannelBlocks);
  for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
  {
    //solutions[chBlock].assign(_nDirections * _nAntennas, 1.0);
    nextSolutions[chBlock].resize(_nDirections * _nAntennas);
    const size_t
      channelIndexStart = chBlock * _nChannels / _nChannelBlocks,
      channelIndexEnd = (chBlock+1) * _nChannels / _nChannelBlocks,
      curChannelBlockSize = channelIndexEnd - channelIndexStart;
    gTimesCs[chBlock].resize(_nAntennas);
    vs[chBlock].resize(_nAntennas);
    
    for(size_t ant=0; ant!=_nAntennas; ++ant)
    {
      // Model matrix [N x D] and visibility vector [N x 1]
      // Also space for the auto correlation is reserved, but they will be set to 0.
      gTimesCs[chBlock][ant] = cx_mat(_nAntennas * nTimes * curChannelBlockSize * 4,
                                      _nDirections, fill::zeros);
      vs[chBlock][ant] = cx_vec(_nAntennas * nTimes * curChannelBlockSize * 4, fill::zeros);
    }
  }
  
  // TODO the data and model data needs to be preweighted.
  // Maybe we can get a non-const pointer from DPPP, that saves copying/allocating
  
  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  double normSum = 0.0, sum = 0.0;
  do {
#pragma omp parallel for
    for(size_t chBlock=0; chBlock<_nChannelBlocks; ++chBlock)
    {
      performScalarIteration(chBlock, gTimesCs[chBlock], vs[chBlock],
                            solutions[chBlock], nextSolutions[chBlock],
                            data, modelData);
    }
      
    // Move the solutions towards nextSolutions
    // (the moved solutions are stored in 'nextSolutions')
    for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
    {
      for(size_t i=0; i!=_nAntennas*_nDirections; ++i)
      {
        if(_phaseOnly)
        {
          double ab = std::abs(nextSolutions[chBlock][i]);
          if(ab != 0.0)
            nextSolutions[chBlock][i] /= ab;
        }
        nextSolutions[chBlock][i] = solutions[chBlock][i]*(1.0-_stepSize) +
          nextSolutions[chBlock][i] * _stepSize;
      }
    }
    
    for(size_t i=0; i!=_constraints.size(); ++i)
    {
      result._results[i] = _constraints[i]->Apply(nextSolutions, time);
    }
    
    //  Calculate the norm of the difference between the old and new solutions
    for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
    {
      for(size_t i=0; i!=_nAntennas*_nDirections; ++i)
      {
        double e = std::norm(nextSolutions[chBlock][i] - solutions[chBlock][i]);
        normSum += e;
        sum += std::abs(solutions[chBlock][i]);
        
        solutions[chBlock][i] = nextSolutions[chBlock][i];
        
        // For debug: output the solutions of the first antenna
        if(i<_nDirections && false)
        {
          std::cout << " |s_" << i << "|=|" << solutions[chBlock][i] << "|="
                    << std::abs(solutions[chBlock][i]);
        }
      }
    }
    normSum /= _nChannelBlocks*_nAntennas*_nDirections;
    sum /= _nChannelBlocks*_nAntennas*_nDirections;
    iteration++;
    
  } while(iteration < _maxIterations && normSum/sum > _accuracy);
  
  if(normSum/sum <= _accuracy)
    result.iterations = iteration;
  else
    result.iterations = _maxIterations+1;
  return result;
}

void MultiDirSolver::performScalarIteration(size_t channelBlockIndex,
                       std::vector<arma::cx_mat>& gTimesCs,
                       std::vector<arma::cx_vec>& vs,
                       const std::vector<DComplex>& solutions,
                       std::vector<DComplex>& nextSolutions,
                       const std::vector<Complex *>& data,
                       const std::vector<std::vector<Complex *> >& modelData) const
{
  for(size_t ant=0; ant!=_nAntennas; ++ant)
  {
    gTimesCs[ant].zeros();
    vs[ant].zeros();
  }
  
  const size_t
    channelIndexStart = channelBlockIndex * _nChannels / _nChannelBlocks,
    channelIndexEnd = (channelBlockIndex+1) * _nChannels / _nChannelBlocks,
    curChannelBlockSize = channelIndexEnd - channelIndexStart,
    nTimes = data.size();
  
  // The following loop fills the matrices for all antennas
  for(size_t timeIndex=0; timeIndex!=nTimes; ++timeIndex)
  {
    std::vector<const Complex*> modelPtrs(_nDirections);
    for(size_t baseline=0; baseline!=_ant1.size(); ++baseline)
    {
      size_t antenna1 = _ant1[baseline];
      size_t antenna2 = _ant2[baseline];
      if(antenna1 != antenna2)
      {
        cx_mat& gTimesC1 = gTimesCs[antenna1];
        cx_vec& v1 = vs[antenna1];
        cx_mat& gTimesC2 = gTimesCs[antenna2];
        cx_vec& v2 = vs[antenna2];
        for(size_t d=0; d!=_nDirections; ++d)
          modelPtrs[d] = modelData[timeIndex][d] + (channelIndexStart + baseline * _nChannels) * 4;
        const Complex* dataPtr = data[timeIndex] + (channelIndexStart + baseline * _nChannels) * 4;
        const size_t p1top2[4] = {0, 2, 1, 3};
        for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
        {
          const size_t
            dataIndex1 = ch-channelIndexStart + (timeIndex + antenna1 * nTimes) * curChannelBlockSize,
            dataIndex2 = ch-channelIndexStart + (timeIndex + antenna2 * nTimes) * curChannelBlockSize;
          for(size_t p1=0; p1!=4; ++p1)
          {
            size_t p2 = p1top2[p1];
            for(size_t d=0; d!=_nDirections; ++d)
            {
              std::complex<double> predicted = *modelPtrs[d];
              
              size_t solIndex1 = antenna1*_nDirections + d;
              size_t solIndex2 = antenna2*_nDirections + d;
              gTimesC2(dataIndex1*4+p1, d) = std::conj(solutions[solIndex1] * predicted); // using a* b* = (ab)*
              gTimesC1(dataIndex2*4+p2, d) = std::conj(solutions[solIndex2]) * predicted;
              
              ++modelPtrs[d]; // Goto the next polarization of this 2x2 matrix.
            }
            v1(dataIndex2*4+p2) = *dataPtr;
            v2(dataIndex1*4+p1) = std::conj(*dataPtr);
            ++dataPtr; // Goto the next polarization of this 2x2 matrix.
          }
        }
      }
    }
  }
  
  // The matrices have been filled; compute the linear solution
  // for each antenna.
  for(size_t ant=0; ant!=_nAntennas; ++ant)
  {
    cx_mat& gTimesC = gTimesCs[ant];
    cx_vec& v = vs[ant];
    // solve [g* C] x  = v
    cx_vec x = solve(gTimesC, v);
    for(size_t d=0; d!=_nDirections; ++d)
      nextSolutions[ant*_nDirections + d] = x(d);
  }
}

MultiDirSolver::SolveResult MultiDirSolver::processFullJones(std::vector<Complex *>& data,
  std::vector<std::vector<Complex *> >& modelData,
  std::vector<std::vector<DComplex> >& solutions, double time) const
{
  // This algorithm is basically the same, but visibility values are
  // extended to 2x2 matrices and concatenated in the matrices
  // equations as block matrices.
  
  // First we pre-apply the left-hand solutions to the model to make JM. Each
  // 2x2 coherence matrix Ji is matrix-multied by the lh solutions, for all
  // directions, and visibilities (times x channels).
  //   JMi = Ji Mi
  // These are stacked in matrix JM :
  //        JM0_d0 JM1_d0 ...
  //   JM = JM0_d1 JM1_d1 
  //        ...          
  // such that JM is a (2D) rows x (2N) col matrix, N=nvis, D=ndir.
  // The solved rh 2D x 2 solution matrix is similarly formed with the rh solution
  // values:
  //       J0
  //   J = J1
  //       ...
  // And the 2N x 2 visibility matrix as well:
  //       V0
  //   V = V1
  //       ...
  // And we solve the equation:
  //   'JM' J* = V
  // Rewritten:
  //   'JM'* J = V*
  // With dimensions:
  //   [ 2N x 2D ] [ 2D x 2 ] = [ 2N x 2 ]

  const size_t nTimes = data.size();
  SolveResult result;
  
  std::vector<std::vector<DComplex> > nextSolutions(_nChannelBlocks);

#ifndef NDEBUG
  if (solutions.size() != _nChannelBlocks) {
    cout << "Error: 'solutions' parameter does not have the right shape" << endl;
    result.iterations = 0;
    return result;
  }
#endif

  result._results.resize(_constraints.size());
  
  // Model matrix ant x [2N x 2D] and visibility matrix ant x [2N x 2],
  // for each channelblock
  // The following loop allocates all structures
  std::vector<std::vector<cx_mat> > gTimesCs(_nChannelBlocks);
  std::vector<std::vector<cx_mat> > vs(_nChannelBlocks);
  for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
  {
    nextSolutions[chBlock].resize(_nDirections * _nAntennas * 4);
    const size_t
      channelIndexStart = chBlock * _nChannels / _nChannelBlocks,
      channelIndexEnd = (chBlock+1) * _nChannels / _nChannelBlocks,
      curChannelBlockSize = channelIndexEnd - channelIndexStart;
    gTimesCs[chBlock].resize(_nAntennas);
    vs[chBlock].resize(_nAntennas);
    
    for(size_t ant=0; ant!=_nAntennas; ++ant)
    {
      // Model matrix [2N x 2D] and visibility matrix [2N x 2]
      // Also space for the auto correlation is reserved, but they will be set to 0.
      gTimesCs[chBlock][ant] = cx_mat(_nAntennas * nTimes * curChannelBlockSize * 2,
                                      _nDirections * 2, fill::zeros);
      vs[chBlock][ant] = cx_mat(_nAntennas * nTimes * curChannelBlockSize * 2, 2, fill::zeros);
    }
  }
  
  // TODO the data and model data needs to be preweighted.
  // Maybe we can get a non-const pointer from DPPP, that saves copying/allocating
  
  ///
  /// Start iterating
  ///
  size_t iteration = 0;
  double normSum = 0.0, sum = 0.0;
  do {
#pragma omp parallel for
    for(size_t chBlock=0; chBlock<_nChannelBlocks; ++chBlock)
    {
      performFullJonesIteration(chBlock, gTimesCs[chBlock], vs[chBlock],
                                solutions[chBlock], nextSolutions[chBlock],
                                data, modelData);
    }
      
    // Move the solutions towards nextSolutions
    // (the moved solutions are stored in 'nextSolutions')
    for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
    {
      for(size_t i=0; i!=_nAntennas*_nDirections*4; ++i)
      {
        if(_phaseOnly)
        {
          double ab = std::abs(nextSolutions[chBlock][i]);
          if(ab != 0.0)
            nextSolutions[chBlock][i] /= ab;
        }
        nextSolutions[chBlock][i] = solutions[chBlock][i]*(1.0-_stepSize) +
          nextSolutions[chBlock][i] * _stepSize;
      }
    }
    
    for(size_t i=0; i!=_constraints.size(); ++i)
    {
      result._results[i] = _constraints[i]->Apply(nextSolutions, time);
    }
    
    //  Calculate the norm of the difference between the old and new solutions
    for(size_t chBlock=0; chBlock!=_nChannelBlocks; ++chBlock)
    {
      for(size_t i=0; i!=_nAntennas*_nDirections*4; ++i)
      {
        double e = std::norm(nextSolutions[chBlock][i] - solutions[chBlock][i]);
        normSum += e;
        sum += std::abs(solutions[chBlock][i]);
        
        solutions[chBlock][i] = nextSolutions[chBlock][i];
      }
    }
    normSum /= _nChannelBlocks*_nAntennas*_nDirections*4;
    sum /= _nChannelBlocks*_nAntennas*_nDirections*4;
    iteration++;
    
  } while(iteration < _maxIterations && normSum/sum > _accuracy);
  
  if(normSum/sum <= _accuracy)
    result.iterations = iteration;
  else
    result.iterations = _maxIterations+1;
  return result;
}

void MultiDirSolver::performFullJonesIteration(size_t channelBlockIndex,
                             std::vector<arma::cx_mat>& gTimesCs,
                             std::vector<arma::cx_mat>& vs,
                             const std::vector<DComplex>& solutions,
                             std::vector<DComplex>& nextSolutions,
                             const std::vector<Complex *>& data,
                             const std::vector<std::vector<Complex *> >& modelData) const
{
  for(size_t ant=0; ant!=_nAntennas; ++ant)
  {
    gTimesCs[ant].zeros();
    vs[ant].zeros();
  }
  
  const size_t
    channelIndexStart = channelBlockIndex * _nChannels / _nChannelBlocks,
    channelIndexEnd = (channelBlockIndex+1) * _nChannels / _nChannelBlocks,
    curChannelBlockSize = channelIndexEnd - channelIndexStart,
    nTimes = data.size();
  
  // The following loop fills the matrices for all antennas
  for(size_t timeIndex=0; timeIndex!=nTimes; ++timeIndex)
  {
    std::vector<const Complex*> modelPtrs(_nDirections);
    for(size_t baseline=0; baseline!=_ant1.size(); ++baseline)
    {
      size_t antenna1 = _ant1[baseline];
      size_t antenna2 = _ant2[baseline];
      if(antenna1 != antenna2)
      {
        cx_mat
          &gTimesC1 = gTimesCs[antenna1],
          &v1 = vs[antenna1],
          &gTimesC2 = gTimesCs[antenna2],
          &v2 = vs[antenna2];
        for(size_t d=0; d!=_nDirections; ++d)
          modelPtrs[d] = modelData[timeIndex][d] + (channelIndexStart + baseline * _nChannels) * 4;
        const Complex* dataPtr = data[timeIndex] + (channelIndexStart + baseline * _nChannels) * 4;
        for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
        {
          const size_t
            dataIndex1 = 2 * (ch-channelIndexStart + (timeIndex + antenna1 * nTimes) * curChannelBlockSize),
            dataIndex2 = 2 * (ch-channelIndexStart + (timeIndex + antenna2 * nTimes) * curChannelBlockSize);
            
          for(size_t d=0; d!=_nDirections; ++d)
          {
            MC2x2
              modelMat(modelPtrs[d]),
              gTimesC1Mat, gTimesC2Mat;
            size_t solIndex1 = antenna1*_nDirections + d;
            size_t solIndex2 = antenna2*_nDirections + d;
            Matrix2x2::HermATimesHermB(gTimesC2Mat.Data(), &solutions[solIndex1*4], modelMat.Data());
            Matrix2x2::HermATimesB(gTimesC1Mat.Data(), &solutions[solIndex2*4], modelMat.Data());
            for(size_t p=0; p!=4; ++p)
            {
              gTimesC2(dataIndex1+(p/2), d+p%2) = gTimesC2Mat[p];
              gTimesC1(dataIndex2+(p/2), d+p%2) = gTimesC1Mat[p];
            }
            
            modelPtrs[d] += 4; // Goto the next 2x2 matrix.
          }
          for(size_t p=0; p!=4; ++p)
          {
            v1(dataIndex2+(p/2), p%2) = *dataPtr;
            v2(dataIndex1+(p/2), p%2) = std::conj(*dataPtr);
            ++dataPtr;  // Goto the next element of 2x2 matrix.
          }
        }
      }
    }
  }
  
  // The matrices have been filled; compute the linear solution
  // for each antenna.
  for(size_t ant=0; ant!=_nAntennas; ++ant)
  {
    cx_mat& gTimesC = gTimesCs[ant];
    cx_mat& v = vs[ant];
    // solve [g* C] x  = v
    cx_mat x = solve(gTimesC, v);
    for(size_t d=0; d!=_nDirections; ++d)
    {
      for(size_t p=0; p!=4; ++p)
        nextSolutions[(ant*_nDirections + d)*4 + p] = x(d+p/2, p%2);
    }
  }
}
