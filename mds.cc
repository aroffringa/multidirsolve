#include "MultiDirSolver.h"
#include "Matrix2x2.h"
#include "TECConstraint.h"

#include "stopwatch.h"

int main(int argc, char* argv[])
{
  std::vector<Stopwatch> watches(2);
  for(size_t fullOrNot=0; fullOrNot!=2; ++fullOrNot)
  {
    watches[fullOrNot].Start();
    typedef std::complex<float> cf;
    MultiDirSolver mds;
    mds.set_max_iterations(1000);
    mds.set_accuracy(1e-7);
    mds.set_step_size(0.5);
    
    mds.set_phase_only(false);
    size_t nPol = 4, nAnt = 200, nDir = 3, nChan = 10, nChanBlocks = 2, nTimes = 1, nBl=nAnt*(nAnt-1)/2;
    
    std::vector<int> ant1s, ant2s;
    for(size_t a1=0; a1!=nAnt; ++a1)
    {
      for(size_t a2=a1+1; a2!=nAnt; ++a2)
      {
        ant1s.push_back(a1);
        ant2s.push_back(a2);
      }
    }
    
    mds.init(nAnt, nDir, nChan, ant1s, ant2s);
    mds.set_channel_blocks(nChanBlocks);
    
    std::vector<double> nu(nChanBlocks);
    for(size_t i=0; i!=nChanBlocks; ++i)
      nu[i] = i+1;
    
    TECConstraint tecConstraint(TECConstraint::TECOnlyMode);
    tecConstraint.initialize(nAnt, nDir, nChanBlocks, nu.data());
    
    cf gain1(0.31415926535, 0.0), gain2(2.0, 1.0), gain3(0.0, 3.0);
    std::vector<cf> inputSolutions(nAnt * nDir);
    for(size_t a=0; a!=nAnt; ++a)
    {
      if(a == 1)
      {
        inputSolutions[a*nDir + 0] = 1.0;
        inputSolutions[a*nDir + 1] = 1.0;
        inputSolutions[a*nDir + 2] = 1.0;
      }
      else {
        inputSolutions[a*nDir + 0] = a;
        inputSolutions[a*nDir + 1] = gain2;
        inputSolutions[a*nDir + 2] = gain3;
      }
    }
    
    //cf gain1(0.5*M_SQRT2, 0.5*M_SQRT2), gain2(0.0, 1.0), gain3(-0.5*M_SQRT2, -0.5*M_SQRT2);
    std::vector<cf*> data;
    std::vector<std::vector<cf*>> modelData;
    for(size_t timestep=0; timestep!=nTimes; ++timestep)
    {
      cf* dataPtr = new cf[nPol * nChan * nBl];
      cf* model1Ptr = new cf[nPol * nChan * nBl];
      cf* model2Ptr = new cf[nPol * nChan * nBl];
      cf* model3Ptr = new cf[nPol * nChan * nBl];
      
      for(size_t i=0; i!=nChan * nBl; ++i)
      {
        float unit = 1.0 * (i+1);
        model1Ptr[i*4 + 0] = unit;
        model1Ptr[i*4 + 1] = 0.0;
        model1Ptr[i*4 + 2] = 0.0;
        model1Ptr[i*4 + 3] = unit;
        model2Ptr[i*4 + 0] = (i%2==0) ? unit : 0.0;
        model2Ptr[i*4 + 1] = 0.0;
        model2Ptr[i*4 + 2] = 0.0;
        model2Ptr[i*4 + 3] = (i%2==0) ? unit : 0.0;
        model3Ptr[i*4 + 0] = (i%3==0) ? unit : 0.0;
        model3Ptr[i*4 + 1] = 0.0;
        model3Ptr[i*4 + 2] = 0.0;
        model3Ptr[i*4 + 3] = (i%3==0) ? unit : 0.0;
      }
      
      size_t baselineIndex = 0;
      for(size_t a1=0; a1!=nAnt; ++a1)
      {
        for(size_t a2=a1+1; a2!=nAnt; ++a2)
        {
          for(size_t j=0; j!=nPol * nChan; ++j)
          {
            cf gain1ant1, gain2ant1, gain3ant1, gain1ant2, gain2ant2, gain3ant2;
            gain1ant1 = inputSolutions[a1*nDir + 0];
            gain2ant1 = inputSolutions[a1*nDir + 1];
            gain3ant1 = inputSolutions[a1*nDir + 2];
            gain1ant2 = inputSolutions[a2*nDir + 0];
            gain2ant2 = inputSolutions[a2*nDir + 1];
            gain3ant2 = inputSolutions[a2*nDir + 2];
            dataPtr[baselineIndex] =
              gain1ant1 * std::conj(gain1ant2) * model1Ptr[baselineIndex] +
              gain2ant1 * std::conj(gain2ant2) * model2Ptr[baselineIndex] +
              gain3ant1 * std::conj(gain3ant2) * model3Ptr[baselineIndex];
            ++baselineIndex;
          }
        }
      }
      
      data.push_back(dataPtr);
      modelData.push_back(std::vector<cf*>{model1Ptr, model2Ptr, model3Ptr});
    }
    
    MultiDirSolver::SolveResult result;
    std::vector<std::vector<std::complex<double> > > solutions(nChanBlocks);
    if(fullOrNot == 0)
    {
      for(auto& vec : solutions)
        vec.assign(nDir * nAnt, 1.0);
      result = mds.processScalar(data, modelData, solutions, 0.0);
      std::cout << '\n';
      for(size_t ch=0; ch!=nChanBlocks; ++ch)
      {
        for(size_t ant=0; ant!=nAnt; ++ant)
        {
          for(size_t d=0; d!=nDir; ++d)
          {
            std::cout << "ch" << ch << ", ant" << ant << ", d" << d << " = " << solutions[ch][d + ant*nDir] << ", ";
            std::complex<double> sol = solutions[ch][d + ant*nDir]/solutions[ch][d + 1*nDir];
            std::complex<double> inp(inputSolutions[ant*nDir+d]);
            std::cout << "inp: " << inp << ' ';
            std::cout << "ref: " <<
            sol << " dist: " << std::abs(sol-inp) << '\n';
          }
        }
      }
    }
    else {
      for(auto& vec : solutions)
      {
        vec.assign(nDir * nAnt * 4, 0.0);
        for(size_t i=0; i!=nDir*nAnt; ++i)
        {
          vec[i * 4 + 0] = 1.0;
          vec[i * 4 + 3] = 1.0;
        }
      }
      result = mds.processFullMatrix(data, modelData, solutions, 0.0);
      std::cout << '\n';
      for(size_t ch=0; ch!=nChanBlocks; ++ch)
      {
        for(size_t ant=0; ant!=nAnt; ++ant)
        {
          for(size_t d=0; d!=nDir; ++d)
          {
            MC2x2 solM(&solutions[ch][(d + ant*nDir) * 4]);
            std::cout << "ch" << ch << ", ant" << ant << ", d" << d << " = " << solM.ToString() << ", ";
            std::complex<double> sol = solM[0]/solutions[ch][(d + 1*nDir)*4];
            std::complex<double> inp(inputSolutions[ant*nDir+d]);
            std::cout << "inp: " << inp << ' ';
            std::cout << "ref: " <<
            sol << " dist: " << std::abs(sol-inp) << '\n';
          }
        }
      }
    }
    std::cout << "Iterations: " << result.iterations << '\n';
    
    while(!data.empty())
    {
      delete[] modelData.back()[0];
      delete[] modelData.back()[1];
      delete[] modelData.back()[2];
      modelData.pop_back();
      delete[] data.back();
      data.pop_back();
    }
    watches[fullOrNot].Pause();
  }
  
  for(const Stopwatch& watch : watches)
    std::cout << watch.ToString() << '\n';
}
