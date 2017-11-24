#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <complex>

#include "PhaseFitter.h"
#include "PieceWisePhaseFitter.h"

void readFile(const std::string& filename)
{
  std::vector<double> frequencies;
  std::vector<double> phases;
  std::ifstream file(filename);
  size_t flagged = 0;
  while(file.good())
  {
    double freq, p;
    file >> freq >> p;
    if(file.good())
    {
      if(p != 0.0)
      {
        frequencies.push_back(freq);
        phases.push_back(p);
      }
      else ++flagged;
    }
  }
  std::cout << "Removed " << flagged << "/" << phases.size() << " values.\n";
  
  PhaseFitter pFitter(phases.size());
  size_t chunkSize = PieceWisePhaseFitter::CalculateChunkSize(frequencies.front(), frequencies.back(), frequencies.size());
  PieceWisePhaseFitter pwFitter(chunkSize);
  
  std::copy(phases.begin(), phases.end(), pFitter.PhaseData());
  std::copy(frequencies.begin(), frequencies.end(), pFitter.FrequencyData());
  pFitter.FitDataToTEC2Model();
  
  std::vector<double> fittedData(phases.size());
  pwFitter.SlidingFit(frequencies.data(), phases, fittedData); 
  std::ofstream outf("fit.txt");
  for(size_t i=0; i!=fittedData.size(); ++i)
  {
    outf << frequencies[i] << '\t' << fmod(fittedData[i]+2.0*M_PI, 2.0*M_PI) << '\t' << fmod(pFitter.PhaseData()[i], 2.0*M_PI) << '\n';
  }
}

void TestWeightedMedian()
{
  using T = std::vector<std::pair<double,double>>;

  T a({ {1,1}, {2,100}, {3,1}, {4,1}, {5,1}});
  std::cout << "2: " << PieceWisePhaseFitter::WeightedMedian(a) << '\n';

  T b({ {16,5}, {8,0}, {8,0} });
  std::cout << "16: " << PieceWisePhaseFitter::WeightedMedian(b) << '\n';
  
  T c({ {1336,0}, {1336.1,1}, {1336.2,0}, {1336.3,0}, {1337,10}, {1340,1} });
  std::cout << "1337: " << PieceWisePhaseFitter::WeightedMedian(c) << '\n';
  
  T d({ {5,0} });
  std::cout << "5: " << PieceWisePhaseFitter::WeightedMedian(d) << '\n';
}

int main(int argc, char* argv[])
{
  //TestWeightedMedian();
  //readFile(argv[1]);
  std::string filenamePrefix = argv[1];
  bool onlyTECFit = false;
  int seed = atoi(argv[2]);
  size_t nIter = atoi(argv[3]);
  //bool useOldStrategy = atoi(argv[4])!=0;
  std::mt19937 rnd(seed);
  std::uniform_real_distribution<double> dUniform(-1.0, 1.0), dStdDev(0.2, 1.0);
  double maxStddev = 3.5;
  double inputAlpha = dUniform(rnd)*7500.e6, inputBeta = dUniform(rnd)*2.0*M_PI;
  double prevAlpha = dUniform(rnd)*7500.e6, prevBeta = dUniform(rnd)*2.0*M_PI;
  const double startFrequency = 120.0e6; // 120 Mhz
  const double bandwidth = 48e6;
  const double stepSizeStart = 0.2, stepSizeEnd = 0.2;
  std::normal_distribution<double> d(0.0, maxStddev);
  PhaseFitter pFitter(600);
  size_t chunkSize = PieceWisePhaseFitter::CalculateChunkSize(startFrequency, startFrequency+bandwidth, pFitter.Size());
  PieceWisePhaseFitter pwFitter(chunkSize);
  std::vector<double> input(pFitter.Size()), solutions(pFitter.Size());
  std::vector<double> noiseless(pFitter.Size());
  std::vector<double> weights(pFitter.Size()), nu(pFitter.Size()), phases(pFitter.Size());
  std::vector<double> scratch(phases.size());
  
  for(size_t ch=0; ch!=pFitter.Size(); ++ch)
  {
    double noise = dStdDev(rnd);
    nu[ch] = (bandwidth/pFitter.Size())*ch + startFrequency;
    input[ch] = fmod(inputAlpha / nu[ch] + inputBeta + d(rnd)*noise, 2.0*M_PI);
		if(input[ch] < 0.0)
			input[ch]+=2.0*M_PI;
    noiseless[ch] = fmod(inputAlpha / nu[ch] + inputBeta, 2.0*M_PI);
		if(noiseless[ch] < 0.0)
			noiseless[ch]+=2.0*M_PI;
    solutions[ch] = fmod(prevAlpha / nu[ch] + prevBeta + d(rnd), 2.0*M_PI);
		if(solutions[ch] < 0.0)
			solutions[ch]+=2.0*M_PI;
    weights[ch] = 1.0/(noise*noise);
    
    if(ch>pFitter.Size()*3/5 && ch<pFitter.Size()*7/10)
    {
      weights[ch] = 0.0;
      input[ch] = dUniform(rnd)*2.0*M_PI;
    }
  }
  std::copy(nu.begin(), nu.end(), pFitter.FrequencyData());
  std::copy(weights.begin(), weights.end(), pFitter.WeightData());
  
  for(size_t i=0; i!=nIter; ++i)
  {
    double stepSize = stepSizeStart*double(nIter-i)/nIter + stepSizeEnd*double(i)/nIter;
    std::ostringstream filename;
    filename << filenamePrefix << i << ".txt";
    std::ofstream outf(filename.str());
    for(size_t i=0; i!=pFitter.Size(); ++i)
    {
      std::complex<double> moveTo = std::polar(1.0, input[i]);
      std::complex<double> moveFrom = std::polar(1.0, solutions[i]);
      std::complex<double> newVal = moveTo*stepSize + moveFrom*(1.0-stepSize);
      solutions[i] = std::arg(newVal);
      phases[i] = fmod(solutions[i] + 2.0*M_PI, 2.0*M_PI);
    }
    
    if(i<nIter*3/4 && !onlyTECFit)
    {
      pwFitter.SlidingFit(nu.data(), phases, weights.data(), scratch);
      std::copy(scratch.begin(), scratch.end(), phases.begin());
    }
    else {
      std::copy(phases.begin(), phases.end(), pFitter.PhaseData());
      pFitter.FitDataToTEC2Model();
      std::copy(pFitter.PhaseData(), pFitter.PhaseData()+pFitter.Size(), phases.begin());
    }
  
    for(size_t i=0; i!=pFitter.Size(); ++i)
    {
      double ph = fmod(phases[i], 2.0*M_PI);
      outf << nu[i] << '\t';
      if(weights[i] > 0.0)
        outf << fmod(input[i]+2.0*M_PI, 2.0*M_PI);
      else
        outf << "?";
      outf << '\t' << fmod(ph+2.0*M_PI, 2.0*M_PI) << '\t' << noiseless[i] << '\n';
    }
  }
}
