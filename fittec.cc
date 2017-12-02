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
  pwFitter.SlidingFit(frequencies.data(), phases.data(), fittedData.data(), phases.size()); 
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
  int seed = atoi(argv[2]);
  size_t nIter = atoi(argv[3]);
  std::string methodStr = argv[4];
  std::string titleMethod = argv[5];
  enum Method { Best, TECOnly, NoBreak, PWOnly, CStep, ComplexNoBreak } method = Best;
  if(methodStr == "best") method = Best;
  else if(methodStr == "teconly") method = TECOnly;
  else if(methodStr == "nobreak") method = NoBreak;
  else if(methodStr == "pwonly") method = PWOnly;
  else if(methodStr == "cstep") method = CStep;
  else if(methodStr == "complexnobreak") method = ComplexNoBreak;
  else std::cerr << "Invalid method!\n";
  bool showTitle = false;
  if(titleMethod == "title")
    showTitle = true;
  //bool useOldStrategy = atoi(argv[4])!=0;
  std::mt19937 rnd(seed);
  std::uniform_real_distribution<double>
    dUniform(-1.0, 1.0), dStdDev(0.2, 1.0), dUnit(0.0, 1.0);
  double maxStddev = dUnit(rnd)*5.0+0.001;
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
  std::vector<double> beforeConstraint;
  
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
  
  for(size_t iter=0; iter!=nIter; ++iter)
  {
    double stepSize = stepSizeStart*double(nIter-iter)/nIter + stepSizeEnd*double(iter)/nIter;
    std::ostringstream filename;
    filename << filenamePrefix << iter << ".txt";
    std::ofstream outf(filename.str());
    
    for(size_t i=0; i!=pFitter.Size(); ++i)
    {
      if(method == CStep || method == ComplexNoBreak)
      {
        std::complex<double> moveTo = std::polar(1.0, input[i]);
        std::complex<double> moveFrom = std::polar(1.0, solutions[i]);
        std::complex<double> newVal = moveTo*stepSize + moveFrom*(1.0-stepSize);
        solutions[i] = std::arg(newVal);
      }
      else {
        double moveTo = fmod(fmod(input[i] + 2.0*M_PI, 2.0*M_PI) + 2.0*M_PI, 2.0*M_PI);
        double moveFrom = fmod(fmod(solutions[i] + 2.0*M_PI, 2.0*M_PI) + 2.0*M_PI, 2.0*M_PI);
        double distance = moveTo - moveFrom;
        if(distance > M_PI) distance = distance - 2.0*M_PI;
        else if(distance < -M_PI) distance = distance + 2.0*M_PI;
        solutions[i] = solutions[i] + stepSize * distance;
      }
      phases[i] = fmod(fmod(solutions[i] + 2.0*M_PI, 2.0*M_PI) + 2.0*M_PI, 2.0*M_PI);
    }
    
    // APPLY CONSTRAINTS
    //
    beforeConstraint = phases;
    if((iter<nIter*3/4 || method==PWOnly) && method!=TECOnly)
    {
      if(method != NoBreak && method != ComplexNoBreak)
      {
        size_t bp =
          pwFitter.SlidingFitWithBreak(nu.data(), phases.data(), weights.data(), scratch.data(), phases.size());
        std::ostringstream bpfilename;
        bpfilename << filenamePrefix << "-bp-" << iter << ".txt";
        std::ofstream bpoutf(bpfilename.str());
        bpoutf << nu[bp] << "\t0\n" << nu[bp] << "\t" << 2*M_PI << "\n";
      }
      else
        pwFitter.SlidingFit(nu.data(), phases.data(), weights.data(), scratch.data(), phases.size());
      phases = scratch;
    }
    else {
      std::copy(phases.begin(), phases.end(), pFitter.PhaseData());
      pFitter.FitDataToTEC2Model();
      std::copy(pFitter.PhaseData(), pFitter.PhaseData()+pFitter.Size(), phases.begin());
    }

    for(double& v : phases)
      v = fmod(fmod(v, 2.0*M_PI)+2.0*M_PI, 2.0*M_PI);
    std::copy(phases.begin(), phases.end(), solutions.begin());
  
    for(size_t i=0; i!=pFitter.Size(); ++i)
    {
      double ph = fmod(phases[i], 2.0*M_PI);
      outf << nu[i] << '\t';
      if(weights[i] > 0.0)
        outf << fmod(input[i]+2.0*M_PI, 2.0*M_PI);
      else
        outf << "?";
      outf 
        << '\t' << ph
        << '\t' << noiseless[i] 
        << '\t' << beforeConstraint[i]
        << '\n';
    }
  }
  
  double sqError = 0.0;
  for(size_t i=0; i!=phases.size(); ++i)
  {
    double distance = fmod(std::fabs(phases[i] - input[i]), 2.0*M_PI);
    if(distance > M_PI)
      distance = 2.0 * M_PI - distance;
    sqError += distance*distance;
  }
  sqError /= phases.size();
  if(showTitle)
  {
    std::cout << "TEC=" << round(PhaseFitter::AlphaToTEC(inputAlpha)*10.0)/10.0 << ", max stddev=" << round(maxStddev*100.0)/100.0 << ", error=" << round(sqError*100.0)/100.0 << '\n';
  }
  else {
    std::cout <<
      PhaseFitter::AlphaToTEC(inputAlpha) << '\t' <<
      maxStddev << '\t' <<
      sqError << '\t' <<
      PhaseFitter::AlphaToTEC(-std::fabs(inputAlpha - prevAlpha)) << '\n';
  }
}
