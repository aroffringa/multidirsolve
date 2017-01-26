#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#ifdef AOPROJECT
#include "phasefitter.h"
#define UPTR std::unique_ptr
#else
#include <DPPP/phasefitter.h>
#define UPTR std::auto_ptr
#endif

#include <complex>
#include <vector>
#include <memory>

class Constraint
{
public:
	typedef std::complex<double> dcomplex;
	
  struct Result { };

  virtual ~Constraint() { };
   
  virtual void init(size_t nAntennas, size_t nDirections, 
                    size_t nChannelBlocks, const double* frequencies) = 0;

  virtual Result Apply(std::vector<std::vector<dcomplex> >& solutions) = 0;
};

class PhaseConstraint : public Constraint
{
public:
  PhaseConstraint() {};

  virtual void init(size_t, size_t, size_t, const double*) {};

  virtual Result Apply(std::vector<std::vector<dcomplex> >& solutions);
};

class TECConstraint : public Constraint
{
public:
  struct TECResult : public Constraint::Result {
    std::vector<std::pair<double, double> > tecvals;
  };

  enum Mode { TECAndCommonScalarMode, TECOnlyMode };
  
  TECConstraint(Mode mode, size_t nAntennas, size_t nDirections, 
                size_t nChannelBlocks, const double* frequencies);
  TECConstraint(Mode mode);

  virtual void init(size_t nAntennas, size_t nDirections, 
                    size_t nChannelBlocks, const double* frequencies);
  
  virtual Result Apply(std::vector<std::vector<dcomplex> >& solutions);
  
private:
  Mode _mode;
  size_t _nAntennas, _nDirections, _nChannelBlocks;
  std::vector<PhaseFitter> _phaseFitters;
};

#endif
