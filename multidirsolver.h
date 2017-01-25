#ifndef MULTI_DIR_SOLVER_H
#define MULTI_DIR_SOLVER_H

#ifdef AOPROJECT
#include "phasefitter.h"
#define UPTR std::unique_ptr
#else
#include <DPPP/phasefitter.h>
#define UPTR std::auto_ptr
#endif

#include <armadillo>

#include <complex>
#include <vector>
#include <memory>

class MultiDirSolver
{
public:
	
	typedef std::complex<double> DComplex;
	typedef std::complex<float> Complex;
	
	class Constraint
	{
	public:
		virtual ~Constraint() { };
		
		virtual void init(size_t nAntennas, size_t nDirections, size_t nChannelBlocks, const double* frequencies) = 0;

		virtual void Apply(std::vector<std::vector<DComplex> >& solutions) = 0;
	};
	
	struct ConstraintResult { };
	
	struct TECConstraintResult : public ConstraintResult
	{
		// This is an example interface how to return things
		// TODO make the TEC constraint do something with this.
		double alpha, beta;
	};
	
	struct SolveResult {
		size_t iterations;
		std::vector<ConstraintResult*> _results;
	};
	
	enum CalibrationMode { CalibrateComplexGain, CalibrateTEC1, CalibrateTEC2, CalibratePhase };
	
  MultiDirSolver(size_t maxIterations, double accuracy, double stepSize);
	
	void init(size_t nAntennas, size_t nDirections, size_t nChannels, const std::vector<int>& ant1, const std::vector<int>& ant2);
	
	// data[i] is een pointer naar de data voor tijdstap i, vanaf die pointer staat het in volgorde als in MS (bl, chan, pol)
	// mdata[i] is een pointer voor tijdstap i naar arrays van ndir model data pointers (elk van die data pointers staat in zelfde volgorde als data)
	// solutions[ch] is een pointer voor channelblock ch naar antenna x directions oplossingen.
	SolveResult process(std::vector<Complex*>& data, std::vector<std::vector<Complex* > >& modelData,
		std::vector<std::vector<DComplex> >& solutions) const;
	
	void set_mode(CalibrationMode mode) { _mode = mode; }
	
	void set_channel_blocks(size_t nChannelBlocks) { _nChannelBlocks = nChannelBlocks; }
	
	void set_max_iterations(size_t maxIterations) { _maxIterations = maxIterations; }
	
	void set_accuracy(double accuracy) { _accuracy = accuracy; }
	
	void set_step_size(double stepSize) { _stepSize = stepSize; }
	
	void add_constraint(Constraint* constraint) { _constraints.push_back(constraint); }
	
private:
	void performSolveIteration(size_t channelBlockIndex,
	                           std::vector<arma::cx_mat>& gTimesCs,
	                           std::vector<arma::cx_vec>& vs,
	                           const std::vector<DComplex>& solutions,
	                           std::vector<DComplex>& nextSolutions,
	                           const std::vector<Complex *>& data,
	                           const std::vector<std::vector<Complex *> >& modelData) const;
	
	size_t _nAntennas, _nDirections, _nChannels, _nChannelBlocks;
	std::vector<int> _ant1, _ant2;
	
	// Calibration setup
	enum CalibrationMode _mode;
	size_t _maxIterations;
	double _accuracy;
	double _stepSize;
	std::vector<Constraint*> _constraints;
};

class PhaseConstraint : public MultiDirSolver::Constraint
{
public:
    PhaseConstraint() {};

    virtual void init(size_t, size_t, size_t, const double*) {};

    virtual void Apply(std::vector<std::vector<MultiDirSolver::DComplex> >& solutions);
};

class TECConstraint : public MultiDirSolver::Constraint
{
public:
	enum Mode { TECAndCommonScalarMode, TECOnlyMode };
	
	TECConstraint(Mode mode, size_t nAntennas, size_t nDirections, size_t nChannelBlocks, const double* frequencies);
	TECConstraint(Mode mode);

	virtual void init(size_t nAntennas, size_t nDirections, size_t nChannelBlocks, const double* frequencies);
	
	virtual void Apply(std::vector<std::vector<MultiDirSolver::DComplex> >& solutions);
	
private:
	Mode _mode;
	size_t _nAntennas, _nDirections, _nChannelBlocks;
	std::vector<PhaseFitter> _phaseFitters;
};

#endif