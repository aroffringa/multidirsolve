#include "multidirsolver.h"

int main(int argc, char* argv[])
{
	typedef std::complex<float> cf;
	MultiDirSolver mds(1000, 1e-6, 0.5);
	size_t nPol = 4, nAnt = 200, nDir = 3, nChan = 10, nChanBlocks = 1, nTimes = 1, nBl=nAnt*(nAnt-1)/2;
	
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
	
	//mds.set_mode(MultiDirSolver::CalibrateTEC2);
	
	std::vector<double> nu(nChanBlocks);
	for(size_t i=0; i!=nChanBlocks; ++i)
		nu[i] = i+1;
	
	TECConstraint tecConstraint(TECConstraint::TECOnlyMode, nAnt, nDir, nChanBlocks, nu.data());
	
	//mds.add_constraint(&tecConstraint);
	
	cf gain1(0.31415926535, 0.0), gain2(2.0, 1.0), gain3(0.0, 3.0);
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
			model1Ptr[i*4 + 0] = 1.0;
			model1Ptr[i*4 + 1] = 0.0;
			model1Ptr[i*4 + 2] = 0.0;
			model1Ptr[i*4 + 3] = 1.0;
			model2Ptr[i*4 + 0] = (i%2==0) ? cf(1.0, 0.0) : 0.0;
			model2Ptr[i*4 + 1] = 0.0;
			model2Ptr[i*4 + 2] = 0.0;
			model2Ptr[i*4 + 3] = (i%2==0) ? cf(1.0, 0.0) : 0.0;
			model3Ptr[i*4 + 0] = (i%3==0) ? cf(1.0, 0.0) : 0.0;
			model3Ptr[i*4 + 1] = 0.0;
			model3Ptr[i*4 + 2] = 0.0;
			model3Ptr[i*4 + 3] = (i%3==0) ? cf(1.0, 0.0) : 0.0;
		}
		
		size_t baselineIndex = 0;
		for(size_t a1=0; a1!=nAnt; ++a1)
		{
			for(size_t a2=a1+1; a2!=nAnt; ++a2)
			{
				for(size_t j=0; j!=nPol * nChan; ++j)
				{
					cf gain1ant1, gain2ant1, gain3ant1, gain1ant2, gain2ant2, gain3ant2;
					if(a1 == 1)
					{
						gain1ant1 = 1.0;
						gain2ant1 = 1.0;
						gain3ant1 = 1.0;
					}
					else {
						gain1ant1 = a1;
						gain2ant1 = gain2;
						gain3ant1 = gain3;
					}
					gain1ant2 = a2;
					gain2ant2 = gain2;
					gain3ant2 = gain3;
					dataPtr[baselineIndex] =
						gain1ant1 * gain1ant2 * model1Ptr[baselineIndex] +
						gain2ant1 * gain2ant2 * model2Ptr[baselineIndex] +
						gain3ant1 * gain3ant2 * model3Ptr[baselineIndex];
					++baselineIndex;
				}
			}
		}
		
		data.push_back(dataPtr);
		modelData.push_back(std::vector<cf*>{model1Ptr, model2Ptr, model3Ptr});
	}
	
	std::vector<std::vector<std::complex<double> > > solutions(nChanBlocks);
	for(auto& vec : solutions)
		vec.assign(nDir * nAnt, 1.0);
	MultiDirSolver::SolveResult result = mds.process(data, modelData, solutions, 0.0);
	std::cout << '\n';
	for(size_t ch=0; ch!=nChanBlocks; ++ch)
	{
		for(size_t ant=0; ant!=nAnt; ++ant)
		{
			for(size_t d=0; d!=nDir; ++d)
				std::cout << "ch" << ch << ", ant" << ant << ", d" << d << " = " << solutions[ch][d + ant*nDir] << '\n';
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
}

