/*Hoi André,

  Na iets meer gedachten hierover zou ik de interface een tikkeltje willen aanpassen, om iets meer aan te sluiten bij de DPPP-filosofie. Ik ben het aan het implementeren, kan zijn dat het niet afkomt tijdens de vlucht. Ik stuur m'n idee nu vast zodat jij niet onnodig werk zit te doen. 

  Stel dat jouw klasse DDESolver heet, dan zou de constructor de parset mee krijgen:*/

class MultiDirSolver {

 public:
  
  MultiDirSolver(const Parset& parset, HDF5bestand*);

  init(size_t nants, size_t ndir, size_t nchan);

  // TODO this should receive weights!
  // Float per vis or per pol x vis?
  process(vector<DComplex*> data, vector<float*> data, vector<vector<DComplex* > > mdata);

  //  -- eventuele opruimacties (bijvoorbeeld wegschrijven van de data)
  finish();

  //  -- hier kun je wat op het scherm dumpen aan statistieken
  showCounts();

};
