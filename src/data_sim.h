// Templated class to simulate a dataset with one or more sample sizes

//#include <Rcpp.h>

//#include "enums.h"
#include "distribution.h"

/*

Tasks to separate:

- Drawing pairs of data points (pre and post) based on different distributions, potentially accepting a vector of reductions to give a vector of post ?? (not yet)
  - Do this using templates as for eggSim
- Keeping track of running mean/variance (for each N & potentially reduction) and also estimating k by ML if/when necessary
  - Use inheritance pattern:  simple_sim tracks only mu/var, ml_sim derives from simple_sim and also estimates k by ML
- Looping over iterations and obtaining BNB / WAAVP / Levecke / Dobson method results
  - Plain function




// Base class handles the data simulation part, derived class changes how much data is stored

template<rdists distribution>
class data_sim
{
protected:
  data_sim()
  {
    
  }
public:
  void draw()
  {
    
  }
};


template<typename t, rdists distribution>
class data_sim_child : public data_sim<distribution>
{
  
};

*/