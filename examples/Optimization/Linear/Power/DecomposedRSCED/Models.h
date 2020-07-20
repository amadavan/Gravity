//
// Created by Avinash Madavan on 6/11/20.
//

#ifndef GRAVITY_EXAMPLES_OPTIMIZATION_LINEAR_POWER_DECOMPOSEDRSCED_MODELS_H_
#define GRAVITY_EXAMPLES_OPTIMIZATION_LINEAR_POWER_DECOMPOSEDRSCED_MODELS_H_

#include <PowerNet.h>
#include <gravity/solver.h>

namespace gravity {

// TODO: This data should be included within the network model.
struct RSCEDdata {
  param<> VoLL = param<>("VoLL");
  param<> ramp_max = param<>("ramp_max");
  double DAL_multiplier = 1.8;
  double STE_multiplier = 1.3;
  double failure_probability;
};

extern Model<> createNominalModel(PowerNet &grid, const RSCEDdata &data, const double risk_aversion);
extern Model<> createContingencyModel(PowerNet &grid,
                                      const RSCEDdata &data,
                                      const double risk_aversion,
                                      const std::pair<std::string, std::pair<Arc *, Gen *>> &contingency);

}

#endif //GRAVITY_EXAMPLES_OPTIMIZATION_LINEAR_POWER_DECOMPOSEDRSCED_MODELS_H_
