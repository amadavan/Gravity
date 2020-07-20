//
// Created by Avinash Madavan on 6/11/20.
//

#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif

#include "Models.h"

int main(int argc, char *argv[]) {
  int output = 0;
  bool use_cplex = false, use_gurobi = false;
  double tol = 1e-6;
  double solver_time_end, total_time_end, solve_time, total_time;
  string mehrotra = "no", log_level = "0";
  string fname = string(prj_dir) + "/data_sets/Power/nesta_case5_pjm.m";

  string path = argv[0];
  string solver_str = "ipopt";

  // TODO: properly include CVaR parameter
  double cvar_param = 0.;

#ifdef USE_OPT_PARSER
  /** Create a OptionParser with options */
  auto options = readOptions(argc, argv);
  options.add_option("f", "file", "Input file name", fname);
  options.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);
  options.add_option("a", "alpha", "CVaR parameter", std::to_string(cvar_param));

  /** Parse the options and verify that all went well. If not, errors and help will be shown */
  bool correct_parsing = options.parse_options(argc, argv);

  if (!correct_parsing) {
    return EXIT_FAILURE;
  }

  output = op::str2int(options["l"]);

  fname = options["f"];
  bool has_help = op::str2bool(options["h"]);
  if (has_help) {
    options.show_help();
    exit(0);
  }
  solver_str = options["s"];
  if (solver_str.compare("gurobi") == 0) {
    use_gurobi = true;
  } else if (solver_str.compare("cplex") == 0) {
    use_cplex = true;
  }
  cvar_param = op::str2double(options["a"]);
#else
  if(argc==2){
    fname=argv[1];
  }
  else{
    fname=string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
  }
#endif

  double total_time_start = get_wall_time();
  PowerNet grid;
  grid.readgrid(fname, true);

  double max_cost = 0, gen_cost;
  for (size_t i = 0; i < grid.get_nb_active_gens(); ++i) {
    grid.c1.get_double_val(i, gen_cost);
    max_cost = std::max(max_cost, gen_cost);
  }

  RSCEDdata data;
  data.VoLL.in(indices(grid.nodes));
  data.VoLL = 10 * max_cost; // TODO: set VoLL
  data.failure_probability = 0.2 / grid.get_nb_active_arcs(); // 20% chance of off-nominal case
  data.ramp_max.in(indices(grid.gens));
  data.ramp_max = grid.pg_max; // TODO: set ramp_max

  Model<> NominalModel = createNominalModel(grid, data, cvar_param);

  /** Handling contingencies */
  std::vector<std::pair<std::string, std::pair<Arc *, Gen *>>> conting_lines;
  for (auto cont : grid.arcs)
    conting_lines.emplace_back("cont_" + to_string(cont->_id), std::make_pair(cont, nullptr));

  std::vector<Model<>> ContingencyModels{};
  ContingencyModels.resize(grid.get_nb_active_arcs());
  for (auto cont : grid.arcs) {
    std::pair<std::string, std::pair<Arc *, Gen *>> contingency = std::make_pair(to_string(cont->_id), std::make_pair(cont, nullptr));
//    Model<> ContingencyModel = createContingencyModel(grid, data, cvar_param, contingency);
    ContingencyModels.emplace_back(createContingencyModel(grid, data, cvar_param, contingency));
  }
}