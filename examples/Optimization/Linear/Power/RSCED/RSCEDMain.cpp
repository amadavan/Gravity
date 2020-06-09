//
// RSCED.cpp
// Gravity
//
// Created by Avinash Madavan on 5/20/20.
//
//

#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif

int main(int argc, char *argv[]) {
  int output = 0;
  bool projected = false, use_cplex = false, use_gurobi = false;
  double tol = 1e-6;
  double solver_time_end, total_time_end, solve_time, total_time;
  string mehrotra = "no", log_level = "0";
  string fname = string(prj_dir) + "/data_sets/Power/nesta_case5_pjm.m";

  string path = argv[0];
  string solver_str = "ipopt";
  string proj_str = "0";

  // TODO: properly include CVaR parameter
  double cvar_param = 0.;

#ifdef USE_OPT_PARSER
  /** Create a OptionParser with options */
  auto options = readOptions(argc, argv);
  options.add_option("f", "file", "Input file name", fname);
  options.add_option("p", "project", "Project the power flow variables", proj_str);
  options.add_option("s", "solver", "Solvers: ipopt/cplex/gurobi, default = ipopt", solver_str);

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
  solver_str = options["p"];
  if (solver_str.compare("1") == 0) {
    projected = true;
  }
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

  /* Grid Stats */
  auto node_pairs = grid.get_node_pairs();
  auto nb_gen = grid.get_nb_active_gens();
  auto nb_lines = grid.get_nb_active_arcs();
  auto nb_buses = grid.get_nb_active_nodes();
  DebugOn("nb active gens = " << nb_gen << endl);
  DebugOn("nb active lines = " << nb_lines << endl);
  DebugOn("nb active buses = " << nb_buses << endl);

  /** Sets */
  auto nodes = indices(grid.nodes);
  auto arcs = indices(grid.arcs);
  auto gens = indices(grid.gens);
  auto gen_nodes = grid.gens_per_node();
  auto out_arcs = grid.out_arcs_per_node();
  auto in_arcs = grid.in_arcs_per_node();

  /* Grid Parameters */
  auto pg_min = grid.pg_min.in(gens);
  auto pg_max = grid.pg_max.in(gens);
  auto pl = grid.pl.in(nodes);
  auto c1 = grid.c1.in(gens);
  auto c2 = grid.c2.in(gens);
  auto c0 = grid.c0.in(gens);
  auto gs = grid.gs.in(nodes);
  auto S_max = grid.S_max.in(arcs);
  auto b = grid.b.in(arcs);
  auto th_min = grid.th_min.in(node_pairs);
  auto th_max = grid.th_max.in(node_pairs);

  // TODO: Include value of lost load (voll)
  param<> voll("VoLL");
  voll.in(nodes);
  voll = 10;

  // TODO: DA/STE limits for line capacities S_max
  double DAL_multiplier = 1.5;
  double STE_multiplier = 1.2;

  // TODO: Line failure probabilities
  double failure_probability = 0.05;

  // TODO: Include CVaR parameter (alpha' = 1/(1 - alpha))
  double alpha_prime = 1. / (1. - cvar_param);

  // TODO: Include generation ramp capacity
  auto ramp_max = pg_max;

  /** Declare model */
  Model<> RSCED("RSCED Model");

  /** Variables */
  /* Power generation variables */
  var<> Pg("Pg", pg_min, pg_max);
  RSCED.add(Pg.in(gens));

  /* Phase angle variables */
  var<> theta("𝛉");
  RSCED.add(theta.in(nodes));

  /* CVaR reformulation variable */
  // y0, z0(?)
  var<> z("z");
  var<> y("y");
  RSCED.add(z);
  RSCED.add(y);
  z.add_lb_only(0);
  y.add_lb_only(0);

  /** Constraints */
  /* CVaR lower bound */
  Constraint<> OBJ_BOUND("OBJ_BOUND");
  OBJ_BOUND = y + z - product(c1, Pg);
  RSCED.add(OBJ_BOUND >= 0);

  /* REF BUS */
  Constraint<> Ref_Bus("Ref_Bus");
  Ref_Bus = theta(grid.ref_bus);
  RSCED.add(Ref_Bus == 0);

  /* Power balance constraint */
  Constraint<> KCL_P("KCL_P");
  KCL_P = b.tr().in(in_arcs) * (theta.from(in_arcs) - theta.to(in_arcs))
      - b.tr().in(out_arcs) * (theta.from(out_arcs) - theta.to(out_arcs)) + pl + gs - sum(Pg, gen_nodes);
  RSCED.add(KCL_P.in(nodes) == 0);

  /* Phase Angle Bounds constraints */
  Constraint<> PAD_UB("PAD_UB");
  PAD_UB = theta.from(node_pairs) - theta.to(node_pairs);
  PAD_UB -= th_max;
  RSCED.add(PAD_UB.in(node_pairs) <= 0);
  Constraint<> PAD_LB("PAD_LB");
  PAD_LB = theta.from(node_pairs) - theta.to(node_pairs);
  PAD_LB -= th_min;
  RSCED.add(PAD_LB.in(node_pairs) >= 0);

  /* Line Limits constraints */
  Constraint<> Thermal_UB("Thermal_UB");
  Thermal_UB = b * (theta.to(arcs) - theta.from(arcs));
  Thermal_UB -= S_max;
  RSCED.add(Thermal_UB.in(arcs) <= 0);
  Constraint<> Thermal_LB("Thermal_LB");
  Thermal_LB = b * (theta.to(arcs) - theta.from(arcs));
  Thermal_LB += S_max;
  RSCED.add(Thermal_LB.in(arcs) >= 0);

  /** Objective */
  auto obj = z + alpha_prime * (1 - nb_lines * failure_probability) * y;

  /** Handling contingencies */
  std::vector<std::pair<std::string, std::pair<Arc *, Gen *>>> conting_lines;
  for (auto cont : grid.arcs)
    conting_lines.emplace_back("cont_" + to_string(cont->_id), std::make_pair(cont, nullptr));

  auto arcs_c = grid.get_conting_arcs(conting_lines);

  indices contingencies("conting");
  for (auto cont : conting_lines)
    contingencies.insert(cont.first);

  auto node_pairs_c = grid.get_node_pairs_cont(conting_lines);
  auto nodes_c = indices(contingencies, nodes);
  auto gens_c = indices(contingencies, gens);

  auto gen_nodes_c1 = grid.gens_per_node_cont1(conting_lines, gens);
  auto gen_nodes_c = grid.gens_per_node_cont(conting_lines, gens_c);
  auto out_arcs_c = grid.out_arcs_per_node_cont(conting_lines, arcs_c);
  auto in_arcs_c = grid.in_arcs_per_node_cont(conting_lines, arcs_c);

  auto matrix_gens_c = grid.get_gens_cont(conting_lines, gens);
  auto matrix_nodes_c = grid.get_nodes_cont(conting_lines, nodes);
  auto matrix_gens_c2 = grid.get_gens_cont2(conting_lines, gens_c);
  auto matrix_nodes_c2 = grid.get_nodes_cont2(conting_lines, nodes_c);

  auto c1_c = grid.c1.in(matrix_gens_c);
  auto voll_c = voll.in(matrix_nodes_c);

  auto S_max_c = grid.S_max.from_ith(1, arcs_c);
  auto b_c = grid.b.from_ith(1, arcs_c);
  auto th_max_c = th_max.from_ith(1, node_pairs_c);
  auto th_min_c = th_min.from_ith(1, node_pairs_c);
  auto pg_max_c = pg_max.from_ith(1, gens_c);
  auto pg_min_c = pg_min.from_ith(1, gens_c);
  auto ramp_max_c = ramp_max.from_ith(1, gens_c);

  auto pl_c = pl.from_ith(1, nodes_c);
  auto gs_c = gs.from_ith(1, nodes_c);

  auto Pg_nom = Pg.from_ith(1, gens_c);

  /** Contingency variable declaration */
  /* Positive recourse generation */
  var<> dPg_p("dPg_p");
  RSCED.add(dPg_p.in(gens_c));

  /* Negative recourse generation */
  var<> dPg_n("dPg_n");
  RSCED.add(dPg_n.in(gens_c));

  /* Power flows */
  var<> Pf_c("Pf_c");
  RSCED.add(Pf_c.in(arcs_c));
  var<> Pf_recourse_c("Pf_recourse_c");
  RSCED.add(Pf_recourse_c.in(arcs_c));

  /* Phase angle variables */
  var<> theta_c("theta_c");
  RSCED.add(theta_c.in(nodes_c));
  var<> theta_recourse_c("theta_recourse_c");
  RSCED.add(theta_recourse_c.in(nodes_c));

  /* Load shed variables */
  var<> dD_c("dD_c");
  RSCED.add(dD_c.in(nodes_c));

  /* CVaR reformulation variable */
  var<> y_c("y_c");
  RSCED.add(y_c.in(contingencies));
  y_c.add_lb_only(0);

  /** Post-failure constraints */
  /* REF BUS */
  Constraint<> Ref_Bus_c("Ref_Bus_c");
  Ref_Bus_c = theta_c.in(grid.ref_bus);
  RSCED.add(Ref_Bus_c == 0);

  /* Power flow constraint */
  Constraint<> Flow_P_c("Flow_P_c");
  Flow_P_c = Pf_c + b_c * (theta_c.from(arcs_c) - theta_c.to(arcs_c));
  RSCED.add(Flow_P_c.in(arcs_c) == 0);

  /* Power balance constraint */
  Constraint<> KCL_P_c("KCL_P_c");
  KCL_P_c = sum(Pf_c, out_arcs_c) - sum(Pf_c, in_arcs_c) + pl_c + gs_c - sum(Pg, gen_nodes_c1);
  RSCED.add(KCL_P_c.in(nodes_c) == 0);

  /* Phase Angle Bounds constraints */
  Constraint<> PAD_UB_c("PAD_UB_c");
  PAD_UB_c = theta_c.from_ith(0, node_pairs_c) - theta_c.in_ignore_ith(1, 1, node_pairs_c);
  PAD_UB_c -= th_max.from_ith(1, node_pairs_c);
  RSCED.add(PAD_UB_c.in(node_pairs_c) <= 0);
  Constraint<> PAD_LB_c("PAD_LB_c");
  PAD_LB_c = theta_c.from_ith(0, node_pairs_c) - theta_c.in_ignore_ith(1, 1, node_pairs_c);
  PAD_LB_c -= th_min.from_ith(1, node_pairs_c);
  RSCED.add(PAD_LB_c.in(node_pairs_c) >= 0);

  /* Line Limits constraints */
  Constraint<> Thermal_UB_c("Thermal_UB_c");
  Thermal_UB_c = b_c * (theta_c.to(arcs_c) - theta_c.from(arcs_c));
  Thermal_UB_c -= S_max_c * DAL_multiplier;
  RSCED.add(Thermal_UB_c.in(arcs_c) <= 0);
  Constraint<> Thermal_LB_c("Thermal_LB_c");
  Thermal_LB_c = b_c * (theta_c.to(arcs_c) - theta_c.from(arcs_c));
  Thermal_LB_c += S_max_c * DAL_multiplier;
  RSCED.add(Thermal_LB_c.in(arcs_c) >= 0);

  /** Post-recourse constraints */
  /* REF BUS */
  Constraint<> Ref_Bus_recourse_c("Ref_Bus_recourse_c");
  Ref_Bus_recourse_c = theta_c.in(grid.ref_bus);
  RSCED.add(Ref_Bus_recourse_c == 0);

  /* CVaR lower bound */
  Constraint<> OBJ_BOUND_c("OBJ_BOUND_c");
  OBJ_BOUND_c =
      y_c.in(contingencies) + z - product(c1_c, Pg.in(matrix_gens_c)) - product(c1_c, dPg_p.in(matrix_gens_c2))
          - product(c1_c, dPg_n.in(matrix_gens_c2)) - product(voll_c, dD_c.in(matrix_nodes_c2));
  RSCED.add(OBJ_BOUND_c.in(contingencies) >= 0);

  /* Power flow constraint */
  Constraint<> Flow_P_recourse_c("Flow_P_recourse_c");
  Flow_P_recourse_c = Pf_recourse_c + b_c * (theta_recourse_c.from(arcs_c) - theta_recourse_c.to(arcs_c));
  RSCED.add(Flow_P_recourse_c.in(arcs_c) == 0);

  /* Power balance constraint */
  Constraint<> KCL_P_recourse_c("KCL_P_recourse_c");
  KCL_P_recourse_c =
      sum(Pf_recourse_c, out_arcs_c) - sum(Pf_recourse_c, in_arcs_c) + pl_c + gs_c - sum(Pg, gen_nodes_c1)
          - sum(dPg_p, gen_nodes_c) + sum(dPg_n, gen_nodes_c) - dD_c;
  RSCED.add(KCL_P_recourse_c.in(nodes_c) == 0);

  /* Phase Angle Bounds constraints */
  Constraint<> PAD_UB_recourse_c("PAD_UB_recourse_c");
  PAD_UB_c = theta_recourse_c.from_ith(0, node_pairs_c) - theta_recourse_c.in_ignore_ith(1, 1, node_pairs_c);
  PAD_UB_c -= th_max.from_ith(1, node_pairs_c);
  RSCED.add(PAD_UB_recourse_c.in(node_pairs_c) <= 0);
  Constraint<> PAD_LB_recourse_c("PAD_LB_recourse_c");
  PAD_LB_c = theta_recourse_c.from_ith(0, node_pairs_c) - theta_recourse_c.in_ignore_ith(1, 1, node_pairs_c);
  PAD_LB_c -= th_min.from_ith(1, node_pairs_c);
  RSCED.add(PAD_LB_recourse_c.in(node_pairs_c) >= 0);

  /* Line Limits constraints */
  Constraint<> Thermal_UB_recourse_c("Thermal_UB_recourse_c");
  Thermal_UB_recourse_c = b_c * (theta_recourse_c.to(arcs_c) - theta_recourse_c.from(arcs_c));
  Thermal_UB_recourse_c -= S_max_c * STE_multiplier;
  RSCED.add(Thermal_UB_recourse_c.in(arcs_c) <= 0);
  Constraint<> Thermal_LB_recourse_c("Thermal_LB_recourse_c");
  Thermal_LB_recourse_c = b_c * (theta_recourse_c.to(arcs_c) - theta_recourse_c.from(arcs_c));
  Thermal_LB_recourse_c += S_max_c * STE_multiplier;
  RSCED.add(Thermal_LB_recourse_c.in(arcs_c) >= 0);

  /* Generation capacity limits */
  Constraint<> Gen_UB_c("Gen_recourse_UB_c");
  Gen_UB_c = Pg_nom + dPg_p - dPg_n;
  Gen_UB_c -= pg_max_c;
  RSCED.add(Gen_UB_c.in(gens_c) <= 0);
  Constraint<> Gen_LB_c("Gen_recourse_LB_c");
  Gen_LB_c = Pg_nom + dPg_p - dPg_n;
  Gen_LB_c += pg_min_c;
  RSCED.add(Gen_LB_c.in(gens_c) >= 0);

  /* Recourse ramp capacity limits */
  Constraint<> Recourse_UB_c("Recourse_UB_c");
  Recourse_UB_c = dPg_p - dPg_n;
  Recourse_UB_c -= ramp_max_c;
  RSCED.add(Recourse_UB_c.in(gens_c) <= 0);
  Constraint<> Recourse_LB_c("Recourse_LB_c");
  Recourse_LB_c = dPg_p - dPg_n;
  Recourse_LB_c += ramp_max_c;
  RSCED.add(Recourse_LB_c.in(gens_c) >= 0);

  /** Objective */
  obj += alpha_prime * failure_probability * sum(y_c, contingencies);

  /* Set objective */
  RSCED.min(obj);

  RSCED.print();

  /** Solve */
  solver<> RSCED_SOLVER(RSCED, ipopt);

  if (use_gurobi)
    RSCED_SOLVER = solver<>(RSCED, gurobi);
  else if (use_cplex)
    RSCED_SOLVER = solver<>(RSCED, cplex);
  else
    RSCED_SOLVER = solver<>(RSCED, ipopt);

  auto solver_time_start = get_wall_time();
  RSCED_SOLVER.run(output = 5, tol = 1e-6);
  solver_time_end = get_wall_time();
  total_time_end = get_wall_time();
  solve_time = solver_time_end - solver_time_start;
  total_time = total_time_end - total_time_start;

  /** Uncomment next line to print expanded model */
  /* DCOPF.print(); */
  string out = "DATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) + ", "
      + to_string(RSCED.get_obj_val()) + ", " + to_string(-numeric_limits<double>::infinity()) + ", "
      + to_string(solve_time) + ", GlobalOptimal, " + to_string(total_time);
  DebugOn(out << endl);
  return 0;
}