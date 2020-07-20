//
// Created by Avinash Madavan on 6/11/20.
//

#include "Models.h"

Model<> gravity::createNominalModel(PowerNet &grid, const RSCEDdata &data, const double risk_aversion) {
  if (risk_aversion < 0. || risk_aversion >= 1.)
    throw std::runtime_error("CVaR risk aversion parameter must be in [0, 1), got " + std::to_string(risk_aversion));
  double alpha_prime = 1./(1. - risk_aversion);

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

  /** Declare model */
  Model<> NominalED("NominalED Model");

  /** Variables */
  /* Power generation variables */
  var<> Pg("Pg", pg_min, pg_max);
  NominalED.add(Pg.in(gens));

  /* Phase angle variables */
  var<> theta("𝛉");
  NominalED.add(theta.in(nodes));

  /* CVaR reformulation variable */
  // y0, z0(?)
  var<> z("z");
  var<> y("y");
  NominalED.add(z);
  NominalED.add(y);
  z.add_lb_only(0);
  y.add_lb_only(0);

  /** Constraints */
  /* CVaR lower bound */
  Constraint<> OBJ_BOUND("OBJ_BOUND");
  OBJ_BOUND = y + z - product(c1, Pg);
  NominalED.add(OBJ_BOUND >= 0);

  /* REF BUS */
  Constraint<> Ref_Bus("Ref_Bus");
  Ref_Bus = theta(grid.ref_bus);
  NominalED.add(Ref_Bus == 0);

  /* Power balance constraint */
  Constraint<> KCL_P("KCL_P");
  KCL_P = b.tr().in(in_arcs) * (theta.from(in_arcs) - theta.to(in_arcs))
      - b.tr().in(out_arcs) * (theta.from(out_arcs) - theta.to(out_arcs)) + pl + gs - sum(Pg, gen_nodes);
  NominalED.add(KCL_P.in(nodes) == 0);

  /* Phase Angle Bounds constraints */
  Constraint<> PAD_UB("PAD_UB");
  PAD_UB = theta.from(node_pairs) - theta.to(node_pairs);
  PAD_UB -= th_max;
  NominalED.add_lazy(PAD_UB.in(node_pairs) <= 0);
  Constraint<> PAD_LB("PAD_LB");
  PAD_LB = theta.from(node_pairs) - theta.to(node_pairs);
  PAD_LB -= th_min;
  NominalED.add_lazy(PAD_LB.in(node_pairs) >= 0);

  /* Line Limits constraints */
  Constraint<> Thermal_UB("Thermal_UB");
  Thermal_UB = b * (theta.to(arcs) - theta.from(arcs));
  Thermal_UB -= S_max;
  NominalED.add_lazy(Thermal_UB.in(arcs) <= 0);
  Constraint<> Thermal_LB("Thermal_LB");
  Thermal_LB = b * (theta.to(arcs) - theta.from(arcs));
  Thermal_LB += S_max;
  NominalED.add_lazy(Thermal_LB.in(arcs) >= 0);

  /** Objective */
  auto obj = z + alpha_prime * (1 - nb_lines * data.failure_probability) * y;

  return NominalED;
}