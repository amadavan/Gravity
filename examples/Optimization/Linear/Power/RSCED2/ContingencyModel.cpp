//
// Created by Avinash Madavan on 6/11/20.
//

#include "Models.h"

func<double> gravity::createContingencyModel(Model<> &model, PowerNet &grid,
                                             const RSCEDdata &data,
                                             const double risk_aversion,
                                             const std::pair<std::string, std::pair<Arc *, Gen *>> &contingency,
                                             int subproblem) {
    if (risk_aversion < 0. || risk_aversion >= 1.)
        throw std::runtime_error(
                "CVaR risk aversion parameter must be in [0, 1), got " + std::to_string(risk_aversion));
    double alpha_prime = 1. / (1. - risk_aversion);

    // Extract pairs to meaningful variable names
    std::string contingency_name = contingency.first;
    Arc *arc_failure = contingency.second.first;
    Gen *gen_failure = contingency.second.second;

    // Set failed components to inactive
    if (arc_failure) arc_failure->_active = false;
    if (gen_failure) gen_failure->_active = false;

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

    std::string subscript = "_c" + contingency_name;

    /** Nominal model variables as parameters */
    param<> Pg("Pg");
    Pg.in(gens);
    Pg = 0;
    param<> z("z");
    z = 0;

    /** Variable definition */
    /* Positive recourse generation */
    var<> dPg_p("dPg_p" + subscript);
    model.add(dPg_p.in(gens));
    dPg_p.add_lb_only(0);
    dPg_p.set_problem(subproblem);

    /* Negative recourse generation */
    var<> dPg_n("dPg_n" + subscript);
    model.add(dPg_n.in(gens));
    dPg_n.add_lb_only(0);
    dPg_n.set_problem(subproblem);

    /* Power flows */
    var<> Pf_c("Pf_c" + subscript);
    model.add(Pf_c.in(arcs));
    Pf_c.set_problem(subproblem);

    var<> Pf_recourse_c("Pf_recourse_c" + subscript);
    model.add(Pf_recourse_c.in(arcs));
    Pf_recourse_c.set_problem(subproblem);

    /* Phase angle variables */
    var<> theta_c("theta_c" + subscript);
    model.add(theta_c.in(nodes));
    theta_c.set_problem(subproblem);

    var<> theta_recourse_c("theta_recourse_c" + subscript);
    model.add(theta_recourse_c.in(nodes));
    theta_recourse_c.set_problem(subproblem);

    /* Load shed variables */
    var<> dD_c("dD_c" + subscript);
    model.add(dD_c.in(nodes));
    dD_c.set_lb(0);
    dD_c.set_problem(subproblem);

    /* CVaR reformulation variable */
    var<> y_c("y_c" + subscript);
    model.add(y_c);
    y_c.add_lb_only(0);
    y_c.set_problem(subproblem);

    /* Slack variables */
    var<> slack_kcl_p_c("slack_kcl_p_c" + subscript);
    model.add(slack_kcl_p_c.in(nodes));
    slack_kcl_p_c.add_lb_only(0);
    slack_kcl_p_c.set_problem(subproblem);

    var<> slack_kcl_n_c("slack_kcl_n_c" + subscript);
    model.add(slack_kcl_n_c.in(nodes));
    slack_kcl_n_c.add_lb_only(0);
    slack_kcl_n_c.set_problem(subproblem);

    var<> slack_pad_c("slack_pad_c" + subscript);
    model.add(slack_pad_c.in(node_pairs));
    slack_pad_c.add_lb_only(0);
    slack_pad_c.set_problem(subproblem);

    var<> slack_pad_recourse_c("slack_pad_recourse_c" + subscript);
    model.add(slack_pad_recourse_c.in(node_pairs));
    slack_pad_recourse_c.add_lb_only(0);
    slack_pad_recourse_c.set_problem(subproblem);

    var<> slack_thermal_c("slack_thermal_c" + subscript);
    model.add(slack_thermal_c.in(arcs));
    slack_thermal_c.add_lb_only(0);
    slack_thermal_c.set_problem(subproblem);

    var<> slack_thermal_recourse_c("slack_thermal_recourse_c" + subscript);
    model.add(slack_thermal_recourse_c.in(arcs));
    slack_thermal_recourse_c.add_lb_only(0);
    slack_thermal_recourse_c.set_problem(subproblem);

    /** Post-failure constraints */
    /* Power flow constraint */
    Constraint<> Flow_P_c("Flow_P_c" + subscript);
    Flow_P_c = Pf_c + b * (theta_c.from(arcs) - theta_c.to(arcs));
    model.add(Flow_P_c.in(arcs) == 0);

    /* Power balance constraint */
    Constraint<> KCL_P_c("KCL_P_c" + subscript);
    KCL_P_c = sum(Pf_c, out_arcs) - sum(Pf_c, in_arcs) + pl + gs - sum(Pg, gen_nodes) + slack_kcl_p_c - slack_kcl_n_c;
    model.add(KCL_P_c.in(nodes) == 0);

    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB_c("PAD_UB_c" + subscript);
    PAD_UB_c = theta_c.from(node_pairs) - theta_c.to(node_pairs);
    PAD_UB_c -= th_max.in(node_pairs) + slack_pad_c;
    model.add_lazy(PAD_UB_c.in(node_pairs) <= 0);
    Constraint<> PAD_LB_c("PAD_LB_c" + subscript);
    PAD_LB_c = theta_c.from(node_pairs) - theta_c.to(node_pairs);
    PAD_LB_c -= th_min.in(node_pairs) - slack_pad_c;
    model.add_lazy(PAD_LB_c.in(node_pairs) >= 0);

    /* Line Limits constraints */
    Constraint<> Thermal_UB_c("Thermal_UB_c" + subscript);
    Thermal_UB_c = Pf_c - S_max * data.DAL_multiplier - slack_thermal_c;
    model.add_lazy(Thermal_UB_c.in(arcs) <= 0);
    Constraint<> Thermal_LB_c("Thermal_LB_c" + subscript);
    Thermal_LB_c = Pf_c + S_max * data.DAL_multiplier + slack_thermal_c;
    model.add_lazy(Thermal_LB_c.in(arcs) >= 0);

    /** Post-recourse constraints */
    /* CVaR lower bound */
    Constraint<> OBJ_BOUND_c("OBJ_BOUND_c" + subscript);
    OBJ_BOUND_c = y_c + z - product(c1, Pg) - product(c1, dPg_p) - product(c1, dPg_n) - product(data.VoLL, dD_c);
    model.add(OBJ_BOUND_c >= 0);

    /* Power flow constraint */
    Constraint<> Flow_P_recourse_c("Flow_P_recourse_c" + subscript);
    Flow_P_recourse_c = Pf_recourse_c + b * (theta_recourse_c.from(arcs) - theta_recourse_c.to(arcs));
    model.add(Flow_P_recourse_c.in(arcs) == 0);

    /* Power balance constraint */
    Constraint<> KCL_P_recourse_c("KCL_P_recourse_c" + subscript);
    KCL_P_recourse_c =
            sum(Pf_recourse_c, out_arcs) - sum(Pf_recourse_c, in_arcs) + pl + gs - sum(Pg, gen_nodes)
            - sum(dPg_p, gen_nodes) + sum(dPg_n, gen_nodes) - dD_c;
    model.add(KCL_P_recourse_c.in(nodes) == 0);

    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB_recourse_c("PAD_UB_recourse_c" + subscript);
    PAD_UB_recourse_c = theta_recourse_c.from(node_pairs) - theta_recourse_c.to(node_pairs);
    PAD_UB_recourse_c -= th_max.in(node_pairs) + slack_pad_recourse_c;
    model.add_lazy(PAD_UB_recourse_c.in(node_pairs) <= 0);
    Constraint<> PAD_LB_recourse_c("PAD_LB_recourse_c" + subscript);
    PAD_LB_recourse_c = theta_recourse_c.from(node_pairs) - theta_recourse_c.to(node_pairs);
    PAD_LB_recourse_c -= th_min.in(node_pairs) - slack_pad_recourse_c;
    model.add_lazy(PAD_LB_recourse_c.in(node_pairs) >= 0);

    /* Line Limits constraints */
    Constraint<> Thermal_UB_recourse_c("Thermal_UB_recourse_c" + subscript);
    Thermal_UB_recourse_c = Pf_recourse_c;
    Thermal_UB_recourse_c -= S_max * data.STE_multiplier + slack_thermal_recourse_c;
    model.add_lazy(Thermal_UB_recourse_c.in(arcs) <= 0);
    Constraint<> Thermal_LB_recourse_c("Thermal_LB_recourse_c" + subscript);
    Thermal_LB_recourse_c = Pf_recourse_c;
    Thermal_LB_recourse_c += S_max * data.STE_multiplier + slack_thermal_recourse_c;
    model.add_lazy(Thermal_LB_recourse_c.in(arcs) >= 0);

    /* Generation capacity limits */
    Constraint<> Gen_UB_c("Gen_recourse_UB_c" + subscript);
    Gen_UB_c = Pg + dPg_p - dPg_n;
    Gen_UB_c -= pg_max;
    model.add_lazy(Gen_UB_c.in(gens) <= 0);
    Constraint<> Gen_LB_c("Gen_recourse_LB_c" + subscript);
    Gen_LB_c = Pg + dPg_p - dPg_n;
    Gen_LB_c += pg_min;
    model.add_lazy(Gen_LB_c.in(gens) >= 0);

    /* Recourse ramp capacity limits */
    Constraint<> Recourse_UB_c("Recourse_UB_c" + subscript);
    Recourse_UB_c = dPg_p - dPg_n;
    Recourse_UB_c -= data.ramp_max;
    model.add_lazy(Recourse_UB_c.in(gens) <= 0);
    Constraint<> Recourse_LB_c("Recourse_LB_c" + subscript);
    Recourse_LB_c = dPg_p - dPg_n;
    Recourse_LB_c += data.ramp_max;
    model.add_lazy(Recourse_LB_c.in(gens) >= 0);

    /** Objective */
    func<double> obj = alpha_prime * data.failure_probability * y_c;
    obj += 9999 * sum(slack_thermal_c, arcs) + 9999 * sum(slack_thermal_recourse_c, arcs)
           + 9999 * sum(slack_pad_c, node_pairs) + 9999 * sum(slack_pad_recourse_c, node_pairs)
           + 9999 * sum(slack_kcl_p_c, nodes) + 9999 * sum(slack_kcl_n_c, nodes);

    // Reset failed components to active
    if (arc_failure) arc_failure->_active = true;
    if (gen_failure) gen_failure->_active = true;

    return obj;
}