 //
//  ACUC.cpp
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
using namespace gravity;

int main (int argc, const char * argv[])
{
    const char* fname;
    double solver_time_end, total_time_end, solve_time, total_time;
    double total_time_start = get_cpu_time();
    int T = 2;
    if (argc >= 2) {
        fname = argv[1];
        T = atoi(argv[2]);
    }
    else {
        // fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        // fname = "/Users/hh/Dropbox/Work/Dev/nesta-0.7.0/opf/nesta_case2383wp_mp.m";
        // fname = "../../data_sets/Power/nesta_case3_lmbd.m";
        // fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
        //fname = "../../data_sets/Power/nesta_case1354_pegase.m";
        // fname = "../../data_sets/Power/nesta_case5_pjm.m";
        //fname = "../../data_sets/Power/nesta_case14_ieee.m";
        //fname = "../../data_sets/Power/nesta_case300_ieee.m";
        fname = "../../data_sets/Power/nesta_case118_ieee.m";
        //string fname = "../../data_sets/Power/anu.m";
    }
    // ACUC
    PowerNet* grid = new PowerNet();
    grid->readgrid(fname);

    // Grid Parameters
    const auto bus_pairs = grid->get_bus_pairs();
    auto nb_bus_pairs = grid->get_nb_active_bus_pairs();
    auto nb_gen = grid->get_nb_active_gens();
    auto nb_lines = grid->get_nb_active_arcs();
    auto nb_buses = grid->get_nb_active_nodes();

    // Schedule
    param<Real> rate_ramp("rate_ramp");
    param<Real> rate_switch("rate_switch");
    param<Real> min_up("min_up");
    param<Real> min_down("min_down");
    param<Real> cost_up("cost_up");
    param<Real> cost_down("cost_down");
    for (auto g: grid->gens) {
        rate_ramp(g->_name) = max(grid->pg_min(g->_name).getvalue(), 1*grid->pg_max(g->_name).getvalue());
        rate_switch(g->_name) = max(grid->pg_min(g->_name).getvalue(), 1*grid->pg_max(g->_name).getvalue());
    }
    min_up = 2;
    min_down = 2;
    cost_up = 50;
    cost_down = 30;

    grid->time_expand(T);
    rate_ramp.time_expand(T);
    rate_switch.time_expand(T);

    // set initial state or read the initial state
    param<Real> Pg_initial("Pg_initial");
    param<bool> On_off_initial("On_off_initial");
    for (auto g: grid->gens){
        Pg_initial(g->_name) = 0;
        On_off_initial(g->_name) = 0;
    }
    
    /** build model */
    Model ACUC("ACUC Model");

    /** Variables */
    // POWER GENERATION
    var<Real> Pg("Pg", grid->pg_min.in(grid->gens, T), grid->pg_max.in(grid->gens, T)); //This changes the lb and rb indices.
    //var<Real> Pg2("Pg2", non_neg_); // new var introduced for the perspective formulation.
    var<Real> Pg2("Pg2", 0, 1000); // new var introduced for the perspective formulation.
    var<Real> Qg ("Qg", grid->qg_min.in(grid->gens, T), grid->qg_max.in(grid->gens, T));
    ACUC.add_var(Pg.in(grid->gens, T));
    ACUC.add_var(Qg.in(grid->gens, T));
    ACUC.add_var(Pg2.in(grid->gens, T));

    //power flow
    var<Real> Pf_from("Pf_from", grid->S_max.in(grid->arcs, T));
    var<Real> Qf_from("Qf_from", grid->S_max.in(grid->arcs, T));
    var<Real> Pf_to("Pf_to", grid->S_max.in(grid->arcs, T));
    var<Real> Qf_to("Qf_to", grid->S_max.in(grid->arcs, T));
    ACUC.add_var(Pf_from.in(grid->arcs, T));
    ACUC.add_var(Qf_from.in(grid->arcs, T));
    ACUC.add_var(Pf_to.in(grid->arcs,  T));
    ACUC.add_var(Qf_to.in(grid->arcs,  T));

    //Lifted variables.
    var<Real>  R_Wij("R_Wij", grid->wr_min.in(bus_pairs, T), grid->wr_max.in(bus_pairs, T)); // real part of Wij
    var<Real>  Im_Wij("Im_Wij", grid->wi_min.in(bus_pairs, T), grid->wi_max.in(bus_pairs, T)); // imaginary part of Wij.
    var<Real>  Wii("Wii", grid->w_min.in(grid->nodes, T), grid->w_max.in(grid->nodes, T));
    ACUC.add_var(Wii.in(grid->nodes, T));
    ACUC.add_var(R_Wij.in(bus_pairs, T));
    ACUC.add_var(Im_Wij.in(bus_pairs, T));
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    // Commitment variables
    var<Real>  On_off("On_off", 0, 1);
    var<Real>  Start_up("Start_up", 0, 1);
    var<Real>  Shut_down("Shut_down", 0, 1);
    ACUC.add_var(On_off.in(grid->gens, -1, T));// -1
    ACUC.add_var(Start_up.in(grid->gens, T));
    ACUC.add_var(Shut_down.in(grid->gens, T));

    /* Construct the objective function*/
    func_ obj;
    for (auto g:grid->gens) {
        if (g->_active) {
            for (int t = 0; t < T; t++) {
                if (g->_active) {
                    string name = g->_name + ","+ to_string(t);
                    obj += grid->c1(name)*Pg(name)+ grid->c2(name)*Pg2(name) + grid->c0(name)*On_off(name);
                    //obj += grid->c1(name)*Pg(name)+ grid->c2(name)*Pg(name)*Pg(name)  + grid->c0(name)*On_off(name);
                    obj += cost_up*Start_up(name)+ cost_down*Shut_down(name);
                }
            }
        }
    }
    ACUC.set_objective(min(obj));

    /** Define constraints */

    /* SOCP constraints */
    Constraint SOC("SOC");
    SOC =  power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to() ;
    ACUC.add_constraint(SOC.in(bus_pairs, T) <= 0);
    
    /* perspective formulation of Pg^2 */
    Constraint Perspective_OnOff_Pg2("Perspective_OnOff_Pg2");
    Perspective_OnOff_Pg2 += power(Pg, 2) - Pg2*On_off;
    ACUC.add_constraint(Perspective_OnOff_Pg2.in(grid->gens, T) <= 0);
//    for (auto g:grid->gens) {
//        if (g->_active) {
//            for (int t = 0; t < T; t++) {
//                if (g->_active) {
//                    string name = g->_name + ","+ to_string(t);
//                    Constraint Perspective_OnOff_Pg2("Perspective_OnOff_Pg2_" + name);
//                    Perspective_OnOff_Pg2 = Pg(name)*Pg(name) - 100*Pg(name);
//                    ACUC.add_constraint(Perspective_OnOff_Pg2 <= 0);
//                }
//            }
//        }
//    }
    
    //KCL
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + grid->pl -sum(Pg.in_gens()) + grid->gs*Wii;
    ACUC.add_constraint(KCL_P.in(grid->nodes, T) == 0);

    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + grid->ql - sum(Qg.in_gens()) - grid->bs*Wii;
    ACUC.add_constraint(KCL_Q.in(grid->nodes, T) == 0);

    //AC Power Flow.
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid->g_ff*Wii.from() + grid->g_ft*R_Wij.in_pairs() + grid->b_ft*Im_Wij.in_pairs());
    ACUC.add_constraint(Flow_P_From.in(grid->arcs, T) == 0);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid->g_tt*Wii.to() + grid->g_tf*R_Wij.in_pairs() - grid->b_tf*Im_Wij.in_pairs());
    ACUC.add_constraint(Flow_P_To.in(grid->arcs, T) == 0);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid->g_ft*Im_Wij.in_pairs() - grid->b_ff*Wii.from() - grid->b_ft*R_Wij.in_pairs());
    ACUC.add_constraint(Flow_Q_From.in(grid->arcs, T) == 0);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (grid->b_tt*Wii.to() + grid->b_tf*R_Wij.in_pairs() + grid->g_tf*Im_Wij.in_pairs());
    ACUC.add_constraint(Flow_Q_To.in(grid->arcs, T) == 0);
    //Phase Angle Bounds constraints
    Constraint PAD_UB("PAD_UB");
    PAD_UB = Im_Wij - grid->tan_th_max*R_Wij;
    ACUC.add_constraint(PAD_UB.in(bus_pairs, T) <= 0);

    Constraint PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij - grid->tan_th_min*R_Wij;
    ACUC.add_constraint(PAD_LB.in(bus_pairs, T) >= 0);
    //* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from,  2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(grid->S_max, 2);
    ACUC.add_constraint(Thermal_Limit_from.in(grid->arcs, T) <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(grid->S_max, 2);
    ACUC.add_constraint(Thermal_Limit_to.in(grid->arcs, T) <= 0);

    // COMMITMENT CONSTRAINTS
    // Inter-temporal constraints 3a, 3d
    for (int t = 0; t < T; t++) {
            Constraint MC1("Inter_temporal_MC1");
            Constraint MC2("Inter_temporal_MC2");
            MC1 = On_off.in_at(grid->gens, t)- On_off.in_at(grid->gens, t-1)-  Start_up.in_at(grid->gens, t);
            MC2 = On_off.in_at(grid->gens, t-1) - On_off.in_at(grid->gens, t) - Shut_down.in_at(grid->gens, t);
            ACUC.add_constraint(MC1 <= 0);
            ACUC.add_constraint(MC2 <= 0);
    }
   // OnOffStartupShutdown_
    for (int t = 0; t < T; t++) {
        Constraint OnOffStartupShutdown("OnOffStartupShutdown_"+ to_string(t));
        OnOffStartupShutdown = On_off.in_at(grid->gens, t) - On_off.in_at(grid->gens, t-1)
                               -Start_up.in_at(grid->gens, t) + Shut_down.in_at(grid->gens, t);
        ACUC.add_constraint(OnOffStartupShutdown == 0);
    }
    // 7b
    for (int t = min_up.getvalue()-1; t < T; t++) {
        Constraint Min_Up("Min_Up_constraint" + to_string(t));
        for (int l = t-min_up.getvalue()+1; l < t+1; l++) {
            Min_Up   += Start_up.in_at(grid->gens, l);
        }
        Min_Up -= On_off.in_at(grid->gens, t);
        ACUC.add_constraint(Min_Up <= 0);
    }
    // 7c
    for (int t = min_down.getvalue()-1; t < T; t++) {
        Constraint Min_Down("Min_Down_constraint" + to_string(t));
        for (int l = t-min_down.getvalue()+1; l < t+1; l++) {
            Min_Down   += Shut_down.in_at(grid->gens, l);
        }
        Min_Down -= 1 - On_off.in_at(grid->gens, t);
        ACUC.add_constraint(Min_Down <= 0);
    }
    //production constraints
    Constraint Production_P_LB("Production_P_LB");
    Constraint Production_P_UB("Production_P_UB");
    Constraint Production_Q_LB("Production_Q_LB");
    Constraint Production_Q_UB("Production_Q_UB");
    // 5A
    Production_P_UB = Pg- grid->pg_max*On_off;
    Production_P_LB = Pg- grid->pg_min*On_off;
    ACUC.add_constraint(Production_P_UB.in(grid->gens, T) <= 0);
    ACUC.add_constraint(Production_P_LB.in(grid->gens, T) >= 0);
    // 5B
    Production_Q_UB = Qg - grid->qg_max*On_off;
    Production_Q_LB = Qg - grid->qg_min*On_off;
    ACUC.add_constraint(Production_Q_UB.in(grid->gens, T) <= 0);
    ACUC.add_constraint(Production_Q_LB.in(grid->gens, T) >= 0);
    // 5C
    for (int t = 1; t < T; t++) {
        Constraint Ramp_up("Ramp_up_constraint" + to_string(t));
        Constraint Ramp_down("Ramp_down_constraint" + to_string(t));
        Ramp_up =  Pg.in_at(grid->gens, t);
        Ramp_up -= Pg.in_at(grid->gens, t-1);
        Ramp_up -= rate_ramp.in_at(grid->gens, t)*On_off.in_at(grid->gens, t-1);
        Ramp_up -= rate_switch.in_at(grid->gens, t)*(1 - On_off.in_at(grid->gens, t-1));

        Ramp_down =  Pg.in_at(grid->gens, t-1);
        Ramp_down -= Pg.in_at(grid->gens, t);
        Ramp_down -= rate_ramp.in_at(grid->gens, t)*On_off.in_at(grid->gens, t);
        Ramp_down -= rate_switch.in_at(grid->gens, t)*(1 - On_off.in_at(grid->gens, t));

        ACUC.add_constraint(Ramp_up <= 0);
        ACUC.add_constraint(Ramp_down <= 0);
    }
    Constraint Ramp_up("Ramp_up_constraint0");
    Ramp_up =  Pg.in_at(grid->gens, 0) - Pg_initial.in(grid->gens) - rate_ramp.in_at(grid->gens, 0)*On_off_initial.in(grid->gens);
    Ramp_up -= rate_switch.in_at(grid->gens, 0)*(1 - On_off_initial.in(grid->gens));
    ACUC.add_constraint(Ramp_up <= 0);
    
    Constraint Ramp_down("Ramp_down_constraint0");
    Ramp_down = -1*Pg.in_at(grid->gens, 0) + Pg_initial.in(grid->gens);
    Ramp_down -= rate_ramp.in_at(grid->gens, 0)*On_off.in_at(grid->gens, 0);
    Ramp_down -= rate_switch.in_at(grid->gens, 0)*(1 - On_off.in_at(grid->gens, 0));
    ACUC.add_constraint(Ramp_down <= 0);

    //set the initial state of generators.
    Constraint gen_initial("gen_initialisation");
    gen_initial += On_off.in_at(grid->gens, -1)- On_off_initial.in(grid->gens);
    ACUC.add_constraint(gen_initial  == 0);

    /* Resolve it! */
    solver OPF(ACUC, cplex);
    //solver OPF(ACUC, ipopt);
    bool relax = true;
    int output = 1;
    double tol = 1e-6;
    OPF.run(output, relax, tol);
    //OPF.run();
    total_time_end = get_cpu_time();
    total_time = total_time_end - total_time_start;
    string out = "DATA_OPF, " + grid->_name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines)
                 +", " + to_string(ACUC._obj_val) + ", " + to_string(-numeric_limits<double>::infinity()) +", CPU time, " + to_string(total_time);

    cout << std::setprecision(4) <<out << endl;

    return 0;
}