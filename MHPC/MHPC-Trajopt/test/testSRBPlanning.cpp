#include "MHPCProblem.h"
#include "QuadReference.h"
#include "MultiPhaseDDP.h"
#include "visualize_quadTraj_lcmt.hpp"
#include <lcm/lcm-cpp.hpp>

int main()
{
    std::string reference_file_path = "../Reference/Data/quad_reference.csv";
    std::shared_ptr<QuadReference> quad_ref;
    quad_ref = std::make_shared<QuadReference>();
    quad_ref->load_top_level_data(reference_file_path, true);

    MHPCProblem<double> problem;
    MHPCProblemData<double> pdata;
    pdata.quad_reference = quad_ref;
    MHPCConfig config;
    std::string mhpc_config_file("../MHPC/MHPC-Trajopt/test/SRBPlanning_config.info");
    loadMHPCConfig(mhpc_config_file, config);
    config.print();

    problem.set_problem_data(&pdata, &config);    
    problem.prepare_initialization();
    problem.pretty_print();
    problem.initialize_parameters();
    problem.initialize_multiPhaseProblem();

    /* Setup solver */
    // Assuming only srb planning phase
    MultiPhaseDDP<double> solver;
    std::deque<shared_ptr<SinglePhaseBase<double>>> multiple_phases;    
    multiple_phases.push_back(pdata.srb_phase);

    /* Print phase information */
    printf("=============MultiplePhase problem============\n");
    for (auto phase : multiple_phases)
    {
        phase->print();
    }

    /* Initial condition */
    VecM<double, SRBM::xs> xinit_srb;
    VecM<double, 3> pos, eul, vel, eulrate;
    pos.setZero();
    eul.setZero();
    vel.setZero();
    eulrate.setZero();
    pos[2] = 0.28;
    xinit_srb << pos, eul, vel, eulrate;

    /* Solve th multiple-phase TO problem */
    HSDDP_OPTION ddp_setting;
    std::string fname_ddp_setting("../MHPC/ddp_setting.info");
    loadHSDDPSetting(fname_ddp_setting, ddp_setting);
    ddp_setting.print();

    
    /* Solve the multi-phase problem */
    solver.set_initial_condition(xinit_srb);
    solver.set_multiPhaseProblem(multiple_phases);
    solver.solve(ddp_setting);

    // Debug
    const std::string logMHPC_folderName = "../MHPC/log/";    
    log_a_trajectory(logMHPC_folderName, pdata.srb_traj );

    /* Publish via LCM */
    lcm::LCM visualize_traj_lcm;
    if (!visualize_traj_lcm.good())
    {
        printf("Failed to initialize LCM \n");
    }
    visualize_quadTraj_lcmt visualization_lcm_data;    
    visualization_lcm_data.WB_plan_dur = config.plan_dur_wb;
    visualization_lcm_data.WB_dt = config.dt_wb;
    visualization_lcm_data.SRB_plan_dur = config.plan_dur_srb;
    visualization_lcm_data.SRB_dt = config.dt_srb;

    std::vector<float> pos_vec(3), eul_vec(3), vWorld_vec(3), eulrate_vec(3);
    std::vector<float> qJ_vec(12), qJd_vec(12), pFoot_vec(12, 0), torque_vec(12), GRF_vec(12);

    const auto tau = pdata.srb_traj;
    tau->measure_dynamics_feasibility();
    for (size_t k = 0; k < tau->size(); k++)
    {
        VecM<float, SRBM::xs> xk = tau->X[k].cast<float>();
        VecM<float, SRBM::us> uk = tau->Ubar[k].cast<float>();

        std::copy(xk.data(), xk.data() + 3, pos_vec.data());
        std::copy(xk.data() + 3, xk.data() + 6, eul_vec.data());        
        std::copy(xk.data() + 6, xk.data() + 9, vWorld_vec.data());
        std::copy(xk.data() + 9, xk.data() + 12, eulrate_vec.data());        
        std::copy(uk.data(), uk.data() + 12, GRF_vec.data());

        visualization_lcm_data.pos.push_back(pos_vec);
        visualization_lcm_data.eul.push_back(eul_vec);
        visualization_lcm_data.vWorld.push_back(vWorld_vec);
        visualization_lcm_data.eulrate.push_back(eulrate_vec);
        visualization_lcm_data.qJ.push_back(qJ_vec);
        visualization_lcm_data.pFoot.push_back(pFoot_vec);
        visualization_lcm_data.torque.push_back(torque_vec);
        visualization_lcm_data.grf.push_back(GRF_vec);
        visualization_lcm_data.feas.push_back(tau->Defect[k].norm());                        
    }

    visualization_lcm_data.len = visualization_lcm_data.pos.size();
    printf("publishing lcm data ...\n");
    visualize_traj_lcm.publish("visualize_mc_motion", &visualization_lcm_data);
}