#include "MHPCProblem.h"
#include "QuadReference.h"
#include "MultiPhaseDDP.h"
#include "visualize_quadTraj_lcmt.hpp"
#include <lcm/lcm-cpp.hpp>

int main()
{
    // std::string reference_file_path = "../Reference/Data/quad_reference.csv";    
    MHPCConfig config;
    std::string mhpc_config_file("../MHPC/mhpc_config.info");
    loadMHPCConfig(mhpc_config_file, config);
    config.print();

    
    std::shared_ptr<QuadReference> quad_ref;
    quad_ref = std::make_shared<QuadReference>();
    std::string quad_reference_file("../Reference/Data/");
    quad_reference_file.append(config.referenceFileName);
    quad_reference_file.append("/quad_reference.csv");
    quad_ref->load_top_level_data(quad_reference_file, true);

    MHPCProblem<double> problem;
    MHPCProblemData<double> pdata;
    pdata.quad_reference = quad_ref;
    problem.set_problem_data(&pdata, &config);
    problem.prepare_initialization();
    problem.pretty_print();
    problem.initialize_parameters();
    problem.initialize_multiPhaseProblem();

    /* Setup solver */
    MultiPhaseDDP<double> solver;
    std::deque<shared_ptr<SinglePhaseBase<double>>> multiple_phases;
    for (auto phase : pdata.wb_phases)
    {
        multiple_phases.push_back(phase);
    }
    if (pdata.srb_phase.get() != nullptr)
    {
        multiple_phases.push_back(pdata.srb_phase);
    }

    /* Print phase information */
    printf("=============MultiplePhase problem============\n");
    for (auto phase : multiple_phases)
    {
        phase->print();
    }

    /* Initial condition */
    VecM<double, WBM::xs> xinit;
    VecM<double, 3> pos, eul, vel, eulrate;
    VecM<double, 12> qJ, qJd;
    pos.setZero();
    eul.setZero();
    vel.setZero();
    eulrate.setZero();
    qJ = Vec3<double>(0, 1.0, -2.2).replicate<4, 1>();
    qJd.setZero();
    pos[2] = 0.22;
    xinit << pos, eul, qJ, vel, eulrate, qJd;

    /* Solve th multiple-phase TO problem */
    HSDDP_OPTION ddp_setting;
    std::string fname_ddp_setting("../MHPC/ddp_setting.info");
    loadHSDDPSetting(fname_ddp_setting, ddp_setting);
    ddp_setting.print();

    /* Initial control guess */
    for (auto &tau_i : pdata.wb_trajs)
    {
        for (size_t k = 0; k < tau_i->size(); k++)
        {
            tau_i->Ubar[k].setConstant(.0);                
        }
    }
    if (pdata.srb_phase.get() != nullptr)
    {
        for (size_t k = 0; k < pdata.srb_traj->size() - 1; k++)
        {
            pdata.srb_traj->Ubar[k].setZero();                    
        }
    }

    /* Solve the multi-phase problem */
    solver.set_initial_condition(xinit);
    solver.set_multiPhaseProblem(multiple_phases);
    solver.solve(ddp_setting);

    // int loop = 0;
    // while (loop <= 5)
    // {
    //     loop++;
    //     problem.update();
    //     multiple_phases.clear();
    //     for (auto phase : pdata.wb_phases)
    //     {
    //         multiple_phases.push_back(phase);
    //     }
    //     if (pdata.srb_phase.get() != nullptr)
    //     {
    //         multiple_phases.push_back(pdata.srb_phase);
    //     }
    //     solver.set_initial_condition(xinit);
    //     solver.set_multiPhaseProblem(multiple_phases);
    //     solver.solve(ddp_setting);
    // }

    // Debug
    const std::string logMHPC_folderName = "../MHPC/log/";
    log_trajectory_sequence(logMHPC_folderName, pdata.wb_trajs);
    
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

    for (const auto tau : pdata.wb_trajs)
    {
        tau->measure_dynamics_feasibility();
        for (size_t k = 0; k < tau->size(); k++)
        {
            VecM<float, WBM::xs> xk = tau->X[k].cast<float>();
            VecM<float, WBM::us> uk = tau->Ubar[k].cast<float>();
            VecM<float, WBM::ys> yk = tau->Y[k].cast<float>();

            std::copy(xk.data(), xk.data() + 3, pos_vec.data());
            std::copy(xk.data() + 3, xk.data() + 6, eul_vec.data());
            std::copy(xk.data() + 6, xk.data() + 18, qJ_vec.data());
            std::copy(xk.data() + 18, xk.data() + 21, vWorld_vec.data());
            std::copy(xk.data() + 21, xk.data() + 24, eulrate_vec.data());
            std::copy(xk.data() + 24, xk.data() + 36, qJd_vec.data());
            std::copy(uk.data(), uk.data() + 12, torque_vec.data());
            std::copy(yk.data(), yk.data() + 12, GRF_vec.data());

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
    }
  
    const auto tau = pdata.srb_traj;                
    if (pdata.srb_traj.get() != nullptr)
    {
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
    }

    
    visualization_lcm_data.len = visualization_lcm_data.pos.size();
    printf("publishing lcm data with size %u ...\n", visualization_lcm_data.len);
    visualize_traj_lcm.publish("visualize_mc_motion", &visualization_lcm_data);
}