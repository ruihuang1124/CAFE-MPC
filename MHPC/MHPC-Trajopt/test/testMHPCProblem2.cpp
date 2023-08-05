#include "MHPCProblem.h"
#include "QuadReference.h"
#include "MultiPhaseDDP.h"
#include "visualize_quadTraj_lcmt.hpp"
#include "visualize_quadState_lcmt.hpp"
#include <lcm/lcm-cpp.hpp>

#include <unistd.h>         // sleep

int main()
{
    // std::string reference_file_path = "../Reference/Data/quad_reference.csv";    
    MHPCConfig config;
    std::string mhpc_config_file("../MHPC/settings/mhpc_config.info");
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
    std::string fname_ddp_setting("../MHPC/settings/ddp_setting.info");
    loadHSDDPSetting(fname_ddp_setting, ddp_setting);
    ddp_setting.print();

    /* Solve the multi-phase problem */
    solver.set_initial_condition(xinit);
    solver.set_multiPhaseProblem(multiple_phases);
    solver.solve(ddp_setting);

    lcm::LCM visualize_traj_lcm;
    if (!visualize_traj_lcm.good())
    {
        printf("Failed to initialize LCM \n");
    }
    visualize_quadState_lcmt state_visual_data;
    problem.pretty_print();
    /* Run in MPC Loop */
    int loop = 0;
    ddp_setting.max_AL_iter = 2;
    ddp_setting.max_DDP_iter = 1;
    while (loop <= 300)
    {
        loop++;
        std::copy(xinit.data(), xinit.data() + 3, state_visual_data.pos);
        std::copy(xinit.data() + 3, xinit.data() + 6, state_visual_data.eul);
        std::copy(xinit.data() + 6, xinit.data() + 18, state_visual_data.qJ);
        std::copy(xinit.data() + 18, xinit.data() + 21, state_visual_data.vWorld);
        std::copy(xinit.data() + 21, xinit.data() + 24, state_visual_data.eulrate);
        std::copy(xinit.data() + 24, xinit.data() + 36, state_visual_data.qJd);
        visualize_traj_lcm.publish("visualize_mc_state", &state_visual_data);
        usleep(1e6 * config.dt_mpc * 5.0);

        problem.update();
        problem.pretty_print();
        multiple_phases.clear();
        for (auto phase : pdata.wb_phases)
        {
            multiple_phases.push_back(phase);
        }
        if (pdata.srb_phase.get() != nullptr)
        {
            multiple_phases.push_back(pdata.srb_phase);
        }
        // if(pdata.wb_trajs.front()->size()>=3)
        // {
        //     xinit = pdata.wb_trajs.front()->Xbar[2];
        // }else
        // {
        //     xinit = pdata.wb_trajs[1]->Xbar[0];
        // }
        
        xinit = pdata.wb_trajs.front()->Xbar[0];                
        solver.set_initial_condition(xinit);
        solver.set_multiPhaseProblem(multiple_phases);
        solver.solve(ddp_setting);
    }    
}