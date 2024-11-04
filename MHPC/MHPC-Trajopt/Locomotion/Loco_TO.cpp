#include <lcm/lcm-cpp.hpp>

#include "LocoProblem.h"
#include "MultiPhaseDDP.h"
#include "MHPC_Command_lcmt.hpp"
#include "wbTraj_lcmt.hpp"

using WBSingleTrajectory_d = Trajectory<double, WBM::xs, WBM::us, WBM::ys>;

void publish_trajectory(const deque<shared_ptr<WBSingleTrajectory_d>>& trajs,
                        const deque<Vec4<int>>& contacts);

void publish_command(const deque<shared_ptr<WBSingleTrajectory_d>>& trajs,
                     const deque<Vec4<int>>& contacts);

int main()
{
    // Problem data
    LocoProblem<double> opt_problem;
    MHPCProblemData<double> pdata;
    MHPCConfig pconfig;
    HSDDP_OPTION ddp_setting;   

   // Load MHPC config
    std::string mhpc_config_file("../MHPC/MHPC-Trajopt/Locomotion/settings/loco_config.info");
    loadMHPCConfig(mhpc_config_file, pconfig);
    pconfig.print();

    // Load DDP setting
    std::string ddp_setting_file("../MHPC/MHPC-Trajopt/Locomotion/settings/loco_ddp_setting.info");
    loadHSDDPSetting(ddp_setting_file, ddp_setting);

    // Load reference trajectory
    printf("Loading quadruped reference ... \n");
    std::string quad_reference_file("../Reference/Data/");
    quad_reference_file.append(pconfig.referenceFileName);
    quad_reference_file.append("/quad_reference.csv");
    pdata.quad_reference = std::make_shared<QuadReference>();
    pdata.quad_reference->load_top_level_data(quad_reference_file, false);

    // Initialize the multi-phase OCP (phase determination, memory allocation)
    opt_problem.set_problem_data(&pdata, &pconfig);
    opt_problem.initialization();

    // Initial condition
    VecM<double, WBM::nu> qJ, qJd;    
    VecM<double, 3> pos, eul, vWorld, eulrate;
    LocoProblem<double>::WBState xinit_wb;
    qJd.setZero();
    eul.setZero(); vWorld.setZero(); eulrate.setZero();
    // qJ = Vec3<double>(0, -0.9, 1.8).replicate<4,1>();
    // pos << 0, 0, 0.2511;
    qJ = Vec3<double>(0, -1.0, 2.0).replicate<4,1>();
    pos << 0, 0, 0.2183;
    xinit_wb << pos, eul, qJ, vWorld, eulrate, qJd;
    
    // Build solve
    MultiPhaseDDP<double> solver;
    deque<shared_ptr<SinglePhaseBase<double>>> multiple_phases;
    for (const auto &phase : pdata.wb_phases)
    {
        multiple_phases.push_back(phase);
    }
    if (pdata.srb_phase.get() != nullptr)
    {
        multiple_phases.push_back(pdata.srb_phase);
    }
    // Print multiphase information
    printf("=============MultiplePhase problem============\n");
    for (auto phase : multiple_phases)
    {
        phase->print();
    }

    solver.set_initial_condition(xinit_wb);
    solver.set_multiPhaseProblem(multiple_phases);
    solver.solve(ddp_setting);

    publish_trajectory(pdata.wb_trajs, pdata.wb_phase_contacts);
    publish_command(pdata.wb_trajs, pdata.wb_phase_contacts);

}

void publish_command(const deque<shared_ptr<WBSingleTrajectory_d>>& trajs,
                     const deque<Vec4<int>>& contacts)
{
    /* Publish wb trajectory for viz via lcm */
    lcm::LCM cmd_lcm;    
    if (!cmd_lcm.good())
    {
        printf("Failed to initialize LCM \n");
    }

    MHPC_Command_lcmt mpc_cmd;

    // Memory allocation
    std::vector<float> torque_k_float(12);
    std::vector<float> eul_k_float(3);
    std::vector<float> pos_k_float(3);
    std::vector<float> qJ_k_float(12);
    std::vector<float> vWorld_k_float(3);
    std::vector<float> eulrate_k_float(3);
    std::vector<float> qJd_k_float(12);
    std::vector<float> GRF_k_float(12);
    std::vector<float> Qu_k_float(12);
    std::vector<float> Quu_k_float(144);
    std::vector<float> Qux_k_float(432);
    std::vector<float> feedback_k_float(432);
    std::vector<int> phase_contact(4);
    std::vector<float> statusDuration_k(4);   

    mpc_cmd.N_mpcsteps = 0;
    float t = 0;
    for (int i(0); i < trajs.size(); i++)
    {        
        const int h = trajs[i]->size()-1;
        const auto& tau = trajs[i];
        for (int k = 0; k < h; k++)
        {                     

            const VecM<float,36>&x_k = tau->Xbar[k].template cast<float>();
            const Vec12<float>&u_k = tau->Ubar[k].template cast<float>();
            const Vec12<float>&GRF_k = tau->Y[k].template cast<float>();
            const Vec12<float>&Qu_k = tau->Qu[k].template cast<float>();
            const MatMN<float, 12, 12>&Quu_k = tau->Quu[k].template cast<float>();
            const MatMN<float, 12, 36>&Qux_k = tau->Qux[k].template cast<float>();
            const MatMN<float, 12, 36>&K_k = tau->K[k].template cast<float>();

            std::copy(u_k.data(), u_k.data() + u_k.size(), torque_k_float.data());
            std::copy(x_k.data(), x_k.data() + 3, pos_k_float.data());
            std::copy(x_k.data() + 3, x_k.data() + 6, eul_k_float.data());
            std::copy(x_k.data() + 6, x_k.data() + 18, qJ_k_float.data());
            std::copy(x_k.data() + 18, x_k.data() + 21, vWorld_k_float.data());
            std::copy(x_k.data() + 21, x_k.data() + 24, eulrate_k_float.data());
            std::copy(x_k.data() + 24, x_k.data() + x_k.size(), qJd_k_float.data());
            std::copy(GRF_k.data(), GRF_k.data() + GRF_k.size(), GRF_k_float.data());
            std::copy(Qu_k.data(), Qu_k.data() + Qu_k.size(), Qu_k_float.data());
    
            std::copy(Quu_k.data(), Quu_k.data()+144, Quu_k_float.data());
            std::copy(Qux_k.data(), Qux_k.data()+432, Qux_k_float.data());
            std::copy(K_k.data(), K_k.data()+432, feedback_k_float.data());
            std::copy(contacts[i].data(), contacts[i].data()+4, phase_contact.data());            
    
            mpc_cmd.pos.push_back(pos_k_float);
            mpc_cmd.eul.push_back(eul_k_float);
            mpc_cmd.qJ.push_back(qJ_k_float);
            mpc_cmd.vWorld.push_back(vWorld_k_float);
            mpc_cmd.eulrate.push_back(eulrate_k_float);
            mpc_cmd.qJd.push_back(qJd_k_float);
            mpc_cmd.torque.push_back(torque_k_float);
            mpc_cmd.GRF.push_back(GRF_k_float);
            mpc_cmd.Qu.push_back(Qu_k_float);
            mpc_cmd.Quu.push_back(Quu_k_float);
            mpc_cmd.Qux.push_back(Qux_k_float);
            mpc_cmd.feedback.push_back(feedback_k_float);
            mpc_cmd.contacts.push_back(phase_contact);
            mpc_cmd.statusTimes.push_back(statusDuration_k);
            mpc_cmd.mpc_times.push_back(t);
            t+=tau->timeStep;            
            mpc_cmd.N_mpcsteps++;
        }        
    }    
    cmd_lcm.publish("MHPC_COMMAND", &mpc_cmd);
    sleep(1);
    printf(GRN);
    printf("published a mpc command message \n");
    printf(RESET);
}                     

void publish_trajectory(const deque<shared_ptr<WBSingleTrajectory_d>>& trajs,
                        const deque<Vec4<int>>& contacts)
{
    /* Publish wb trajectory for viz via lcm */
    lcm::LCM viz_lcm(getLcmUrl(255));    
    if (!viz_lcm.good())
    {
        printf("Failed to initialize LCM \n");
    }

    wbTraj_lcmt wbtraj_lcmt;
    std::vector<double> pos_vec(3), eul_vec(3), vWorld_vec(3), eulrate_vec(3);
    std::vector<double> qJ_vec(12), qJd_vec(12), pFoot_vec(12), torque_vec(12);
    std::vector<double> hg_vec(3), dhg_vec(3);
    std::vector<int> phase_contact(4);
    wbtraj_lcmt.sz = 0;
    wbtraj_lcmt.wb_sz = 0;
    double t = 0;
    for (int i(0); i < trajs.size(); i++)
    {        
        const int h = trajs[i]->size() - 1;
        const auto& tau = trajs[i];
        std::copy(contacts[i].data(), contacts[i].data()+4, phase_contact.data());
        for (int k = 0; k < h; k++)
        {
            wbtraj_lcmt.sz ++;
            wbtraj_lcmt.wb_sz ++;

            const auto& xk = tau->Xsim[k];
            const auto& uk = tau->U[k];
            const auto& defect_k = tau->Defect[k].lpNorm<Eigen::Infinity>();
    
            

            std::copy(xk.data(), xk.data() + 3, pos_vec.data());
            std::copy(xk.data() + 3, xk.data() + 6, eul_vec.data());
            std::copy(xk.data() + 6, xk.data() + 18, qJ_vec.data());
            std::copy(xk.data() + 18, xk.data() + 21, vWorld_vec.data());
            std::copy(xk.data() + 21, xk.data() + 24, eulrate_vec.data());
            std::copy(xk.data() + 24, xk.data() + 36, qJd_vec.data());
            std::copy(uk.data(), uk.data() + 12, torque_vec.data());

            wbtraj_lcmt.pos.push_back(pos_vec);
            wbtraj_lcmt.vWorld.push_back(vWorld_vec);
            wbtraj_lcmt.eul.push_back(eul_vec);
            wbtraj_lcmt.eulrate.push_back(eulrate_vec);
            wbtraj_lcmt.qJ.push_back(qJ_vec);
            wbtraj_lcmt.qJd.push_back(qJd_vec); 
            wbtraj_lcmt.torque.push_back(torque_vec);
            wbtraj_lcmt.contact.push_back(phase_contact);
            wbtraj_lcmt.defect.push_back(defect_k);
            wbtraj_lcmt.hg.push_back(hg_vec);
            wbtraj_lcmt.dhg.push_back(dhg_vec);
            wbtraj_lcmt.time.push_back(t);
            if (k<h)
            {
                t+=tau->timeStep;            
            }                        
        }        
    }

    viz_lcm.publish("visualize_wb_traj", &wbtraj_lcmt);
}