/*!
 * @file opt_main.cpp
 * @brief Main Function for the standalone DDP trajectory optimizer
 *
 * The main function initilizes the DDP solver updates upon the request
 * of the main Imitation Controller
 */

#include "MHPCLocomotion.h"


#ifdef TIME_BENCHMARK
#include <chrono>
using namespace std::chrono;
using duration_ms = std::chrono::duration<float, std::chrono::milliseconds::period>;
#endif // TIME_BENCHMARK

template <typename T>
void MHPCLocomotion<T>::initialize()
{
    // Clear problem data (used for Monte Carlo sim)
    opt_problem_data.clear();   

    // Load MHPC config    
    std::string mhpc_config_file("../MHPC/mhpc_config.info");    
    loadMHPCConfig(mhpc_config_file, mpc_config);
    mpc_config.print();
    
    // Load DDP setting        
    std::string ddp_setting_file("../MHPC/ddp_setting.info");
    loadHSDDPSetting(ddp_setting_file, ddp_setting);        

    // Load reference trajectory
    printf("Loading quadruped reference ... \n");    
    std::string quad_reference_file("../Reference/Data/");
    quad_reference_file.append(mpc_config.referenceFileName);
    quad_reference_file.append("/quad_reference.csv");
    opt_problem_data.quad_reference =  std::make_shared<QuadReference>();
    opt_problem_data.quad_reference->load_top_level_data(quad_reference_file, true);     
    
    // Initialize the multi-phase OCP (phase determination, memory allocation)
    opt_problem.set_problem_data(&opt_problem_data, &mpc_config);
    opt_problem.initialization();

#ifdef DEBUG_MODE
    opt_problem.pretty_print();
#endif    

    mpc_time = 0;
    mpc_time_prev = 0;    
    
    // Assemble the multi-phase probelm and solve it with MSDDP
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (const auto& phase : opt_problem_data.wb_phases)
    {
        multiple_phases.push_back(phase);
    }
    if (opt_problem_data.srb_phase.get() != nullptr)
    {
        multiple_phases.push_back(opt_problem_data.srb_phase);
    }
    // Print multiphase information
    printf("=============MultiplePhase problem============\n");
    for (auto phase : multiple_phases)
    {
        phase->print();
    }

    // set the initial condition
    eul.setZero();
    pos << 0, 0, 0.2486;
    eulrate.setZero();
    vWorld.setZero();
    qJ = Vec3<T>(0, 0.8, -1.6).template replicate<4,1>();
    qJd.setZero();   
    x_init_wb << pos, eul, qJ, vWorld, eulrate, qJd;
    solver.set_initial_condition(x_init_wb);    

#ifdef TIME_BENCHMARK
    auto start = high_resolution_clock::now();
#endif
    solver.set_multiPhaseProblem(multiple_phases);
    solver.solve(ddp_setting);
#ifdef TIME_BENCHMARK
    auto stop = high_resolution_clock::now();
    auto duration = duration_ms(stop - start);
    solve_time = duration.count();
#endif
   
    mpc_iter = 0;

    printf("MHPC solver is initialized successfully \n\n");

    publish_mpc_cmd();
    
}

template <typename T>
void MHPCLocomotion<T>::update()
{
    mpc_mutex.lock(); // lock mpc to prevent updating while the previous hasn't finished

    // run-time DDP setting when re-solving DDP in MPC 
    ddp_setting.max_AL_iter = 2;
    ddp_setting.max_DDP_iter = 1;
    mpc_iter++;

    printf("************************************* \n");
    printf("************Resolving MPC************ \n");
    printf("********MPC Iteration = %d*********** \n", mpc_iter);
    printf("************************************* \n");

    /* update the problem */
    opt_problem.update();

#ifdef DEBUG_MODE    
    opt_problem.pretty_print();
#endif

    /* update current state*/
    Eigen::Map<VecM<float,3>> pos_map(mpc_data.pos);
    Eigen::Map<VecM<float,3>> eul_map(mpc_data.eul);
    Eigen::Map<VecM<float,12>> qJ_map(mpc_data.qJ);
    Eigen::Map<VecM<float,3>> vWorld_map(mpc_data.vWorld);
    Eigen::Map<VecM<float,3>> eulrate_map(mpc_data.eulrate);
    Eigen::Map<VecM<float,12>> qJd_map(mpc_data.qJd);     
    x_init_wb << pos_map.cast<T>(), eul_map.cast<T>(), qJ_map.cast<T>(),
                 vWorld_map.cast<T>(), eulrate_map.cast<T>(), qJd_map.cast<T>();

    /* build solver and solve the TO problem */
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (const auto& phase : opt_problem_data.wb_phases)
    {
        multiple_phases.push_back(phase);
    }
    if (opt_problem_data.srb_phase.get() != nullptr)
    {
        multiple_phases.push_back(opt_problem_data.srb_phase);
    }
    solver.set_multiPhaseProblem(multiple_phases);
    solver.set_initial_condition(x_init_wb);

#ifdef TIME_BENCHMARK
    auto start = high_resolution_clock::now();
#endif
    solver.solve(ddp_setting);
#ifdef TIME_BENCHMARK
    auto stop = high_resolution_clock::now();
    auto duration = duration_ms(stop - start);
    solve_time = duration.count();
    printf("solve time = %f \n", solve_time);
#endif

    publish_mpc_cmd();
#ifdef DEBUG_MODE
    // publish_debugfoot();
    // opt_problem.lcm_publish();
#endif
    mpc_mutex.unlock();
}

template <typename T>
void MHPCLocomotion<T>::mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                                          const MHPC_Data_lcmt *msg)
{
    (void) (rbuf);
    (void) (chan);
    mpc_mutex.lock();
    printf(GRN);
    printf("Received resolving request\n");
    printf(RESET);

    if (msg->reset_mpc)
    {
        ddp_setting.MS = msg->MS;
        mpc_mutex.unlock();
        initialize();
        return;
    }

    std::memcpy(&mpc_data, msg, sizeof(mpc_data));
    mpc_time_prev = mpc_time;
    mpc_time = mpc_data.mpctime;

    
    mpc_mutex.unlock();
    std::thread solve_mpc_thread(&MHPCLocomotion::update, this);
    solve_mpc_thread.detach(); // detach the thread from the main thread. The thread would exit once it completes
}


template <typename T>
void MHPCLocomotion<T>::publish_mpc_cmd()
{
    
    int nControlSteps = opt_problem.get_num_control_steps();

    nControlSteps = 3; // use 3 more controls than control duration to account for delay

    mpc_cmd.N_mpcsteps = nControlSteps;

    const auto& trajs = opt_problem_data.wb_trajs;
    const auto& ctacts = opt_problem_data.wb_phase_contacts;    
    const auto& statusDurations = opt_problem_data.wb_contact_durations;
    
    int pidx(0), k_rel(0); // phase index and relative time instant w.r.t. each phase
    for (int k = 0; k < nControlSteps; k++)
    {        
        opt_problem_data.get_index(k, pidx, k_rel);

        const VecM<float, WBM::xs>& x_float = trajs[pidx]->Xbar[k_rel].template cast<float>();
        const VecM<float, WBM::us>& u_float = trajs[pidx]->Ubar[k_rel].template cast<float>();

        std::copy(x_float.begin(), x_float.begin()+3, mpc_cmd.pos[k]);

        std::copy(x_float.data()+3, x_float.data()+6, mpc_cmd.eul[k]);

        std::copy(x_float.data()+6, x_float.data()+18, mpc_cmd.qJ[k]);

        std::copy(x_float.data()+18, x_float.data()+21, mpc_cmd.vWorld[k]);

        std::copy(x_float.data()+21, x_float.data()+24, mpc_cmd.eulrate[k]);        

        std::copy(x_float.data()+24, x_float.data()+36,mpc_cmd.qJd[k]);

        std::copy(u_float.begin(), u_float.end(), mpc_cmd.torque[k]);        

        std::copy(ctacts[pidx].begin(), ctacts[pidx].end(), mpc_cmd.contacts[k]);

        std::copy(statusDurations[pidx].begin(), statusDurations[pidx].end(), mpc_cmd.statusTimes[k]);

        mpc_cmd.mpc_times[k] = mpc_time + k*mpc_config.dt_mpc;
    }            
    mpc_cmd.solve_time = solve_time;
    mpc_lcm.publish("MHPC_COMMAND", &mpc_cmd);

    printf(GRN);
    printf("published a mpc command message \n");
    printf(RESET);
}



template class MHPCLocomotion<double>;