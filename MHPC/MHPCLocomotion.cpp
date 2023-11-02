/*!
 * @file opt_main.cpp
 * @brief Main Function for the standalone DDP trajectory optimizer
 *
 * The main function initilizes the DDP solver updates upon the request
 * of the main Imitation Controller
 */

#include "MHPCLocomotion.h"

#include <chrono>
using namespace std::chrono;
using duration_ms = std::chrono::duration<float, std::chrono::milliseconds::period>;

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
    std::string mhpc_config_file("../MHPC/settings/mhpc_config.info");
    loadMHPCConfig(mhpc_config_file, mpc_config);
    mpc_config.print();

    // Load DDP setting
    std::string ddp_setting_file("../MHPC/settings/ddp_setting.info");
    loadHSDDPSetting(ddp_setting_file, ddp_setting);

    // Load reference trajectory
    printf("Loading quadruped reference ... \n");
    std::string quad_reference_file("../Reference/Data/");
    quad_reference_file.append(mpc_config.referenceFileName);
    quad_reference_file.append("/quad_reference.csv");
    opt_problem_data.quad_reference = std::make_shared<QuadReference>();
    opt_problem_data.quad_reference->load_top_level_data(quad_reference_file, false);

    // Initialize the multi-phase OCP (phase determination, memory allocation)
    opt_problem.set_problem_data(&opt_problem_data, &mpc_config);
    opt_problem.initialization();

#ifdef DEBUG_MODE
    opt_problem.pretty_print();
#endif

    // mpc_time = 0;
    // mpc_time_prev = 0;

    // Assemble the multi-phase probelm and solve it with MSDDP
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (const auto &phase : opt_problem_data.wb_phases)
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

#ifdef DEBUG_MODE
    mhpc_viz.publishWBTrajectory(&opt_problem_data, mpc_config);
#endif

    // run-time DDP setting when re-solving DDP in MPC
    ddp_setting.max_AL_iter = ddp_setting.max_AL_iter_runtime;
    ddp_setting.max_DDP_iter = ddp_setting.max_DDP_iter_runtime;
}

template <typename T>
void MHPCLocomotion<T>::update()
{
    mpc_mutex.lock(); // lock mpc to prevent updating while the previous hasn't finished
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

    /* build solver and solve the TO problem */
    MultiPhaseDDP<T> solver;
    deque<shared_ptr<SinglePhaseBase<T>>> multiple_phases;
    for (const auto &phase : opt_problem_data.wb_phases)
    {
        multiple_phases.push_back(phase);
    }
    if (opt_problem_data.srb_phase.get() != nullptr)
    {
        multiple_phases.push_back(opt_problem_data.srb_phase);
    }
    solver.set_multiPhaseProblem(multiple_phases);
    solver.set_initial_condition(x_init_wb);

// #ifdef TIME_BENCHMARK
    static float solve_time_acc = 0;
    static int solve_count = 0;
    auto solve_start = high_resolution_clock::now();
// #endif
    solver.solve(ddp_setting, mpc_config.dt_mpc*1000);
// #ifdef TIME_BENCHMARK
    auto solve_stop = high_resolution_clock::now();
    auto solve_duration = duration_ms(solve_stop - solve_start);
    solve_time = solve_duration.count();
    solve_time_acc += solve_time;
    solve_count ++;
    printf("average solve time = %f ms \n", solve_time_acc/solve_count);
// #endif

#ifdef TIME_BENCHMARK
    static float time_publish_mpc_acc = 0;
    static int count_publish_mpc = 0;
    auto time_pub_mpc_start = high_resolution_clock::now();
#endif
    publish_mpc_cmd();
#ifdef TIME_BENCHMARK
    auto time_pub_mpc_end = high_resolution_clock::now();
    auto duration_pub_mpc = duration_ms(time_pub_mpc_end - time_pub_mpc_start);
    time_publish_mpc_acc += duration_pub_mpc.count();
    count_publish_mpc ++;
    std::cout << "average publishing time = " << time_publish_mpc_acc/count_publish_mpc << " ms \n";
#endif

#ifdef DEBUG_MODE
    mhpc_viz.publishWBTrajectory(&opt_problem_data, mpc_config);
#endif

    mpc_mutex.unlock();
}

template <typename T>
void MHPCLocomotion<T>::mpcdata_lcm_handler(const lcm::ReceiveBuffer *rbuf, const std::string &chan,
                                            const MHPC_Data_lcmt *msg)
{
    (void)(rbuf);
    (void)(chan);
    mpc_mutex.lock();
    printf(GRN);
    printf("Received resolving request\n");
    printf(RESET);

    /* update current state*/
    Eigen::Map<const VecM<float, 3>> pos_map(msg->pos);
    Eigen::Map<const VecM<float, 3>> eul_map(msg->eul);
    Eigen::Map<const VecM<float, 12>> qJ_map(msg->qJ);
    Eigen::Map<const VecM<float, 3>> vWorld_map(msg->vWorld);
    Eigen::Map<const VecM<float, 3>> eulrate_map(msg->eulrate);
    Eigen::Map<const VecM<float, 12>> qJd_map(msg->qJd);
    x_init_wb << pos_map.cast<T>(), eul_map.cast<T>(), qJ_map.cast<T>(),
        vWorld_map.cast<T>(), eulrate_map.cast<T>(), qJd_map.cast<T>();
    mpc_time_prev = mpc_time;
    mpc_time = msg->mpctime;
    
    mpc_mutex.unlock();

    if (is_first_mpc)
    {        
        initialize();
        is_first_mpc = false;
    }else
    {        
        std::thread solve_mpc_thread(&MHPCLocomotion::update, this);
        solve_mpc_thread.detach(); // detach the thread from the main thread. The thread would exit once it completes    
    }
    
}

template <typename T>
void MHPCLocomotion<T>::publish_mpc_cmd()
{

    int nControlSteps = opt_problem.get_num_control_steps();

    nControlSteps = 5; // use 4 controls than control duration to account for delay

    mpc_cmd.N_mpcsteps = nControlSteps;

    const auto &trajs = opt_problem_data.wb_trajs;
    const auto &ctacts = opt_problem_data.wb_phase_contacts;
    const auto &statusDurations = opt_problem_data.wb_contact_durations;

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
    std::vector<int> contact_k(4);
    std::vector<float> statusDuration_k(4);

    // Clear old data
    mpc_cmd.pos.clear();
    mpc_cmd.eul.clear();
    mpc_cmd.qJ.clear();
    mpc_cmd.vWorld.clear();
    mpc_cmd.eulrate.clear();
    mpc_cmd.qJd.clear();
    mpc_cmd.torque.clear();
    mpc_cmd.GRF.clear();
    mpc_cmd.Qu.clear();
    mpc_cmd.Quu.clear();
    mpc_cmd.Qux.clear();
    mpc_cmd.feedback.clear();
    mpc_cmd.contacts.clear();
    mpc_cmd.statusTimes.clear();
    mpc_cmd.mpc_times.clear();

    int pidx(0), k_rel(0); // phase index and relative time instant w.r.t. each phase
    for (int k = 0; k < nControlSteps; k++)
    {
        opt_problem_data.get_index(k, pidx, k_rel);

        const VecM<float,36>&x_k = trajs[pidx]->Xbar[k_rel].template cast<float>();
        const Vec12<float>&u_k = trajs[pidx]->Ubar[k_rel].template cast<float>();
        const Vec12<float>&GRF_k = trajs[pidx]->Y[k_rel].template cast<float>();
        const Vec12<float>&Qu_k = trajs[pidx]->Qu[k_rel].template cast<float>();
        const MatMN<float, 12, 12>&Quu_k = trajs[pidx]->Quu[k_rel].template cast<float>();
        const MatMN<float, 12, 36>&Qux_k = trajs[pidx]->Qux[k_rel].template cast<float>();
        const MatMN<float, 12, 36>&K_k = trajs[pidx]->K[k_rel].template cast<float>();        
        
        std::copy(u_k.begin(), u_k.end(), torque_k_float.data());
        std::copy(x_k.begin(), x_k.begin() + 3, pos_k_float.data());
        std::copy(x_k.begin() + 3, x_k.begin() + 6, eul_k_float.data());
        std::copy(x_k.begin() + 6, x_k.begin() + 18, qJ_k_float.data());
        std::copy(x_k.begin() + 18, x_k.begin() + 21, vWorld_k_float.data());
        std::copy(x_k.begin() + 21, x_k.begin() + 24, eulrate_k_float.data());
        std::copy(x_k.begin() + 24, x_k.end(), qJd_k_float.data());
        std::copy(GRF_k.begin(), GRF_k.end(), GRF_k_float.data());
        std::copy(Qu_k.begin(), Qu_k.end(), Qu_k_float.data());

        std::copy(Quu_k.data(), Quu_k.data()+144, Quu_k_float.data());
        std::copy(Qux_k.data(), Qux_k.data()+432, Qux_k_float.data());
        std::copy(K_k.data(), K_k.data()+432, feedback_k_float.data());
        std::copy(ctacts[pidx].begin(), ctacts[pidx].end(), contact_k.data());
        std::copy(statusDurations[pidx].begin(), statusDurations[pidx].end(), statusDuration_k.data());

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
        mpc_cmd.contacts.push_back(contact_k);
        mpc_cmd.statusTimes.push_back(statusDuration_k);
        mpc_cmd.mpc_times.push_back(mpc_time + k * mpc_config.dt_wb);
    }
    mpc_cmd.solve_time = solve_time;
    mpc_lcm.publish("MHPC_COMMAND", &mpc_cmd);

    printf(GRN);
    printf("published a mpc command message \n");
    printf(RESET);
}

template class MHPCLocomotion<double>;