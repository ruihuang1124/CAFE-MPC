
#include "PinocchioInteface.h"
#include "WBM.h"
#include "MHPCReset.h"
#include "BarrelRollConstraints.h"
#include "SinglePhaseInterface.h"
#include "SinglePhase.h"
#include "MultiPhaseDDP.h"

#include <functional>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iterator>

#include <lcm/lcm-cpp.hpp>
#include "wbTraj_lcmt.hpp"

using vectord = std::vector<double>;
using WBM_d = WBM::Model<double>;
using WBSinglePhase_d = SinglePhase<double, WBM::xs, WBM::us, WBM::ys>;
using SinglePhase_d = SinglePhaseBase<double>;
using WBSingleTrajectory_d = Trajectory<double, WBM::xs, WBM::us, WBM::ys>;
using WBTrackCost_d = QuadraticTrackingCost<double, WBM::xs, WBM::us, WBM::ys>;
using Vec3d = Vec3<double>;
using Vec4d = Vec4<double>;
using Vec12d = VecM<double, 12>;
using Vec36d = VecM<double, WBM::xs>;
using REB_Paramd = REB_Param_Struct<double>;
using AL_Paramd = AL_Param_Struct<double>;

namespace pc = std::placeholders;

template <typename T1, typename T2>
void lerp_eigen_vectors(const Eigen::MatrixBase<T1> &m0,
                        const Eigen::MatrixBase<T1> &m1,
                        float dur,
                        float t,
                        Eigen::MatrixBase<T2> const &mt)
{

    const_cast<Eigen::MatrixBase<T2> &>(mt) =
        m0 + (m1 - m0) * (t / dur);
}

void load_desired_final_states(vector<Vec36d> &x_des);

void load_cost_weights(vector<vectord> &weights, const string &fname);

void load_reb_params(REB_Paramd &params, const string &fname, std::string c_type);

void load_al_params(AL_Paramd &params, const string &fname, std::string c_type);

void publish_trajectory(const deque<shared_ptr<WBSingleTrajectory_d>>& trajs,
                        const vector<Vec4<int>>& contacts);

int main()
{
    /* Swiching times: Full stance -> Right foot stance -> Flight -> Landing  */
    const int num_phases = 4;
    double dt = 0.01;
    vectord switching_times{0.0, 0.15, 0.3, 0.7, 1.2};
    vector<int> horizons(num_phases);
    deque<shared_ptr<SinglePhase_d>> phases(num_phases);
    deque<shared_ptr<WBSingleTrajectory_d>> trajs(num_phases);
    vector<Vec4<int>> contacts(num_phases);
    contacts[0].setOnes();
    contacts[1] << 0, 1, 0, 1; // FL, FR, HL, HR
    contacts[2].setZero();
    contacts[3].setOnes();

    /* Create a whole-body model */
    const std::string urdf_filename = "../urdf/mini_cheetah_simple_correctedInertia.urdf";
    pinocchio::ModelTpl<double> pin_model;
    double BG_alpha = 10.0; // First-order Baumgard stabilizatoin
    buildPinModelFromURDF(urdf_filename, pin_model);
    shared_ptr<WBM_d> wbm_ptr;
    wbm_ptr = std::make_shared<WBM_d>(pin_model, BG_alpha);

    /* Reset map */
    MHPCReset<double> reset_ptr(wbm_ptr.get(), nullptr);

    /* Initial condition */
    Vec3d pos, eul, vWorld, euld;
    Vec12d qJ, qJd;
    Vec36d xinit;
    pos.setZero();
    eul.setZero();
    vWorld.setZero();
    euld.setZero();
    pos[2] = 0.1464;
    qJ = Vec3d(0, 1.2, -2.4).replicate<4, 1>();
    qJd.setZero();
    xinit << pos, eul, qJ, vWorld, euld, qJd;

    /* Desired final states of each phase */
    vector<Vec36d> xf_des(num_phases); // first one is initial condition
    load_desired_final_states(xf_des);

    /* Configuration files */
    string fpath = "../MHPC/MHPC-Trajopt/BarrelRoll/setting/";
    string ddp_setting_fname = fpath + "br_ddp_setting.info";
    string cost_weights_fname = fpath + "br_cost_weights.JSON";
    string constraint_params_fname = fpath + "br_constraint_params.info";

    /* Cost weights */
    vector<vectord> cost_weights;
    load_cost_weights(cost_weights, cost_weights_fname);

    for (size_t i = 0; i < num_phases; i++)
    {
        horizons[i] = static_cast<int>(std::round((switching_times[i + 1] - switching_times[i]) / dt));
        shared_ptr<WBSinglePhase_d> phase;
        phase = make_shared<WBSinglePhase_d>();
        trajs[i] = make_shared<WBSingleTrajectory_d>(dt, horizons[i]);
        phase->set_time_offset(switching_times[i]);

        std::cout << "Trajectory size = " << trajs[i]->size() << "\n";
        // Initialize state trajectory using linear interpolation
        float t = 0.0;
        float t_dur = switching_times[i + 1] - switching_times[i];
        for (int k = 0; k <= horizons[i]; k++)
        {
            if (i == 0)
            {
                lerp_eigen_vectors(xinit, xf_des[i], t_dur, t, trajs[i]->Xbar[k]);
            }
            else
            {
                lerp_eigen_vectors(xf_des[i - 1], xf_des[i], t_dur, t, trajs[i]->Xbar[k]);
            }
            t += dt;
        }

        // assign trajectory
        phase->set_trajectory(trajs[i]);

        // assign dynamics
        auto dynamics_callback = bind(&WBM_d::dynamics, wbm_ptr,
                                      pc::_1, pc::_2, pc::_3, pc::_4, pc::_5,
                                      contacts[i], dt);
        phase->set_dynamics(dynamics_callback);

        // assign dynamics partial
        auto dynamics_partial_callback =
            bind(&WBM_d::dynamics_partial, wbm_ptr,
                 pc::_1, pc::_2, pc::_3, pc::_4, pc::_5, pc::_6, pc::_7,
                 contacts[i], dt);
        phase->set_dynamics_partial(dynamics_partial_callback);

        // assign reset map except for the last phase
        if (i < num_phases - 1)
        {
            auto resetmap_callback = bind(&MHPCReset<double>::reset_map, reset_ptr,
                                          pc::_1, pc::_2,
                                          contacts[i], contacts[i + 1], ModelType::WB, ModelType::WB);
            auto resetmap_partial_callback = bind(&MHPCReset<double>::reset_map_partial, reset_ptr,
                                                  pc::_1, pc::_2,
                                                  contacts[i], contacts[i + 1], ModelType::WB, ModelType::WB);

            phase->set_resetmap(resetmap_callback);
            phase->set_resetmap_partial(resetmap_partial_callback);
        }

        // tracking cost
        auto tracking_cost = make_shared<WBTrackCost_d>();
        tracking_cost->set_reference_state(xf_des[i], Vec12d::Zero());
        tracking_cost->update_weighting_matrix(cost_weights[i]);
        phase->add_cost(tracking_cost);

        /* Torque limit */
        REB_Paramd torque_reb_param;
        load_reb_params(torque_reb_param, constraint_params_fname, "Torque");
        shared_ptr<BarrelRoll::TorqueLimit<double>> torqueLimit;
        torqueLimit = std::make_shared<BarrelRoll::TorqueLimit<double>>();
        torqueLimit->update_horizon_len(horizons[i]);
        torqueLimit->create_data();
        torqueLimit->initialize_params(torque_reb_param);
        phase->add_pathConstraint(torqueLimit);

        /* Joint speed limit */
        REB_Paramd jointvel_reb_param;
        load_reb_params(jointvel_reb_param, constraint_params_fname, "JointVel");
        shared_ptr<BarrelRoll::JointSpeedLimit<double>> jointSpeedLimit;
        jointSpeedLimit = std::make_shared<BarrelRoll::JointSpeedLimit<double>>();
        jointSpeedLimit->update_horizon_len(horizons[i]);
        jointSpeedLimit->create_data();
        jointSpeedLimit->initialize_params(jointvel_reb_param);
        phase->add_pathConstraint(jointSpeedLimit);

        /* Minimum Height constraint */
        REB_Paramd minheight_reb_param;
        load_reb_params(minheight_reb_param, constraint_params_fname, "MinHeight");
        shared_ptr<BarrelRoll::MinimumHeight<double>> wbMinHeightConstraint;
        wbMinHeightConstraint = std::make_shared<BarrelRoll::MinimumHeight<double>>();
        wbMinHeightConstraint->update_horizon_len(horizons[i]);
        wbMinHeightConstraint->create_data();
        wbMinHeightConstraint->initialize_params(minheight_reb_param);
        phase->add_pathConstraint(wbMinHeightConstraint);

        /* Set GRF constraints if any*/
        if (contacts[i].cwiseEqual(1).any())
        {
            REB_Paramd grf_reb_param;
            load_reb_params(grf_reb_param, constraint_params_fname, "GRF");
            shared_ptr<BarrelRoll::GRF<double>> grfConstraint;
            grfConstraint = std::make_shared<BarrelRoll::GRF<double>>(contacts[i]);
            grfConstraint->update_horizon_len(horizons[i]);
            grfConstraint->create_data();
            grfConstraint->initialize_params(grf_reb_param);
            phase->add_pathConstraint(grfConstraint);
        }

        // Touchdown constraint for landing phase
        if (i == 2)
        {
            AL_Paramd td_al_param;
            load_al_params(td_al_param, constraint_params_fname, "TD");
            shared_ptr<BarrelRoll::TouchDown<double>> tdConstraint;
            tdConstraint = std::make_shared<BarrelRoll::TouchDown<double>>(Vec4<int>::Ones(), wbm_ptr.get());
            tdConstraint->create_data();
            tdConstraint->initialize_params(td_al_param);
            phase->add_terminalConstraint(tdConstraint);
        }

        phase->initialization();
        phase->update_SS_config(horizons[i] + 1);
        phases[i] = phase;
    }
    MultiPhaseDDP<double> solver;
    HSDDP_OPTION ddp_setting;
    loadHSDDPSetting(ddp_setting_fname, ddp_setting);
    solver.set_initial_condition(xinit);
    solver.set_multiPhaseProblem(phases);    
    solver.solve(ddp_setting);

    publish_trajectory(trajs, contacts);
}

void load_desired_final_states(vector<Vec36d> &x_des)
{
    Vec3d pos, eul, vWorld, euld;
    Vec12d qJ, qJd;
    Vec36d xdes_phase_i;
    pos.setZero();
    eul.setZero();
    vWorld.setZero();
    euld.setZero();
    qJ = Vec3d(0, 1.2, -2.4).replicate<4, 1>();
    qJd.setZero();

    // Desired final state for the first phase  (stance) 
    pos[2] = 0.16; 
    eul[2] = M_PI/3;
    euld[2] = 1.5*M_PI;
    vWorld[2] = 1.0;
    xdes_phase_i << pos, eul, qJ, vWorld, euld, qJd;
    x_des[0] = xdes_phase_i;

    // Desired final state for the second phase (right stance)
    pos[2] = 0.18; 
    eul[2] = 0.7*M_PI;
    euld[2] = 3.5*M_PI;
    vWorld[2] = 2.5;
    xdes_phase_i << pos, eul, qJ, vWorld, euld, qJd;
    x_des[1] = xdes_phase_i;

    // Desired final state for the third phase (air stance)
    pos[2] = 0.18; 
    eul[2] = 2*M_PI;
    euld[2] = 3.5*M_PI;
    vWorld[2] = -2.5;
    // qJ[0] = 0.2; qJ[3] = -0.2; qJ[7] = 0.2; qJ[10] = -0.2;
    xdes_phase_i << pos, eul, qJ, vWorld, euld, qJd;
    x_des[2] = xdes_phase_i;

    // Desired final state for the last phase (stance stance)
    pos[2] = 0.18; 
    eul[2] = 2*M_PI;
    euld[2] = 0;
    vWorld[2] = 0;
    xdes_phase_i << pos, eul, qJ, vWorld, euld, qJd;
    x_des[3] = xdes_phase_i;
}

void load_cost_weights(vector<vectord> &weights, const string &fname)
{
    // weights for phase 1
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(fname, pt);

    for (int pidx = 0; pidx < 4; pidx++)
    {
        string cost_phase = "cost_phase_" + to_string(pidx+1) + ".";

        vectord qw;
        for (auto &item : pt.get_child(cost_phase+"qw_qB"))
        {
            qw.push_back(item.second.get_value<double>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto &item : pt.get_child(cost_phase+"qw_qJ"))
            {
                qw.push_back(item.second.get_value<double>());
            }
        }
        for (auto &item : pt.get_child(cost_phase+"qw_vB"))
        {
            qw.push_back(item.second.get_value<double>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto &item : pt.get_child(cost_phase+"qw_vJ"))
            {
                qw.push_back(item.second.get_value<double>());
            }
        }
        vectord rw;
        for (size_t i = 0; i < 12; i++)
        {
            rw.push_back(pt.get<double>(cost_phase+"rw"));
        }
        vectord qfw;
        for (auto &item : pt.get_child(cost_phase+"qfw_qB"))
        {
            qfw.push_back(item.second.get_value<double>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto &item : pt.get_child(cost_phase+"qfw_qJ"))
            {
                qfw.push_back(item.second.get_value<double>());
            }
        }
        for (auto &item : pt.get_child(cost_phase+"qfw_vB"))
        {
            qfw.push_back(item.second.get_value<double>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto &item : pt.get_child(cost_phase+"qfw_vJ"))
            {
                qfw.push_back(item.second.get_value<double>());
            }
        }

        vectord weight_phase;
        weight_phase.insert(weight_phase.end(),
                            std::make_move_iterator(qw.begin()),
                            std::make_move_iterator(qw.end()));
        weight_phase.insert(weight_phase.end(),
                            std::make_move_iterator(rw.begin()),
                            std::make_move_iterator(rw.end()));
        weight_phase.insert(weight_phase.end(),
                            std::make_move_iterator(qfw.begin()),
                            std::make_move_iterator(qfw.end()));
        weights.push_back(weight_phase);
    }
}

void load_reb_params(REB_Paramd &params, const string &fname, std::string constr_type)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(fname, pt);

    string reb_name = constr_type+"_ReB";
    params.delta = pt.get<double>(reb_name+".delta");
	params.delta_min = pt.get<double>(reb_name+".delta_min");
	params.eps = pt.get<double>(reb_name+".eps");
    
}

void load_al_params(AL_Paramd &params, const string &fname, std::string constr_type)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(fname, pt);

    string td_name = constr_type+"_AL";
    params.sigma = pt.get<double>(td_name+".sigma");
	params.lambda = pt.get<double>(td_name+".lambda");	
    params.sigma_max = pt.get<double>(td_name+".sigma_max");
}

void publish_trajectory(const deque<shared_ptr<WBSingleTrajectory_d>>& trajs,
                        const vector<Vec4<int>>& contacts)
{
    /* Publish wb trajectory for viz via lcm */
    lcm::LCM viz_lcm;    
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
            t+=tau->timeStep;            
        }        
    }

    viz_lcm.publish("visualize_wb_traj", &wbtraj_lcmt);
}