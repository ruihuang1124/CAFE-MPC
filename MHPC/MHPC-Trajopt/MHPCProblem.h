#pragma once
#ifndef MHPC_PROBLEM_H
#define MHPC_PROBLEM_H

#include "SRBM.h"
#include "WBM.h"
#include "MHPCReset.h"
#include "PinocchioInteface.h"
#include "QuadReference.h"
#include "MHPCReference.h"
#include "SinglePhase.h"
#include "TrajectoryManagement.h"
#include "MHPCFootStep.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

struct MHPCConfig
{
    // WB plan duration (in seconds)
    double plan_dur_wb;

    // SRB plan duration (in seconds)
    double plan_dur_srb;    

    // WB simulation timestep (in seconds)
    double dt_wb;

    // SRB simulation timestep (in seconds)
    double dt_srb;

    // how often MPC is updated
    float dt_mpc;       

    // Bamguart parameter (velocity feedback gain)
    float BG_alpha;   

    std::string referenceFileName;

    void print(){
        std::cout << "===================== MHPC Config =================== \n";
        std::cout << "WB plan duration = \t" << plan_dur_wb << "\n";
        std::cout << "WB simulation timestep = \t" << dt_wb << "\n";
        std::cout << "SRB plan duration = \t" << plan_dur_srb << "\n";        
        std::cout << "SRB simulation timestep = \t" << dt_srb << "\n";
        std::cout << "MPC updates every " << dt_mpc << " seconds" << "\n";
        std::cout << "Baumgart alpha (vel) " << BG_alpha << "\n";        
        std::cout << "Reference " << referenceFileName << "\n";
    }
};

inline void loadMHPCConfig(const std::string filename, MHPCConfig& config)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);

    std::cout << "********* loading MHPC configuration from file *********\n" << filename << "\n\n";
    config.plan_dur_wb = pt.get<double>("config.plan_dur_wb");
    config.plan_dur_srb = pt.get<double>("config.plan_dur_srb");
    config.dt_mpc = pt.get<double>("config.dt_mpc");
    config.dt_wb = pt.get<double>("config.dt_wb");
    config.dt_srb = pt.get<double>("config.dt_srb");
    config.BG_alpha = pt.get<double>("config.BG_alpha");    
    config.referenceFileName = pt.get<std::string>("config.referenceFile");
}

template<typename T>
inline void loadConstrintParameters(const std::string&fileName, 
                                    REB_Param_Struct<T>& GRF_reb_param,
                                    REB_Param_Struct<T>& Torque_reb_param,
                                    AL_Param_Struct<T>& TD_al_param)
{
	boost::property_tree::ptree pt;
    boost::property_tree::read_info(fileName, pt);
	std::cout << "********* loading MHPC Constraint Parameter from file **********\n" << fileName << "\n\n";

	GRF_reb_param.delta = pt.get<T>("GRF_ReB.delta");
	GRF_reb_param.delta_min = pt.get<T>("GRF_ReB.delta_min");
	GRF_reb_param.eps = pt.get<T>("GRF_ReB.eps");

    Torque_reb_param.delta = pt.get<T>("Torque_ReB.delta");
	Torque_reb_param.delta_min = pt.get<T>("Torque_ReB.delta_min");
	Torque_reb_param.eps = pt.get<T>("Torque_ReB.eps");	

    TD_al_param.sigma = pt.get<T>("TD_AL.sigma");
	TD_al_param.lambda = pt.get<T>("TD_AL.lambda");	
    TD_al_param.sigma_max = pt.get<T>("TD_AL.sigma_max");
}



template <typename T>
struct MHPCProblemData
{
    std::deque<shared_ptr<Trajectory<T, WBM::xs, WBM::us, WBM::ys>>> wb_trajs;
    std::deque<shared_ptr<SinglePhase<T, WBM::xs, WBM::us, WBM::ys>>> wb_phases;

    shared_ptr<Trajectory<T, SRBM::xs, SRBM::us, SRBM::ys>> srb_traj = nullptr;
    shared_ptr<SinglePhase<T, SRBM::xs, SRBM::us, SRBM::ys>> srb_phase = nullptr;

    std::shared_ptr<QuadReference> quad_reference = nullptr;

    std::deque<int> wb_phase_horizons;
    std::deque<bool> wb_is_phase_reach_end;
    std::deque<float> wb_phase_start_times;
    std::deque<float> wb_phase_end_times;
    std::deque<VecM<double, 4>> wb_contact_durations;
    std::deque<VecM<int, 4>> wb_phase_contacts;
    int n_wb_phases = 0;

    float srb_start_time = 0;
    float srb_end_time = 0;
    int   srb_phase_horizon = 0;
    int n_srb_phase = 0;

    // Return where I am  (the closest) along and within the phase
    void get_index(const int& k_cur, int& pidx, int& k_pidx){
        pidx = 0; k_pidx = 0;                                
        int s_i(0);
        for (int i(0); i < n_wb_phases; ++i)
        {
            int h = wb_phase_horizons[i];
            if (k_cur >= s_i && k_cur < s_i + h)
            {
                pidx = i;
                k_pidx = k_cur - s_i;
                break;
            }
            s_i += h;            
        }        
    }

    void clear(){
        wb_trajs.clear();
        wb_phases.clear();
        srb_traj = nullptr;
        srb_phase = nullptr;
        quad_reference = nullptr;

        wb_phase_horizons.clear();
        wb_is_phase_reach_end.clear();
        wb_phase_start_times.clear();
        wb_phase_end_times.clear();
        wb_contact_durations.clear();
        wb_phase_contacts.clear();
        n_wb_phases = 0;

        srb_start_time = 0;
        srb_end_time = 0;
        srb_phase_horizon = 0;
        n_srb_phase = 0;
    }

    void pop_front_phase(){
        if (!wb_trajs.empty())
        {
            wb_trajs.pop_front();
            wb_phases.pop_front();
            wb_phase_horizons.pop_front();
            wb_is_phase_reach_end.pop_front();
            wb_phase_start_times.pop_front();
            wb_phase_end_times.pop_front();
            wb_contact_durations.pop_front();
            wb_phase_contacts.pop_front();
            n_wb_phases--;
        }
        
    }
    
};



template <typename T>
class MHPCProblem
{
public:
    MHPCProblem():pdata(nullptr),
                  quad_reference(nullptr),
                  pconfig(nullptr),
                  plan_dur_all(0.0),
                  wb_nsteps_between_mpc(0),
                  srb_nsteps_between_mpc(0) {}

    void set_problem_data(MHPCProblemData<T>* pdata_in, 
                          const MHPCConfig* pconfig_in)
    {
        pdata = pdata_in;
        pconfig = pconfig_in;        

        quad_reference = pdata->quad_reference;
        plan_dur_all = pconfig->plan_dur_wb + pconfig->plan_dur_srb;
        wb_nsteps_between_mpc = (int) round(pconfig->dt_mpc/pconfig->dt_wb);
        srb_nsteps_between_mpc = (int) round(pconfig->dt_mpc/pconfig->dt_srb);
    }

    void clear_problem_data() { 
        pdata->clear();
        plan_dur_all = 0;
        wb_nsteps_between_mpc = 0;
        srb_nsteps_between_mpc = 0;

        reset_ptr.reset();
        wbm_ptr.reset();
        srbm_ptr.reset();        
    }
    
    void initialization();

    void initialize_multiPhaseProblem();    

    void prepare_initialization();

    void initialize_parameters();

    void update();

    void update_WB_plan();

    void update_SRB_plan();

    void create_problem_one_phase(shared_ptr<SinglePhase<T, WBM::xs, WBM::us, WBM::ys>>, int phase_idx);

    void create_problem_one_phase(shared_ptr<SinglePhase<T, SRBM::xs, SRBM::us, SRBM::ys>>);

    void update_resetmap(shared_ptr<SinglePhase<T, WBM::xs, WBM::us, WBM::ys>>, int phase_idx);

    void add_tconstr_one_phase(shared_ptr<SinglePhase<T, WBM::xs, WBM::us, WBM::ys>>, int phase_idx);

    int get_num_control_steps() {return (int) round(pconfig->dt_mpc / pconfig->dt_wb);}    

    void pretty_print();   

private:
    MHPCProblemData<T>* pdata;
    const MHPCConfig* pconfig;
    std::shared_ptr<QuadReference> quad_reference;
    
    std::shared_ptr<WBM::Model<T>> wbm_ptr;
    std::shared_ptr<SRBM::Model<T>> srbm_ptr;
    std::shared_ptr<MHPCReset<T>> reset_ptr;
    
    pinocchio::ModelTpl<T> pin_model;

    WBReference wb_reference;
    SRBReference srb_reference;
    std::shared_ptr<MHPCFootStep<T>> footStepPlanner;


    // Overall planning horizon in seconds
    float plan_dur_all;

    // Number of WB simulation timesteps beween two consecutive MPC updates                       
    int wb_nsteps_between_mpc;          

    // Number of SRB simulation timesteps between two consecutie MPC udpates
    int srb_nsteps_between_mpc;

    // Define constraint parameters (for WB plan)
    REB_Param_Struct<T> grf_reb_param;
    REB_Param_Struct<T> torque_reb_param;
    AL_Param_Struct<T> td_al_param;  
};

#endif