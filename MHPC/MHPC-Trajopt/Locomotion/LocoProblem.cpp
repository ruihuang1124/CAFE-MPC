#include "LocoProblem.h"
#include "MHPCCostUtil.h"

namespace pc = std::placeholders;

template <typename T>
void LocoProblem<T>::initialize_parameters()
{
    /* Initialize REB and AL parameters for ineq and eq constraints */    
    const std::string& constraint_setting_fname= "../"+ this->pconfig->constraintParamFileName;
    load_reb_params(this->grf_reb_param, constraint_setting_fname, "GRF");
    load_reb_params(this->torque_reb_param, constraint_setting_fname, "Torque");
    load_al_params(this->td_al_param, constraint_setting_fname, "TD");    

    if (!this->pconfig->costFileName.empty())
    {
        const std::string& cost_weights_fname= "../"+ this->pconfig->costFileName;
        loadCostWeights(cost_weights_fname, 
                        this->wb_track_cost,
                        this->wb_foot_reg,
                        this->swing_pos_tracking,
                        this->swing_vel_tracking,
                        this->srb_track_cost);    
    }
        
}

template <typename T>
void LocoProblem<T>::create_problem_one_phase(shared_ptr<WBPhase_T>phase, int idx)
{
    /* specialize dynamics and resetmap  */
    const auto &phase_contact = this->pdata->wb_phase_contacts[idx];

    auto dynamics_callback = bind(&WBM::Model<T>::dynamics, this->wbm_ptr[0],
                                  pc::_1, pc::_2, pc::_3, pc::_4, pc::_5,
                                  phase_contact, this->pconfig->dt_wb);

    vector<function<void(WBStateMap&, WBContrlMap&, WBOutputMap&, 
                         WBDirectMap&, WBState&, WBContrl&, T)>> dynamics_partial_callback(this->pconfig->num_threads);

    for (int i = 0; i < this->pconfig->num_threads; i++)
    {        
        dynamics_partial_callback[i] =
            bind(&WBM::Model<T>::dynamics_partial, this->wbm_ptr[i],
             pc::_1, pc::_2, pc::_3, pc::_4, pc::_5, pc:: _6, pc::_7,
             phase_contact, this->pconfig->dt_wb);
    }
        

    /* set dynamics */
    phase->set_dynamics(dynamics_callback);
    phase->set_dynamics_partial(dynamics_partial_callback);

    /* set tracking cost */    
    phase->add_cost(this->wb_track_cost);

    /*Set foot regularization */
    phase->add_cost(this->wb_foot_reg);

    /* Swing foot posiiton tracking */   
    phase->add_cost(this->swing_pos_tracking);

    /* Swing foot velocity tracking */
    phase->add_cost(this->swing_vel_tracking);

    /* Torque limit */    
    shared_ptr<MHPCConstraints::TorqueLimit<T>> torqueLimit;
    torqueLimit = std::make_shared<MHPCConstraints::TorqueLimit<T>>();
    torqueLimit->update_horizon_len(this->pdata->wb_phase_horizons[idx]);
    torqueLimit->create_data();
    torqueLimit->initialize_params(this->torque_reb_param);
    phase->add_pathConstraint(torqueLimit);    


    /* Set GRF constraints if any*/
    if (phase_contact.cwiseEqual(1).any())
    {
        shared_ptr<MHPCConstraints::WBGRF<T>> grfConstraint;
        grfConstraint = std::make_shared<MHPCConstraints::WBGRF<T>>(phase_contact);
        grfConstraint->update_horizon_len(this->pdata->wb_phase_horizons[idx]);
        grfConstraint->create_data();
        grfConstraint->initialize_params(this->grf_reb_param);
        phase->add_pathConstraint(grfConstraint);
    }    
}

template class LocoProblem<double>;
