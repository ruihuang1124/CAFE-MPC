#include <functional> // std::bind, and std::placeholder
#include <tabulate/table.hpp>
#include "MHPCProblem.h"
#include "MHPCConstraint.h"
#include "TrajectoryManagement.h"
#include "HSDDP_Utils.h"
#include "utilities.h"
#include <tabulate/table.hpp>   // pretty print of tables


namespace pc = std::placeholders;

template <typename T>
void MHPCProblem<T>::initialization()
{
    printf("\n============Initializing MHPCProblem ============ \n");

    prepare_initialization();
    printf("Preprared for initialization \n");

    initialize_parameters();
    printf("Parameters are initialized \n");

    initialize_multiPhaseProblem();
    printf("Multiphase problem is initialized \n");
}

/*
@brief: Go through all the contacts along the pre-defined schedule.
        Assign a hybrid event with time whenever contact change is detected
*/
template <typename T>
void MHPCProblem<T>::prepare_initialization()
{    
    /* Load pinocchio model */
    const std::string urdf_filename = "../urdf/mini_cheetah_simple_correctedInertia.urdf";    
    buildPinModelFromURDF(urdf_filename, pin_model);

    /* Initialize the reference */
    quad_reference->initialize(plan_dur_all);

    /* Initialize the wb_reference interface */
    wb_reference.set_quadruped_reference(quad_reference);

    /* Initialize the srb refernce interface */
    srb_reference.set_quadruped_reference(quad_reference);

    wbm_ptr.clear();
    srbm_ptr.clear();
    footStepPlanner.clear();
    /* Build whole-body and srb models for each thread */
    for (int i = 0; i < pconfig->num_threads; i++)
    {
        // Whole-body model
        wbm_ptr.push_back(std::make_shared<WBM::Model<T>>(pin_model, pconfig->BG_alpha));

        // SRB model
        srbm_ptr.push_back(std::make_shared<SRBM::Model<T>>());

        // Create the foot step planner
        // It doesn't matter to use which wbm for the footstep planner
        footStepPlanner.push_back(std::make_shared<MHPCFootStep<T>>(wbm_ptr[i].get(),quad_reference.get())); 

        srbm_ptr[i]->set_footStepPlanner(footStepPlanner[i].get());        
    }             

    // Reset model
    reset_ptr = std::make_shared<MHPCReset<T>>(wbm_ptr[0].get(), footStepPlanner[0].get());    

    /* Build cost model */
    wb_track_cost = std::make_shared<WBTrackingCost<T>>();
    wb_track_cost->set_reference(&wb_reference);

    wb_foot_reg = std::make_shared<WBFootPlaceReg<T>>(quad_reference, wbm_ptr[0]);

    swing_pos_tracking = std::make_shared<SwingFootPosTracking<T>>(quad_reference, wbm_ptr[0]);

    swing_vel_tracking = std::make_shared<SwingFootVelTracking<T>>(quad_reference, wbm_ptr[0]);

    srb_track_cost = std::make_shared<SRBTrackingCost<T>>();
    srb_track_cost->set_reference(&srb_reference);

    if (approx_leq_scalar(plan_dur_all, .0))
    {
        printf("Warning: total plannining horizon cannot be zero \n");
        return;
    }

    /* If WB planning horizon is > 0, initialize the WB plan */
    if (pconfig->plan_dur_wb > 1e-5)
    {
        float wb_phase_start_time(0.0);
        float wb_phase_end_time(0.0);
        int wb_phase_horizon(0);
        int n_wb_phases(0);
        float t(0);

        VecM<int, 4> contact_prev, contact_cur;
        VecM<double, 4> contact_duration;
        contact_prev.setZero();
        contact_cur.setZero();
        contact_duration.setZero();

        quad_reference->get_contact_at_t(contact_prev, t);
        quad_reference->get_contact_duration_at_t(contact_duration, t);

        /* Iterate through all points along the reference */
        while (approx_leq_scalar(t, pconfig->plan_dur_wb))
        {
            // Get a contact at current time
            quad_reference->get_contact_at_t(contact_cur, t);

            // Determine whether to construct one phase: change of contact or reach to the WB planning horizon
            if ((contact_cur.cwiseNotEqual(contact_prev)).any() ||
                approx_eq_scalar(t, pconfig->plan_dur_wb))
            {
                wb_phase_end_time = t;
                n_wb_phases++;
                wb_phase_horizon = (int)round((wb_phase_end_time - wb_phase_start_time) / pconfig->dt_wb);

                pdata->wb_phase_start_times.push_back(wb_phase_start_time);
                pdata->wb_phase_end_times.push_back(wb_phase_end_time);
                pdata->wb_phase_horizons.push_back(wb_phase_horizon);
                pdata->wb_phase_contacts.push_back(contact_prev);
                pdata->wb_contact_durations.push_back(contact_duration);

                // Check whether a phase reach to a switching surface (contact change)
                bool wb_phase_reach_to_end = (contact_prev.cwiseNotEqual(contact_prev)).any();
                pdata->wb_is_phase_reach_end.push_back(wb_phase_reach_to_end);

                contact_prev = contact_cur;
                quad_reference->get_contact_duration_at_t(contact_duration, t);
                wb_phase_start_time = wb_phase_end_time;
            }            
            t += pconfig->dt_wb;
        }
        pdata->n_wb_phases = n_wb_phases;
    }
    /* If SRB planning horizon > 0, initialize and append the srb plan */    
    if (pconfig->plan_dur_srb > 1e-5)
    {
        pdata->srb_start_time = pconfig->plan_dur_wb;
        pdata->srb_end_time = plan_dur_all;
        pdata->srb_phase_horizon = (int) round(pconfig->plan_dur_srb/pconfig->dt_srb);
        pdata->n_srb_phase = 1;
    }
}

template <typename T>
void MHPCProblem<T>::initialize_parameters()
{
    /* Initialize REB and AL parameters for ineq and eq constraints */    
    const std::string& constraint_setting_fname= "../"+pconfig->constraintParamFileName;
    load_reb_params(grf_reb_param, constraint_setting_fname, "GRF");
    load_reb_params(torque_reb_param, constraint_setting_fname, "Torque");
    load_reb_params(jointspeed_reb_param, constraint_setting_fname, "JointSpeed");
    load_reb_params(joint_reb_param, constraint_setting_fname, "Joint");
    load_reb_params(minheight_reb_param, constraint_setting_fname, "MinHeight");
    load_al_params(td_al_param, constraint_setting_fname, "TD");    

    if (!pconfig->costFileName.empty())
    {
        const std::string& cost_weights_fname= "../"+pconfig->costFileName;
        loadCostWeights(cost_weights_fname, 
                        wb_track_cost,
                        wb_foot_reg,
                        swing_pos_tracking,
                        swing_vel_tracking,
                        srb_track_cost);    
    }
        
}

template <typename T>
void MHPCProblem<T>::initialize_multiPhaseProblem()
{       
    /* Create multi-phase WB planning problem */
    for (int i = 0; i < pdata->n_wb_phases; i++)
    {
        shared_ptr<WBPhase_T> phase;
        phase = make_shared<WBPhase_T>(pconfig->num_threads);

        shared_ptr<Trajectory<T, WBM::xs, WBM::us, WBM::ys>> traj;
        traj = make_shared<Trajectory<T, WBM::xs, WBM::us, WBM::ys>>(pconfig->dt_wb, pdata->wb_phase_horizons[i]);

        /* Initialize the state trajectory using the WB reference */
        VecM<double, WBM::xs> xr_k;
        xr_k.setZero();
        for (int k(0); k <= pdata->wb_phase_horizons[i]; k++)
        {            
            wb_reference.get_reference_at_t(xr_k, pdata->wb_phase_start_times[i] + k * pconfig->dt_wb);
            traj->X.at(k) = xr_k.cast<T>();
            traj->Xbar.at(k) = xr_k.cast<T>();
        }

        // Add trajectory to the phase (phase_idx)
        phase->set_trajectory(traj);

        create_problem_one_phase(phase, i);

        update_resetmap(phase, i);

        add_tconstr_one_phase(phase, i);

        phase->set_time_offset(pdata->wb_phase_start_times[i] - pdata->wb_phase_start_times[0]);

        phase->initialization();

        // Configure the set of shooting state
        phase->update_SS_config(pdata->wb_phase_horizons[i]+1);

        pdata->wb_trajs.push_back(traj);
        pdata->wb_phases.push_back(phase);        
    }
    
    /* Create single-phase SRB planning problem */
    if (pdata->srb_phase_horizon > 0)
    {
        shared_ptr<SRBPhase_T> phase;
        phase = make_shared<SRBPhase_T>(pconfig->num_threads);

        shared_ptr<Trajectory<T, SRBM::xs, SRBM::us, SRBM::ys>> traj;
        traj = make_shared<Trajectory<T, SRBM::xs, SRBM::us, SRBM::ys>>(pconfig->dt_srb, pdata->srb_phase_horizon);

        /* Initialize the state trajectory using the SRB reference */
        VecM<double, SRBM::xs> xr_k;
        xr_k.setZero();
        for (int k(0); k <= pdata->srb_phase_horizon; k++)
        {            
            srb_reference.get_reference_at_t(xr_k, pdata->srb_start_time + k * pconfig->dt_srb);            
            traj->Xbar.at(k) = xr_k.cast<T>();
        }

        // Add trajectory to the phase (phase_idx)
        phase->set_trajectory(traj);

        create_problem_one_phase(phase);

        phase->set_time_offset(pdata->srb_start_time);

        phase->initialization();

         // Configure the set of shooting state
        phase->update_SS_config(pdata->srb_phase_horizon+1);

        pdata->srb_phase = phase;
        pdata->srb_traj = traj;
    }    
    
}

template<typename T>
void MHPCProblem<T>::update()
{
    /* Shift the reference forward by a mpc step */ 
    quad_reference->step(pconfig->dt_mpc);

    /* Update the WB planning problem if necessary */
    if (pconfig->plan_dur_wb > 0)
    {
        printf("updating WB plan\n");
        update_WB_plan();
    }

    /* Update the SRB planning problem if necessary */
    if (pconfig->plan_dur_srb > 0)
    {
        printf("updating SRB plan\n");
        update_SRB_plan();
    }       
   
}

template <typename T>
void MHPCProblem<T>::update_WB_plan()
{        
    int nsteps = (int) round(pconfig->dt_mpc / pconfig->dt_wb);    
    float new_start_time = quad_reference->get_start_time();

    /* update the front end of WB phases */
    for (size_t j = 0; j < nsteps; j++)
    {                
        float first_phase_start_time = pdata->wb_phase_start_times.front() + pconfig->dt_wb;        
        
        // check whether the first phase shrinks to a point
        if (approx_eq_scalar(pdata->wb_phase_end_times.front(), first_phase_start_time))
        {
            // If yes, remove the front phase
            pdata->pop_front_phase();
        }
        else
        {
            // else, pop_front one time step from the front phase and update phase horizon
            pdata->wb_phases.front()->pop_front();
            pdata->wb_phase_horizons.front()--;
            pdata->wb_phase_start_times.front() = first_phase_start_time;
        }
    }

    /* update the back end of WB phases*/        
    for (size_t j = 0; j < nsteps; j++)    
    {   
        VecM<int, 4> new_contact;
        float new_end_time = pdata->wb_phase_end_times.back() + pconfig->dt_wb;                                        
        quad_reference->get_contact_at_t(new_contact, new_end_time - new_start_time);        
        bool contact_change = (new_contact.cwiseNotEqual(pdata->wb_phase_contacts.back())).any();        
        
        // If there is contact change at next timestep, and the last phase  of current paln reaches to the end, grow the current problem by a new phase
        if (contact_change && pdata->wb_is_phase_reach_end.back())
        {                        
            float new_phase_start_time = pdata->wb_phase_end_times.back();
            float new_phase_end_time = new_end_time;
            int new_phase_horizon = 1;

            VecM<double, 4> new_contact_duration;                        
            quad_reference->get_contact_duration_at_t(new_contact_duration, new_end_time - new_start_time);
            pdata->wb_phase_start_times.push_back(new_phase_start_time);
            pdata->wb_phase_end_times.push_back(new_phase_end_time);            
            pdata->wb_phase_horizons.push_back(new_phase_horizon);
            pdata->wb_is_phase_reach_end.push_back(false);
            pdata->wb_phase_contacts.push_back(new_contact);
            pdata->wb_contact_durations.push_back(new_contact_duration);
            pdata->n_wb_phases++;

            shared_ptr<Trajectory<T, WBM::xs, WBM::us, WBM::ys>> traj_to_add;
            traj_to_add = make_shared<Trajectory<T, WBM::xs, WBM::us, WBM::ys>>(pconfig->dt_wb, new_phase_horizon);         

            shared_ptr<WBPhase_T> phase_to_add;
            phase_to_add = make_shared<WBPhase_T>(pconfig->num_threads);

            phase_to_add->set_trajectory(traj_to_add);

            create_problem_one_phase(phase_to_add, pdata->n_wb_phases - 1);            

            phase_to_add->initialization();

            pdata->wb_trajs.push_back(traj_to_add);

            pdata->wb_phases.push_back(phase_to_add);            
            
        }
        else
        // Else, grow the last phase by one time step
        {
            pdata->wb_phase_end_times.back() = new_end_time;
            pdata->wb_phase_horizons.back()++;            
            if (contact_change)
            {
                pdata->wb_is_phase_reach_end.back() = true;
                add_tconstr_one_phase(pdata->wb_phases.back(), pdata->n_wb_phases - 1);
            }           
            pdata->wb_phases.back()->push_back_default();            
        }
                
    }

    /* Other updates */
    for (int i = 0; i < pdata->n_wb_phases; i++)
    {
        update_resetmap(pdata->wb_phases[i], i);
        
        /* Update the time offset of each phase */
        pdata->wb_phases[i]->set_time_offset(pdata->wb_phase_start_times[i] - pdata->wb_phase_start_times[0]);

        pdata->wb_phases[i]->reset_params();          

         if (i < pdata->n_wb_phases - 1 || (i == pdata->n_wb_phases - 1 && pdata->wb_phase_horizons[i] > nsteps))
        {
            pdata->wb_phases[i]->update_SS_config(pdata->wb_phase_horizons[i]+1);
        }

    }    
}

template <typename T>
void MHPCProblem<T>::update_SRB_plan()
{
    float new_start_time = quad_reference->get_start_time() + pconfig->plan_dur_wb;
    float new_end_time = new_start_time + pconfig->plan_dur_srb;
    int nsteps = (int) floor(pconfig->dt_mpc / pconfig->dt_srb + 1e-6);

    // If dt_mpc >= dt_srb, grow the SRB_phase by nsteps
    for (size_t j = 0; j < nsteps; j++)
    {
        // pop_front one time step from the front phase and update phase horizon
        pdata->srb_phase->pop_front();
        pdata->srb_phase_horizon--;        

        pdata->srb_phase_horizon++;            
        pdata->srb_phase->push_back_default();
    }
    
    // Otherwise, keep the data, and only update the times
    pdata->srb_start_time = new_start_time;
    pdata->srb_end_time = new_end_time;
    pdata->srb_phase->reset_params();    
    
}

/*  
    @brief: Create a single-phase WB planning problem 
*/
template <typename T>
void MHPCProblem<T>::create_problem_one_phase(shared_ptr<WBPhase_T>phase, int idx)
{
    /* specialize dynamics and resetmap  */
    const auto &phase_contact = pdata->wb_phase_contacts[idx];

    auto dynamics_callback = bind(&WBM::Model<T>::dynamics, wbm_ptr[0],
                                  pc::_1, pc::_2, pc::_3, pc::_4, pc::_5,
                                  phase_contact, pconfig->dt_wb);

    vector<function<void(WBStateMap&, WBContrlMap&, WBOutputMap&, 
                         WBDirectMap&, WBState&, WBContrl&, T)>> dynamics_partial_callback(pconfig->num_threads);

    for (int i = 0; i < pconfig->num_threads; i++)
    {        
        dynamics_partial_callback[i] =
            bind(&WBM::Model<T>::dynamics_partial, wbm_ptr[i],
             pc::_1, pc::_2, pc::_3, pc::_4, pc::_5, pc:: _6, pc::_7,
             phase_contact, pconfig->dt_wb);
    }
        

    /* set dynamics */
    phase->set_dynamics(dynamics_callback);
    phase->set_dynamics_partial(dynamics_partial_callback);

    /* set tracking cost */    
    phase->add_cost(wb_track_cost);

    /*Set foot regularization */
    phase->add_cost(wb_foot_reg);

    /* Swing foot posiiton tracking */   
    phase->add_cost(swing_pos_tracking);

    /* Swing foot velocity tracking */
    phase->add_cost(swing_vel_tracking);

    /* Torque limit */    
    shared_ptr<MHPCConstraints::TorqueLimit<T>> torqueLimit;
    torqueLimit = std::make_shared<MHPCConstraints::TorqueLimit<T>>();
    torqueLimit->update_horizon_len(pdata->wb_phase_horizons[idx]);
    torqueLimit->create_data();
    torqueLimit->initialize_params(torque_reb_param);
    phase->add_pathConstraint(torqueLimit);

    /* Joint speed limit */
    shared_ptr<MHPCConstraints::JointSpeedLimit<T>> jointSpeedLimit;
    jointSpeedLimit = std::make_shared<MHPCConstraints::JointSpeedLimit<T>>();
    jointSpeedLimit->update_horizon_len(pdata->wb_phase_horizons[idx]);
    jointSpeedLimit->create_data();
    jointSpeedLimit->initialize_params(jointspeed_reb_param);
    phase->add_pathConstraint(jointSpeedLimit);

    /* Joint limit */
    shared_ptr<MHPCConstraints::JointLimit<T>> jointLimit;
    jointLimit = std::make_shared<MHPCConstraints::JointLimit<T>>();
    jointLimit->update_horizon_len(pdata->wb_phase_horizons[idx]);
    jointLimit->create_data();
    jointLimit->initialize_params(joint_reb_param);
    phase->add_pathConstraint(jointLimit);

    /* Minimum Height constraint */
    shared_ptr<MHPCConstraints::WBMinimumHeight<T>> wbMinHeightConstraint;
    wbMinHeightConstraint = std::make_shared<MHPCConstraints::WBMinimumHeight<T>>();
    wbMinHeightConstraint->update_horizon_len(pdata->wb_phase_horizons[idx]);
    wbMinHeightConstraint->create_data();
    wbMinHeightConstraint->initialize_params(minheight_reb_param);
    phase->add_pathConstraint(wbMinHeightConstraint);

    /* Set GRF constraints if any*/
    if (phase_contact.cwiseEqual(1).any())
    {
        shared_ptr<MHPCConstraints::WBGRF<T>> grfConstraint;
        grfConstraint = std::make_shared<MHPCConstraints::WBGRF<T>>(phase_contact);
        grfConstraint->update_horizon_len(pdata->wb_phase_horizons[idx]);
        grfConstraint->create_data();
        grfConstraint->initialize_params(grf_reb_param);
        phase->add_pathConstraint(grfConstraint);
    }    
}


/*  
    @brief: Create a single-phase SRB planning problem
*/
template <typename T>
void MHPCProblem<T>::create_problem_one_phase(shared_ptr<SRBPhase_T> phase)
{
    
    auto dynamics_callback = bind(&SRBM::Model<T>::dynamics, srbm_ptr[0],
                                  pc::_1, pc::_2, pc::_3, pc::_4, pc::_5,
                                  pconfig->dt_srb);

    vector<function<void(SRBMStateMap&, SRBMContrlMap&, SRBMOutputMap&, 
                         SRBMDirectMap&, SRBMState&, SRBMContrl&, T)>> dynamics_partial_callback;

    for (int i = 0; i < pconfig->num_threads; i++)
    {        
        dynamics_partial_callback.push_back(
            bind(&SRBM::Model<T>::dynamics_partial, srbm_ptr[i],
             pc::_1, pc::_2, pc::_3, pc::_4, pc::_5, pc:: _6, pc::_7,
             pconfig->dt_srb));
    }
        
    phase->set_dynamics(dynamics_callback);
    phase->set_dynamics_partial(dynamics_partial_callback);

    /* set tracking cost */    
    phase->add_cost(srb_track_cost);   

     /* Minimum Height constraint */
    shared_ptr<MHPCConstraints::SRBMMinimumHeight<T>> srbMinHeightConstraint;
    srbMinHeightConstraint = std::make_shared<MHPCConstraints::SRBMMinimumHeight<T>>();
    srbMinHeightConstraint->update_horizon_len(pdata->srb_phase_horizon);
    srbMinHeightConstraint->create_data();
    srbMinHeightConstraint->initialize_params(minheight_reb_param);
    phase->add_pathConstraint(srbMinHeightConstraint);
        
}

template <typename T>
void MHPCProblem<T>::update_resetmap(shared_ptr<WBPhase_T>phase, int idx)
{
    /* Determine the touchdown status */
    const VecM<int, 4> &phase_contact_cur = pdata->wb_phase_contacts[idx];    
    VecM<int, 4> phase_contact_next;    

    if (idx < pdata->n_wb_phases - 1) // If it is an intermediate phase
    {
        phase_contact_next = pdata->wb_phase_contacts[idx + 1];
    }
    else // else the last phase
    {
        quad_reference->get_contact_at_t(phase_contact_next, pconfig->plan_dur_wb + pconfig->dt_mpc);
    }

    ModelType mtype_current, mtype_next;
    mtype_current = ModelType::WB;
    mtype_next = ModelType::WB;

    // If this is the last whole-body phase
    if (idx == pdata->n_wb_phases-1 && pdata->n_srb_phase > 0) 
    {
        mtype_next = ModelType::SRB;
    }
    
    /* set resetmap */
    auto resetmap_callback = bind(&MHPCReset<T>::reset_map, reset_ptr,
                                  pc::_1, pc::_2, phase_contact_cur, phase_contact_next, mtype_current, mtype_next);
    auto resetmap_partial_callback = bind(&MHPCReset<T>::reset_map_partial, reset_ptr,
                                          pc::_1, pc::_2, phase_contact_cur, phase_contact_next, mtype_current, mtype_next);

    phase->set_resetmap(resetmap_callback);
    phase->set_resetmap_partial(resetmap_partial_callback);                                          
}

/*  
    @brief: Add terminal constraints and set reset map for a single phase of WB planning
*/
template <typename T>
void MHPCProblem<T>::add_tconstr_one_phase(shared_ptr<WBPhase_T>phase, int idx)
{
    /* Determine the touchdown status */
    const VecM<int, 4> &phase_contact_cur = pdata->wb_phase_contacts[idx];
    VecM<int, 4> touchdown_status;
    VecM<int, 4> phase_contact_next;
    touchdown_status.setZero();

    if (idx < pdata->n_wb_phases - 1) // If it is an intermediate phase
    {
        phase_contact_next = pdata->wb_phase_contacts[idx + 1];
    }
    else // else the last phase
    {
        quad_reference->get_contact_at_t(phase_contact_next, pconfig->plan_dur_wb + pconfig->dt_mpc);
    }

    for (int leg = 0; leg < 4; leg++)
    {
        if (phase_contact_cur[leg] == 0 && phase_contact_next[leg] == 1)
        {
            touchdown_status[leg] = 1;
        }
    }    

    /* Set touchdown constraint and cost if any touchdown happens */
    if (find_eigen(touchdown_status, 1).size() > 0)
    {
        shared_ptr<MHPCConstraints::WBTouchDown<T>> tdConstraint;
        tdConstraint = std::make_shared<MHPCConstraints::WBTouchDown<T>>(touchdown_status, wbm_ptr[0].get());
        tdConstraint->create_data();
        tdConstraint->initialize_params(td_al_param);
        phase->add_terminalConstraint(tdConstraint);

        shared_ptr<TDVelocityPenalty<T>> tdVelPenalty;
        tdVelPenalty = std::make_shared<TDVelocityPenalty<T>>(touchdown_status, wbm_ptr[0]);
        phase->add_cost(tdVelPenalty);               
    }
}

template <typename T>
void MHPCProblem<T>::pretty_print()
{
    /* Print MHPC Problem Configuration */
    printf("\n");
    printf("**************MHPC Problem Config*******************\n");
    tabulate::Table config_print;
    config_print.add_row({"plan_dur", "wb_plan", "srb_plan", "dt_wb", "dt_srb", "dt_mpc"});
    config_print.add_row({std::to_string(plan_dur_all),
                          std::to_string(pconfig->plan_dur_wb),
                          std::to_string(pconfig->plan_dur_srb),
                          std::to_string(pconfig->dt_wb),
                          std::to_string(pconfig->dt_srb),
                          std::to_string(pconfig->dt_mpc)});

    config_print.column(1).format() 
        .font_align(tabulate::FontAlign::center);

    for (size_t i = 0; i < config_print[0].size(); i++)
    {
        config_print[0][i].format()
        .font_color(tabulate::Color::yellow)
        .font_align(tabulate::FontAlign::center)
        .font_style({tabulate::FontStyle::bold});
    }

    std::cout << config_print << "\n";

    /* Print WB plan details  */
    printf("************Whole-Body plan details***************\n");
    if (pdata->n_wb_phases > 0)
    {
        tabulate::Table wbPlan_print;
        wbPlan_print.add_row({"Phase Index", "Horizon", "Start Time", "End Time", "Contact", "Contact Duration"});    
        for (int i = 0; i < pdata->n_wb_phases; i++)
        {
            wbPlan_print.add_row({std::to_string(i), 
                            std::to_string(pdata->wb_phase_horizons[i]),
                            std::to_string(pdata->wb_phase_start_times[i]), 
                            std::to_string(pdata->wb_phase_end_times[i]),
                            eigenToString(pdata->wb_phase_contacts[i].transpose()),
                            eigenToString(pdata->wb_contact_durations[i].transpose())});                            
        }    
        // center-align and color header cells
        for (size_t i = 0; i < wbPlan_print[0].size(); ++i) {
            wbPlan_print[0][i].format()
            .font_color(tabulate::Color::yellow)
            .font_align(tabulate::FontAlign::center)
            .font_style({tabulate::FontStyle::bold});

            wbPlan_print.column(i).format()
                .font_align(tabulate::FontAlign::center);       
        }    
        std::cout << wbPlan_print << "\n";
    }else{
        printf("Whole-body plan is empty \n");
    }
    printf("\n");

    /* Print SRB plan details */
    printf("********Single-Rigeid-Body plan details*********\n");
    if (pdata->n_srb_phase > 0)
    {
        tabulate::Table srbPlan_print;
        srbPlan_print.add_row({"Phase Index", "Horizon", "Start Time", "End Time"});    
        srbPlan_print.add_row({std::to_string(pdata->n_srb_phase-1), 
                            std::to_string(pdata->srb_phase_horizon),
                            std::to_string(pdata->srb_start_time), 
                            std::to_string(pdata->srb_end_time)});   
        // center-align and color header cells
        for (size_t i = 0; i < srbPlan_print[0].size(); ++i) {
            srbPlan_print[0][i].format()
            .font_color(tabulate::Color::yellow)
            .font_align(tabulate::FontAlign::center)
            .font_style({tabulate::FontStyle::bold});

            srbPlan_print.column(i).format()
                .font_align(tabulate::FontAlign::center);       
        }    
        std::cout << srbPlan_print << "\n";
    }else{
        printf("SRB plan is empty \n");
    }    
    printf("\n");
}

template class MHPCProblem<double>;

