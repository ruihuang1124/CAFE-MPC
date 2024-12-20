#include <algorithm>
#include <iomanip>
#include "MultiPhaseDDP.h"
#include "HSDDP_CompoundTypes.h"
#include "HSDDP_Utils.h"

/*
    @brief: perform linear rollout to compute search direction for shooting state
            use before running hybrid rollout
*/
template <typename T>
void MultiPhaseDDP<T>::linear_rollout(T eps, HSDDP_OPTION &option)
{
    DVec<T> dx_init, dx_end;
    DVec<T> xend, xsim_end;
    DMat<T> Px; // resetmap partial
    dx0.setZero();
    dx_init = dx0;
    dV_1 = 0;
    dV_2 = 0;

    for (int i = 0; i < n_phases; i++)
    {
        T dV_1_i(0), dV_2_i(0);
        if (i > 0)
        {
            phases[i - 1]->get_terminal_state_dx(dx_end);
            phases[i - 1]->get_terminal_state(xend, xsim_end);
            phases[i - 1]->resetmap_partial(Px, xend);
            dx_init = Px * dx_end;
        }

        phases[i]->set_initial_condition_dx(dx_init);

        phases[i]->linear_rollout(eps, option);

        phases[i]->get_exp_cost_change(dV_1_i, dV_2_i);

        dV_1 += dV_1_i;
        dV_2 += dV_2_i;
    }
}

/*
    @brief: perform hybrid rollout
            nominal system rollout is obtained with eps = 0
*/
template <typename T>
bool MultiPhaseDDP<T>::hybrid_rollout(T eps, HSDDP_OPTION &option)
{
    actual_cost = 0;
    max_pconstr = 0;
    max_tconstr = 0;
    DVec<T> xinit = x0; // initial condition for each phase
    DVec<T> xsim_init = x0;
    DVec<T> xend;
    DVec<T> xsim_end;
    bool is_last_phase = false;
    bool success = true;

    run_before_forward_sweep(); // currently not used

    for (int i = 0; i < n_phases; i++)
    {
        if (!option.MS)
        {
            phases[i]->update_SS_config(0);
        }
        
        if (i > 0)
        {
            // If not the first phase, run resetmap at the end of previous phase
            phases[i - 1]->get_terminal_state(xend, xsim_end);
            xinit = phases[i - 1]->resetmap(xend);
            xsim_init = phases[i - 1]->resetmap(xsim_end);
        }

        phases[i]->set_initial_condition(xinit, xsim_init);    // Set initial condition of current phase
        
        if (!phases[i]->hybrid_rollout(eps, option, is_last_phase)) // run hybrid rollout for each phase)
        {            
            success = false;
            break;
        }
        

        // update the maximum constraint violations
        max_pconstr = std::min(max_pconstr, phases[i]->get_max_pconstrs()); // should have non-positive value
        max_tconstr = std::max(max_tconstr, phases[i]->get_max_tconstrs()); // should have non-negative value
    }
    return success;
}

template <typename T>
std::pair<bool, int> MultiPhaseDDP<T>::line_search(HSDDP_OPTION &option)
{
    T exp_cost_change = 0;
    T exp_merit_change = 0;
    T eps = 1;
    T cost_prev = actual_cost;
    T merit_prev = merit;
    T feas_prev = feas;

    bool success = false;
    bool rollout_success = true;
    int iter = 0;

    while (eps > 1e-3)
    {        
        iter ++;
        rollout_success = hybrid_rollout(eps, option);
        compute_cost(option);
        feas = measure_dynamics_feasibility();
        merit = actual_cost + merit_rho * feas;

        exp_cost_change = eps * dV_1 + 0.5 * eps * eps * dV_2;
        exp_merit_change = exp_cost_change - eps * merit_rho* feas_prev;

#ifdef DEBUG_MODE
        printf("\t eps=%.3e \t actual cost change =%.3e \t expeced cost change =%.3e \t feas =%.3e\n",
               eps, actual_cost - cost_prev, option.gamma * exp_cost_change, feas);
#endif
        if ((merit <= merit_prev + option.gamma * exp_merit_change) && rollout_success)
        {
            success = true;
            break;
        }        
            
        eps *= option.alpha;

    }
    return std::make_pair<bool, int>(std::move(success), std::move(iter));
}

template <typename T>
std::pair<bool, int> MultiPhaseDDP<T>::backward_sweep_regularized(T &regularization, HSDDP_OPTION &option)
{
    bool success = false;
    int iter = 0;

    while (!success)
    {
#ifdef DEBUG_MODE
        printf("\t regularization = %f\n", regularization);
#endif
        iter++;
        success = backward_sweep(regularization);
        if (success)
            break;

        regularization = std::max(regularization * option.update_regularization, T(1e-03));        

        if (regularization > 1e2)
        {
            printf("Regularization term exceeds maximum value! \n");
            printf("Optimization terminates! \n");
            break;
        }
    }
    
    regularization = regularization / 20;
    if (regularization < 1e-06)
        regularization = 0;
    return std::make_pair<bool, int>(std::move(success), std::move(iter));
}

/*
  @brief: perform backward sweep for the multi-phase problem
  @params:
          regularization: regularization parameter
  @return: success of backward sweep
*/
template <typename T>
bool MultiPhaseDDP<T>::backward_sweep(T regularization)
{
    bool success = true;
    // T dVprime = 0;
    DVec<T> Gprime;
    DMat<T> Hprime;
    DVec<T> xend;
    DMat<T> Px;
    size_t xs(0), xs_next(0);   
    dV_1 = 0;
    dV_2 = 0;
    T dV_1_i(0), dV_2_i(0);

    for (int i = n_phases - 1; i >= 0; i--)
    {
        xs = phases[i]->get_state_dim();
        xs_next = (i < n_phases - 1) ? phases[i + 1]->get_state_dim() : xs;
        phases[i]->get_terminal_state(xend);
        Gprime.setZero(xs_next);
        Hprime.setZero(xs_next, xs_next);
        Px.setZero(xs, xs_next);
        // If not the last phase, udpate approximated value function back propagated from next phase
        if (i <= n_phases - 2)
        {
            phases[i]->resetmap_partial(Px, xend);
            phases[i + 1]->get_value_approx(Gprime, Hprime); // get the value infomation at the begining of phase i+1
            impact_aware_step(Gprime, Hprime, Px);
        }

        // If backward sweep of current phase fails, break and return false
        success = phases[i]->backward_sweep(regularization, Gprime, Hprime);
        if (!success)
            return success;       
        
        phases[i]->get_exp_cost_change(dV_1_i, dV_2_i);
        dV_1 += dV_1_i;
        dV_2 += dV_2_i;
    }
    return success;
}

template <typename T>
void MultiPhaseDDP<T>::solve(HSDDP_OPTION& option, const float& max_cputime)
{
    iter_ = 0;
    ls_iter_total_ = 0;
    int iter_ou = 0;
    int iter_in = 0;

    T cost_prev(0), merit_prev(0);
    bool success = true; // currently defined only for backward sweep
    int print_precision = 8;

    /* clear all buffer */    
    cost_buffer.clear();
    dyn_feas_buffer.clear();
    eqn_feas_buffer.clear();
    ineq_feas_buffer.clear();
    
    bool max_cputime_reached = false;
    auto solve_start = high_resolution_clock::now();
    auto check_point = high_resolution_clock::now();
    auto solve_time_elapse = duration_ms(check_point - solve_start);

    hybrid_rollout(0, option);
    update_nominal_trajectory();
    compute_cost(option);
    feas = measure_dynamics_feasibility();
    

    printf("Initial total cost = %f \n", actual_cost);
    printf("Initial dynamics feasibility = %f \n", feas);
    if (option.AL_active)
    {
        printf("Initial terminal constraint = %f \n", max_tconstr);
    }
    if (option.ReB_active)
    {
        printf("Initial path constraint = %f \n", max_pconstr);
    }
    
    

    /* buffer the initial information */    
    cost_buffer.push_back(actual_cost);
    dyn_feas_buffer.push_back(feas);
    eqn_feas_buffer.push_back(max_tconstr);
    ineq_feas_buffer.push_back(max_pconstr);

    T regularization = 0;
    while (iter_ou < option.max_AL_iter && 
           !max_cputime_reached)
    {
        iter_ou++;
        max_tconstr_prev = max_tconstr;
        max_pconstr_prev = max_pconstr;

#ifdef DEBUG_MODE
        printf("outer loop iteration %d \n", iter_ou);
#endif
                       
        regularization = 0;
        iter_in = 0;
        while (iter_in < option.max_DDP_iter && 
           !max_cputime_reached)
        {
            compute_cost(option);
            feas = measure_dynamics_feasibility();            

            // printf("total cost = %f, dynamics infeasibility = %f \n", actual_cost, feas);
            iter_in++;
            iter_++;

#ifdef DEBUG_MODE
            printf("\t inner loop iteration %d \n", iter_in);
#endif                  
      
            // Checkpoint 1
            check_point = high_resolution_clock::now();
            solve_time_elapse = duration_ms(check_point - solve_start);            
            printf("Checkpoint 1 time elapsed %f\n", solve_time_elapse.count());
            if (approx_geq_scalar(solve_time_elapse.count(), max_cputime))
            {
                max_cputime_reached = true;                
                break;
            }
            

            LQ_approximation(option);

             // Checkpoint 2
            check_point = high_resolution_clock::now();
            solve_time_elapse = duration_ms(check_point - solve_start);            
            printf("checkpoint 2 time elapsed %f\n", solve_time_elapse.count());
            if (approx_geq_scalar(solve_time_elapse.count(), max_cputime))
            {
                max_cputime_reached = true;                
                break;
            }

            int reg_iter(0);
            std::tie<bool, int>(success, reg_iter) = backward_sweep_regularized(regularization, option);
            reg_iter_total_ += reg_iter;
            if (!success)
            {
                goto bad_solve;
            }

            // Checkpoint 3
            check_point = high_resolution_clock::now();
            solve_time_elapse = duration_ms(check_point - solve_start);            
            printf("checkpoint 3 time elapsed %f\n", solve_time_elapse.count());
            if (approx_geq_scalar(solve_time_elapse.count(), max_cputime))
            {
                max_cputime_reached = true;                
                break;
            }

            if (option.MS)
            {
                linear_rollout(1.0, option);                
            }
            
            T dV_abs = fabs(dV_1 + 0.5*dV_2);
        
            merit_rho = (feas > option.dynamics_feas_thresh) ? dV_abs/((1-option.merit_scale)*feas) + option.merit_offset : 0; 

            merit = actual_cost + merit_rho * feas;
            cost_prev = actual_cost;
            merit_prev = merit;

            // Early terminates before running line search
            if ((dV_abs < option.cost_thresh) && (feas <= option.dynamics_feas_thresh))
            {                
                break;
            }            
          

            bool ls_success(false);
            int ls_iter(0);
            std::tie<bool, int>(ls_success, ls_iter) = line_search(option);
            ls_iter_total_ += ls_iter;
            if (ls_success)
            {
                // if line search succeeds, accept the step
                update_nominal_trajectory();
            }
            else
            {
                // else do not update
                actual_cost = cost_prev;
                merit = merit_prev;
            }

            // Later terminates if actual cost change very small
            if ((fabs((cost_prev - actual_cost)/cost_prev) < option.cost_thresh) && (feas <= option.dynamics_feas_thresh))
                break;         

            // Checkpoint 4
            check_point = high_resolution_clock::now();
            solve_time_elapse = duration_ms(check_point - solve_start);       
            printf("checkpoint 4 time elapsed %f\n", solve_time_elapse.count());     
            if (approx_geq_scalar(solve_time_elapse.count(), max_cputime))
            {
                max_cputime_reached = true;                
                break;
            }                           
     
            cost_buffer.push_back(actual_cost);
            dyn_feas_buffer.push_back(feas);
            eqn_feas_buffer.push_back(max_tconstr);
            ineq_feas_buffer.push_back(max_pconstr);            

        }       

#ifdef DEBUG_MODE
        printf("terminal constraint violation = %f \n", max_tconstr);
        printf("path constraint violation = %f \n", fabs(max_pconstr));
#endif

        if (max_tconstr < option.tconstr_thresh && fabs(max_pconstr) < option.pconstr_thresh && feas <= option.dynamics_feas_thresh)
        {
            printf("all constraints satisfied \n");
            printf("optimization finished \n");
            break;
        }
        if (fabs(max_tconstr - max_tconstr_prev) < 0.0001 && fabs(max_pconstr - max_pconstr_prev) < 0.0001 && feas <= option.dynamics_feas_thresh)
        {
            printf("constraints stop decreasing \n");
            printf("optimization terminates \n");
            break;
        }
        
        if (max_cputime_reached)
        {
            printf("maximum cputime reached \n");
            break;
        }
                           
        if (option.AL_active)
        {
            update_AL_params(option);
        }
        if (option.ReB_active)
        {
            update_REB_params(option);
        }

        if (iter_ou >= option.max_AL_iter)
        {
            printf("maximum iteration reached \n");
            break;
        }      
    }
    

    printf("total cost = %f \n", actual_cost);
    printf("dynamics infeasibility = %f \n", feas);
    printf("terminal constraint violation = %f \n", max_tconstr);    
    printf("path constraint violation = %f \n", max_pconstr);
    

bad_solve:
    if (!success)
    {
        printf(RED);
        printf("Failed to solve the optimization due to too large regularization \n");
        printf(RESET);
    }

    auto stop = high_resolution_clock::now();
    solve_time_elapse = duration_ms(stop - solve_start);         
    solve_time_ = solve_time_elapse.count();
}

template <typename T>
void MultiPhaseDDP<T>::compute_cost(const HSDDP_OPTION &option)
{
    actual_cost = 0;
    for (int i = 0; i < n_phases; i++)
    {
        phases[i]->compute_cost(option);
        actual_cost += phases[i]->get_actual_cost();
    }
}

template <typename T>
void MultiPhaseDDP<T>::LQ_approximation(HSDDP_OPTION &option)
{
    for (int i = 0; i < n_phases; i++)
    {
        phases[i]->LQ_approximation(option);
    }
}
/*
    @brief: Set a method to initialize the dynamics everytime berfore running forward sweep
    @Note:  In future version, this function need modified for general purpose. Currently it supports foothold
            initialization.
*/
template <typename T>
void MultiPhaseDDP<T>::set_dynamics_init_callback(function<void(DVec<T>)> dynamics_init_callback_)
{
    dynamics_init_callback = dynamics_init_callback_;
}

/*
    @brief: Always run before forward sweep to initialize the dyanmics model
*/
template <typename T>
void MultiPhaseDDP<T>::run_before_forward_sweep()
{
    if (nullptr != dynamics_init_callback)
    {
        dynamics_init_callback(x0);
    }
}

/*
  @brief: perform impact-aware backward DDP step
  @params:
          Px: reset map
          G, H: gradient and hessian of value function propaged from next phase through phase transition
*/

template <typename T>
void MultiPhaseDDP<T>::impact_aware_step(DVec<T> &G, DMat<T> &H, const DMat<T> &Px)
{
    G = Px.transpose() * G;
    H = Px.transpose() * H * Px;
}

template <typename T>
void MultiPhaseDDP<T>::update_REB_params(HSDDP_OPTION &option)
{
    for (auto &phase : phases)
    {
        phase->update_REB_params(option);
    }
}

template <typename T>
void MultiPhaseDDP<T>::update_AL_params(HSDDP_OPTION &option)
{
    for (auto &phase : phases)
    {
        phase->update_AL_params(option);
    }
}

template <typename T>
void MultiPhaseDDP<T>::update_nominal_trajectory()
{
    for (auto &phase : phases)
    {
        phase->update_nominal_trajectory();
    }
}

template <typename T>
T MultiPhaseDDP<T>::measure_dynamics_feasibility(int norm_id)
{
    T feasibility = 0;

    for (auto &phase : phases)
    {
        feasibility += phase->measure_dynamics_feasibility(norm_id);
    }

    if (norm_id == 2)
    {
        feasibility = sqrt(feasibility);
    }

    return feasibility;
}

template <typename T>
void MultiPhaseDDP<T>::get_solver_info( std::vector<float>& cost_out,
                                        std::vector<float>& dyn_feas_out,
                                        std::vector<float>& eqn_feas_out, 
                                        std::vector<float>& ineq_feas_out)
{
    cost_out = cost_buffer;
    dyn_feas_out = dyn_feas_buffer;
    eqn_feas_out = eqn_feas_buffer;
    ineq_feas_out = ineq_feas_buffer;
}

template class MultiPhaseDDP<double>;
