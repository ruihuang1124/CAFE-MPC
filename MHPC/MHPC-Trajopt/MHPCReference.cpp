#include "MHPCReference.h"
#include <cassert>

/* 
    @brief: Get a quadruped reference at time t, and converts to the WBM state reference
            WB state rerefence: [pos, eul, qJ, vWorld, eulrate, qJd]
    @params: t: time relative to the start of the reference
*/

void WBReference::get_reference_at_t(
    VecM<double, WBM::xs>& xt, VecM<double, WBM::us>& ut, VecM<double, WBM::ys>& yt, float t)
{            
    get_reference_at_t(xt, t);

    ut = quad_state_t_ptr->torque;       // torque      
    yt = quad_state_t_ptr->grf;
}

/* 
    @brief: Get a quadruped reference at time t, and converts to the WBM state reference
            WB state rerefence: [pos, eul, qJ, vWorld, eulrate, qJd]
    @params: t: time relative to the start of the reference
*/

void WBReference::get_reference_at_t(VecM<double, WBM::xs>& xt, float t)
{
    
#ifdef DEBUG_MHPC
    assert(("quad_ref_ cannot be nullptr", nullptr != quad_ref_));
#endif    

    quad_state_t_ptr = quad_ref_->get_a_reference_ptr_at_t(t);

#ifdef DEBUG_MHPC
    assert(("quad_state_ref_ptr is nullptr \n", nullptr != quad_state_t_ptr));
#endif

    xt.head<6>() = quad_state_t_ptr->body_state.head<6>();          //body state: pos, euler, vel, eulrate
    xt.segment<12>(6) = quad_state_t_ptr->qJ;                       //joint angle
    xt.segment<6>(18) = quad_state_t_ptr->body_state.tail<6>();
    xt.tail<12>() = quad_state_t_ptr->qJd;                          //joint velocity
}

/* 
    @brief: Get a quadruped reference state and control at time t, and converts to the SRBM state refrence
    @params: t: time relative to the start of the reference
*/

void SRBReference::get_reference_at_t(
    VecM<double, SRBM::xs>& xt, VecM<double, SRBM::us>& ut, VecM<double, SRBM::ys>& yt, float t)
{            
    (void) (yt);
    get_reference_at_t(xt, t);

    ut = quad_state_t_ptr->grf;       // GRF      
}

/* 
    @brief: Get a quadruped reference state at time t, and converts to the SRBM state refrence
    @params: t: time relative to the start of the reference
*/

void SRBReference::get_reference_at_t(VecM<double, SRBM::xs>& xt, float t)
{
#ifdef DEBUG_MHPC
    assert(("quad_ref_ cannot be null \n", nullptr != quad_ref));
#endif

quad_state_t_ptr = quad_ref_->get_a_reference_ptr_at_t(t);

#ifdef DEBUG_MHPC
    assert(("quad_state_ref_ptr cannot be nullptr \n", nullptr != quad_state_t_ptr));
#endif

    xt = quad_state_t_ptr->body_state;   //body state: pos, euler, vel, ang
}

