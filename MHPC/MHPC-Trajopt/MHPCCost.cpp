#include "MHPCCost.h"

template <typename T>
void WBFootPlaceReg<T>::running_cost(RCost &rcost, const State &x, const Contrl &u, const Output &y, T dt, float t)
{
    (void)(y);
    (void)(u);
    quad_astate = quad_reference->get_a_reference_ptr_at_t(t);

    pCoM = x.template head<3>();
    pCoM_des = quad_astate->body_state.head<3>();

    wbm_ptr->get_footPositions(pFoot, x);
    rcost.l = 0;
    for (size_t l = 0; l < 4; l++)
    {
        if (contact[l] > 0)
        {
            prel = pFoot[l] - pCoM;
            prel_des = quad_astate->foot_placements.segment<3>(l * 3) - pCoM_des;
            d_prel = prel - prel_des;
            T footReg = .5 * d_prel.transpose() * QFoot * d_prel;
            footReg *= dt;
            rcost.l += footReg;
        }
    }
}

template <typename T>
void WBFootPlaceReg<T>::running_cost_par(RCost &rcost, const State &x, const Contrl &u, const Output &y, T dt, float t)
{
    (void)(y);
    (void)(u);
    quad_astate = quad_reference->get_a_reference_ptr_at_t(t);

    pCoM = x.template head<3>();
    pCoM_des = quad_astate->body_state.head<3>();

    wbm_ptr->get_footPositions(pFoot, x);
    wbm_ptr->get_footJacobians(J_foot, x);

    rcost.lx.setZero();
    rcost.lxx.setZero();
    for (size_t l = 0; l < 4; l++)
    {
        if (contact[l] > 0)
        {
            prel = pFoot[l] - pCoM;
            prel_des = quad_astate->foot_placements.segment<3>(l * 3) - pCoM_des;
            d_prel = prel - prel_des;

            J_foot[l].template leftCols<3>().setZero();
            lq = (J_foot[l].transpose() * QFoot * d_prel) * dt;
            lqq = (J_foot[l].transpose() * QFoot * J_foot[l]) * dt;

            rcost.lx.template head<WBM::nq>() += lq;
            rcost.lxx.template topLeftCorner<WBM::nq, WBM::nq>() += lqq;
        }
    }
}

template <typename T>
void WBFootPlaceReg<T>::terminal_cost(TCost &tcost, const State &x, float tend)
{
    quad_astate = quad_reference->get_a_reference_ptr_at_t(tend);

    pCoM = x.template head<3>();
    pCoM_des = quad_astate->body_state.head<3>();

    wbm_ptr->get_footPositions(pFoot, x);
    tcost.Phi = 0;
    for (size_t l = 0; l < 4; l++)
    {
        if (contact[l] > 0)
        {
            prel = pFoot[l] - pCoM;
            prel_des = quad_astate->foot_placements.segment<3>(l * 3) - pCoM_des;
            d_prel = prel - prel_des;
            T footReg = .5 * d_prel.transpose() * QFoot * d_prel;
            tcost.Phi += 10*footReg;
        }
    }
}

template <typename T>
void WBFootPlaceReg<T>::terminal_cost_par(TCost &tcost, const State &x, float tend)
{
    quad_astate = quad_reference->get_a_reference_ptr_at_t(tend);

    pCoM = x.template head<3>();
    pCoM_des = quad_astate->body_state.head<3>();
    
    wbm_ptr->get_footPositions(pFoot, x);
    wbm_ptr->get_footJacobians(J_foot, x);

    tcost.Phix.setZero();
    tcost.Phixx.setZero();
    for (size_t l = 0; l < 4; l++)
    {
        if (contact[l] > 0)
        {
            prel = pFoot[l] - pCoM;
            prel_des = quad_astate->foot_placements.segment<3>(l * 3) - pCoM_des;
            d_prel = prel - prel_des;

            J_foot[l].template leftCols<3>().setZero();
            lq = J_foot[l].transpose() * QFoot * d_prel;
            lqq = J_foot[l].transpose() * QFoot * J_foot[l];

            tcost.Phix.template head<WBM::nq>() += 10*lq;
            tcost.Phixx.template topLeftCorner<WBM::nq, WBM::nq>() += 15*lqq;
        }
    }
}

template <typename T>
void SwingFootClearance<T>::running_cost(RCost&rcost, const State& x, const Contrl& u, const Output& y, T dt, float t)
{
    (void)(u);
    (void)(y);
    (void)(t);

    rcost.l = 0;
    for (int l = 0; l < 4; l++)
    {
        if (contact[l] == 0)
        {
            Vec3<T> pF_local;
            Vec3<T> ql = x.template segment<3>(6+3*l);
            computeLegPosition(pF_local, ql, l);
            T dh = pF_local[2] - des_height;
            rcost.l += 0.5 * weight * dh*dh *dt;
        }        
    }
    
}

template <typename T>
void SwingFootClearance<T>::running_cost_par(RCost&rcost, const State& x, const Contrl& u, const Output& y, T dt, float t)
{
    (void)(u);
    (void)(y);
    (void)(t);

    rcost.lx.setZero();
    rcost.lxx.setZero();
    for (int l = 0; l < 4; l++)
    {
        if (contact[l] == 0)
        {
            Vec3<T> pF_local;
            Mat3<T> J_local;
            Vec3<T> ql = x.template segment<3>(6+3*l);
            computeLegPosition(pF_local, ql, l);
            computeLegJacobian(J_local, ql, l);
            T dh = pF_local[2] - des_height;
            const auto& Jz = J_local.template bottomRows<1>();
            lq = dh * Jz.transpose() * weight *dt;
            lqq = Jz.transpose() * Jz * weight * dt;

            rcost.lx.template segment<3>(6+3*l) += lq;
            rcost.lx.template block<3,3>(6+3*l, 6+3*l) += lqq;
        }            
    }   
}

template <typename T>
void ImpulseTerminalCost<T>::terminal_cost(TCost& tcost, const State& x, float tend)
{
    (void) (tend);
    wbm_ptr->impact(xpostImp, impulse, x, contact_cur, contact_next);
    tcost.Phi = 0;
    for (int l = 0; l < 4; l++)
    {
        if (TDStatus(l) > 0)
        {
            const auto& impulse_normal_l = impulse(3*l + 2);
            tcost.Phi += 0.5 * weight * impulse_normal_l * impulse_normal_l;
        }        
    }
    
}

template <typename T>
void ImpulseTerminalCost<T>::terminal_cost_par(TCost& tcost, const State& x, float tend)
{
    (void) (tend);

    wbm_ptr->impact(xpostImp, impulse, x, contact_cur, contact_next);
    wbm_ptr->impact_partial(dxpost_dx, dImp_dx, x, contact_cur, contact_next);
    tcost.Phix.setZero();
    tcost.Phixx.setZero();
    for(int l(0); l < 4; l++)
    {
        if (TDStatus[l] > 0)
        {
            const auto& impulse_normal_l = impulse(3*l + 2);
            const auto& impNormal_partial = dImp_dx.row(3*l + 2);
            tcost.Phix += impNormal_partial * weight * impulse_normal_l;
            tcost.Phixx += weight * impNormal_partial.transpose() * impNormal_partial;
        }        
    }
}

template class WBFootPlaceReg<double>;
template class SwingFootClearance<double>;
template class ImpulseTerminalCost<double>;
