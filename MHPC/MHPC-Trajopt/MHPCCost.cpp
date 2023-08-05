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
            tcost.Phi += footReg;
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

            tcost.Phix.template head<WBM::nq>() += 2*lq;
            tcost.Phixx.template topLeftCorner<WBM::nq, WBM::nq>() += 2*lqq;
        }
    }
}

template <typename T>
void SwingFootPosTracking<T>::running_cost(RCost &rcost, const State &x, const Contrl &u, const Output &y, T dt, float t)
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
        if (contact[l] == 0)
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
void SwingFootPosTracking<T>::running_cost_par(RCost &rcost, const State &x, const Contrl &u, const Output &y, T dt, float t)
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
        if (contact[l] == 0)
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
void SwingFootVelTracking<T>::running_cost(RCost &rcost, const State &x, const Contrl &u, const Output &y, T dt, float t)
{
    (void)(y);
    (void)(u);
    quad_astate = quad_reference->get_a_reference_ptr_at_t(t);

    wbm_ptr->get_footVelocities(vFoot, x);
    rcost.l = 0;
    for (size_t l = 0; l < 4; l++)
    {
        if (contact[l] == 0)
        {            
            Vec3<T> dvFoot = vFoot[l] - quad_astate->foot_velocities.segment<3>(l * 3);           
            T footReg = .5 * dvFoot.transpose() * QFoot * dvFoot;
            footReg *= dt;
            rcost.l += footReg;
        }
    }
}

template <typename T>
void SwingFootVelTracking<T>::running_cost_par(RCost &rcost, const State &x, const Contrl &u, const Output &y, T dt, float t)
{
    (void)(y);
    (void)(u);
    quad_astate = quad_reference->get_a_reference_ptr_at_t(t);

    wbm_ptr->get_footVelocities(vFoot, x);
    wbm_ptr->get_footVelDerivatives(Jv_foot, x);
    wbm_ptr->get_footJacobians(J_foot, x);
    
    rcost.lx.setZero();
    rcost.lxx.setZero();
    for (size_t l = 0; l < 4; l++)
    {
        if (contact[l] == 0)
        {
            Vec3<T> dvFoot = vFoot[l] - quad_astate->foot_velocities.segment<3>(l * 3);           
            MatMN<T, 3, WBM::xs> J;
            J.template leftCols<WBM::nq>() = Jv_foot[l];
            J.template rightCols<WBM::nv>() = J_foot[l];
            rcost.lx += (J.transpose() * QFoot * dvFoot) * dt;
            rcost.lxx += (J.transpose() * QFoot * J) * dt;            
        }
    }
}

template <typename T>
void TDVelocityPenalty<T>::terminal_cost(TCost&tcost , const State& x, float tend)
{    
    wbm_ptr->get_footVelocities(vFoot, x);
    tcost.Phi = 0;
    for (size_t l = 0; l < 4; l++)
    {
        if (TDStatus[l] == 1)
        {            
            T dvFoot = vFoot[l][2];
            T footReg = .5 * dvFoot * qFoot * dvFoot;            
            tcost.Phi += footReg;
        }
    }
}

template <typename T>
void TDVelocityPenalty<T>::terminal_cost_par(TCost&tcost , const State& x, float tend)
{        
    wbm_ptr->get_footVelocities(vFoot, x);
    wbm_ptr->get_footVelDerivatives(Jv_foot, x);
    wbm_ptr->get_footJacobians(J_foot, x);
    
    tcost.Phix.setZero();
    tcost.Phixx.setZero();
    for (size_t l = 0; l < 4; l++)
    {
        if (TDStatus[l] == 1)
        {
            T dvFoot = vFoot[l][2];
            MatMN<T, 1, WBM::xs> J;
            J << Jv_foot[l].template bottomRows<1>(), J_foot[l].template bottomRows<1>(); 
           
            tcost.Phix += J.transpose() * qFoot * dvFoot ;
            tcost.Phixx += J.transpose() * qFoot * J;            
        }
    }
}



template class WBFootPlaceReg<double>;
template class SwingFootPosTracking<double>;
template class SwingFootVelTracking<double>;
template class TDVelocityPenalty<double>;

