#include "WBM.h"
#include "casadiGen_MHPC.h"

#include <pinocchio/fwd.hpp>
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/rnea.hpp"


/*
    @brief: discrete-time dynamics
*/
template <typename T>
void WBM::Model<T>::dynamics(StateType &xnext, OutputType &y,
                             StateType &x, ContrlType &u, float t,
                             CtactStatusType &contact, const T &dt)
{
    (void)(t);
    dynamics_continuousTime(xnext, y, x, u, contact);

    // Update next state with Forward Euler
    xnext.template head<nq>() = q + v * dt;
    xnext.template tail<nv>() = v + qdd * dt;   
}

/*
    @brief: continuous-time dynamics
*/
template <typename T>
void WBM::Model<T>::dynamics_continuousTime(StateType &xdot, OutputType &y,
                                            StateType &x, ContrlType &u, CtactStatusType &contact)
{
    q << x.template head<nq>(); // x,y,z, yaw, pitch roll, joint angles
    v << x.template tail<nv>();
    tau = SelectionMat * u;

    // Update contacts
    updateContactFrameIds(contact);

    // Run KKT dynamics
    KKTContactDynamics();

    // Update xdot
    xdot.template head<nq>() = v;
    xdot.template tail<nv>() = qdd;

    // Update contact force
    y = GRF;
}

template <typename T>
void WBM::Model<T>::dynamics_partial(StateMapType &A, ContrlMapType &B, OutputMapType &C, DirectMapType &D,
                                     StateType &x, ContrlType &u, float t,
                                     CtactStatusType &contact, const T &dt)
{
    (void)(t);

    dynamics_partial_continuousTime(Ac_, B, C, D, x, u, contact);

    // Forward Euler
    A = Ac_ * dt + StateMapType::Identity();    
    B *= dt;

}

template <typename T>
void WBM::Model<T>::dynamics_partial_continuousTime(StateMapType &Ac, ContrlMapType &Bc, OutputMapType &C, DirectMapType &D,
                                                    StateType &x, ContrlType &u, CtactStatusType &contact)
{
    // Update q, v
    q << x.template head<nq>(); // x,y,z, yaw, pitch roll, joint angles
    v << x.template tail<nv>();
    tau = SelectionMat * u;

    // Update contacts
    updateContactFrameIds(contact);

    // Compute partials
    KKTContactDynamicsDerivatives();

    Ac.setZero();
    Ac.template topRightCorner<nq, nv>() = MatMN<T, nv, nv>::Identity();
    Ac.template bottomLeftCorner<nv, nq>() = dqdd_dq;
    Ac.template bottomRightCorner<nv, nv>() = dqdd_dv;

    Bc.setZero();
    Bc.template bottomRows<nv>() = dqdd_dtau * SelectionMat;

    C.setZero();
    D.setZero();
    const auto dGRF_du = dGRF_dtau * SelectionMat;
    for (size_t i = 0; i < N_contacts; i++)
    {
        C.template block<3, nv>(3 * contactFootIds_cur[i], 0) = dGRF_dq.block(3 * i, 0, 3, nv);
        C.template block<3, nv>(3 * contactFootIds_cur[i], nv) = dGRF_dv.block(3 * i, 0, 3, nv);
        D.template block<3, nu>(3 * contactFootIds_cur[i], 0) = dGRF_du.block(3 * i, 0, 3, nu);
    }
}

template <typename T>
void WBM::Model<T>::impact(DVec<T> &xnext, DVec<T> &y, const DVec<T> &x,
                           const CtactStatusType &contact_cur,
                           const CtactStatusType &contact_next)
{

    impact(xnext, x, contact_cur, contact_next);

    // Update impulse
    y = impulse;
}

template <typename T>
void WBM::Model<T>::impact(DVec<T> &xnext, const DVec<T> &x,
                           const CtactStatusType &contact_cur,
                           const CtactStatusType &contact_next)
{
    // Update q, v
    q << x.template head<nq>(); // x,y,z, yaw, pitch roll, joint angles
    v << x.template tail<nv>();

    CtactStatusType impact_status;
    impact_status.setZero();
    for (size_t foot = 0; foot < 4; ++foot)
    {
        if (contact_cur[foot] == 0 && contact_next[foot] == 1)
        {
            impact_status[foot] = 1;
        }
    }

    // Update contacts
    updateContactFrameIds(impact_status);

    // Run KKT impact dynamics
    KKTImpact();

    // Update post-impact state
    xnext.setZero(xs);
    xnext.template head(18) = q;
    xnext.template tail(18) = v_post_imp;
}

template <typename T>
void WBM::Model<T>::impact_partial(DMat<T> &dP_dx, DMat<T> &dy_dx,
                                   const DVec<T> &x,
                                   const CtactStatusType &contact_cur,
                                   const CtactStatusType &contact_next)
{
    impact_partial(dP_dx, x, contact_cur, contact_next);

    dy_dx.setZero(ys, xs);
    for (size_t i = 0; i < N_contacts; i++)
    {
        dy_dx.template block<3, nq>(3 * contactFootIds_cur[i], 0) = dImp_dq.template block(i, 0, 3, nq);
        dy_dx.template block<3, nv>(3 * contactFootIds_cur[i], nq) = dImp_dv.template block(i, 0, 3, nv);
    }
}

template <typename T>
void WBM::Model<T>::impact_partial(DMat<T> &dP_dx,
                                   const DVec<T> &x,
                                   const CtactStatusType &contact_cur,
                                   const CtactStatusType &contact_next)
{
    // Update q, v
    q << x.template head<nq>(); // x,y,z, yaw, pitch roll, joint angles
    v << x.template tail<nv>(); // vx, vy, vz, yaw rate, pitch rate, roll rate, joint vel

    CtactStatusType impact_status;
    impact_status.setZero();
    for (size_t foot = 0; foot < 4; ++foot)
    {
        if (contact_cur[foot] == 0 && contact_next[foot] == 1)
        {
            impact_status[foot] = 1;
        }
    }

    // Update contacts
    updateContactFrameIds(impact_status);

    // Run KKT Impact derivatives
    KKTImpactDerivatives();

    dP_dx.setZero(xs, xs);
    dP_dx.template topLeftCorner(nq, nq) = MatMN<T, nq, nq>::Identity();
    dP_dx.template bottomLeftCorner(nv, nq) = dvpost_dq;
    dP_dx.template bottomRightCorner(nv, nv) = dvpost_dv;
}

/*
    @brief: Get foot heights only for contact foot
*/
template <typename T>
void WBM::Model<T>::get_contactFootHeights(DVec<T> &EE_heights, const StateType &x, const CtactStatusType &contact)
{
    if (contact.sum() > 0)
    {
        updateContactFrameIds(contact);
        EE_heights.setZero(N_contacts);
        pinocchio::forwardKinematics(pin_model, pin_data, x.template head<nq>());
        for (size_t i = 0; i < N_contacts; i++)
        {
            pinocchio::updateFramePlacement(pin_model, pin_data, contactFrameIds_cur[i]);
            EE_heights[i] = pin_data.oMf[contactFrameIds_cur[i]].translation()[2];
        }
    }
}

/*
    @brief: Get XY positions for all foot
*/
template <typename T>
void WBM::Model<T>::get_footXYPositions(VecM<T, 8> &EE_XY, const StateType &x)
{    
    EE_XY.setZero();
    pinocchio::forwardKinematics(pin_model, pin_data, x.template head<nq>());
    for (size_t i = 0; i < 4; i++)
    {
        pinocchio::updateFramePlacement(pin_model, pin_data, contactFootFrameIds[i]);
        EE_XY.template segment<2>(2 * i) = pin_data.oMf[contactFootFrameIds[i]].translation().template head<2>();
    }
}

/*
    @brief: Get all foot positions in world frame
*/
template <typename T>
void WBM::Model<T>::get_footPositions(VecM<T,3> EE_pos[4], const StateType &x)
{            
    pinocchio::forwardKinematics(pin_model, pin_data, x.template head<nq>());
    for (size_t i = 0; i < 4; i++)
    {
        pinocchio::updateFramePlacement(pin_model, pin_data, contactFootFrameIds[i]);
        EE_pos[i] = pin_data.oMf[contactFootFrameIds[i]].translation().template head<3>();
    }
}


/*
    @brief: Get all foot velocities in world frame
*/
template <typename T>
void WBM::Model<T>::get_footVelocities(VecM<T,3> EE_vel[4], const StateType &x)
{            
    const auto& q_temp = x.template head<nq>();
    const auto& v_temp = x.template tail<nv>();
    pinocchio::forwardKinematics(pin_model, pin_data, q_temp, v_temp);
    for (size_t i = 0; i < 4; i++)
    {
        pinocchio::updateFramePlacement(pin_model, pin_data, contactFootFrameIds[i]);        

        EE_vel[i] = pinocchio::getFrameVelocity(pin_model, pin_data, contactFootFrameIds[i], pinocchio::LOCAL_WORLD_ALIGNED).toVector().head(3);        
    }
}

/*
    @brief: Get foot Jacobian only for contact foot
*/
template <typename T>
void WBM::Model<T>::get_contactFootJacobians(DMat<T> &J_EE, const StateType &x, const CtactStatusType &contact)
{
    if (contact.sum() > 0)
    {
        updateContactFrameIds(contact);
        pinocchio::computeJointJacobians(pin_model, pin_data, x.template head<nq>());
        J_EE.setZero(3 * N_contacts, nv);
        for (size_t i(0); i < N_contacts; ++i)
        {
            // Get spatial Jacobian
            J_6D.setZero();
            pinocchio::getFrameJacobian(pin_model, pin_data, contactFrameIds_cur[i], pinocchio::LOCAL_WORLD_ALIGNED, J_6D);

            // The top three rows of spatial Jacian corresponds to translation
            J_EE.block(3 * i, 0, 3, pin_model.nv) = J_6D.topRows(3);
        }
    }
}

/*
    @brief: Get all end-effector Jacobians
*/
template <typename T>
void WBM::Model<T>::get_footJacobians(MatMN<T, 3, nq> J_EE[4], const StateType& x)
{    
    pinocchio::computeJointJacobians(pin_model, pin_data, x.template head<nq>());
    
    for (size_t i(0); i < 4; ++i)
    {
        pinocchio::updateFramePlacement(pin_model, pin_data, contactFootFrameIds[i]);

        // Get spatial Jacobian
        J_6D.setZero();
        pinocchio::getFrameJacobian(pin_model, pin_data, contactFootFrameIds[i], pinocchio::LOCAL_WORLD_ALIGNED, J_6D);

        // The top three rows of spatial Jacian corresponds to translation
        J_EE[i] = J_6D.topRows(3);
    }
}


template <typename T>
void WBM::Model<T>::KKTContactDynamics()
{
    GRF.setZero();
    J_activeEE.setZero(3 * N_contacts, nv);
    if (N_contacts > 0)
    {
        // Run forward kinematics, and compute the drift term Jdot * qdot
        pinocchio::forwardKinematics(pin_model, pin_data, q, v, VecM<T, nv>::Zero());

        // Compute contact Jacobian for each contact point
        pinocchio::computeJointJacobians(pin_model, pin_data);        
        for (size_t i(0); i < N_contacts; ++i)
        {
            // Update contact foot placement
            pinocchio::updateFramePlacement(pin_model, pin_data, contactFrameIds_cur[i]);

            // Get spatial Jacobian
            J_6D.setZero();
            pinocchio::getFrameJacobian(pin_model, pin_data, contactFrameIds_cur[i], pinocchio::LOCAL_WORLD_ALIGNED, J_6D);

            // The top three rows of spatial Jacian corresponds to translation
            J_activeEE.block(3 * i, 0, 3, pin_model.nv) = J_6D.topRows(3);
        }
        
        gamma_activeEE.setZero(3 * N_contacts);
        for (size_t i(0); i < N_contacts; ++i)
        {
            // Get spatial acceleration of contact frame            

            // Linear component of spatial acceleration
            gamma_activeEE.segment(3 * i, 3) = pinocchio::getFrameAcceleration(pin_model, pin_data, contactFrameIds_cur[i], pinocchio::LOCAL_WORLD_ALIGNED).linear();

            // Convert to conventional acceleration
            vc_6D = pinocchio::getFrameVelocity(pin_model, pin_data, contactFrameIds_cur[i], pinocchio::LOCAL_WORLD_ALIGNED).toVector();
            vv_c = vc_6D.head(3);
            vw_c = vc_6D.tail(3);
            gamma_activeEE.segment(3 * i, 3) += vw_c.cross(vv_c);

            // Bamguart stabilization                        
            gamma_activeEE.segment(3 * i, 3) += 2 * BG_alpha*vv_c;            
        }
        pinocchio::crba(pin_model, pin_data, q);
        pinocchio::nonLinearEffects(pin_model, pin_data, q, v);        
        qdd = pinocchio::forwardDynamics(pin_model, pin_data, tau, J_activeEE, gamma_activeEE, 1e-12);

        // Map contact force
        for (size_t i = 0; i < N_contacts; i++)
        {
            GRF.template segment<3>(3 * contactFootIds_cur[i]) = pin_data.lambda_c.segment(3 * i, 3);
        }
    }
    else
    {
        pinocchio::aba(pin_model, pin_data, q, v, tau);
        qdd = pin_data.ddq;
    }
}

template <typename T>
void WBM::Model<T>::KKTImpact()
{
    // Run forward kinematics
    pinocchio::forwardKinematics(pin_model, pin_data, q);

    // Compute contact Jacobian for each contact point
    pinocchio::computeJointJacobians(pin_model, pin_data);
    impulse.setZero();
    J_activeEE.setZero(3 * N_contacts, nv);
    for (size_t i(0); i < N_contacts; ++i)
    {   
        // Update contact foot placement
        pinocchio::updateFramePlacement(pin_model, pin_data, contactFrameIds_cur[i]);

        // Get spatial Jacobian
        J_6D.setZero();
        pinocchio::getFrameJacobian(pin_model, pin_data, contactFrameIds_cur[i], pinocchio::LOCAL_WORLD_ALIGNED, J_6D);

        // The top three rows of spatial Jacian corresponds to translation
        J_activeEE.block(3 * i, 0, 3, pin_model.nv) = J_6D.topRows(3);
    }
    pinocchio::impulseDynamics(pin_model, pin_data, q, v, J_activeEE);

    v_post_imp = pin_data.dq_after;

    for (size_t i = 0; i < N_contacts; i++)
    {
        impulse.template segment<3>(3 * contactFootIds_cur[i]) = pin_data.impulse_c.template segment<3>(i);
    }    
}

template <typename T>
void WBM::Model<T>::KKTContactDynamicsDerivatives()
{
    // Improvement needed here to avoid re-calling the contact dynamics
    // This function will compute the Jacobians meeded for KKT inverse, qdd needed for RENA
    KKTContactDynamics();

    // Compute the inverse of the KKT matrix
    Kinv.setZero(nv + 3 * N_contacts, nv + 3 * N_contacts);
    pinocchio::computeKKTContactDynamicMatrixInverse(pin_model, pin_data, q, J_activeEE, Kinv);

    // Partial w.r.t fully actuated joint torque, need to right-multiply the selection matrix later
    dqdd_dtau = Kinv.topRows(nv).leftCols(nv);
    dGRF_dtau = -Kinv.bottomRows(3 * N_contacts).leftCols(nv);

    // Derivatives of inverse dynamics (assuming no contacts)
    pinocchio::computeRNEADerivatives(pin_model, pin_data, q, v, qdd);
    auto dtau_dv = pin_data.dtau_dv;
    auto dtau_dq = pin_data.dtau_dq;

    // Account for contacts in inverse dynamics
    computeContactForceDerivatives(q, GRF);
    dtau_dq -= dJTF_dq;

    /* KKT dynamics partial w.r.t. q and v */
    // Compute the partials of foot accerleration
    computeFootAccDerivatives();

    // Bamguart stabilization
    computeFootVelDerivatives(q, v);
    da_dq_activeEE += 2 * BG_alpha * dv_dq_activeEE;    
    da_dv_activeEE += 2 * BG_alpha * J_activeEE;

    // Compute dqdd_dq, dGRF_dq
    dqdd_dq = -Kinv.topRows(nv).leftCols(nv) * dtau_dq;
    dqdd_dq -= Kinv.topRows(nv).rightCols(3 * N_contacts) * da_dq_activeEE;

    dGRF_dq = Kinv.bottomRows(3 * N_contacts).leftCols(nv) * dtau_dq;
    dGRF_dq += Kinv.bottomRows(3 * N_contacts).rightCols(3 * N_contacts) * da_dq_activeEE;

    // Compute dqdd_dv, dGRF_dv
    dqdd_dv = -Kinv.topRows(nv).leftCols(nv) * dtau_dv;
    dqdd_dv -= Kinv.topRows(nv).rightCols(3 * N_contacts) * da_dv_activeEE;

    dGRF_dv = Kinv.bottomRows(3 * N_contacts).leftCols(nv) * dtau_dv;
    dGRF_dv += Kinv.bottomRows(3 * N_contacts).rightCols(3 * N_contacts) * da_dv_activeEE;    

}

template <typename T>
void WBM::Model<T>::KKTImpactDerivatives()
{
    // Compute contact Jacobians needed for KKT matrix, and v_pose, and impulse needed for RNEA
    KKTImpact(); // Improvement needed here to avoid re-computing the impulse

    // Derivatives of ID, assuming no external forces, and zero velocith
    pinocchio::computeRNEADerivatives(pin_model, pin_data, q, VecM<T, 18>::Zero(), v_post_imp - v);
    auto dtau_dq = pin_data.dtau_dq;

    MatMN<T, nv, nv> dgrav_dq;
    dgrav_dq.setZero();
    // Generalized gravity is not active in impulse dynamics, so we would like to take this out
    pinocchio::computeGeneralizedGravityDerivatives(pin_model, pin_data, q, dgrav_dq);

    // Correct dtau_dq. Removing the gravity effect
    dtau_dq -= dgrav_dq;

    // Account for contacts in inverse dynamics
    computeContactForceDerivatives(q, impulse);
    dtau_dq -= dJTF_dq;

    // Get the KKT Matix inverse
    Kinv.setZero(nv + 3 * N_contacts, nv + 3 * N_contacts);
    pinocchio::getKKTContactDynamicMatrixInverse(pin_model, pin_data, J_activeEE, Kinv);

    computeFootVelDerivatives(q, v_post_imp);

    dvpost_dq = -Kinv.topLeftCorner(nv, nv) * dtau_dq;
    dvpost_dq -= Kinv.topRightCorner(nv, 3 * N_contacts) * dv_dq_activeEE;
    pin_data.M.transpose().template triangularView<Eigen::Upper>() = pin_data.M.template triangularView<Eigen::Upper>();
    dvpost_dv = Kinv.topLeftCorner(nv, nv) * pin_data.M;

    dImp_dq = Kinv.bottomLeftCorner(3 * N_contacts, nv) * dtau_dq;
    dImp_dq += Kinv.bottomRightCorner(3 * N_contacts, 3 * N_contacts) * dv_dq_activeEE;
    dImp_dv = Kinv.bottomRightCorner(3 * N_contacts, 3 * N_contacts) * J_activeEE; // Different from what croccodyl has
}

template <typename T>
void WBM::Model<T>::updateContactFrameIds(const CtactStatusType &contacts)
{
    contactFrameIds_cur.clear();
    contactFootIds_cur.clear();
    for (size_t foot(0); foot < 4; ++foot)
    {
        if (contacts[foot] > 0)
        {
            contactFrameIds_cur.push_back(WBM::contactFootFrameIds[foot]);
            contactFootIds_cur.push_back(foot);
        }
    }
    N_contacts = contactFrameIds_cur.size();
}

/*
    @brief: Compute the Partial of foot velocity w.r.t. q for all foot
*/
template <typename T>
void WBM::Model<T>::computeFootVelDerivatives(VecM<T,nq>&q_, VecM<T,nv>&v_)
{

    std::vector<T *> arg = {q_.data(), v_.data()};

    std::vector<T *> res = {dv_dq_EE[0].data(),
                            dv_dq_EE[1].data(),
                            dv_dq_EE[2].data(),
                            dv_dq_EE[3].data()};

    casadi_interface(arg, res, dv_dq_EE[0].size(), footVelPartialDq,
                     footVelPartialDq_sparsity_out,
                     footVelPartialDq_work);

    // Map to active contacts
    dv_dq_activeEE.setZero(3 * N_contacts, nv);
    for (size_t i = 0; i < N_contacts; i++)
    {
        dv_dq_activeEE.block(3 * i, 0, 3, nv) = dv_dq_EE[contactFootIds_cur[i]];
    }
}

/*
    @brief: Compute the Partial of foot velocity w.r.t. q for all foot
*/
template <typename T>
void WBM::Model<T>::get_footVelDerivatives(MatMN<T, 3, nv> J_EE[4], const StateType& x)
{
    VecM<T, nq> q_temp= x.template head<nq>();
    VecM<T, nv> v_temp= x.template tail<nv>();

    std::vector<T *> arg = {q_temp.data(), v_temp.data()};

    for (int i = 0; i < 4; i++)
    {
        J_EE[i].setZero();
    }
    
    std::vector<T *> res = {J_EE[0].data(),
                            J_EE[1].data(),
                            J_EE[2].data(),
                            J_EE[3].data()};

    casadi_interface(arg, res, J_EE[0].size(), footVelPartialDq,
                     footVelPartialDq_sparsity_out,
                     footVelPartialDq_work);    
}

/*
    @brief: Compute the Partial of foot acceleration w.r.t. q for all foot
*/
template <typename T>
void WBM::Model<T>::computeFootAccDerivatives()
{

    std::vector<T *> arg = {q.data(), v.data(), qdd.data()};

    std::vector<T *> res = {da_dq_EE[0].data(),
                            da_dq_EE[1].data(),
                            da_dq_EE[2].data(),
                            da_dq_EE[3].data()};

    casadi_interface(arg, res, da_dq_EE[0].size(), footAccPartialDq,
                     footAccPartialDq_sparsity_out,
                     footAccPartialDq_work);

    res.clear();

    for (auto &da_dv : da_dv_EE)
    {
        res.push_back(da_dv.data());
    }

    casadi_interface(arg, res, da_dv_EE[0].size(), footAccPartialDv,
                     footAccPartialDv_sparsity_out,
                     footAccPartialDv_work);

    // Map to active contacts
    da_dq_activeEE.setZero(3 * N_contacts, nv);
    da_dv_activeEE.setZero(3 * N_contacts, nv);
    for (size_t i = 0; i < N_contacts; i++)
    {
        da_dq_activeEE.block(3 * i, 0, 3, nv) = da_dq_EE[contactFootIds_cur[i]];
        da_dv_activeEE.block(3 * i, 0, 3, nv) = da_dv_EE[contactFootIds_cur[i]];
    }
}

/*
    @Brief: Compute Partial of JTF w.r.t. q for all foot
*/
template <typename T>
void WBM::Model<T>::computeContactForceDerivatives(VecM<T,nq>&q_, VecM<T,12>&GRF_)
{
    std::vector<T *> arg = {q_.data(), GRF_.data()};

    std::vector<T *> res = {dJTF_dq_EE[0].data(),
                            dJTF_dq_EE[1].data(),
                            dJTF_dq_EE[2].data(),
                            dJTF_dq_EE[3].data()};

    casadi_interface(arg, res, dJTF_dq_EE[0].size(), footForcePartialDq,
                     footForcePartialDq_sparsity_out,
                     footForcePartialDq_work);

    // Use only the active contacts
    dJTF_dq.setZero();
    for (size_t i = 0; i < N_contacts; i++)
    {
        dJTF_dq += dJTF_dq_EE[contactFootIds_cur[i]];
    }
}

template <typename T>
void WBM::Model<T>::printModelInfo()
{
    // Print model information
    std::cout << "--------------- Model Information -----------------\n";
    std::cout << "Robot name: " << pin_model.name << " \n";
    std::cout << "Number of joints: " << pin_model.njoints << " \n";
    std::cout << "Number of bodies: " << pin_model.nbodies << " \n";
    std::cout << "Number of position variables nq: " << pin_model.nq << " \n";
    std::cout << "Number of velocity variables nv: " << pin_model.nv << " \n"
              << " \n";

    std::cout << "--------------- Joint Information -----------------\n";
    for (int j(0); j < pin_model.njoints; ++j)
    {
        std::cout << "Joint " << j << " name: " << pin_model.names[j] << std::endl;
    }

    // std::cout << "--------------- Frame Information -----------------\n";
    // for (const auto &frame : pin_model.frames)
    // {
    //     std::cout << "Frame id: " << pin_model.getFrameId(frame.name) << " Frame name: " << frame.name << std::endl;
    // }

    // // // Print contact information
    // std::cout << "--------------- Current Contact Information -----------------\n";
    // std::cout << "Numner of contacts: " << N_contacts << "\n";
    // for (const auto & f : contactFrameIds_cur)
    // {
    //     std::cout << "contactFrameId = " << f << " frameName: " << pin_model.frames[f].name << " \n";
    // }
}

template class WBM::Model<double>;

template <typename T>
void computeLegPosition(Vec3<T>& p, const Vec3<T>& qleg, int leg)
{
    T l1 = 0.062;
    T l2 = 0.209;
    T l3 = 0.195;
    T l4 = 0.004;
    T sideSigns[] = {1,-1,1,-1};
    T sideSign = sideSigns[leg];

    T s1 = std::sin(qleg(0));
    T s2 = std::sin(-qleg(1));
    T s3 = std::sin(-qleg(2));

    T c1 = std::cos(qleg(0));
    T c2 = std::cos(-qleg(1));
    T c3 = std::cos(-qleg(2));

    T c23 = c2 * c3 - s2 * s3;
    T s23 = s2 * c3 + c2 * s3;

    p[0] = l3 * s23 + l2 * s2;
    p[1] = (l1+l4) * sideSign * c1 + l3 * (s1 * c23) + l2 * c2 * s1;
    p[2] = (l1+l4) * sideSign * s1 - l3 * (c1 * c23) - l2 * c1 * c2;
}

template <typename T>
void computeLegJacobian(MatMN<T, 3, 3>& J, const Vec3<T>& qleg, int leg)
{
    
    T l1 = 0.062;
    T l2 = 0.209;
    T l3 = 0.195;
    T l4 = 0.004;
    T sideSigns[] = {1,-1,1,-1};
    T sideSign = sideSigns[leg];

    T s1 = std::sin(qleg(0));
    T s2 = std::sin(-qleg(1));
    T s3 = std::sin(-qleg(2));

    T c1 = std::cos(qleg(0));
    T c2 = std::cos(-qleg(1));
    T c3 = std::cos(-qleg(2));

    T c23 = c2 * c3 - s2 * s3;
    T s23 = s2 * c3 + c2 * s3;

    J(0, 0) = 0;
    J(0, 1) = l3 * c23 + l2 * c2;
    J(0, 2) = l3 * c23;
    J(1, 0) = l3 * c1 * c23 + l2 * c1 * c2 - (l1+l4) * sideSign * s1;
    J(1, 1) = -l3 * s1 * s23 - l2 * s1 * s2;
    J(1, 2) = -l3 * s1 * s23;
    J(2, 0) = l3 * s1 * c23 + l2 * c2 * s1 + (l1+l4) * sideSign * c1;
    J(2, 1) = l3 * c1 * s23 + l2 * c1 * s2;
    J(2, 2) = l3 * c1 * s23;

    J.template rightCols<2>() *= -1;
}
template void computeLegPosition<double>(Vec3<double>& p, const Vec3<double>& qleg, int leg);
template void computeLegJacobian<double>(MatMN<double, 3, 3>& J, const Vec3<double>& qleg, int leg);