#pragma once
#ifndef WBM_H
#define WBM_H

#include <vector>
#include <memory>
#include "HSDDP_CPPTypes.h"
#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/algorithm/contact-dynamics.hpp"

namespace WBM
{
    // Whole-Body state: eul, pos, omega_body, vel, qJ, qdJ
    const size_t xs = 36;
    const size_t us = 12;
    const size_t ys = 12;
    const size_t nu = 12; // number of actuated joints
    const size_t nq = 18; // number of generalized joints
    const size_t nv = 18; // number of generalized velocity

    const std::vector<size_t> contactFootFrameIds ={11, 19, 27, 35}; //FL, FR, HL, HR(obtained from loaded urdf model)

    template <typename T>
    class Model
    {
    public:
        typedef VecM<T, xs> StateType;
        typedef VecM<T, us> ContrlType;
        typedef VecM<T, ys> OutputType;

        typedef MatMN<T, xs, xs> StateMapType;
        typedef MatMN<T, xs, us> ContrlMapType;
        typedef MatMN<T, ys, xs> OutputMapType;
        typedef MatMN<T, ys, us> DirectMapType;

        typedef VecM<int, 4> CtactStatusType;

        MatMN<T, nv, nu> SelectionMat;  // selection matrix

    public:
        Model(const pinocchio::ModelTpl<T>& pin_model_in, T BG_alpha_in=10, T BG_beta_in=0) :
        pin_model(pin_model_in),
        pin_data(pinocchio::DataTpl<T>(pin_model_in)),
        BG_alpha(BG_alpha_in),
        BG_beta(BG_beta_in)
        {                      
            SelectionMat.setZero();            
            SelectionMat.diagonal(-(nv - nu)).setOnes();

            contactFrameIds_cur.clear();
            N_contacts = 0;

            dv_dq_EE = std::vector<MatMN<T, 3, nv>>(4, MatMN<T, 3, nv>::Zero());
            da_dq_EE = std::vector<MatMN<T, 3, nv>>(4, MatMN<T, 3, nv>::Zero());
            da_dv_EE = std::vector<MatMN<T, 3, nv>>(4, MatMN<T, 3, nv>::Zero());
            dJTF_dq_EE = std::vector<MatMN<T, nv, nv>>(4, MatMN<T, nv, nv>::Zero());            
        }

        void dynamics(StateType &xnext, OutputType &y,
                      StateType &x, ContrlType &u, float t,
                      CtactStatusType &contact, const T &dt);

        void dynamics_continuousTime(StateType &xdot, OutputType &y,
                      StateType &x, ContrlType &u, CtactStatusType &contact);                      

        void dynamics_partial(StateMapType &A, ContrlMapType &B, OutputMapType &C, DirectMapType &D,
                              StateType &x, ContrlType &u, float t,
                              CtactStatusType &contact, const T &dt); // To Be Implemented

        void dynamics_partial_continuousTime(StateMapType &Ac, ContrlMapType &Bc, OutputMapType &C, DirectMapType &D,
                              StateType &x, ContrlType &u,
                              CtactStatusType &contact);                              

        void impact(DVec<T> &xnext, DVec<T> &y, const DVec<T> &x,
                    const CtactStatusType &contact_cur,
                    const CtactStatusType &contact_next);

        void impact(DVec<T> &xnext, const DVec<T> &x,
                    const CtactStatusType &contact_cur,
                    const CtactStatusType &contact_next);                     

        void impact_partial(DMat<T> &df_dx, DMat<T> &dy_dx,
                            const DVec<T>& x,
                            const CtactStatusType &contact_cur,
                            const CtactStatusType &contact_next);   

        void impact_partial(DMat<T>& df_dx,
                            const DVec<T>& x,
                            const CtactStatusType& contact_cur,
                            const CtactStatusType& contact_next);                                  

        void get_contactFootJacobians(DMat<T>& J_EE, 
                                     const StateType& x, 
                                     const CtactStatusType& contact);

        void get_contactFootHeights(DVec<T>& EE_pos, 
                                   const StateType& x,
                                   const CtactStatusType& contact);

        void get_footXYPositions(VecM<T,8>& EE_XY,
                                    const StateType& x);                                   

        void get_footPositions(VecM<T,12>& EE_pos, 
                               const StateType& x);

        void printModelInfo();                                                

    private:

        void KKTContactDynamics();                            

        void KKTContactDynamicsDerivatives();

        void KKTImpact();

        void KKTImpactDerivatives();

        void updateContactFrameIds(const CtactStatusType &contact);

        void computeContactForceDerivatives(VecM<T,nq>&q_, VecM<T,12>&GRF_);

        void computeFootAccDerivatives();

        void computeFootVelDerivatives(VecM<T,nq>&q_, VecM<T,nv>&v_);

    private:        
        VecM<T, nq> q;              // generalized joint
        VecM<T, nv> v;              // generalized veloctiy        
        VecM<T, nv> v_post_imp;     // post-impact generalized velocity
        VecM<T, nv> tau;            // actuated joint torque (assume fully actuated)
        VecM<T, nv> qdd;            // generalized acceleration
        VecM<T, 12> GRF;            // GRF for all foot (output)
        VecM<T, 12> impulse;        // Impulse for all foot (output of touchdown)

        DVec<T> GRF_activeEE;        // GRF for active foot
        DVec<T> impulse_acticeEE;    // Contact impulse (3D) of all active contacts

        MatMN<T, nv, nv> dqdd_dtau;     // Assume fully actuated
        MatMN<T, nv, nq> dqdd_dq;
        MatMN<T, nv, nv> dqdd_dv;
        MatMN<T, nv, nv> dvpost_dq;
        MatMN<T, nv, nv> dvpost_dv;

        DMat<T> dGRF_dtau;          // Assume fully actuated
        DMat<T> dGRF_dq;            // Dynamics size depending on foot contacts
        DMat<T> dGRF_dv;
        DMat<T> dImp_dq;
        DMat<T> dImp_dv;

        DMat<T> K;          // KKT matrix (dynamic size because contact size may vary)
        DMat<T> Kinv;       // KKT matrix inverse (same as above)


        DMat<T> J_activeEE;             // 3D contact Jacobian for all active contact feet
        MatMN<T, 6, nv> J_6D;           // spatial contact Jacobian for one foot        

        VecM<T, 6> vc_6D;               // spatial velocity of one contact
        VecM<T, 6> ac_6D;               // spatial acceleration of one contact
        VecM<T, 3> vv_c;                // linear component of v_c
        VecM<T, 3> vw_c;                // angular componnet of v_c
        DVec<T> gamma_activeEE;         // 3D drift (Jd*qd) terms of active contact points

        std::vector<size_t> contactFrameIds_cur;    // contact frame ids for current contacts
        std::vector<size_t> contactFootIds_cur;     // [FL, FR],     [0, 1]
                                                    // [HL, HR] ->   [2, 3]
        size_t N_contacts;

        // Bamguart stablization
        T BG_alpha;
        T BG_beta;
    
    private: 
        // pinocchio model and data
        pinocchio::ModelTpl<T> pin_model;
        pinocchio::DataTpl<T> pin_data;
        
        // Kinematics derivatives (for all foot)
        std::vector<MatMN<T,3,nv>> dv_dq_EE;
        std::vector<MatMN<T,3,nv>> da_dq_EE;
        std::vector<MatMN<T,3,nv>> da_dv_EE;
        DMat<T> dv_dq_activeEE;
        DMat<T> da_dq_activeEE;
        DMat<T> da_dv_activeEE;

        // Partial of J'F in world frame (for all foot)
        std::vector<MatMN<T, nv, nv>>  dJTF_dq_EE;       
        MatMN<T, nv, nv> dJTF_dq; 
    };

}



#endif