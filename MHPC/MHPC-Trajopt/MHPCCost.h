#ifndef MHPCCOST_H
#define MHPCCOST_H

#include "WBM.h"    // WBM.h has to be called before HSDDP_CompoundTypes.h due to pinocchio setup
#include "SRBM.h"
#include "SinglePhaseInterface.h"

template <typename T>
class WBTrackingCost : public QuadraticTrackingCost<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    WBTrackingCost(const VecM<int, 4>& contact) : 
    QuadraticTrackingCost<T, WBM::xs, WBM::us, WBM::ys>()
    {   
        /* Intermediate state weighting matrix */
        VecM<T, 3> q_pos(0.0, 0, 50);
        VecM<T, 3> q_eul(2.0, 10, 5);
        VecM<T, 3> q_v(2.0, 4.0, 4.0);
        VecM<T, 3> q_eulrate(1.0, 2.0, 2.0);        
        VecM<T, 12> q_qJ, q_qJd;                    
        
        q_qJ.setConstant(1.0);        
        q_qJd.setConstant(.01);
      
        this->Q.setZero();
        this->Q.diagonal() << q_pos, q_eul, q_qJ, q_v, q_eulrate, q_qJd;
        
        /* Terminal state weighting matrix */
        this->Qf = this->Q;
        this->Qf.template block<12,12>(6,6).diagonal().setConstant(0.5);
        this->Qf.template bottomRightCorner<12,12>().diagonal().setConstant(0.01);        

        /* Control (torque) weighting matricx */
        VecM<T, WBM::us> r;        
        r.setConstant(.1);      
        this->R = r.asDiagonal();              
    }       
};

template <typename T>
class WBFootPlaceReg : public CostBase<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::State;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Contrl;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Output;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::RCost;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::TCost;    

    WBFootPlaceReg(const VecM<int, 4>& contact_in,
                   std::shared_ptr<QuadReference> quad_reference_in,
                   std::shared_ptr<WBM::Model<T>> wbm_ptr_in) :
                   contact(contact_in), 
                   quad_astate(nullptr),
                   CostBase<T, WBM::xs, WBM::us, WBM::ys>("Foot regularization")
    {
        VecM<T, 3> qFoot(10.0,10.0,1.0);
        QFoot = qFoot.asDiagonal();        
        
        quad_reference = quad_reference_in;
        wbm_ptr = wbm_ptr_in;        
    }

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void terminal_cost(TCost&, const State& x, float tend=0) override;

    void terminal_cost_par(TCost&, const State&x, float tend=0) override;    

private:
    Vec4<int> contact;
    MatMN<T, 3, 3> QFoot;
    VecM<double, 3> pCoM_des;
    VecM<double, 3> prel_des;
    
    Vec3<T> pFoot[4];
    Vec3<T> prel;
    Vec3<T> pCoM;        
    Vec3<T> d_prel;
    VecM<T, WBM::nq> lq;
    MatMN<T, WBM::nq, WBM::nq> lqq;
    MatMN<T, 3, WBM::nq> J_foot[4];

    std::shared_ptr<QuadReference> quad_reference;
    std::shared_ptr<WBM::Model<T>> wbm_ptr;
    QuadAugmentedState* quad_astate;
};


template <typename T>
class SwingFootPosTracking : public CostBase<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::State;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Contrl;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Output;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::RCost;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::TCost;    

    SwingFootPosTracking(const VecM<int, 4>& contact_in, 
                      std::shared_ptr<QuadReference> quad_reference_in,
                      std::shared_ptr<WBM::Model<T>> wbm_ptr_in) :
                      contact(contact_in),                    
                      quad_reference(quad_reference_in),
                      wbm_ptr(wbm_ptr_in),
                      CostBase<T, WBM::xs, WBM::us, WBM::ys>("Swing_Pos_Tracking")
    {
        VecM<T, 3> qFoot(10,10,40);
        QFoot = qFoot.asDiagonal();                 
    }

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void terminal_cost(TCost&, const State& x, float tend=0) override {}

    void terminal_cost_par(TCost&, const State&x, float tend=0) override {}    

private:
    Vec4<int> contact;
    MatMN<T, 3, 3> QFoot;
    VecM<double, 3> pCoM_des;
    VecM<double, 3> prel_des;
    
    Vec3<T> pFoot[4];
    Vec3<T> prel;
    Vec3<T> pCoM;        
    Vec3<T> d_prel;
    VecM<T, WBM::nq> lq;
    MatMN<T, WBM::nq, WBM::nq> lqq;
    MatMN<T, 3, WBM::nq> J_foot[4];

    std::shared_ptr<QuadReference> quad_reference;
    QuadAugmentedState* quad_astate;
    std::shared_ptr<WBM::Model<T>> wbm_ptr;
};

template <typename T>
class SwingFootVelTracking : public CostBase<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::State;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Contrl;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Output;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::RCost;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::TCost;    

    SwingFootVelTracking(const VecM<int, 4>& contact_in, 
                      std::shared_ptr<QuadReference> quad_reference_in,
                      std::shared_ptr<WBM::Model<T>> wbm_ptr_in) :
                      contact(contact_in),                    
                      quad_reference(quad_reference_in),
                      wbm_ptr(wbm_ptr_in),
                      CostBase<T, WBM::xs, WBM::us, WBM::ys>("Swing_Vel_Tracking")
    {
        VecM<T, 3> qFoot(2.0,2.0,4.0);
        QFoot = qFoot.asDiagonal();                 
    }

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void terminal_cost(TCost&, const State& x, float tend=0) override {}

    void terminal_cost_par(TCost&, const State&x, float tend=0) override {}    

private:
    Vec4<int> contact;
    MatMN<T, 3, 3> QFoot;
   
    Vec3<T> vFoot[4];
    MatMN<T, 3, WBM::nq> J_foot[4];     // Foot Jacobian: d_EEvel/d_v
    MatMN<T, 3, WBM::nq> Jv_foot[4];    // Partial of foot vel w.r.. v: d_EEvel/dq

    std::shared_ptr<QuadReference> quad_reference;
    QuadAugmentedState* quad_astate;
    std::shared_ptr<WBM::Model<T>> wbm_ptr;
};

template<typename T>
class TDVelocityPenalty : public CostBase<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::State;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Contrl;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Output;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::RCost;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::TCost;    

    TDVelocityPenalty(const VecM<int, 4>& TDStatus_in,                       
                      std::shared_ptr<WBM::Model<T>> wbm_ptr_in) :
                      TDStatus(TDStatus_in),                                          
                      wbm_ptr(wbm_ptr_in),
                      CostBase<T, WBM::xs, WBM::us, WBM::ys>("TDVelPenalty") {}

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override {}

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override {}

    void terminal_cost(TCost&, const State& x, float tend=0) override;

    void terminal_cost_par(TCost&, const State&x, float tend=0) override;    

private:
    Vec4<int> TDStatus;
    T qFoot{1.0};
   
    Vec3<T> vFoot[4];
    MatMN<T, 3, WBM::nq> J_foot[4];     // Foot Jacobian: d_EEvel/d_v
    MatMN<T, 3, WBM::nq> Jv_foot[4];    // Partial of foot vel w.r.. v: d_EEvel/dq
    
    std::shared_ptr<WBM::Model<T>> wbm_ptr;
};



template <typename T>
class SRBTrackingCost : public QuadraticTrackingCost<T, SRBM::xs, SRBM::us, SRBM::ys>
{
public:    
    SRBTrackingCost() : QuadraticTrackingCost<T, SRBM::xs, SRBM::us, SRBM::ys>()
    {   
        /* Intermediate state weighting matrix */
        VecM<T, 3> q_pos(0, 0, 50);
        VecM<T, 3> q_eul(0, 10, 5);       
        VecM<T, 3> q_v(2.0, 3.0, 3.0); 
        VecM<T, 3> q_eulrate(.5, .5, .5);        
        this->Q.setZero();
        this->Q.diagonal() << q_pos, q_eul, q_v, q_eulrate;
        
        /* Terminal state weighting matrix */        
        this->Qf = .5*this->Q;                 

        /* Control (torque) weighting matricx */
        VecM<T, SRBM::us> r;        
        r.setConstant(.01);      
        this->R = r.asDiagonal();       
        
    }       
};


#endif