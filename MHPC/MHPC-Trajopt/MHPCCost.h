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
        VecM<T, 3> q_pos(0.1, 0.1, 20);
        VecM<T, 3> q_eul(1, 2, 2);
        VecM<T, 3> q_v(.5, .5, .5);
        VecM<T, 3> q_eulrate(.2, .2, .2);        
        VecM<T, 12> q_qJ, q_qJd;                    
        
        q_qJ.setConstant(.05);        
        q_qJd.setConstant(.01);
      
        this->Q.setZero();
        this->Q.diagonal() << q_pos, q_eul, q_qJ, q_v, q_eulrate, q_qJd;

        VecM<T, 12> q_scale;
        q_scale << -4.5*contact[0]+5,-4.5*contact[0]+5,-4.5*contact[0]+5,
                   -4.5*contact[1]+5,-4.5*contact[1]+5,-4.5*contact[1]+5,
                   -4.5*contact[2]+5,-4.5*contact[2]+5,-4.5*contact[2]+5,
                   -4.5*contact[3]+5,-4.5*contact[3]+5,-4.5*contact[3]+5;
        this->Q.template block<12,12>(6,6) *= q_scale.asDiagonal();
        
        /* Terminal state weighting matrix */
        this->Qf = this->Q;
        this->Qf.template block<12,12>(6,6).diagonal().setConstant(0.01);
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
        VecM<T, 3> qFoot(5,5,1);
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
class SwingFootClearance : public CostBase<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::State;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Contrl;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Output;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::RCost;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::TCost;    

    SwingFootClearance(const VecM<int, 4>& contact_in) :
                   contact(contact_in),                    
                   CostBase<T, WBM::xs, WBM::us, WBM::ys>("Swing Clearance")
    {
        weight = 20;              
        des_height = 0.1;            
    }

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t=0) override;

    void terminal_cost(TCost&, const State& x, float tend=0) override {}

    void terminal_cost_par(TCost&, const State&x, float tend=0) override {}    

private:
    Vec4<int> contact;
    T weight;     
    T des_height;
                   
    Vec3<T> lq;
    Mat3<T> lqq;
};

template <typename T>
class ImpulseTerminalCost : public CostBase<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::State;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Contrl;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::Output;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::RCost;
    using typename CostBase<T,WBM::xs, WBM::us, WBM::ys>::TCost;    

    ImpulseTerminalCost(const VecM<int, 4>& contact_cur_in,
                        const VecM<int, 4>& contact_next_in,
                        std::shared_ptr<WBM::Model<T>> wbm_ptr_in):
    contact_cur(contact_cur_in),
    contact_next(contact_cur_in),                    
    wbm_ptr(wbm_ptr_in),
    CostBase<T, WBM::xs, WBM::us, WBM::ys>("ImpulsePenalty")
    {
        TDStatus.setZero();
        weight = 10;
        for (int leg = 0; leg < 4; leg++)
        {
            if (contact_cur[leg] == 0 && contact_next[leg] == 1)
            {
                TDStatus[leg] = 1;
            }
        }  
        xpostImp.setZero(WBM::xs);
        impulse.setZero(WBM::ys);
        dImp_dx.setZero(WBM::ys, WBM::xs);
        dxpost_dx.setZero(WBM::xs, WBM::xs);
    }

    void running_cost(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t) override {}

    void running_cost_par(RCost&, const State& x, const Contrl& u, const Output& y, T dt, float t) override{}

    void terminal_cost(TCost&, const State& x, float tend=0) override;

    void terminal_cost_par(TCost&, const State&x, float tend=0) override;

private:
    Vec4<int> contact_cur;
    Vec4<int> contact_next;    
    Vec4<int> TDStatus;
    T weight;
    DVec<T> xpostImp;
    DVec<T> impulse;
    DMat<T> dImp_dx;    
    DMat<T> dxpost_dx;
    std::shared_ptr<WBM::Model<T>> wbm_ptr;
};

template <typename T>
class SRBTrackingCost : public QuadraticTrackingCost<T, SRBM::xs, SRBM::us, SRBM::ys>
{
public:    
    SRBTrackingCost() : QuadraticTrackingCost<T, SRBM::xs, SRBM::us, SRBM::ys>()
    {   
        /* Intermediate state weighting matrix */
        VecM<T, 3> q_pos(.1, .1, 10);
        VecM<T, 3> q_eul(1, 3, 3);       
        VecM<T, 3> q_v(1, 1, 3); 
        VecM<T, 3> q_eulrate(.2, .5, .5);        
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