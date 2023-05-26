#ifndef MHPCCOST_H
#define MHPCCOST_H

#include "WBM.h"    // WBM.h has to be called before HSDDP_CompoundTypes.h due to pinocchio setup
#include "SRBM.h"
#include "SinglePhaseInterface.h"

template <typename T>
class WBTrackingCost : public QuadraticTrackingCost<T, WBM::xs, WBM::us, WBM::ys>
{
public:
    WBTrackingCost(const VecM<int, 4>& contact) : QuadraticTrackingCost<T, WBM::xs, WBM::us, WBM::ys>()
    {   
        /* Intermediate state weighting matrix */
        VecM<T, 3> q_pos(.1, .1, 10);
        VecM<T, 3> q_eul(1, 1, 1);
        VecM<T, 3> q_v(1, 1, 4);
        VecM<T, 3> q_eulrate(.5, .5, .5);        
        VecM<T, 12> q_qJ, q_qJd;                    
        VecM<T, 12> q_scale;
        
        q_qJ.setConstant(5);        
        q_qJd.setConstant(.1);
      

        this->Q.setZero();
        this->Q.diagonal() << q_pos, q_eul, q_qJ, q_v, q_eulrate, q_qJd;
        
        /* Terminal state weighting matrix */
        this->Qf = this->Q;
        this->Qf.template block<12,12>(6,6).diagonal().setConstant(0.1);
        this->Qf.template bottomRightCorner<12,12>().diagonal().setConstant(0.1);

        /* Control (torque) weighting matricx */
        VecM<T, WBM::us> r;        
        r.setConstant(.1);      
        this->R = r.asDiagonal();       

        /* Output weighting matrix */
        this->S.setZero();
    }       
};

template <typename T>
class SRBTrackingCost : public QuadraticTrackingCost<T, SRBM::xs, SRBM::us, SRBM::ys>
{
public:    
    SRBTrackingCost() : QuadraticTrackingCost<T, SRBM::xs, SRBM::us, SRBM::ys>()
    {   
        /* Intermediate state weighting matrix */
        VecM<T, 3> q_pos(.1, .1, 20);
        VecM<T, 3> q_eul(1, 1, 1);       
        VecM<T, 3> q_v(1, 1, 4); 
        VecM<T, 3> q_eulrate(.5, .5, .5);        
        this->Q.setZero();
        this->Q.diagonal() << q_pos, q_eul, q_v, q_eulrate;
        this->Q *= 1;
        
        /* Terminal state weighting matrix */        
        this->Qf = this->Q;

        /* Control (torque) weighting matricx */
        VecM<T, SRBM::us> r;        
        r.setConstant(.001);      
        this->R = r.asDiagonal();       
        
    }       
};


#endif