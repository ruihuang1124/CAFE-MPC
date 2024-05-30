#ifndef SRBM_H
#define SRBM_H

#include "HSDDP_CPPTypes.h"
#include <vector>
#include <cassert>
#include "casadiGen_MHPC.h"
#include "QuadReference.h"
#include "MHPCFootStep.h"

namespace SRBM
{
    const size_t xs = 12;
    const size_t us = 12;
    const size_t ys = 0;

    template <typename T>
    class Model
    {
    public:        
        typedef VecM<T, xs> StateType;
        typedef VecM<T, us> ContrlType;
        typedef VecM<T, ys> OutputType;

        typedef MatMN<T, xs, xs> StateMapType;
        typedef MatMN<T, xs, us> contrlMapType;
        typedef MatMN<T, ys, xs> OutputMapType;
        typedef MatMN<T, ys, us> DirectMapType;

        typedef VecM<int, 4> CtactStatusType;

    public:
        Model(){            
            contact_status_T.setOnes();
            foothold_location_T.setZero();            
        }             

        void set_footStepPlanner(MHPCFootStep<T>* footStepPlanner_in)
        {
            footStepPlanner = footStepPlanner_in;
        }        

        void dynamics(StateType &xnext, OutputType &y,
                      StateType &x, ContrlType &u, float t, const T &dt)
        {
            (void)(y);                        
            dynamics_continuousTime(xdot_, x, u, t);
            xnext = x + xdot_ * dt;                             
        }

        void dynamics_continuousTime(StateType& xdot, StateType& x, ContrlType&u, float t)
        {            
            footStepPlanner->updateFootPositions(t);            
            footStepPlanner->getContactStatus(contact_status);
            footStepPlanner->getFootPositions(foothold_location_T);

            xdot.setZero();           
            contact_status_T = contact_status.cast<T>();
            std::vector<T *> arg = {x.data(), u.data(), foothold_location_T.data(), contact_status_T.data()};
            std::vector<T *> res = {xdot.data()};
            casadi_interface(arg, res, xdot.size(), SRBDynamics,
                             SRBDynamics_sparsity_out,
                             SRBDynamics_work);
        }        

        void dynamics_partial(StateMapType &A, contrlMapType &B, OutputMapType &C, DirectMapType &D,
                              StateType &x, ContrlType &u, float t, const T &dt)
        {
            (void)(C);
            (void)(D);
            dynamics_partial_continuousTime(Ac_, Bc_, x, u, t);
                        
            A = StateMapType::Identity() + Ac_*dt;
            B = Bc_ * dt;                             
        }

        void dynamics_partial_continuousTime(StateMapType &Ac, contrlMapType &Bc,
                              StateType &x, ContrlType &u, float t)
        {
            Ac.setZero();
            Bc.setZero();            

            footStepPlanner->updateFootPositions(t);
            footStepPlanner->getContactStatus(contact_status);
            footStepPlanner->getFootPositions(foothold_location_T);

            contact_status_T = contact_status.cast<T>();                
            std::vector<T *> arg = {x.data(), u.data(), foothold_location_T.data(), contact_status_T.data()};
            std::vector<T *> res = {Ac.data(), Bc.data()};
            casadi_interface(arg, res, Ac.size(), SRBDynamicsDerivatives,
                             SRBDynamicsDerivatives_sparsity_out,
                             SRBDynamicsDerivatives_work);
        }

    private:      
        MHPCFootStep<T>* footStepPlanner = nullptr;
        Vec4<T> contact_status_T;
        CtactStatusType contact_status;
        VecM<T, 12> foothold_location_T;     
        StateType xdot_;
        StateMapType Ac_;
        contrlMapType Bc_;
    };

}

#endif // SRBM_H