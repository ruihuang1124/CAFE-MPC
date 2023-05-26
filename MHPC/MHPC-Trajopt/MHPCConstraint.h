#ifndef MHPC_CONSTRAINTS_H
#define MHPC_CONSTRAINTS_H

#include "HSDDP_CPPTypes.h"
#include "HSDDP_Utils.h"
#include "WBM.h"
#include "SRBM.h"
#include "ConstraintsBase.h"

namespace MHPCConstraints
{
    template <typename T>
    class WBGRF : public PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>
    {
    private:
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::State;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Contrl;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Output;

        T mu_fric = .7; // Default Columb friction coefficient
        DMat<T> A;
        DVec<T> b;

    public:
        WBGRF(const VecM<int, 4> &ctact_);

        void set_friction_coefficient(T mu_fric_in)
        {
            mu_fric = mu_fric_in;
        }

        void compute_violation(const State &, const Contrl &, const Output &, int k) override;

        void compute_partial(const State &, const Contrl &, const Output &, int k) override;

    public:
        VecM<int, 4> ctact_status;
    };

    template <typename T>
    class TorqueLimit : public PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>
    {
    private:
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::State;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Contrl;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Output;

        MatMN<T, 2*WBM::us, WBM::us> C;
        VecM<T, 2*WBM::us> b;
        const T torque_limit;

    public:
        TorqueLimit();
        
        void compute_violation(const State &, const Contrl &, const Output &, int k) override;

        void compute_partial(const State &, const Contrl &, const Output &, int k) override;
    };

    template <typename T>
    class JointLimit : public PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>
    {
    private:
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::State;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Contrl;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Output;

        MatMN<T, 2*WBM::nu, WBM::nu> C;
        VecM<T, 2*WBM::nu> b;

    public:
        JointLimit();
        
        void compute_violation(const State &, const Contrl &, const Output &, int k) override;

        void compute_partial(const State &, const Contrl &, const Output &, int k) override;
    };

    template <typename T>
    class WBTouchDown : public TerminalConstraintBase<T, WBM::xs>
    {
    public:
        typedef VecM<int, 4> CtactStatusType;        

    private:
        using typename TerminalConstraintBase<T, WBM::xs>::State;                
        
    public:
        WBTouchDown(const CtactStatusType & TDStatus, WBM::Model<T>*wbm_ptr_in);

        void update_ground_height(T ground_height_in)
        {
            ground_height = ground_height_in;
        }

        void compute_violation(const State &);

        void compute_partial(const State &);

    private:        
        T ground_height;       
        CtactStatusType TDStatus;
        WBM::Model<T>* wbm_ptr;

        DVec<T> TD_EE_heights;
        DMat<T> TD_EE_Jacobians;                
    };
}


#endif // MHPC_CONSTRAINT_H