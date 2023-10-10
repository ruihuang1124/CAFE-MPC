#ifndef BARRELROLL_CONSTRAINTS_H
#define BARRELROLL_CONSTRAINTS_H

#include "HSDDP_CPPTypes.h"
#include "HSDDP_Utils.h"
#include "WBM.h"
#include "ConstraintsBase.h"

namespace BarrelRoll
{
    template <typename T>
    class GRF : public PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>
    {
    private:
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::State;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Contrl;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Output;

        T mu_fric; // Default Columb friction coefficient
        DMat<T> A;
        DVec<T> b;

    public:
        GRF(const VecM<int, 4> &ctact_);

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
    class JointSpeedLimit : public PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>
    {
    private:
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::State;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Contrl;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Output;

        MatMN<T, 2*WBM::nu, WBM::nu> C;
        VecM<T, 2*WBM::nu> b;

        // Default lower and upper bound
        T lb_{-20.0};
        T ub_{20.0};

    public:
        JointSpeedLimit();

        void update_joint_limit(const T lb, const T ub) {lb_ = lb; ub_ = ub;}
        
        void compute_violation(const State &, const Contrl &, const Output &, int k) override;

        void compute_partial(const State &, const Contrl &, const Output &, int k) override;
    };

    template <typename T>
    class TouchDown : public TerminalConstraintBase<T, WBM::xs>
    {
    public:
        typedef VecM<int, 4> CtactStatusType;        

    private:
        using typename TerminalConstraintBase<T, WBM::xs>::State;                
        
    public:
        TouchDown(const CtactStatusType & TDStatus, WBM::Model<T>*wbm_ptr_in);

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

    template <typename T>
    class MinimumHeight : public PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>
    {
    private:
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::State;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Contrl;
        using typename PathConstraintBase<T, WBM::xs, WBM::us, WBM::ys>::Output;

        MatMN<T, 1, WBM::xs> C;        

        // Default lower and upper bound
        T h_min_{0.13};

    public:
        MinimumHeight(){
            this->update_constraint_size(1);           
        }

        void set_min_height(double h_min) {h_min_ = h_min;}
        
        void compute_violation(const State &, const Contrl &, const Output &, int k) override;

        void compute_partial(const State &, const Contrl &, const Output &, int k) override;
    };
}
#endif // BARRELROLL_CONSTRAINTS_H