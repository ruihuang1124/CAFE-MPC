#pragma once
#ifndef MHPCRESET_H
#define MHPCRESET_H

#include "WBM.h"
#include "SRBM.h"
#include "MHPCFootStep.h"

enum class ModelType {WB, SRB };

template<typename T>
class MHPCReset
{
public:
    typedef VecM<int, 4> CtactStatusType;

    MHPCReset(WBM::Model<T>* wbm_ptr_in,
              MHPCFootStep<T>* footStepPlanner_ptr_in) :
        wbm_ptr(wbm_ptr_in),
        footStepPlanner_ptr(footStepPlanner_ptr_in) {
            I_wb.setIdentity();
            I_srb.setIdentity();
            StateProjection.setZero();
            StateProjection.template topLeftCorner<6, 6>().setIdentity();
            StateProjection.template block<6,6>(6,18).setIdentity();
        }
    // void reset_map(DVec<T>& x_next, DVec<T>& y, 
    //                DVec<T>& x,
    //                const VecM<int, 4>& contact_cur,
    //                const VecM<int, 4>& contact_next,
    //                const ModelType& mtype_cur,
    //                const ModelType& mtype_next){                   
    //     // To be implemented
    //     printf("This function is not yet implemented \n");
    //     printf("Use the function that not producing the output for now \n");
    //     }

    void reset_map(DVec<T>& x_next, 
                   DVec<T>& x,
                   const VecM<int, 4>& contact_cur,
                   const VecM<int, 4>& contact_next,
                   const ModelType& mtype_cur,
                   const ModelType& mtype_next);                   

    // void reset_map_partial(DMat<T>& df_dx, DMat<T>& dy_dx,
    //                        DVec<T>& x,
    //                        const VecM<int, 4>& contact_cur,
    //                        const VecM<int, 4>& contact_next,
    //                        const ModelType& mtype_cur,
    //                        const ModelType& mtype_next){                   
    //     // To be implemented
    //     printf("This function is not yet implemented \n");
    //     printf("Use the function that not producing the output for now \n");
    //     }

    void reset_map_partial(DMat<T>& df_dx, 
                           DVec<T>& x,
                           const VecM<int, 4>& contact_cur,
                           const VecM<int, 4>& contact_next,
                           const ModelType& mtype_cur,
                           const ModelType& mtype_next);                           

private:
    WBM::Model<T>* wbm_ptr;
    MHPCFootStep<T>* footStepPlanner_ptr;
    MatMN<T, WBM::xs, WBM::xs> I_wb;
    MatMN<T, SRBM::xs, SRBM::xs> I_srb;
    MatMN<T, SRBM::xs, WBM::xs> StateProjection;
};



#endif //MHPCRESET_H