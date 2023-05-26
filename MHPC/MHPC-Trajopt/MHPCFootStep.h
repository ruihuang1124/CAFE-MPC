#pragma once
#ifndef MHPCFOOTSTEP_H
#define MHPCFOOTSTEP_H

#include "WBM.h"
#include "QuadReference.h"
#include <cassert>

template <typename T>
class MHPCFootStep
{
public:
    typedef Vec4<int> ContactStatusType;

public:
    MHPCFootStep(WBM::Model<T> *wbm_ptr_in,
                 QuadReference *quad_reference_in):
                 wbm_ptr(wbm_ptr_in),
                 quad_reference_ptr(quad_reference_in),
                 ground_height(0) {}

    /*
        @brief: Update the foot step locations on the XY plane at transition
                Called once at the model transition time
    */
    void updateFootPosAtTransition(const VecM<T, WBM::xs> &x)
    {
        wbm_ptr->get_footXYPositions(footXYPositions_transition, x);
        useFootPosAtTransition.setOnes();
    }

    /*
        @brief: Update the foot step locations on the XY plan for SRB model
                Determined by whole-body kinematics for contact foot at transition
                Determined by reference foot placements for future contact foot
    */
    void updateFootPositions(float t)
    {        
        quad_astate = quad_reference_ptr->get_a_reference_ptr_at_t(t);                

        for (size_t i = 0; i < 4; i++)
        {
            footPositions[i][2] = ground_height;            

            if (useFootPosAtTransition[i] &&
                quad_astate->contact[i] > 0)
            {
                useFootPosAtTransition[i] = true;
                footPositions[i].template head<2>() = footXYPositions_transition.template segment<2>(2 * i);
            }
            else
            {
                useFootPosAtTransition[i] = false;
                footPositions[i].template head<2>() = quad_astate->foot_placements.segment<2>(3 * i).cast<T>();
            }
        }
    }

    void getFootPositions(VecM<T, 12>& EE_pos)
    {
        // EE_pos << footPositions[0], footPositions[1], footPositions[2], footPositions[3];

        // For now, we consider foot placements are given by the reference
        EE_pos = quad_astate->foot_placements.cast<T>();
    }

    void getContactStatus(ContactStatusType& contactStatus)
    {
        contactStatus = quad_astate->contact;
    }

private:
    WBM::Model<T> *wbm_ptr = nullptr;
    QuadReference *quad_reference_ptr = nullptr;
    QuadAugmentedState *quad_astate = nullptr;

    VecM<T, 8> footXYPositions_transition;
    VecM<T, 3> footPositions[4];

    ContactStatusType cStatus_cur;
    Vec4<bool> useFootPosAtTransition;

    float ground_height;
};

#endif // MHPCFOOTSTEP_H