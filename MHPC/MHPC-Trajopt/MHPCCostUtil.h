#ifndef MHPCCOST_UTILS_H
#define MHPCCOST_UTILS_H

#include <string>
#include <boost/property_tree/json_parser.hpp>

#include "MHPCCost.h"

template<typename T>
inline void loadCostWeights(const std::string&fileName,
                            shared_ptr<WBTrackingCost<T>> wb_tracking_cost,
                            shared_ptr<WBFootPlaceReg<T>> wb_foot_reg,
                            shared_ptr<SwingFootPosTracking<T>> swing_pos_tracking,
                            shared_ptr<SwingFootVelTracking<T>> swing_vel_tracking,
                            shared_ptr<SRBTrackingCost<T>> srb_tracking_cost)
{
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(fileName, pt);
	std::cout << "********* loading Cost weights from file **********\n" << fileName << "\n\n";
    
    /* Load weights for WB tracking */
    {
        std::vector<T> qw;
        for (auto& item : pt.get_child("WB_Tracking_Cost.qw_qB"))
        {            
            qw.push_back(item.second.get_value<T>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto& item : pt.get_child("WB_Tracking_Cost.qw_qJ"))
            {
                qw.push_back(item.second.get_value<T>());
            }
        }        
        for (auto& item : pt.get_child("WB_Tracking_Cost.qw_vB"))
        {
            qw.push_back(item.second.get_value<T>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto& item : pt.get_child("WB_Tracking_Cost.qw_vJ"))
            {
                qw.push_back(item.second.get_value<T>());
            }
        }    
        std::vector<T> rw;
        for (size_t i = 0; i < 12; i++)
        {
            rw.push_back(pt.get<T>("WB_Tracking_Cost.rw"));
        }
        std::vector<T> qfw;
        for (auto& item : pt.get_child("WB_Tracking_Cost.qfw_qB"))
        {
            qfw.push_back(item.second.get_value<T>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto& item : pt.get_child("WB_Tracking_Cost.qfw_qJ"))
            {
                qfw.push_back(item.second.get_value<T>());
            }
        }     
        for (auto& item : pt.get_child("WB_Tracking_Cost.qfw_vB"))
        {
            qfw.push_back(item.second.get_value<T>());
        }
        for (size_t i = 0; i < 4; i++)
        {
            for (auto& item : pt.get_child("WB_Tracking_Cost.qfw_vJ"))
            {
                qfw.push_back(item.second.get_value<T>());
            }
        }    
       
        std::vector<T> weights;
        weights.insert(weights.end(), std::move_iterator(qw.begin()), std::move_iterator(qw.end()));
        weights.insert(weights.end(), std::move_iterator(rw.begin()), std::move_iterator(rw.end()));
        weights.insert(weights.end(), std::move_iterator(qfw.begin()), std::move_iterator(qfw.end()));
        wb_tracking_cost->update_weighting_matrix(weights);        
    }             

    /* Load weights for SRB tracking */
    {
        std::vector<T> qw;
        for (auto& item : pt.get_child("SRB_Tracking_Cost.qw_qB"))
        {
            qw.push_back(item.second.get_value<T>());
        }       
        for (auto& item : pt.get_child("SRB_Tracking_Cost.qw_vB"))
        {
            qw.push_back(item.second.get_value<T>());
        }        
        std::vector<T> rw;
        for (size_t i = 0; i < 12; i++)
        {
            rw.push_back(pt.get<T>("SRB_Tracking_Cost.rw"));
        }
        std::vector<T> qfw;
        for (auto& item : pt.get_child("SRB_Tracking_Cost.qfw_qB"))
        {
            qfw.push_back(item.second.get_value<T>());
        }       
        for (auto& item : pt.get_child("SRB_Tracking_Cost.qfw_vB"))
        {
            qfw.push_back(item.second.get_value<T>());
        }   
        std::vector<T> weights;
        weights.insert(weights.end(), std::move_iterator(qw.begin()), std::move_iterator(qw.end()));
        weights.insert(weights.end(), std::move_iterator(rw.begin()), std::move_iterator(rw.end()));
        weights.insert(weights.end(), std::move_iterator(qfw.begin()), std::move_iterator(qfw.end()));
        srb_tracking_cost->update_weighting_matrix(weights);        
    }       

    /* Load weigths for foot placement regularization */
    {
        std::vector<T> qw;
        for (auto& item : pt.get_child("WB_FootPlace_Reg.qw_per_foot"))
        {
            qw.push_back(item.second.get_value<T>());
        }  
        wb_foot_reg->update_weighting_matrix(qw);        
    }

    /* Load weigths for swing foot position tracking */
    {
        std::vector<T> qw;
        for (auto& item : pt.get_child("Swing_Pos_Tracking.qw_per_foot"))
        {
            qw.push_back(item.second.get_value<T>());
        }  
        swing_pos_tracking->update_weighting_matrix(qw);        
    }

    /* Load weigths for swing foot velocity tracking */
    {
        std::vector<T> qw;
        for (auto& item : pt.get_child("Swing_Vel_Tracking.qw_per_foot"))
        {
            qw.push_back(item.second.get_value<T>());
        }  
        swing_vel_tracking->update_weighting_matrix(qw);        
    }
}

#endif