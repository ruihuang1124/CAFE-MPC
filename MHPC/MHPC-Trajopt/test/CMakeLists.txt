include_directories("../"
                    "../CasadiGen/header"                                                        
                    "../../../Reference")

# # Test of KKT dynamics
# add_executable(testKKTDynamics testKKTDynamics.cpp 
#                                ../PinocchioInteface.cpp
#                                ../WBM.cpp)

# target_link_libraries(testKKTDynamics PUBLIC 
#                                       pinocchio::pinocchio 
#                                       casadigen_mhpc
#                                       hsddp)                                

# # Test of SRB dynamics
# add_executable(testSRBM testSRBDynamics.cpp 
#                         ../../../Reference/QuadReference.cpp
#                         ../PinocchioInteface.cpp
#                         ../WBM.cpp)


# target_link_libraries(testSRBM PUBLIC 
#                                casadigen_mhpc
#                                hsddp
#                                pinocchio::pinocchio)      
                               
# # Test of MHPCConstraint
# add_executable(testConstraint testMHPCConstraint.cpp
#                               ../PinocchioInteface.cpp                              
#                               ../WBM.cpp
#                               ../MHPCConstraint.cpp)

# target_link_libraries(testConstraint PUBLIC 
#                       pinocchio::pinocchio 
#                       casadigen_mhpc
#                       hsddp
#                       lcm)                               

# # Test of MHPCReset
# add_executable(testReset testMHPCReset.cpp
#                         ../MHPCReset.cpp
#                         ../PinocchioInteface.cpp                              
#                         ../WBM.cpp
#                         ../../../Reference/QuadReference.cpp)    

# target_link_libraries(testReset PUBLIC pinocchio::pinocchio casadigen_mhpc hsddp)                                          


# # Test of MHPCReference
# add_executable(testReference testMHPCReference.cpp  
#                         ../MHPCReference.cpp                                           
#                         ../../../Reference/QuadReference.cpp)   

# target_link_libraries(testReference PUBLIC hsddp pinocchio::pinocchio)      

# # Test of MHPCCost
# add_executable(testCost testMHPCCost.cpp
#                         ../MHPCReference.cpp                                           
#                         ../../../Reference/QuadReference.cpp)

# target_link_libraries(testCost PUBLIC hsddp pinocchio::pinocchio) 

# Test of SRBPlanning Problem
file(GLOB sources ../*.cpp)
add_executable(testSRBPlanning testSRBPlanning.cpp ${sources}
                ../../../Reference/QuadReference.cpp)
target_link_libraries(testSRBPlanning PUBLIC hsddp pinocchio::pinocchio casadigen_mhpc)   

# Test of MHPCProblem
file(GLOB sources ../*.cpp)
add_executable(testMHPCProblem testMHPCProblem.cpp ${sources}
                ../../../Reference/QuadReference.cpp)
target_link_libraries(testMHPCProblem PUBLIC hsddp pinocchio::pinocchio casadigen_mhpc)                

# Test of MHPCProblem2
file(GLOB sources ../*.cpp)
add_executable(testTrajOptInLoop testTrajOptInLoop.cpp ${sources}
                ../../../Reference/QuadReference.cpp)
target_link_libraries(testTrajOptInLoop PUBLIC hsddp pinocchio::pinocchio casadigen_mhpc)                

# Test of trajectory optimization
file(GLOB sources ../*.cpp)
add_executable(testBarrelRollTO testBarrelRollTO.cpp ${sources}
                ../../../Reference/QuadReference.cpp)
target_link_libraries(testBarrelRollTO PUBLIC hsddp pinocchio::pinocchio casadigen_mhpc)  