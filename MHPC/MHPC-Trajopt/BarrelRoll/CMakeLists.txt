add_executable(br_to BarrelRollTO.cpp BarrelRollConstraints.cpp)

target_link_libraries(br_to PUBLIC mhpc-trajopt)      

target_include_directories(br_to PUBLIC ../)

