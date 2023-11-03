add_executable(loco_to Loco_TO.cpp LocoProblem.cpp)

target_link_libraries(loco_to PUBLIC mhpc-trajopt)      

target_include_directories(loco_to PUBLIC ../)

