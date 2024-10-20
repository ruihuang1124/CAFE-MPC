file(GLOB_RECURSE sources "source/*.cpp" "common/*.cpp")

add_library(hsddp SHARED ${sources})

# Optional debug mode 
option(DEBUG_DDP OFF)

find_package(Boost REQUIRED)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(hsddp PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(hsddp PUBLIC
								 "common"
                                 "header"
                                 "/usr/local/include/eigen3")   
target_link_libraries(hsddp PUBLIC lcm)                                 

# If in debug mode, set precompiler to DEBUG_MODE
if (DEBUG_DDP)
    message("DEBUG_MODE is set for HSDDP")
    target_compile_definitions(hsddp PUBLIC DEBUG_MODE)
endif()

