option(DEBUG_DDP OFF)
option(DEBUG_MHPC OFF)
option(TIME_BENCHMARK_DDP OFF)

include_directories("../Reference")

# Add subdirectories
add_subdirectory(MHPC-Trajopt)

# Add *.cpp files to the executable file
add_executable(mhpc_run mpc_main.cpp MHPCLocomotion.cpp)

# Add target specific include include_directories
target_include_directories(mhpc_run PUBLIC "MHPC-Trajopt")

# Link executable file to other libraries
target_link_libraries(mhpc_run mhpc-trajopt pthread lcm)

if (DEBUG_MHPC)
    message("DEBUG_MODE is set for MHPC")
    target_compile_definitions(mhpc_run PUBLIC DEBUG_MODE)
endif()