
find_package(pinocchio REQUIRED)

add_subdirectory(CasadiGen)

file(GLOB sources "*.cpp"
                  "../../Reference/*.cpp")                           

add_library(mhpc-trajopt SHARED ${sources})

target_include_directories(mhpc-trajopt PUBLIC ../../Reference)

target_link_libraries(mhpc-trajopt PUBLIC casadigen_mhpc hsddp pinocchio::pinocchio)

add_subdirectory(test)