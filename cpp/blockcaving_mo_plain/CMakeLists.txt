include_directories("../plain_ga/src")

add_executable(blockcaving_mo_plain max_npv.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(HDF5)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    set(_hdf5_libs hdf5 hdf5_cpp)
endif()

find_package(Boost REQUIRED COMPONENTS system program_options filesystem )
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries(blockcaving_mo_plain ${Boost_LIBRARIES})
endif()


target_link_libraries (blockcaving_mo_plain plain_ga blockcaving ${_hdf5_libs} jsoncpp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g ")

