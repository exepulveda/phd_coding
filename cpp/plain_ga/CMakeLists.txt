add_library(plain_ga 
        src/Random.cpp 
        src/Population.cpp 
    )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")

target_include_directories (plain_ga PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
