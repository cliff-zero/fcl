cmake_minimum_required(VERSION 3.12)
project(qlat.x)

set(CMAKE_CXX_STANDARD 17)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
endif()


    #set(eigen_path /thfs1/home/fengxu/zywang/local/include/eigen3)
    #set(mpi_prefix /thfs1/home/fengxu/zywang/software/mpich-3.4.1)
    set(qlat_prefix /thfs1/home/fengxu/Yusheng_Gao/software/Qlat_0_58)
    set(QLAT_INCLUDE ${qlat_prefix}/include)
    set(QLAT_LIB ${qlat_prefix}/lib)
    set(EIGEN_INCLUDE ${qlat_prefix}/eigen3)
    #set(CMAKE_CXX_COMPILER /vol8/home/fengxu1/zywang/software/gcc10/bin/g++)
    set(CMAKE_CXX_COMPILER mpic++)
    set(CMAKE_CXX_FLAGS "-fopenmp -lfftw3 -lfftw3f -lmpi -Wall -Wno-unused-variable")
    #add_definitions(-DNO_ALIGNED_ALLOC)
    include_directories(
            ${QLAT_INCLUDE}
            ${EIGEN_INCLUDE}
	    #${QLAT_INCLUDE}/eigen3
	    #${QLAT_INCLUDE}/qutils
	    #${mpi_prefix}/include
    )
    link_directories(${QLAT_LIB}) #/vol8/home/fengxu1/zywang/software/gcc10/lib64

add_executable(qlat.x main.cpp)

target_link_libraries(qlat.x ${QLAT_LIB}/libqlat.so ${QLAT_LIB}/libqutils.so)
