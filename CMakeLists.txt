cmake_minimum_required(VERSION 3.2)
project(LJ_MC_Multithreaded)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")
set(Boost_INCLUDE_DIR /usr/local/include)
set(Boost_LIBRARY_DIR /usr/local/lib)
find_package(Boost COMPONENTS system filesystem thread REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

set(SOURCE_FILES main.cpp)
add_executable(LJ_MC_Multithreaded ${SOURCE_FILES})
target_link_libraries(LJ_MC_Multithreaded ${Boost_LIBRARIES} )
