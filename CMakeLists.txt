cmake_minimum_required(VERSION 3.26)
project(optimisation)

set(CMAKE_CXX_STANDARD 20)

add_executable(optimisation main.cpp
        matrix.h
        matrix.cpp)
