cmake_minimum_required(VERSION 3.26)
project(Projeto2)

set(CMAKE_CXX_STANDARD 17)

add_executable(Projeto2 main.cpp
        Graph.h
        Logic/LoadingFunctions.cpp
        Interface/UI.h
        Interface/UI.cpp
        Logic/Algorithms.cpp
        Logic/Algorithms.h
        Logic/Clock.h
)
