cmake_minimum_required(VERSION 3.30)
project(PG_2)

set(CMAKE_CXX_STANDARD 20)

add_compile_options(-pedantic -Wall -Wno-long-long -O2 -c -Werror)

add_executable(PG_2 main_time.cpp)
#add_executable(PG_2 main.cpp)
