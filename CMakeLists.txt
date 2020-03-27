project(fem LANGUAGES CXX)
cmake_minimum_required(VERSION 3.0)

add_library(femlib fem.cpp)
target_include_directories(femlib PUBLIC include)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE femlib)