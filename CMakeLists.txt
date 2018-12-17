cmake_minimum_required(VERSION 2.6)
project(muscle)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/shared/cmake)
include("${CMAKE_CURRENT_SOURCE_DIR}/shared/cmake/CMakeLists.txt")

add_executable(info "generate_info.cpp" include/move_patch.h include/PCA_elipse.h)
target_link_libraries(info igl::core igl::opengl igl::opengl_glfw)