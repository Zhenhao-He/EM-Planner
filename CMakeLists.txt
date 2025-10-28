cmake_minimum_required(VERSION 3.10)
project(planner_experiment)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

# === Python3 with NumPy ===
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# === Eigen3 ===
find_package(Eigen3 REQUIRED)

# === OSQP ===
find_package(osqp REQUIRED)
set(MATPLOTLIBCPP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../matplotlib-cpp")
include_directories(${MATPLOTLIBCPP_INCLUDE_DIR})
# === matplotlib-cpp ===
# 你可以通过 apt 或手动安装:
# sudo apt install libmatplotlib-cpp-dev
# 如果没安装，就把 include 路径指向系统的 matplotlibcpp.hpp 所在目录
find_path(MATPLOTLIBCPP_INCLUDE_DIR matplotlibcpp.h)
if (NOT MATPLOTLIBCPP_INCLUDE_DIR)
  message(WARNING "matplotlibcpp.h not found, please install or specify path manually")
endif()

# === Include directories ===
include_directories(
  include
  scenario
  emplanner/include
  ${EIGEN3_INCLUDE_DIRS}
  ${MATPLOTLIBCPP_INCLUDE_DIR}
  ${Python3_INCLUDE_DIRS}
)

# === Source files ===
add_executable(${PROJECT_NAME}
  main.cpp
  scenario/src/sce3.cpp
  common/cartesian_frenet_conversion.cc
  common/discrete_points_math.cc
  emplanner/src/em_planner.cpp
  emplanner/src/sl_dp.cpp
  emplanner/src/sl_qp.cpp
  emplanner/src/st_dp.cpp
  emplanner/src/st_qp.cpp
)

# === Link libraries ===
target_link_libraries(${PROJECT_NAME}
  PRIVATE
  Python3::Python
  Python3::NumPy
  osqp::osqp
)

# === Optional message ===
message(STATUS "Python3 include dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "Eigen3 include dir: ${EIGEN3_INCLUDE_DIRS}")
message(STATUS "matplotlibcpp include dir: ${MATPLOTLIBCPP_INCLUDE_DIR}")
