project(Ohtsuka)
cmake_minimum_required(VERSION 3.16)

set(Eigen_DIR "/usr/include/eigen3")
set(Matplotlib_DIR "/usr/include/python3.8")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCE_FILES
  # robot_viewer.cpp
  # eigen_test.cpp
  main.cpp
  )

foreach(TARGET_FILE ${SOURCE_FILES})
  # common
  get_filename_component(TARGET_NAME ${TARGET_FILE} NAME_WE)
  add_executable(${TARGET_NAME} ${TARGET_FILE} traj_opt.cpp)
  target_compile_options(${TARGET_NAME} PRIVATE -O2 -fdiagnostics-color)
  target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(${TARGET_NAME} PRIVATE -lpython3.8 ${TORCH_LIBRARIES}) 
  
  # Eigen
  target_include_directories(${TARGET_NAME} PRIVATE ${Eigen_DIR})
  # Matplotlib 
  target_include_directories(${TARGET_NAME} PRIVATE ${Matplotlib_DIR})

endforeach()