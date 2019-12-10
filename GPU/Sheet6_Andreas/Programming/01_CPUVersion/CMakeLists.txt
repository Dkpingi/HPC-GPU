cmake_minimum_required(VERSION 3.8)

project(ReductionCPU)

include_directories($PROJECT_SOURCE_DIR)
add_executable(${PROJECT_NAME} "main.cpp" ${SOURCES})

# Disable in-source builds to avoid cluttering the source directory
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(${PROJECT_NAME} LANGUAGES CXX)



# executable and properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17 CXX_STANDARD_REQUIRED OFF CXX_EXTENSIONS OFF)

target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} )
