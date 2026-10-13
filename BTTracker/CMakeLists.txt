cmake_minimum_required(VERSION 3.16)

project (code)

SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -g2 -ggdb")

SET(CMAKE_CXX_STANDARD 20)

find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})
include(${ROOT_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_library(DTOFBT STATIC ${sources} ${headers})
target_link_libraries(DTOFBT ${ROOT_LIBRARIES})

add_executable(ReadTracker "ReadTracker.C")
target_link_libraries(ReadTracker ${ROOT_LIBRARIES} DTOFBT)

file(COPY ${PROJECT_SOURCE_DIR}/config DESTINATION ${PROJECT_BINARY_DIR})
# Launch run.sh/run_batch.sh from the source tree; they resolve scripts/ there.


