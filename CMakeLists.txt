CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(WBMatrix)

include_directories(include)

aux_source_directory(src/WBMatrix WBM_SRC)
aux_source_directory(src WBM_SRC)

add_library(WB_LIB ${WBM_SRC})

add_executable(LMC test/main.c)

target_link_libraries(LMC WB_LIB)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2")