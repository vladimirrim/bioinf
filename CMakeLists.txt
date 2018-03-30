cmake_minimum_required(VERSION 3.10)
project(bioinf)

set(CMAKE_CXX_STANDARD 17)

add_executable(bioinf src/main.cpp include/xmlParser.h src/xmlParser.cpp include/tsvParser.h include/massGeneration.h src/massGeneration.cpp src/tsvParser.cpp)