cmake_minimum_required(VERSION 3.10)
project(bioinf)

set(CMAKE_CXX_STANDARD 17)

set(OpenCV_DIR "C:\\opencv\\mingw-build\\install")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(bioinf src/main.cpp include/xmlParser.h src/xmlParser.cpp include/tsvParser.h include/massGeneration.h src/massGeneration.cpp src/tsvParser.cpp include/spectrumAnalyzer.h src/spectrumAnalyzer.cpp)
set(OpenCV_LIBS opencv_core opencv_imgproc opencv_highgui opencv_imgcodecs)
target_link_libraries( bioinf ${OpenCV_LIBS} )