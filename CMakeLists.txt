cmake_minimum_required(VERSION 2.8.3)
project(acoustic_locator)
set (CMAKE_CXX_STANDARD 11)
add_compile_options("-fpermissive")

#find_package(Qt4 REQUIRED COMPONENTS QtCore QtGui)
#include(${QT_USE_FILE})

include_directories(${CMAKE_CURRENT_BINARY_DIR})

#QT4_WRAP_UI(QT_FORMS_HPP ${QT_FORMS})
#QT4_WRAP_CPP(QT_MOC_HPP ${QT_MOC})

include_directories(include)
include_directories(
	${QT_INCLUDE_DIRS}
)
link_directories(/usr/local/lib/libliquid.so)
link_directories(/usr/local/lib/libbcm2835.a)
link_directories(/usr/local/lib)
include_directories(dsp/include)

add_library(dsp_lib dsp/include/dsp.h dsp/src/dsp.cpp )
add_executable(simulate dsp/src/dsp.cpp)
target_link_libraries(simulate liquid bcm2835 dsp_lib )
