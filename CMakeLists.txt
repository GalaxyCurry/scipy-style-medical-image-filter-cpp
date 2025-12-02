cmake_minimum_required(VERSION 3.13)
project(imageFilter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wextra")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()


set(DCMTK_ROOT "E:/vscode/itk/itk-prefix/")
find_package(ITK REQUIRED)  


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include 
    ${DCMTK_INCLUDE_DIRS}              
)


set(SOURCE_FILES src/main.cpp src/filterfuns.cpp)
add_executable(filterFuns ${SOURCE_FILES})


target_link_libraries(filterFuns
    ${ITK_LIBRARIES}  
)

if(MINGW)
    target_link_libraries(filterFuns stdc++fs)
endif()



set_target_properties(filterFuns PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${ITK_LIBRARY_DIRS}"
    BUILD_RPATH "${ITK_LIBRARY_DIRS}"
)
