cmake_minimum_required(VERSION 3.13)
project(imageFilter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wextra")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()


set(DCMTK_ROOT "D:/tasks/cuda_learn/scipy-style-medical-image-filter-cpp/lib/dcmtk-install/")
find_package(DCMTK REQUIRED)  


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include  
    ${DCMTK_INCLUDE_DIRS}              
)


set(SOURCE_FILES src/main.cpp src/filterfuns.cpp)
add_executable(filterFuns ${SOURCE_FILES})


target_link_libraries(filterFuns
    ${DCMTK_LIBRARIES}  
)

target_link_libraries(filterFuns
    "D:/tasks/cuda_learn/scipy-style-medical-image-filter-cpp/lib/dcmtk-install/lib/libdcmimgle.a"
    "D:/tasks/cuda_learn/scipy-style-medical-image-filter-cpp/lib/libdcmdata.a"
    "D:/tasks/cuda_learn/scipy-style-medical-image-filter-cpp/lib/liboflog.a"
    "D:/tasks/cuda_learn/scipy-style-medical-image-filter-cpp/lib/libofstd.a"
    stdc++fs 
)

if(MINGW)
    target_link_libraries(filterFuns stdc++fs)
endif()



set_target_properties(filterFuns PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${DCMTK_LIBRARY_DIRS}"
    BUILD_RPATH "${DCMTK_LIBRARY_DIRS}"
)
