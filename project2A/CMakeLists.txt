cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project2A)
SET(VTK_DIR /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project2A MACOSX_BUNDLE project2A)

if(VTK_LIBRARIES)
  target_link_libraries(project2A ${VTK_LIBRARIES})
else()
  target_link_libraries(project2A vtkHybrid)
endif()

