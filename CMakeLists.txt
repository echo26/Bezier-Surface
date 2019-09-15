cake_minimum_required(VERSION 2.8.12.1)

PROJECT(bezier)
SET(VTK_DIR /Users/edward/VTK-8.1.2)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(bezier MACOSX_BUNDLE bezier)

if(VTK_LIBRARIES)
  target_link_libraries(bezier ${VTK_LIBRARIES})
else()
  target_link_libraries(bezier vtkHybrid)
endif()

