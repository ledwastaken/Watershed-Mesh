cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(watershed-mesh)

find_package(VTK COMPONENTS 
  IOGeometry
  CommonColor
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersGeneral
  IOXML
  InteractionStyle
  InteractionWidgets
  RenderingAnnotation
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2

  CommonColor
  CommonComputationalGeometry
  CommonCore
  CommonDataModel
  CommonTransforms
  FiltersCore
  FiltersGeneral
  FiltersModeling
  FiltersSources
  InteractionStyle
  InteractionWidgets
  RenderingAnnotation
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "WatershedMesh: Unable to find the VTK build folder.")
endif()

add_executable(watershed-mesh
 src/main.cpp
)

target_link_libraries(watershed-mesh PRIVATE
  ${VTK_LIBRARIES}
)

target_compile_options(watershed-mesh PRIVATE
 -Wall -Wextra -Wold-style-cast -pedantic -std=c++20
)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
 target_compile_options(watershed-mesh PRIVATE
 -g -fsanitize=leak
 )
 target_link_options(watershed-mesh PRIVATE
 -fsanitize=leak
 )
endif()

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS watershed-mesh
  MODULES ${VTK_LIBRARIES}
)
