cmake_minimum_required(VERSION 2.6)

Project(ImageColorSpaceInspector)

# Where to copy executables when 'make install' is run
SET( INSTALL_DIR ${CMAKE_INSTALL_PREFIX} )

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

QT4_WRAP_UI(UISrcs MainWindow.ui)
QT4_WRAP_CPP(MOCSrcs MainWindow.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE( ${USE_VTK_FILE} )

FIND_PACKAGE(ITK REQUIRED)
INCLUDE( ${USE_ITK_FILE} )

ADD_EXECUTABLE(ImageColorSpaceInspector main.cpp MainWindow.cpp DisplayPoints.cpp Conversions.cpp RegionSelectionWidget.cpp Layer.cpp
Helpers.cpp ${MOCSrcs} ${UISrcs})
TARGET_LINK_LIBRARIES(ImageColorSpaceInspector ${VTK_LIBRARIES} QVTK ${QT_LIBRARIES} ${ITK_LIBRARIES})
INSTALL( TARGETS ImageColorSpaceInspector RUNTIME DESTINATION ${INSTALL_DIR} )

# ADD_EXECUTABLE(Test Test.cpp Conversions.cpp)
# TARGET_LINK_LIBRARIES(Test ${VTK_LIBRARIES})
