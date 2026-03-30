# Sources
set(SOURCES
    main.cpp
    geomAlgoLib/io.cpp
    geomAlgoLib/measures.cpp
    geomAlgoLib/myAlgo.cpp
    geomAlgoLib/myOldAlgo.cpp
    geomAlgoLib/segmentation.cpp
    geomAlgoLib/example_archive_uncommented.cpp
    geomAlgoLib/batch_renderer.hpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

# Links GLFW + libigl
target_link_libraries(${PROJECT_NAME} 
    CGAL::CGAL
    Eigen3::Eigen 
    glfw 
    OpenGL::GL 
    Boost::headers
    "-framework Cocoa -framework IOKit -framework CoreVideo"
)

target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE 
    "${CMAKE_SOURCE_DIR}/external/libigl/include"
    "${CMAKE_SOURCE_DIR}/external/stb"
)


# libigl
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/external/libigl/include)
