add_library(geomAlgoLib
    #example.hpp
    #example.cpp
    myAlgo.hpp
    myAlgo.cpp
    io.hpp
    io.cpp
    measures.hpp
    measures.cpp
    types.hpp
)

# Link the library with CGAL
target_link_libraries(geomAlgoLib PUBLIC CGAL::CGAL)

target_include_directories(geomAlgoLib
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR} 
)