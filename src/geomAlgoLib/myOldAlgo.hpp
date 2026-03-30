/**
 * @file example.hpp
 * @brief Mesh segmentation and analysis library header
 * 
 * Public interface for geometric computation on CGAL polyhedron meshes:
 *  - Surface properties (area, genus, normals)
 *  - Face segmentation by geometric criteria
 *  - Colored OFF export for visualization
 *
 * @author ISIMA ZZ3 - TP3 Segmentation
 * @date 2026-02-16
 */

#pragma once

#include "types.hpp"

#include <array>
#include <map>
#include <string>

namespace geomAlgoLib
{

//==============================================================================
// SECTION 1: DATA STRUCTURES
//==============================================================================

/**
 * @brief RGB color representation with utility operations.
 */
struct Color {
    double r;  ///< Red component [0, 1]
    double g;  ///< Green component [0, 1]
    double b;  ///< Blue component [0, 1]

    /// Convert to std::array for serialization
    std::array<double, 3> to_array() const { return {r, g, b}; }
    
    /// Scalar multiplication for brightness adjustment
    std::array<double, 3> operator*(const double scalar) const { 
        return {r * scalar, g * scalar, b * scalar}; 
    }
    
    /// Color addition for blending
    const Color operator+(std::array<double, 3> other) const { 
        return {r + other[0], g + other[1], b + other[2]}; 
    }
};

/**
 * @brief Container for multiple boolean tags and associated colors per face.
 * 
 * Allows a single face to carry multiple classification tags simultaneously
 * (e.g., "Grande" AND "Plane"), each with its own visualization color.
 */
struct TagMaps
{
    std::map<std::string, bool> tags;    ///< Tag name -> presence flag
    std::map<std::string, Color> colors; ///< Tag name -> display color
};

//==============================================================================
// SECTION 2: COLOR CONSTANTS
//==============================================================================

inline Color RED     = {1., 0., 0.};
inline Color GREEN   = {0., 1., 0.};
inline Color BLUE    = {0., 0., 1.};
inline Color YELLOW  = {1., 1., 0.};
inline Color CYAN    = {0., 1., 1.};
inline Color MAGENTA = {1., 0., 1.};
inline Color WHITE   = {1., 1., 1.};
inline Color BLACK   = {0., 0., 0.};

/// String-to-color lookup for tag color assignment
inline std::map<std::string, Color> colorDictionary = {
    {"RED",     RED},
    {"GREEN",   GREEN},
    {"BLUE",    BLUE},
    {"YELLOW",  YELLOW},
    {"CYAN",    CYAN},
    {"MAGENTA", MAGENTA},
    {"WHITE",   WHITE},
    {"BLACK",   BLACK}
};

//==============================================================================
// SECTION 3: GEOMETRIC UTILITIES - Angles and Normals
//==============================================================================

/**
 * @brief Compute angle in degrees between two 3D vectors.
 * @param vec First vector
 * @param axis Second vector (reference axis)
 * @return Angle in degrees [0, 180]
 */
double angle_degrees(const CGAL::Vector_3<Kernel>& vec, 
                     const CGAL::Vector_3<Kernel>& axis);

/**
 * @brief Compute angle between triangle normal and reference vector.
 * @param p First triangle vertex
 * @param q Second triangle vertex
 * @param r Third triangle vertex
 * @param reference_vector Axis to measure angle against
 * @return Angle in degrees
 */
double get_normal_angle(const CGAL::Point_3<Kernel> & p,
                        const CGAL::Point_3<Kernel> & q,
                        const CGAL::Point_3<Kernel> & r,
                        const CGAL::Vector_3<Kernel> & reference_vector);

/**
 * @brief Compute normal vector of a triangle.
 * @param p First triangle vertex
 * @param q Second triangle vertex
 * @param r Third triangle vertex
 * @return Unit normal vector
 */
const CGAL::Vector_3<Kernel> get_normal_vector(const CGAL::Point_3<Kernel> & p,
                                                const CGAL::Point_3<Kernel> & q,
                                                const CGAL::Point_3<Kernel> & r);

/**
 * @brief Compute angles between face normal and cardinal axes (X, Y, Z).
 * @param i Face iterator
 * @return Array of three angles [angle_X, angle_Y, angle_Z] in degrees
 */
std::array<double, 3> get_facet_normals_angle(const FacetCstIt &i);

/**
 * @brief Extract normal vector of a triangular face.
 * @param i Face iterator
 * @return Face normal vector
 */
CGAL::Vector_3<Kernel> get_facet_normals_vector(const FacetCstIt &i);

//==============================================================================
// SECTION 4: AREA COMPUTATION
//==============================================================================

/**
 * @brief Compute area of a triangular face using CGAL primitives.
 * @param i Face iterator (must point to a triangle)
 * @return Face area in mesh units squared
 */
double compute_area_triangle(FacetCstIt &i);

/**
 * @brief Compute area of a quadrilateral face by triangulation.
 * @param i Face iterator (must point to a quad)
 * @return Face area (sum of two triangles)
 * @note Assumes convex quad; splits along diagonal p0-p2
 */
double compute_area_quad(FacetCstIt &i);

/**
 * @brief Compute total area of a purely triangular mesh.
 * @param mesh Input mesh (assumes all faces are triangles)
 * @return Total mesh surface area
 */
double compute_Area(const Mesh &mesh);

/**
 * @brief Compute total area of a mesh with mixed triangle/quad faces.
 * @param mesh Input mesh
 * @return Total mesh surface area
 */
double compute_Area_Quad_Mesh(const Mesh &mesh);

/**
 * @brief Build map associating each face with its area.
 * @param mesh Input mesh (triangle/quad faces supported)
 * @return Map: face_handle -> area
 */
std::map<Mesh::Facet_const_handle, double> compute_Area_Map(const Mesh &mesh);

//==============================================================================
// SECTION 5: MESH PROPERTIES
//==============================================================================

/**
 * @brief Compute genus of a closed mesh using Euler characteristic.
 * @param mesh Input mesh (must be closed, single-component)
 * @return Genus g, where χ = V - E + F = 2 - 2g
 * @note Prints vertex/edge/face counts to stdout for verification
 */
int computeGenus(const Mesh &mesh);

/**
 * @brief Colorize mesh by dominant normal component and export normal map.
 * @param mesh Input mesh
 * @param file_name Output colored OFF file path
 * @return Map: face_handle -> normal_vector (for further processing)
 * @note Produces visual representation of face orientation
 */
std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel>> 
colorize_by_greatest_normal(const Mesh &mesh, std::string file_name);

//==============================================================================
// SECTION 6: DISTANCE COMPUTATION
//==============================================================================

/**
 * @brief Compute Euclidean distance from face centroid to reference point.
 * @param mesh Parent mesh
 * @param i Face iterator
 * @param p Reference point (e.g., robot position)
 * @return Array {dx, dy, dz, d_total} where d_total = sqrt(dx²+dy²+dz²)
 * @note Returns signed component distances for axis-specific analysis
 */
std::array<double, 4> compute_distance_face_point(const Mesh &mesh, 
                                                   const FacetCstIt &i, 
                                                   const geomAlgoLib::Point_3 &p);

//==============================================================================
// SECTION 7: SEGMENTATION
//==============================================================================

/**
 * @brief Segment mesh faces based on geometric criteria.
 * 
 * Assigns multiple boolean tags to each face:
 *  - Grande: area > 5.0
 *  - OrienteHaut: normal dominated by Z component
 *  - Plane: angle with vertical < 15° (horizontal surface)
 *  - VerySteep: angle with vertical > 30° (steep/vertical surface)
 *  - ZoneMarche: plane surface (walkable for robot)
 *  - Obstacle: steep surface with non-negligible area
 *  - ZoneAPortee: plane surface within horizontal reach and height tolerance
 * 
 * @param areas Map of face -> area
 * @param normal_angles Map of face -> normal vector
 * @param mesh Parent mesh (for distance computations)
 * @param reference_point Robot/avatar position for reachability analysis
 * @return Map: face_handle -> TagMaps (tags + associated colors)
 */
std::map<Mesh::Facet_const_handle, TagMaps> segmentationMap(
    std::map<Mesh::Facet_const_handle, double> areas,
    std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel>> normal_angles,
    const Mesh &mesh,
    const geomAlgoLib::Point_3 &reference_point
);

//==============================================================================
// SECTION 8: I/O - Export Functions
//==============================================================================

/**
 * @brief Export single face with RGB color (utility function).
 * @param mesh Parent mesh
 * @param F Face handle
 * @param color RGB array [0,1] range
 * @param file_name Output COFF file path
 */
void save_Face_Color(const Mesh &mesh, 
                     Mesh::Facet_const_handle F, 
                     std::array<double, 3> color, 
                     std::string file_name);

/**
 * @brief Export mesh with per-face RGB colors (array format).
 * @param mesh Input mesh
 * @param m Map: face_handle -> RGB array
 * @param file_name Output COFF file path
 */
void save_Mesh_Color(const Mesh &mesh,
                     std::map<Mesh::Facet_const_handle, std::array<double, 3>> m,
                     std::string file_name);

/**
 * @brief Export mesh with normalized scalar-to-grayscale mapping.
 * @param mesh Input mesh
 * @param m Map: face_handle -> scalar value
 * @param file_name Output COFF file path
 * @note Normalizes values to [0,1] then maps to grayscale: white=min, black=max
 */
void save_Mesh_Color(const Mesh &mesh,
                     std::map<Mesh::Facet_const_handle, double> m,
                     std::string file_name);

//==============================================================================
// SECTION 9: HIGH-LEVEL PIPELINE
//==============================================================================

/**
 * @brief Complete segmentation pipeline for a single tag visualization.
 * 
 * Workflow:
 *  1. Compute area map
 *  2. Compute normal map
 *  3. Run segmentation (assign all tags)
 *  4. Extract color for requested tag
 *  5. Export colored COFF file
 * 
 * @param Tag Tag name to visualize ("Plane", "Obstacle", "ZoneAPortee", etc.)
 * @param mesh Input mesh
 * @param file_name Output COFF file path
 * @param reference_point Robot position for ZoneAPortee computation
 */
void ColorizedTag(std::string Tag,
                  const Mesh &mesh,
                  std::string file_name,
                  const geomAlgoLib::Point_3 &reference_point);

//==============================================================================
// SECTION 10: VISUALIZATION UTILITIES
//==============================================================================

/**
 * @brief Generate cube marker for robot position visualization.
 * @param robot_pos 3D position of robot/avatar
 * @param file_name Output COFF file (red cube)
 * @param size Cube edge length in mesh units
 * @note Load this marker alongside mesh in MeshLab for visual reference
 */
void save_robot_marker(const geomAlgoLib::Point_3 &robot_pos,
                       std::string file_name,
                       double size);

//==============================================================================
// SECTION 11: DEBUG UTILITIES
//==============================================================================

/**
 * @brief Print map contents to stdout (debugging helper).
 * @param comment Prefix label for output
 * @param m Map to display
 */
void print_map(std::string_view comment, 
               std::map<geomAlgoLib::Mesh::Facet_const_handle, double>& m);

} // namespace geomAlgoLib
