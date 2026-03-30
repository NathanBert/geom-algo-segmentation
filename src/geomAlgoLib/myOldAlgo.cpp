/**
 * @file example.cpp
 * @brief Mesh segmentation and analysis library for geometric algorithms
 * 
 * Provides computation of mesh properties (area, genus, normals),
 * face segmentation based on geometric criteria, and colored OFF export
 * for visualization in MeshLab or similar tools.
 *
 * @author ISIMA ZZ3 - TP3 Segmentation
 * @date 2026-02-16
 */

#include "myAlgo.hpp"

#include <iostream>
#include <math.h>

namespace geomAlgoLib
{

//==============================================================================
// SECTION 1: GEOMETRIC UTILITIES - Angles and Normals
//==============================================================================

/**
 * @brief Compute angle in degrees between two 3D vectors.
 * @param vec First vector
 * @param axis Second vector (reference axis)
 * @return Angle in degrees [0, 180]
 */
double angle_degrees(const CGAL::Vector_3<Kernel>& vec, const CGAL::Vector_3<Kernel>& axis) {
    double cos_theta = CGAL::scalar_product(vec, axis) /
                       std::sqrt(CGAL::squared_length(vec) * CGAL::squared_length(axis));
    return std::acos(std::clamp(cos_theta, -1.0, 1.0)) * 180.0 / M_PI;
}

/**
 * @brief Compute angle between triangle normal and reference vector.
 * @param p First triangle vertex
 * @param q Second triangle vertex
 * @param r Third triangle vertex
 * @param reference_vector Axis to measure angle against
 * @return Angle in degrees
 */
double get_normal_angle(const CGAL::Point_3< Kernel > &     p,
                        const CGAL::Point_3< Kernel > &   q,
                        const CGAL::Point_3< Kernel > &   r,
                        const CGAL::Vector_3<Kernel> & reference_vector)
{
    const CGAL::Vector_3<Kernel> normal = CGAL::normal(p,q,r);
    double angle = angle_degrees(normal, reference_vector);
    return angle;
}

/**
 * @brief Compute normal vector of a triangle.
 * @param p First triangle vertex
 * @param q Second triangle vertex
 * @param r Third triangle vertex
 * @return Unit normal vector
 */
const CGAL::Vector_3<Kernel> get_normal_vector(const CGAL::Point_3< Kernel > &  p,
                                                const CGAL::Point_3< Kernel > &   q,
                                                const CGAL::Point_3< Kernel > &   r)
{
    const CGAL::Vector_3<Kernel> normal = CGAL::normal(p,q,r);
    return normal;
}

/**
 * @brief Compute angles between face normal and cardinal axes (X, Y, Z).
 * @param i Face iterator
 * @return Array of three angles [angle_X, angle_Y, angle_Z] in degrees
 */
std::array<double,3>  get_facet_normals_angle(const FacetCstIt &i)
{
    CGAL::Vector_3<Kernel> x_vector = CGAL::Vector_3<Kernel>(1.,0.,0.);
    CGAL::Vector_3<Kernel> y_vector = CGAL::Vector_3<Kernel>(0.,1.,0.);
    CGAL::Vector_3<Kernel> z_vector = CGAL::Vector_3<Kernel>(0.,0.,1.);

    HalfedgeCstIt halfEdges[3];
    halfEdges[0] = i->halfedge();
    halfEdges[1] = i->halfedge()->next();
    halfEdges[2] = halfEdges[1]->next();

    auto p0 = halfEdges[0]->vertex_begin ()->vertex ()->point ();
    auto p1 = halfEdges[1]->vertex_begin ()->vertex ()->point ();
    auto p2 = halfEdges[2]->vertex_begin ()->vertex ()->point ();

    std::array<double,3 > arr = {
        get_normal_angle(p0,p1,p2,x_vector),
        get_normal_angle(p0,p1,p2,y_vector),
        get_normal_angle(p0,p1,p2,z_vector),
    };

    return arr;
}

/**
 * @brief Extract normal vector of a triangular face.
 * @param i Face iterator
 * @return Face normal vector
 */
CGAL::Vector_3<Kernel> get_facet_normals_vector(const FacetCstIt &i)
{
    HalfedgeCstIt halfEdges[3];
    halfEdges[0] = i->halfedge();
    halfEdges[1] = i->halfedge()->next();
    halfEdges[2] = halfEdges[1]->next();

    auto p0 = halfEdges[0]->vertex_begin ()->vertex ()->point ();
    auto p1 = halfEdges[1]->vertex_begin ()->vertex ()->point ();
    auto p2 = halfEdges[2]->vertex_begin ()->vertex ()->point ();

    CGAL::Vector_3<Kernel> arr = get_normal_vector(p0,p1,p2);
    return arr;
}

//==============================================================================
// SECTION 2: AREA COMPUTATION - Individual Faces
//==============================================================================

/**
 * @brief Compute area of a triangular face using CGAL primitives.
 * @param i Face iterator (must point to a triangle)
 * @return Face area in mesh units squared
 */
double compute_area_triangle(FacetCstIt &i)
{
    double this_area = 0;
    HalfedgeCstIt halfEdges[3];

    halfEdges[0] = i->halfedge();
    halfEdges[1] = i->halfedge()->next();
    halfEdges[2] = halfEdges[1]->next();
    this_area = 0;

    auto p0 = halfEdges[0]->vertex_begin ()->vertex ()->point ();
    auto p1 = halfEdges[1]->vertex_begin ()->vertex ()->point ();
    auto p2 = halfEdges[2]->vertex_begin ()->vertex ()->point ();

    this_area = std::sqrt(CGAL::squared_area(p0,p1,p2));
    return this_area;
}

/**
 * @brief Compute area of a quadrilateral face by triangulation.
 * @param i Face iterator (must point to a quad)
 * @return Face area (sum of two triangles)
 * @note Assumes convex quad; splits along diagonal p0-p2
 */
double compute_area_quad(FacetCstIt &i)
{
    double this_area = 0;
    HalfedgeCstIt halfEdges[4];

    halfEdges[0] = i->halfedge();
    halfEdges[1] = i->halfedge()->next();
    halfEdges[2] = halfEdges[1]->next();
    halfEdges[3] = halfEdges[2]->next();
    this_area = 0;

    auto p0 = halfEdges[0]->vertex_begin ()->vertex ()->point ();
    auto p1 = halfEdges[1]->vertex_begin ()->vertex ()->point ();
    auto p2 = halfEdges[2]->vertex_begin ()->vertex ()->point ();
    auto p3 = halfEdges[3]->vertex_begin ()->vertex ()->point ();

    this_area += std::sqrt(CGAL::squared_area(p0,p1,p2));
    this_area += std::sqrt(CGAL::squared_area(p2,p3,p0));
    return this_area;
}

//==============================================================================
// SECTION 3: MESH PROPERTY AGGREGATION
//==============================================================================

/**
 * @brief Compute total area of a purely triangular mesh.
 * @param mesh Input mesh (assumes all faces are triangles)
 * @return Total mesh surface area
 */
double compute_Area(const Mesh &mesh)
{
    double aera_face = 0;

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        aera_face += compute_area_triangle(i);
    }

    return aera_face;
}

/**
 * @brief Compute total area of a mesh with mixed triangle/quad faces.
 * @param mesh Input mesh
 * @return Total mesh surface area
 */
double compute_Area_Quad_Mesh(const Mesh &mesh)
{
    std::map<Mesh::Facet_const_handle, double> m;
    double aera_face = 0;
    HalfedgeCstIt halfEdges[4];

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        if (i->is_triangle())
        {
            aera_face += compute_area_triangle(i);
        }
        else if (i->is_quad())
        {
            aera_face += compute_area_quad(i);
        }
    }

    return aera_face;
}

/**
 * @brief Build map associating each face with its area.
 * @param mesh Input mesh (triangle/quad faces supported)
 * @return Map: face_handle -> area
 */
std::map<Mesh::Facet_const_handle, double> compute_Area_Map(const Mesh &mesh)
{
    std::map<Mesh::Facet_const_handle, double> m;
    double aera_face;
    HalfedgeCstIt halfEdges[4];

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        aera_face = 0;

        if (i->is_triangle())
        {
            aera_face += compute_area_triangle(i);
        }
        else if (i->is_quad())
        {
            aera_face += compute_area_quad(i);
        }
        m[i] = aera_face;
    }

    return m;
}

/**
 * @brief Colorize mesh by dominant normal component and export normal map.
 * @param mesh Input mesh
 * @param file_name Output colored OFF file path
 * @return Map: face_handle -> normal_vector (for further processing)
 * @note Produces visual representation of face orientation
 */
/*std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel> > colorize_by_greatest_normal(const Mesh &mesh, std::string file_name)
{
    std::vector<std::array<double,3>> normals;
    CGAL::Vector_3<Kernel> x_vector = CGAL::Vector_3<Kernel>(1.,0.,0.);
    CGAL::Vector_3<Kernel> y_vector = CGAL::Vector_3<Kernel>(0.,1.,0.);
    CGAL::Vector_3<Kernel> z_vector = CGAL::Vector_3<Kernel>(0.,0.,1.);
    std::map<Mesh::Facet_const_handle, std::array<double,3>> m;
    std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel> > normal_map;

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        std::array<double,3 > f = get_facet_normals_angle(i);
        CGAL::Vector_3<Kernel> v = get_facet_normals_vector(i);
        normal_map[i] = v;

        double x = v.x(); double y = v.y(); double z = v.z();
        x = abs(v[0])*cos(f[0]);
        y = abs(v[1])*cos(f[1]);
        z = abs(v[2])*cos(f[2]);

        double max_axis = std::max(std::max(fabs(x), fabs(y)), fabs(z));
        x = ((int)(x/max_axis)) * 255;
        y = ((int)(y/max_axis)) * 255;
        z = ((int)(z/max_axis)) * 255;

        std::array<double,3> rgb = {x,y,z};
        m[i] = rgb;
    }

    save_Mesh_Color(mesh,m,file_name);
    return(normal_map);
}*/

std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel>> colorize_by_greatest_normal(const Mesh &mesh, std::string file_name)
{
    std::map<Mesh::Facet_const_handle, std::array<double,3>> m;
    std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel>> normal_map;

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        CGAL::Vector_3<Kernel> v = get_facet_normals_vector(i);
        normal_map[i] = v;

        double abs_x = fabs(v.x());
        double abs_y = fabs(v.y());
        double abs_z = fabs(v.z());

        double max_component = std::max({abs_x, abs_y, abs_z});

        std::array<double, 3> rgb;
        
        if (abs_x == max_component) {
            rgb = {255.0, 0.0, 0.0};
        } else if (abs_y == max_component) {
            rgb = {0.0, 255.0, 0.0};
        } else {
            rgb = {0.0, 0.0, 255.0};
        }

        m[i] = rgb;
    }

    save_Mesh_Color(mesh, m, file_name);
    return normal_map;
}


//==============================================================================
// SECTION 4: TOPOLOGICAL PROPERTIES
//==============================================================================

/**
 * @brief Compute genus of a closed mesh using Euler characteristic.
 * @param mesh Input mesh (must be closed, single-component)
 * @return Genus g, where χ = V - E + F = 2 - 2g
 * @note Prints vertex/edge/face counts to stdout for verification
 */
int computeGenus(const Mesh &mesh)
{
    unsigned int nbVerts = 0;
    for (VertexCstIt i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
    {
        ++nbVerts;
    }
    std::cout << "# Vertices : " << nbVerts << std::endl;

    unsigned int nbEdges = 0;
    for (HalfedgeCstIt i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
    {
        ++nbEdges;
    }
    nbEdges /= 2;
    std::cout << "# Edges: " << nbEdges << std::endl;

    unsigned int nbFaces = 0;
    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        ++nbFaces;
    }
    std::cout << "# Faces: " << nbFaces << std::endl;

    unsigned int euler = nbVerts - nbEdges + nbFaces;
    unsigned int genus = (2 - euler) / 2;
    return genus;
}

//==============================================================================
// SECTION 5: DISTANCE COMPUTATIONS
//==============================================================================

/**
 * @brief Compute Euclidean distance from face centroid to reference point.
 * @param mesh Parent mesh
 * @param i Face iterator
 * @param p Reference point (e.g., robot position)
 * @return Array {dx, dy, dz, d_total} where d_total = sqrt(dx²+dy²+dz²)
 * @note Returns signed component distances for axis-specific analysis
 */
std::array<double, 4> compute_distance_face_point(const Mesh &mesh, const FacetCstIt &i, const geomAlgoLib::Point_3 &p)
{
    double sx = 0.0, sy = 0.0, sz = 0.0;
    unsigned int nb_sommets = 0;

    HalfedgeFacetCstCirc circ = i->facet_begin();
    HalfedgeFacetCstCirc start = circ;
    do {
        const auto& pt = circ->vertex()->point();
        sx += CGAL::to_double(pt.x());
        sy += CGAL::to_double(pt.y());
        sz += CGAL::to_double(pt.z());
        ++nb_sommets;
        ++circ;
    } while (circ != start);

    if (nb_sommets == 0) return {0.0, 0.0, 0.0, 0.0};

    geomAlgoLib::Point_3 center(sx / nb_sommets, sy / nb_sommets, sz / nb_sommets);

    double dx = p.x() - center.x();
    double dy = p.y() - center.y();
    double dz = p.z() - center.z();
    double d_total = std::sqrt(dx*dx + dy*dy + dz*dz);

    return {dx, dy, dz, d_total};
}

//==============================================================================
// SECTION 6: SEGMENTATION - Tag Assignment
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
 *  - ZoneAPortee: plane surface within horizontal reach and height tolerance of reference point
 * 
 * @param areas Map of face -> area
 * @param normal_angles Map of face -> normal vector
 * @param mesh Parent mesh (for distance computations)
 * @param reference_point Robot/avatar position for reachability analysis
 * @return Map: face_handle -> TagMaps (tags + associated colors)
 * 
 * @note Uses immediately-invoked lambda expressions for clear encapsulation of each tag logic
 */
std::map<Mesh::Facet_const_handle, TagMaps> segmentationMap(
    std::map<Mesh::Facet_const_handle, double> areas,
    std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel>> normal_angles,
    const Mesh &mesh,
    const geomAlgoLib::Point_3 &reference_point = geomAlgoLib::Point_3(0,0,-1000)
)
{
    std::map<Mesh::Facet_const_handle, TagMaps> M;
    double reach_radius = 2.0;


    for (auto const &pair : areas)
    {
        Mesh::Facet_const_handle key = pair.first;
        TagMaps T;


        auto distances = compute_distance_face_point(mesh, key, reference_point);
        double dx = distances[0];
        double dy = distances[1];
        double dz = distances[2];
        double d_total = distances[3];
        double distance_xy = std::sqrt(dx*dx + dy*dy);

        // Tag: Grande (large face by area threshold)
        auto grande = [&] () -> int { if (areas[key] > 5.)
        {
            T.tags["Grande"] = true;
            T.colors["Grande"] = colorDictionary["RED"];
        }
        else
        {
            T.tags["Grande"] = false;
            T.colors["Grande"] = colorDictionary["BLACK"];
        }return 0;
        }();

        // Tag: OrienteHaut (upward-facing: Z component dominant)
        auto orienteHaut = [&] () -> int {
            auto face_normal_vector = normal_angles[key];
            double max_axis = std::max(fabs(face_normal_vector.x()), std::max(fabs(face_normal_vector.z()), fabs(face_normal_vector.y())));
            if (fabs(normal_angles[key].z()) == max_axis)
            {
                T.tags["OrienteHaut"] = true;
                T.colors["OrienteHaut"] = colorDictionary["BLUE"];
            }
            else
            {
                T.tags["OrienteHaut"] = false;
                T.colors["OrienteHaut"] = colorDictionary["WHITE"];
            }return 0;
        }();

        // Tag: Plane (horizontal surface, angle < 15° with Z axis)
        auto plane = [&]() -> int {
            auto face_normal_vector = normal_angles[key];
            double nz = fabs(face_normal_vector.z());
            double length = std::sqrt(CGAL::squared_length(face_normal_vector));
            double cos_angle = nz / length;

            if (cos_angle > 0.966) {
                T.tags["Plane"] = true;
                T.colors["Plane"] = colorDictionary["GREEN"];
            } else {
                T.tags["Plane"] = false;
                T.colors["Plane"] = colorDictionary["BLACK"];
            }
            return 0;
        }();

        // Tag: VerySteep (inclined surface, angle > 30° from horizontal)
        auto verySteep = [&]() -> int {
            auto face_normal_vector = normal_angles[key];
            double nz = fabs(face_normal_vector.z());
            double length = std::sqrt(CGAL::squared_length(face_normal_vector));
            double cos_angle = nz / length;

            if (cos_angle < 0.866) {
                T.tags["VerySteep"] = true;
                T.colors["VerySteep"] = colorDictionary["YELLOW"];
            } else {
                T.tags["VerySteep"] = false;
                T.colors["VerySteep"] = colorDictionary["BLACK"];
            }
            return 0;
        }();

        // Tag: ZoneMarche (walkable: horizontal surface)
        auto zoneMarche = [&]() -> int {
            bool is_plane = T.tags["Plane"];

            if (is_plane) {
                T.tags["ZoneMarche"] = true;
                T.colors["ZoneMarche"] = colorDictionary["GREEN"];
            } else {
                T.tags["ZoneMarche"] = false;
                T.colors["ZoneMarche"] = colorDictionary["BLACK"];
            }
            return 0;
        }();

        // Tag: ZoneAPortee (reachable zone: plane surface near robot in XY, similar height)
        auto zoneAPortee = [&]() -> int {
            double max_height_diff = 0.3;

            if (T.tags["Plane"] && distance_xy < reach_radius && fabs(dz) < max_height_diff)
            {
                T.tags["ZoneAPortee"] = true;
                T.colors["ZoneAPortee"] = colorDictionary["MAGENTA"];
            }
            else
            {
                T.tags["ZoneAPortee"] = false;
                T.colors["ZoneAPortee"] = colorDictionary["WHITE"];
            }
            return 0;
        }();


        // Tag: Obstacle (vertical/steep surface blocking navigation)
        auto obstacle = [&]() -> int {
            bool is_steep = T.tags["VerySteep"];
            bool is_significant = areas[key] > 0.0;

            if (is_steep && is_significant && distance_xy < reach_radius) {
                T.tags["Obstacle"] = true;
                T.colors["Obstacle"] = colorDictionary["RED"];
            } else {
                T.tags["Obstacle"] = false;
                T.colors["Obstacle"] = colorDictionary["BLACK"];
            }
            return 0;
        }();

        // Suppress unused variable warnings
        (void)grande;
        (void)orienteHaut;
        (void)plane;
        (void)verySteep;
        (void)obstacle;
        (void)zoneMarche;
        (void)zoneAPortee;

        M[key] = T;
    }

    return M;
}

//==============================================================================
// SECTION 7: I/O - Export Functions
//==============================================================================

/**
 * @brief Export single face with RGB color (utility function, rarely used alone).
 * @param mesh Parent mesh
 * @param F Face handle
 * @param color RGB array [0,1] range
 * @param file_name Output COFF file path
 */
void save_Face_Color(const Mesh &mesh, Mesh::Facet_const_handle F, std::array<double, 3> color, std::string file_name)
{
    std::ofstream in_myfile;
    in_myfile.open(file_name);
    CGAL::set_ascii_mode(in_myfile);
    in_myfile << "COFF" << std::endl;

    HalfedgeFacetCstCirc j = F->facet_begin();
    CGAL_assertion(CGAL::circulator_size(j) >= 3);
    in_myfile << CGAL::circulator_size(j) << ' ';
    do
    {
        in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
    } while (++j != F->facet_begin());

    for (int n = 0; n < 3; n++)
    {
        in_myfile << ' ' << (int)(color[n])*255;
    }
    in_myfile << '\n';

    in_myfile.close();
}

/**
 * @brief Export mesh with per-face RGB colors (array format).
 * @param mesh Input mesh
 * @param m Map: face_handle -> RGB array
 * @param file_name Output COFF file path
 */
void save_Mesh_Color(const Mesh &mesh,
                     std::map<Mesh::Facet_const_handle, std::array<double,3>> m,
                     std::string file_name)
{
    std::ofstream in_myfile;
    in_myfile.open(file_name);
    CGAL::set_ascii_mode(in_myfile);

    in_myfile << "COFF" << std::endl
              << mesh.size_of_vertices() << ' '
              << mesh.size_of_facets() << " 0" << std::endl;

    std::copy(mesh.points_begin(), mesh.points_end(),
              std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        HalfedgeFacetCstCirc j = i->facet_begin();
        CGAL_assertion(CGAL::circulator_size(j) >= 3);
        
        in_myfile << CGAL::circulator_size(j) << ' ';
        do
        {
            in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
        } while (++j != i->facet_begin());

        for (int n = 0; n < 3; n++)  
        {
            int color_val = static_cast<int>(m.at(i)[n] * 255.0);
            color_val = std::clamp(color_val, 0, 255);
            in_myfile << ' ' << color_val;
        }
        in_myfile << '\n';
    }

    in_myfile.close();
    std::cout << "Successfully exported at path: " << file_name << " !" << std::endl;
}


/**
 * @brief Export mesh with normalized scalar-to-grayscale mapping.
 * @param mesh Input mesh
 * @param m Map: face_handle -> scalar value
 * @param file_name Output COFF file path
 * @note Normalizes values to [0,1] then maps to grayscale: white=min, black=max
 */
void save_Mesh_Color(const Mesh &mesh,
                     std::map<Mesh::Facet_const_handle, double> m,
                     std::string file_name)
{
    double min_val = std::numeric_limits<double>::max();
    double max_val = std::numeric_limits<double>::lowest();

    for (const auto& [facet, value] : m) {
        min_val = std::min(min_val, value);
        max_val = std::max(max_val, value);
    }

    std::cout << "Normalizing colors: min=" << min_val
              << ", max=" << max_val << std::endl;

    double range = max_val - min_val;
    if (range < 1e-10) range = 1.0;

    std::ofstream in_myfile;
    in_myfile.open(file_name);
    CGAL::set_ascii_mode(in_myfile);

    in_myfile << "COFF" << std::endl
              << mesh.size_of_vertices() << ' '
              << mesh.size_of_facets() << " 0" << std::endl;

    std::copy(mesh.points_begin(), mesh.points_end(),
              std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        HalfedgeFacetCstCirc j = i->facet_begin();
        CGAL_assertion(CGAL::circulator_size(j) >= 3);

        in_myfile << CGAL::circulator_size(j);
        do {
            in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
        } while (++j != i->facet_begin());

        double t = (m.at(i) - min_val) / range;
        int gray = static_cast<int>((1.0 - t) * 255.0);
        gray = std::clamp(gray, 0, 255);

        in_myfile << ' ' << gray << ' ' << gray << ' ' << gray << " 255";
        in_myfile << '\n';
    }

    in_myfile.close();
    std::cout << "Successfully exported at path: " << file_name << " !" << std::endl;
}

//==============================================================================
// SECTION 8: HIGH-LEVEL PIPELINE - Segmentation + Export
//==============================================================================

/**
 * @brief Complete segmentation pipeline for a single tag.
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
                  const geomAlgoLib::Point_3 &reference_point = geomAlgoLib::Point_3(0,0,0))
{
    auto area_map = compute_Area_Map(mesh);
    auto normal_angle_map = colorize_by_greatest_normal(mesh, "temp_normals.off");
    auto M = segmentationMap(area_map, normal_angle_map, mesh, reference_point);

    std::map<Mesh::Facet_const_handle, std::array<double,3>> color_map;
    for (auto const &pair : M) {
        color_map[pair.first] = pair.second.colors.at(Tag).to_array();
    }

    save_Mesh_Color(mesh, color_map, file_name);
}

//==============================================================================
// SECTION 9: VISUALIZATION UTILITIES - Robot Marker
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
                       double size)
{
    std::ofstream out(file_name);
    CGAL::set_ascii_mode(out);

    double x = robot_pos.x();
    double y = robot_pos.y();
    double z = robot_pos.z();
    double s = size / 2.0;

    out << "COFF" << std::endl;
    out << "8 6 0" << std::endl;

    out << (x - s) << " " << (y - s) << " " << (z - s) << std::endl;
    out << (x + s) << " " << (y - s) << " " << (z - s) << std::endl;
    out << (x + s) << " " << (y + s) << " " << (z - s) << std::endl;
    out << (x - s) << " " << (y + s) << " " << (z - s) << std::endl;
    out << (x - s) << " " << (y - s) << " " << (z + s) << std::endl;
    out << (x + s) << " " << (y - s) << " " << (z + s) << std::endl;
    out << (x + s) << " " << (y + s) << " " << (z + s) << std::endl;
    out << (x - s) << " " << (y + s) << " " << (z + s) << std::endl;

    out << "4 0 1 2 3 255 0 0 255" << std::endl;
    out << "4 4 5 6 7 255 0 0 255" << std::endl;
    out << "4 0 1 5 4 255 0 0 255" << std::endl;
    out << "4 2 3 7 6 255 0 0 255" << std::endl;
    out << "4 0 3 7 4 255 0 0 255" << std::endl;
    out << "4 1 2 6 5 255 0 0 255" << std::endl;

    out.close();
    std::cout << "Robot marker (cube) saved at: " << file_name << std::endl;
}

//==============================================================================
// SECTION 10: DEBUG UTILITIES
//==============================================================================

/**
 * @brief Print map contents to stdout.
 * @param comment Prefix label for output
 * @param m Map to display
 */
void print_map(std::string_view comment, std::map<geomAlgoLib::Mesh::Facet_const_handle, double>& m)
{
    std::cout << comment;
    for (const auto& [key, value] : m)
        std::cout << '[' << key << "] = " << value << "; ";
    std::cout << '\n';
}

} // namespace geomAlgoLib
