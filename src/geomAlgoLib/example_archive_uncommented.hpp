#pragma once

#include "types.hpp"
#include <map>
#include <math.h>
#include <stdlib.h>

namespace geomAlgoLib
{


struct Color {
	double r; 
	double g; 
	double b; 

	std::array<double, 3> to_array() const { return {r, g, b}; }
	std::array<double, 3> operator*(const double scalar) const { return {r * scalar, g * scalar, b * scalar}; }
	const Color operator+ (std::array<double, 3> other) const { return {r + other[0], g + other[1], b + other[2]}; }

}typedef Color;

struct TagMaps
{
	std::map<std::string, bool> tags;
	std::map<std::string,Color> colors; 
}typedef TagMaps;


inline Color RED      = {1.,0.,0.}; 
inline Color GREEN    = {0.,1.,0.}; 
inline Color BLUE     = {0.,0.,1.}; 
inline Color YELLOW   = {1.,1.,0.}; 
inline Color CYAN     = {0.,1.,1.}; 
inline Color MAGENTA  = {1.,0.,1.};
inline Color WHITE    = {1.,1.,1.}; 
inline Color BLACK    = {0.,0.,0.};

inline std::map<std::string, Color> colorDictionary = {
    {"RED", RED},
    {"GREEN", GREEN},
    {"BLUE", BLUE},
    {"YELLOW", YELLOW},
    {"CYAN", CYAN},
    {"MAGENTA", MAGENTA},
    {"WHITE", WHITE},
    {"BLACK", BLACK}
};








int computeGenus(const Mesh &mesh);
double compute_Area(const Mesh &mesh);
double compute_area_triangle(FacetCstIt &i);
double compute_area_quad(FacetCstIt &i);
double compute_Area_Quad_Mesh(const Mesh &mesh);
std::map<Mesh::Facet_const_handle, double> compute_Area_Map(const Mesh &mesh);
void save_Mesh_Color(const Mesh &mesh, 
					std::map<Mesh::Facet_const_handle, std::array<double, 3>> m, 
					std::string file_name);
void save_Mesh_Color(const Mesh &mesh, 
					std::map<Mesh::Facet_const_handle, double> m, 
					std::string file_name);
double get_normal_angle(const CGAL::Point_3< Kernel > & 	p,
						const CGAL::Point_3< Kernel > & 	q,
						const CGAL::Point_3< Kernel > & 	r,
						const CGAL::Vector_3<Kernel> & reference_vector
					);
double angle_degrees(const CGAL::Vector_3<Kernel>& vec, const CGAL::Vector_3<Kernel>& axis);
std::array<double,3 >  get_facet_normals_angle(const FacetCstIt &i);
const CGAL::Vector_3<Kernel> get_normal_vector(const CGAL::Point_3< Kernel > & 	p,
										  const CGAL::Point_3< Kernel > & 	q,
										  const CGAL::Point_3< Kernel > & 	r
										);
CGAL::Vector_3<Kernel>  get_facet_normals_vector(const FacetCstIt &i);
std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel> > colorize_by_greatest_normal(const Mesh &mesh, std::string file_name);
std::array<double, 4> compute_distance_face_point(const Mesh &mesh, const FacetCstIt &i, const geomAlgoLib::Point_3 &p);
std::map<Mesh::Facet_const_handle, TagMaps> segmentationMap(
    std::map<Mesh::Facet_const_handle, double> areas, 
    std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel>> normal_angles,
    const Mesh &mesh, 
    const geomAlgoLib::Point_3 &reference_point
);
void save_robot_marker(const geomAlgoLib::Point_3 &robot_pos, 
                       std::string file_name, 
                       double size);

void ColorizedTag(std::string Tag, 
                  const Mesh &mesh, 
                  std::string file_name,
                  const geomAlgoLib::Point_3 &reference_point);
void save_Face_Color(const Mesh &mesh, Mesh::Facet_const_handle F, std::array<double, 3> color, std::string file_name);



}
void print_map(std::string_view comment, std::map<geomAlgoLib::Mesh::Facet_const_handle, double>& m);
