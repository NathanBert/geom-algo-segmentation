/*#include "example.hpp"

#include <iostream>
#include <math.h>

namespace geomAlgoLib
{
    
double compute_Area(const Mesh &mesh)
{
    double aera_face = 0;



	for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
		aera_face += compute_area_triangle(i);
    }   

	return aera_face;
}

std::map<Mesh::Facet_const_handle, TagMaps> segmentationMap(
    std::map<Mesh::Facet_const_handle, double> areas, 
    std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel>> normal_angles,
    const Mesh &mesh, 
    const geomAlgoLib::Point_3 &reference_point = geomAlgoLib::Point_3(0,0,-1000)
)
{

	std::map<Mesh::Facet_const_handle, TagMaps> M;

	for (auto const &pair : areas)
    {
		Mesh::Facet_const_handle key = pair.first;
		TagMaps T;

		// Traitement Tag Grande
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

		// Traitement Tag OrienteHaut 
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

		// Traitment Tag Plane 
		auto plane = [&]() -> int {
            auto face_normal_vector = normal_angles[key];
            double nz = fabs(face_normal_vector.z());
            double length = std::sqrt(CGAL::squared_length(face_normal_vector));
            double cos_angle = nz / length;  
            
            // Plane si angle < 15 degrés avec la verticale
            if (cos_angle > 0.966) {
                T.tags["Plane"] = true;
                T.colors["Plane"] = colorDictionary["GREEN"];
            } else {
                T.tags["Plane"] = false;
                T.colors["Plane"] = colorDictionary["BLACK"];
            }
            return 0;
        }();

		// Traitement Tag VerySteep 
		auto verySteep = [&]() -> int {
			auto face_normal_vector = normal_angles[key];
			double nz = fabs(face_normal_vector.z()); 
			double length = std::sqrt(CGAL::squared_length(face_normal_vector));
			double cos_angle = nz / length;
			
			if (cos_angle < 0.866) {  // cos(30°) = 0.866
				T.tags["VerySteep"] = true;
				T.colors["VerySteep"] = colorDictionary["YELLOW"];
			} else {
				T.tags["VerySteep"] = false;
				T.colors["VerySteep"] = colorDictionary["BLACK"];
			}
			return 0;
		}();


		// Traitement Tag ZoneMarhe (plane + grande) 
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

		// Traitement Tag Obstacle (inclinee + grande)
        auto obstacle = [&]() -> int {
            bool is_steep = T.tags["VerySteep"];
            bool is_significant = areas[key] > 0.0;  //ajustable en fonction du robot
            
            if (is_steep && is_significant) {
                T.tags["Obstacle"] = true;
                T.colors["Obstacle"] = colorDictionary["RED"];
            } else {
                T.tags["Obstacle"] = false;
                T.colors["Obstacle"] = colorDictionary["BLACK"];
            }
            return 0;
        }();

		auto zoneAPortee = [&]() -> int {
			auto distances = compute_distance_face_point(mesh, key, reference_point);
			double dx = distances[0];
			double dy = distances[1];
			double dz = distances[2];
			double d_total = distances[3];
			
			double reach_radius = 2.0;
			double max_height_diff = 0.2; 
			
			double distance_xy = std::sqrt(dx*dx + dy*dy);
			

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

		// AUTRES TRAITEMENT A AJOUTER ICI


		// Regroupement de tous les calculs de tags dans des lambdas pour une meilleur encapsulation et lisibilité du code. 
		// Chaque lambda traite un tag spécifique et met à jour les tags et couleurs correspondants dans la structure TagMaps T.

		// Les lambdas sont appelées immédiatement pour effectuer les calculs et mises à jour des tags pour chaque face du mesh.
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

void ColorizedTag(std::string Tag, 
                  const Mesh &mesh, 
                  std::string file_name,
                  const geomAlgoLib::Point_3 &reference_point = geomAlgoLib::Point_3(0,0,0))
{
    auto area_map = compute_Area_Map(mesh);
    auto normal_angle_map = colorize_by_greatest_normal(mesh, "temp_normals.off");
    auto M = segmentationMap(area_map, normal_angle_map, mesh, reference_point);
    
    // Construire map face -> couleur RGB pour le tag demandé
    std::map<Mesh::Facet_const_handle, std::array<double,3>> color_map;
    for (auto const &pair : M) {
        color_map[pair.first] = pair.second.colors.at(Tag).to_array();
    }
    
    save_Mesh_Color(mesh, color_map, file_name);
}


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
		//i->halfedge()->vertex_begin()->vertex()->point()
		in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

	} while (++j != F->facet_begin());

	for (int n = 0; n < 3; n++)
	{
		in_myfile << ' ' << (int)(color[n])*255;
		
	} 
	in_myfile << '\n';
	
	in_myfile.close();



}

void save_Mesh_Color(const Mesh &mesh, 
					std::map<Mesh::Facet_const_handle, std::array<double, 3>> m, 
					std::string file_name)
{
	std::ofstream in_myfile;
	in_myfile.open(file_name);
	int facet_size;

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
		facet_size = 0;
		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
			facet_size++;

		} while (++j != i->facet_begin());

		for (int n = 0; n < facet_size; n++)
		{
			in_myfile << ' ' << (int)(m[i][n]*1000)%255;
			
		} 
		in_myfile << '\n';


	}

	in_myfile.close();

	std::cout << "Successfully exported at path: " << file_name << " !" << std::endl;
}



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


double angle_degrees(const CGAL::Vector_3<Kernel>& vec, const CGAL::Vector_3<Kernel>& axis) {
    double cos_theta = CGAL::scalar_product(vec, axis) /
                       std::sqrt(CGAL::squared_length(vec) * CGAL::squared_length(axis));
    return std::acos(std::clamp(cos_theta, -1.0, 1.0)) * 180.0 / M_PI;
}

double get_normal_angle(const CGAL::Point_3< Kernel > & 	p,
										  const CGAL::Point_3< Kernel > & 	q,
										  const CGAL::Point_3< Kernel > & 	r,
										  const CGAL::Vector_3<Kernel> & reference_vector
										)
{
	const CGAL::Vector_3<Kernel> normal = CGAL::normal(p,q,r);

	double angle = angle_degrees(normal, reference_vector);

	return angle;
}

const CGAL::Vector_3<Kernel> get_normal_vector(const CGAL::Point_3< Kernel > & 	p,
										  const CGAL::Point_3< Kernel > & 	q,
										  const CGAL::Point_3< Kernel > & 	r
										)
{
	const CGAL::Vector_3<Kernel> normal = CGAL::normal(p,q,r);

	return normal;
}

std::array<double,3 >  get_facet_normals_angle(const FacetCstIt &i)
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

std::map<Mesh::Facet_const_handle, CGAL::Vector_3<Kernel> > colorize_by_greatest_normal(const Mesh &mesh, std::string file_name)
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
		y = ((int)(x/max_axis)) * 255;
		z = ((int)(x/max_axis)) * 255;


		
		std::array<double,3> rgb = {x,y,z};
		m[i] = rgb;
	}

	save_Mesh_Color(mesh,m,file_name);

	return(normal_map);

}

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
    
    // Distances selon chaque axe
    double dx = p.x() - center.x();
    double dy = p.y() - center.y();
    double dz = p.z() - center.z();
    double d_total = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    return {dx, dy, dz, d_total};  // [0]=dx, [1]=dy, [2]=dz, [3]=distance euclidienne totale
}


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
    
    out << (x - s) << " " << (y - s) << " " << (z - s) << std::endl;  // 0
    out << (x + s) << " " << (y - s) << " " << (z - s) << std::endl;  // 1
    out << (x + s) << " " << (y + s) << " " << (z - s) << std::endl;  // 2
    out << (x - s) << " " << (y + s) << " " << (z - s) << std::endl;  // 3
    out << (x - s) << " " << (y - s) << " " << (z + s) << std::endl;  // 4
    out << (x + s) << " " << (y - s) << " " << (z + s) << std::endl;  // 5
    out << (x + s) << " " << (y + s) << " " << (z + s) << std::endl;  // 6
    out << (x - s) << " " << (y + s) << " " << (z + s) << std::endl;  // 7
    
    out << "4 0 1 2 3 255 0 0 255" << std::endl;  // Face bas
    out << "4 4 5 6 7 255 0 0 255" << std::endl;  // Face haut
    out << "4 0 1 5 4 255 0 0 255" << std::endl;  // Face avant
    out << "4 2 3 7 6 255 0 0 255" << std::endl;  // Face arrière
    out << "4 0 3 7 4 255 0 0 255" << std::endl;  // Face gauche
    out << "4 1 2 6 5 255 0 0 255" << std::endl;  // Face droite
    
    out.close();
    std::cout << "Robot marker (cube) saved at: " << file_name << std::endl;
}



}

void print_map(std::string_view comment, std::map<geomAlgoLib::Mesh::Facet_const_handle, double>& m)
{
    std::cout << comment;
    // Iterate using C++17 facilities
    for (const auto& [key, value] : m)
        std::cout << '[' << key << "] = " << value << "; ";
 

    std::cout << '\n';
}*/



