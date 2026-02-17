#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>


#include <map>

namespace geomAlgoLib{





using K = CGAL::Simple_cartesian<double>;
using Point_3 = K::Point_3;
using Traits = CGAL::Search_traits_3<K>;
using K_neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
using Tree = K_neighbor_search::Tree;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Polyhedron_3<Kernel>;

using FacetCstIt = Mesh::Facet_const_iterator;
using VertexCstIt = Mesh::Vertex_const_iterator;
using HalfedgeCstIt = Mesh::Halfedge_const_iterator;
using HalfedgeFacetCstCirc = Mesh::Halfedge_around_facet_const_circulator;


using FacetDoubleMap = std::map<Mesh::Facet_const_handle, double>;
using FacetIntMap = std::map<Mesh::Facet_const_handle, int>;

}