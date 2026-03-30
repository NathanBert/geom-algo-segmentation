#include "myAlgo.hpp"
#include <iostream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <CGAL/boost/graph/iterator.h>  

namespace geomAlgoLib
{


int MeshUtils::getNumberOfNeighbours(VertexCstIt vertex)
{
    int count = 0;
    auto circulator = vertex->vertex_begin();
    auto end        = circulator;
    if (circulator != nullptr) {
        do { ++count; ++circulator; } while (circulator != end);
    }
    return count;
}

Vector3 MeshUtils::getCentroid(VertexCstIt vertex, const Mesh& mesh)
{
    int nb = getNumberOfNeighbours(vertex);
    if (nb == 0) return Vector3(0, 0, 0);

    double x = 0, y = 0, z = 0;

    Mesh::Halfedge_around_vertex_const_circulator h   = vertex->vertex_begin();
    Mesh::Halfedge_around_vertex_const_circulator end = h;
    do {
        auto p = h->opposite()->vertex()->point();
        x += CGAL::to_double(p.x());
        y += CGAL::to_double(p.y());
        z += CGAL::to_double(p.z());
        ++h;
    } while (h != end);

    auto vp = vertex->point();
    return Vector3(
        x / nb - CGAL::to_double(vp.x()),
        y / nb - CGAL::to_double(vp.y()),
        z / nb - CGAL::to_double(vp.z())
    );
}

Point3 MeshUtils::applyTransformOn(VertexCstIt vertex, const std::map<VertexCstIt,Vector3>& trans)
{
    
    auto it = trans.find(vertex);
    auto p = vertex->point();

    auto new_point = Point3(
        CGAL::to_double(p.x()) + it->second.x(),
        CGAL::to_double(p.y()) + it->second.y(),
        CGAL::to_double(p.z()) + it->second.z()
    );

    return new_point;
}

Lissage::Lissage(std::string meshPath)
{
    geomAlgoLib::readOFF(meshPath, mesh);  
    mesh_map["original"] = mesh;
}

void Lissage::LissageLaplacien(int iterations, int depth, std::string name)
{
    Mesh working_mesh = mesh;  

    for (int iter = 0; iter < iterations; ++iter)
    {
        auto centroids = MeshUtils::applyOnAllPoints(working_mesh,
            [&](VertexCstIt v) { return MeshUtils::getCentroid(v, working_mesh); });

        for (auto vit = working_mesh.vertices_begin(); vit != working_mesh.vertices_end(); ++vit) {
            auto it = centroids.find(vit);
            if (it == centroids.end()) continue;
            auto p = vit->point();
            vit->point() = Point3(
                CGAL::to_double(p.x()) + CGAL::to_double(it->second.x()),
                CGAL::to_double(p.y()) + CGAL::to_double(it->second.y()),
                CGAL::to_double(p.z()) + CGAL::to_double(it->second.z())
            );
        }
    }

    mesh_map[name] = working_mesh;
}

void Lissage::LissageFacteurDiffusion(int iterations, int depth, std::string name, double diffusion_factor)
{   
    Mesh working_mesh = mesh;  

    for (int iter = 0; iter < iterations; ++iter)
    {
        auto centroids = MeshUtils::applyOnAllPoints(working_mesh,
            [&](VertexCstIt v) { return MeshUtils::getCentroid(v, working_mesh); });
        
        for (auto vit = working_mesh.vertices_begin(); vit != working_mesh.vertices_end(); ++vit) {
            auto it = centroids.find(vit);
            if (it == centroids.end()) continue;
            auto p = vit->point();
            vit->point() = Point3(
                CGAL::to_double(p.x()) + CGAL::to_double(it->second.x()) * diffusion_factor,
                CGAL::to_double(p.y()) + CGAL::to_double(it->second.y()) * diffusion_factor,
                CGAL::to_double(p.z()) + CGAL::to_double(it->second.z()) * diffusion_factor
            );
        }

    }

    mesh_map[name] = working_mesh;
}

void Lissage::LissageTaubin(int iterations, int depth, std::string name, double labda, double mu)
{   
    assert(mu < 0);
    assert(std::abs(mu) > labda);

    Mesh working_mesh = mesh;  

    for (int iter = 0; iter < iterations; ++iter)
    {
        auto centroids_lambda = MeshUtils::applyOnAllPoints(working_mesh,
            [&](VertexCstIt v) { return MeshUtils::getCentroid(v, working_mesh); });

        for (auto vit = working_mesh.vertices_begin(); vit != working_mesh.vertices_end(); ++vit) {
            auto it = centroids_lambda.find(vit);
            if (it == centroids_lambda.end()) continue;
            auto p = vit->point();
            vit->point() = Point3(
                CGAL::to_double(p.x()) + labda * CGAL::to_double(it->second.x()),
                CGAL::to_double(p.y()) + labda * CGAL::to_double(it->second.y()),
                CGAL::to_double(p.z()) + labda * CGAL::to_double(it->second.z())
            );
        }


        auto centroids_mu = MeshUtils::applyOnAllPoints(working_mesh,
            [&](VertexCstIt v) { return MeshUtils::getCentroid(v, working_mesh); });

        for (auto vit = working_mesh.vertices_begin(); vit != working_mesh.vertices_end(); ++vit) {
            auto it = centroids_mu.find(vit);
            if (it == centroids_mu.end()) continue;
            auto p = vit->point();
            vit->point() = Point3(
                CGAL::to_double(p.x()) + mu * CGAL::to_double(it->second.x()),
                CGAL::to_double(p.y()) + mu * CGAL::to_double(it->second.y()),
                CGAL::to_double(p.z()) + mu * CGAL::to_double(it->second.z())
            );
        }
    }

    mesh_map[name] = working_mesh;
}

const Mesh& Lissage::getMesh() const
{
    return mesh;
}

Mesh& Lissage::getMeshMap(std::string key)
{
    return mesh_map[key];
}

void Lissage::setMeshMap(const std::map<std::string, Mesh>& map)
{
    mesh_map = map;
}

}



// Source - https://stackoverflow.com/a/46931770
// Posted by Arafat Hasan, modified by community. See post 'Timeline' for change history
// Retrieved 2026-02-26, License - CC BY-SA 4.0
// pour séparer les strings
std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}
