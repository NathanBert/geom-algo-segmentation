#pragma once

#include "types.hpp"
#include "io.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>
#include <functional>
#include "myOldAlgo.hpp"

#define PI 3.14159265358979323846

namespace geomAlgoLib
{

using Epick   = CGAL::Epick;
using Vector3 = CGAL::Vector_3<Epick>;
using Point3  = CGAL::Point_3<Epick>;

class MeshUtils
{
    public:

        static int getNumberOfNeighbours(VertexCstIt vertex);
        static Vector3 getCentroid(VertexCstIt vertex, const Mesh& mesh);
        static Point3 applyTransformOn(VertexCstIt vertex, const std::map<VertexCstIt,Vector3>& trans);

        template<typename Func, typename... Args>
        static auto applyOnAllPoints(const Mesh& mesh, Func&& func, Args&&... args)
        -> std::map<VertexCstIt, decltype(func(std::declval<VertexCstIt&>(), std::declval<Args>()...))>
        {
            using RetT = decltype(func(std::declval<VertexCstIt&>(), std::declval<Args>()...));
            std::map<VertexCstIt, RetT> results;
            for (auto vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
                results[vit] = func(vit, std::forward<Args>(args)...);
            }
            return results;
        }

        template<typename DistFunc>
        static std::map<VertexCstIt, double> computeWeightCentroid(VertexCstIt vertex, DistFunc dist_func)
        {
            std::map<VertexCstIt, double> weights;

            Mesh::Halfedge_around_vertex_const_circulator h   = vertex->vertex_begin();
            Mesh::Halfedge_around_vertex_const_circulator end = h;
            if (h == nullptr) return weights;

            double sum_weights = 0.0;

            do {
                VertexCstIt neighbor = h->opposite()->vertex();
                double w = dist_func(vertex, neighbor);
                weights[neighbor] = w;
                sum_weights += w;
                ++h;
            } while (h != end);

            for (auto& [vit, w] : weights) {
                if (!std::isfinite(sum_weights) || sum_weights == 0.0) {
                    w = 1.0 / weights.size();
                } else {
                    w /= sum_weights;
                }
                w = std::clamp(w, 0.0, 1.0);  
            }


            return weights;
        }

};

class Lissage
{
    private:
        Mesh mesh;                             
        std::map<std::string, Mesh> mesh_map; 

    public:
        Lissage(std::string meshPath);

        void LissageLaplacien(int iterations, int depth, std::string name);
        void LissageFacteurDiffusion(int iterations, int depth, std::string name, double diffusion_factor);
        void LissageTaubin(int iterations, int depth, std::string name, double labda, double mu);

        const Mesh& getMesh() const;
        Mesh& getMeshMap(std::string key);
        const std::map<std::string, Mesh>& getMeshMap() const;

        void setMeshMap(const std::map<std::string, Mesh>& map);
        void setMesh(const Mesh& m);


        template<typename DistFunc>
        void LissageLaplacienPondere(int iterations, int depth, std::string name, DistFunc dist_func)
        {
            Mesh working_mesh = mesh;

            for (int iter = 0; iter < iterations; ++iter)
            {
                std::map<VertexCstIt, Vector3> weighted_centroids;

                for (auto vit = working_mesh.vertices_begin(); vit != working_mesh.vertices_end(); ++vit) {
                    auto weights = MeshUtils::computeWeightCentroid(vit, dist_func);

                    double x = 0, y = 0, z = 0;
                    for (auto& [neighbor, w] : weights) {
                        auto p = neighbor->point();
                        x += w * CGAL::to_double(p.x());
                        y += w * CGAL::to_double(p.y());
                        z += w * CGAL::to_double(p.z());
                    }

                    auto vp = vit->point();
                    weighted_centroids[vit] = Vector3(
                        x - CGAL::to_double(vp.x()),
                        y - CGAL::to_double(vp.y()),
                        z - CGAL::to_double(vp.z())
                    );
                }

                for (auto vit = working_mesh.vertices_begin(); vit != working_mesh.vertices_end(); ++vit) {
                    auto it = weighted_centroids.find(vit);
                    if (it == weighted_centroids.end()) continue;
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
        
};

class DistFuncs
{
    public:
        static double inverseDist(VertexCstIt v, VertexCstIt neighbor) {
            double d2 = CGAL::to_double(CGAL::squared_distance(v->point(), neighbor->point()));
            if (d2 < 1e-12) return 0.0;  
            double w = 1.0 / std::sqrt(d2);
            return std::min(w, 100.0); 
        }

        static double uniforme(VertexCstIt v, VertexCstIt neighbor)
        {
            return 1.0;
        }

        static double cotangentes(VertexCstIt v, VertexCstIt neighbor)
        {
            using SCD    = CGAL::Simple_cartesian<double>;
            using Pt3d   = SCD::Point_3;
            using Vec3d  = CGAL::Vector_3<SCD>;

            auto toPt = [](auto p) -> Pt3d {
                return Pt3d(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
            };

            auto circ     = v->vertex_begin();
            auto end_circ = circ;
            Mesh::Halfedge_const_handle h_edge;
            bool found = false;

            do {
                if (circ->opposite()->vertex() == neighbor) {
                    h_edge = &*circ;         
                    found = true;
                    break;
                }
                ++circ;
            } while (circ != end_circ);

            if (!found || h_edge->is_border() || h_edge->opposite()->is_border())
                return 0.0;

            Pt3d pv  = toPt(v->point());
            Pt3d pn  = toPt(neighbor->point());


            Pt3d palpha = toPt(h_edge->next()->vertex()->point());


            Pt3d pbeta  = toPt(h_edge->opposite()->next()->vertex()->point());

            auto cotan = [](Pt3d& apex, Pt3d& a, Pt3d& b) -> double {
                Vec3d u(apex, a), w(apex, b);
                double dot   = CGAL::scalar_product(u, w);
                double norme = std::sqrt(CGAL::cross_product(u, w).squared_length());
                return (norme < 1e-10) ? 0.0 : dot / norme;
            };

            return 0.5 * (cotan(palpha, pv, pn) + cotan(pbeta, pv, pn));
        }

};


} 


std::vector<std::string> split(std::string s, std::string delimiter);