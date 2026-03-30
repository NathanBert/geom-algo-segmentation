#include "geomAlgoLib/io.hpp"
//#include <example.hpp>
#include <iostream>
#include "geomAlgoLib/myAlgo.hpp"
#include "geomAlgoLib/batch_renderer.hpp"

#include <iostream>
#include <string>
#include <string>
//#include <types.hpp>
#include <map>

struct Config {
    std::string meshPath;
    std::string outputDir;
    std::string base_path;
    
    Config(int argc, char** argv) {
        meshPath = "../data/cube-bruit-2.off";
        outputDir = "../output/";

        base_path = split(meshPath, "/").back();
        base_path = split(base_path, ".").front(); 

        if (argc >= 2)
        {
            meshPath = argv[1];
            base_path = split(meshPath, "/").back();
            base_path = split(base_path, ".").front(); 
        } 
        std::cout << "Base path: " << base_path << std::endl;

    }
    
};

int main(int argc, char** argv) {

    Config config(argc, argv);

    geomAlgoLib::Lissage L = geomAlgoLib::Lissage(config.meshPath);

    auto nb = geomAlgoLib::MeshUtils::applyOnAllPoints(L.getMesh(), [](geomAlgoLib::VertexCstIt& v){ return geomAlgoLib::MeshUtils::getNumberOfNeighbours(v); });

    std::vector<int> iterations_number = {1, 2, 3, 4, 5, 6, 7};//, 40, 100, 200};
    std::vector<double> lambda_values = {0.33};
    std::vector<double> mu_values = {-0.34};

    int iterations;
    double lambda, mu;

    for (int i = 0; i < iterations_number.size(); ++i) {

        iterations = iterations_number[i];
        lambda = lambda_values[0];
        mu = mu_values[0];

        std::string name_laplacien = config.base_path + "_LissageLaplacien_" + std::to_string(iterations);
        L.LissageLaplacien(iterations,1, name_laplacien);
        geomAlgoLib::writeOFF(L.getMeshMap(name_laplacien), config.outputDir + name_laplacien + ".off");

        std::string name_diffusion = config.base_path + "_LissageFacteurDiffusion_" + std::to_string(iterations) + "_" + std::to_string(lambda);
        L.LissageFacteurDiffusion(iterations,1, name_diffusion, lambda);
        geomAlgoLib::writeOFF(L.getMeshMap(name_diffusion), config.outputDir + name_diffusion + ".off");

        std::string name_taubin = config.base_path + "_LissageTaubin_" + std::to_string(iterations) + "_" + std::to_string(lambda) + "_" + std::to_string(mu);
        L.LissageTaubin(iterations,1, name_taubin, lambda, mu);
        geomAlgoLib::writeOFF(L.getMeshMap(name_taubin), config.outputDir + name_taubin + ".off");




        std::vector<std::pair<std::string, std::function<double(geomAlgoLib::VertexCstIt,geomAlgoLib::VertexCstIt)>>> dists = {
            {"uniforme", geomAlgoLib::DistFuncs::uniforme},
            {"inverseDist", geomAlgoLib::DistFuncs::inverseDist},
            {"cotangentes", geomAlgoLib::DistFuncs::cotangentes}
        };

        std::string name_laplacien_pondere = config.outputDir + config.base_path + "_LissageLaplacienPondere_" + std::to_string(iterations);

        for (auto& [func_name, func] : dists) { 
                auto fname = name_laplacien_pondere+ "_" + func_name;
                L.LissageLaplacienPondere(iterations, 0, fname, func);
                geomAlgoLib::writeOFF(L.getMeshMap(fname), config.outputDir + fname + ".off");
        }
        
        //L.LissageLaplacienPondere(iterations,1, name_laplacien_pondere);

    }
    std::cout << "Lissage vers: " << config.outputDir << std::endl;
    std::cout << "Render depuis: " << config.outputDir << std::endl;
    std::cout << "CWD: " << fs::current_path() << std::endl;
    batch_render_off(config.outputDir, config.outputDir + "renders/");
    batch_render_off("../data/", config.outputDir + "renders/");


    return 0;
}