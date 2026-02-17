#include <io.hpp>
//#include <example.hpp>
#include <iostream>
#include "myAlgo.hpp"
#include <iostream>
#include <string>
#include <string>
//#include <types.hpp>
#include <map>

/**

int main(int argc, char *argv[]){


    if(argc < 2){
        throw std::invalid_argument("This program expects at least 1 argument (path to a mesh).");
    }

    const std::string meshPath = std::string{argv[1]};
    
    geomAlgoLib::Mesh myMesh;
    
    geomAlgoLib::readOFF(meshPath, myMesh);
    

    auto genus = geomAlgoLib::computeGenus(myMesh);
    std::cout << "The Genus of [" << meshPath << "] is = " << std::to_string(genus) << std::endl;

    geomAlgoLib::writeOFF(myMesh,"../data/cube_save.off");

    std::cout << "The end..." << std::endl;

    auto area_face = geomAlgoLib::compute_Area(myMesh);
    std::cout << "Area Triangle Mesh = "<< area_face << std::endl;

    area_face = geomAlgoLib::compute_Area_Quad_Mesh(myMesh);
    std::cout << "Area Quad Mesh = "<< area_face << std::endl;



	std::map<geomAlgoLib::Mesh::Facet_const_handle, double> m = geomAlgoLib::compute_Area_Map(myMesh);

    print_map("map", m);


    //geomAlgoLib::save_Mesh_Color(myMesh, m, meshPath);


    //geomAlgoLib::colorize_by_greatest_normal(myMesh,meshPath);


    geomAlgoLib::FacetCstIt f = myMesh.facets_begin();
    geomAlgoLib::Point_3 test_p(0,0,0);
    double d = geomAlgoLib::compute_distance_face_point(myMesh, f, test_p);
    std::cout << "Dist origine->face0: " << d << std::endl;



    geomAlgoLib::ColorizedTag("OrienteHaut", myMesh, meshPath);
    std::cout << "ColorizedTag OrienteHaut done." << std::endl;


    return 0;
}*/


// Version pour calcul d'aires et visualisation de taille d'aire par face
/*
int main(int argc, char *argv[]){
    
    if(argc < 2){
        throw std::invalid_argument("Usage: ./mainApp <mesh.off>");
    }
    
    const std::string meshPath = std::string{argv[1]};
    geomAlgoLib::Mesh myMesh;
    geomAlgoLib::readOFF(meshPath, myMesh);
    
    auto area_map = geomAlgoLib::compute_Area_Map(myMesh);
    
    std::string output_path = "../data/output_aires.off";
    geomAlgoLib::save_Mesh_Color(myMesh, area_map, output_path);
    
    std::cout << "Aire visualisation saved to: " << output_path << std::endl;
    
    return 0;
}

*/

// Version pour segmentation et visualisation des zones de marche/portée/obstacles
/*
int main(int argc, char *argv[]){
    
    if(argc < 2){
        throw std::invalid_argument("Usage: ./mainApp <mesh.off> [x y z]");
    }
    
    const std::string meshPath = std::string{argv[1]};
    geomAlgoLib::Mesh myMesh;
    geomAlgoLib::readOFF(meshPath, myMesh);
    
    // Position du robot
    geomAlgoLib::Point_3 robot_position(2, 2, 0);
    if (argc >= 5) {
        robot_position = geomAlgoLib::Point_3(
            std::stod(argv[2]),
            std::stod(argv[3]),
            std::stod(argv[4])
        );
    }
    
    std::cout << "Robot position: (" << robot_position.x() << ", " 
              << robot_position.y() << ", " << robot_position.z() << ")" << std::endl;
    
    // Sauver le marqueur (cube rouge de 30cm de côté)
    geomAlgoLib::save_robot_marker(robot_position, "../data/robot_marker.off", 0.3);
    
    // Segmentations
    geomAlgoLib::ColorizedTag("ZoneAPortee", myMesh, "../data/zone_portee.off", robot_position);
    geomAlgoLib::ColorizedTag("ZoneMarche", myMesh, "../data/zone_marche.off", robot_position);
    geomAlgoLib::ColorizedTag("Obstacle", myMesh, "../data/obstacles.off", robot_position);
    
    std::cout << "Done!" << std::endl;
    std::cout << "Open in MeshLab: zone_portee.off + robot_marker.off" << std::endl;
    
    return 0;
}
*/



/**
 * @file main.cpp
 * @brief TP3 Segmentation - Programme principal organisé par questions
 * 
 * Structure:
 *  - Partie I  : Calcul de propriétés locales
 *  - Partie II : Segmentation par seuillage et application
 * 
 * @author ISIMA ZZ3 - TP3 Segmentation
 * @date 2026-02-17
 */



using namespace geomAlgoLib;

//==============================================================================
// CONFIGURATION GLOBALE
//==============================================================================

struct Config {
    std::string meshPath;
    std::string outputDir;
    Point_3 robotPosition;
    
    Config(int argc, char** argv) {
        meshPath = "../data/improved_home.off";
        outputDir = "../output/";
        robotPosition = Point_3(0.0, 0.0, 0.0);
        
        if (argc >= 2) meshPath = argv[1];
        if (argc >= 4) {
            double x = std::stod(argv[2]);
            double y = std::stod(argv[3]);
            double z = (argc >= 5) ? std::stod(argv[4]) : 0.0;
            robotPosition = Point_3(x, y, z);
        }
    }
    
    void print() const {
        std::cout << "Configuration du TP3\n";
        std::cout << "========================================\n";
        std::cout << "Mesh d'entrée : " << meshPath << "\n";
        std::cout << "Dossier sortie : " << outputDir << "\n";
        std::cout << "Position robot : (" << robotPosition.x() << ", " 
                  << robotPosition.y() << ", " << robotPosition.z() << ")\n";
        std::cout << "========================================\n\n";
    }
};

//==============================================================================
// PARTIE I.0 : PRISE EN MAIN - Genre du maillage
//==============================================================================

void partie_I0_Genre(const Mesh& mesh) {
    
    int genus = computeGenus(mesh);
    
    std::cout << "Genre du maillage : " << genus << "\n";
    std::cout << "Interprétation : ";
    if (genus == 0) {
        std::cout << "Topologie sphérique (surface fermée sans trou)\n";
    } else if (genus == 1) {
        std::cout << "Topologie torique (surface avec 1 trou)\n";
    } else {
        std::cout << "Surface avec " << genus << " trous\n";
    }
    std::cout << "\n";
}

//==============================================================================
// PARTIE I.1 : CALCUL D'AIRES
//==============================================================================

void partie_I1_Aires(const Mesh& mesh, const std::string& outputDir) {
    
    double totalAreaTri = compute_Area(mesh);
    double totalAreaMixed = compute_Area_Quad_Mesh(mesh);
    
    std::cout << "Aire totale (méthode triangles) : " << totalAreaTri << " unités²\n";
    std::cout << "Aire totale (méthode mixte)     : " << totalAreaMixed << " unités²\n";
    
    auto areaMap = compute_Area_Map(mesh);
    std::cout << "Nombre de faces traitées : " << areaMap.size() << "\n";
    
    double minArea = std::numeric_limits<double>::max();
    double maxArea = std::numeric_limits<double>::lowest();
    for (const auto& [facet, area] : areaMap) {
        minArea = std::min(minArea, area);
        maxArea = std::max(maxArea, area);
    }
    std::cout << "Aire min : " << minArea << " | Aire max : " << maxArea << "\n";
    std::cout << "\n";
}

//==============================================================================
// PARTIE I.2 : EXPORT MESH COLORÉ
//==============================================================================

void partie_I2_ExportColore(const Mesh& mesh, const std::string& outputDir, const std::string& savedMeshName = "mesh_colored_by_area.off") {
    
    // I.2.a : Sauvegarder mesh avec couleurs basées sur aires
    auto areaMap = compute_Area_Map(mesh);
    std::string outputFile = outputDir + savedMeshName;
    
    save_Mesh_Color(mesh, areaMap, outputFile);
    
    std::cout << "Maillage coloré (dégradé blanc-noir selon aire) exporté\n";
    std::cout << "Fichier : " << outputFile << "\n";
    std::cout << "Visualisation : ouvrir avec MeshLab (Render > Color > Per Face)\n";
    std::cout << "\n";
}

//==============================================================================
// PARTIE I.4 : AUTRES MESURES LOCALES
//==============================================================================

void partie_I4_AutresMesures(const Mesh& mesh, const std::string& outputDir, const std::string& savedMeshName = "mesh_colored_by_normal.off") {
    
    std::cout << "Calcul des orientations des faces (angle normal vs axes)...\n";
    
    auto normalMap = colorize_by_greatest_normal(mesh, 
                                                  outputDir + savedMeshName);
    
    std::cout << "Maillage coloré selon la composante dominante de la normale exporté\n";
    std::cout << "Rouge : faces orientées selon X\n";
    std::cout << "Vert  : faces orientées selon Y\n";
    std::cout << "Bleu  : faces orientées selon Z\n";
    std::cout << "\n";
}

//==============================================================================
// PARTIE II.1 : SEGMENTATION PAR SEUILLAGE
//==============================================================================

void partie_II1_Segmentation(const Mesh& mesh, 
                              const std::string& outputDir,
                              const Point_3& robotPos) {
    
    auto areaMap = compute_Area_Map(mesh);
    auto normalMap = colorize_by_greatest_normal(mesh, "temp_normals.off");
    auto segmentMap = segmentationMap(areaMap, normalMap, mesh, robotPos);
    
    std::cout << "Segmentation calculée avec les tags suivants :\n";
    std::cout << "  - Grande       : aire > 5.0\n";
    std::cout << "  - OrienteHaut  : composante Z dominante\n";
    std::cout << "  - Plane        : angle < 15° avec horizontal\n";
    std::cout << "  - VerySteep    : angle > 30° avec horizontal\n";
    std::cout << "  - ZoneMarche   : surface plane praticable\n";
    std::cout << "  - Obstacle     : surface inclinée bloquante\n";
    std::cout << "  - ZoneAPortee  : surface accessible depuis robot\n";
    std::cout << "\n";
    
    std::vector<std::string> tags = {
        "Grande", "OrienteHaut", "Plane", "VerySteep"
    };
    
    std::cout << "Export des maillages colorés par tag...\n";
    for (const auto& tag : tags) {
        std::string filename = outputDir + "mesh_tag_" + tag + ".off";
        ColorizedTag(tag, mesh, filename, robotPos);
        std::cout << "  ✓ " << tag << " -> " << filename << "\n";
    }
    std::cout << "\n";
}

//==============================================================================
// PARTIE II.2 : APPLICATION ROBOTIQUE
//==============================================================================

void partie_II2_ApplicationRobot(const Mesh& mesh, 
                                  const std::string& outputDir,
                                  const Point_3& robotPos) {
    
    std::cout << "Analyse de l'environnement intérieur pour navigation autonome\n";
    std::cout << "Position du robot : (" << robotPos.x() << ", " 
              << robotPos.y() << ", " << robotPos.z() << ")\n\n";
        
    std::cout << "Identification des zones de marche\n";
    ColorizedTag("ZoneMarche", mesh, outputDir + "scenario_ZoneMarche.off", robotPos);
    std::cout << "Surfaces planes (< 15° inclinaison) en VERT\n\n";
    
    std::cout << "Détection des obstacles\n";
    ColorizedTag("Obstacle", mesh, outputDir + "scenario_Obstacles.off", robotPos);
    std::cout << "Surfaces inclinées (> 30°) en ROUGE\n\n";
    
    std::cout << "Calcul des zones à portée d'interaction\n";
    ColorizedTag("ZoneAPortee", mesh, outputDir + "scenario_ZoneAPortee.off", robotPos);
    std::cout << "Zones accessibles en MAGENTA\n";
    
    std::cout << "Export du marqueur de position robot\n";
    save_robot_marker(robotPos, outputDir + "robot_marker.off", 0.3);
    std::cout << " Cube rouge de 0.3m exporté\n";
    
    std::cout << "Analyse complète. Fichiers de visualisation générés dans " 
              << outputDir << "\n";
    std::cout << "\n";
}

//==============================================================================
// STATISTIQUES RÉCAPITULATIVES
//==============================================================================

void afficherStatistiques(const Mesh& mesh, const Point_3& robotPos) {
    std::cout << "=== STATISTIQUES FINALES ===\n";
    
    auto areaMap = compute_Area_Map(mesh);
    auto normalMap = colorize_by_greatest_normal(mesh, "temp_normals.off");
    auto segmentMap = segmentationMap(areaMap, normalMap, mesh, robotPos);
    
    // Comptage par tag
    std::map<std::string, int> tagCounts;
    for (const auto& [facet, tagMaps] : segmentMap) {
        for (const auto& [tagName, isPresent] : tagMaps.tags) {
            if (isPresent) {
                tagCounts[tagName]++;
            }
        }
    }
    
    std::cout << "Nombre de faces par classe :\n";
    for (const auto& [tag, count] : tagCounts) {
        double percentage = 100.0 * count / mesh.size_of_facets();
        std::cout << "  " << tag << " : " << count 
                  << " faces (" << std::fixed << std::setprecision(1) 
                  << percentage << "%)\n";
    }
    
    std::cout << "\n";
}

//==============================================================================
// FONCTION PRINCIPALE
//==============================================================================

int main(int argc, char** argv) {
    std::cout << "\n";
    std::cout << "████████╗██████╗ ██████╗     ███████╗███████╗ ██████╗ \n";
    std::cout << "╚══██╔══╝██╔══██╗╚════██╗    ██╔════╝██╔════╝██╔════╝ \n";
    std::cout << "   ██║   ██████╔╝ █████╔╝    ███████╗█████╗  ██║  ███╗\n";
    std::cout << "   ██║   ██╔═══╝  ╚═══██╗    ╚════██║██╔══╝  ██║   ██║\n";
    std::cout << "   ██║   ██║     ██████╔╝    ███████║███████╗╚██████╔╝\n";
    std::cout << "   ╚═╝   ╚═╝     ╚═════╝     ╚══════╝╚══════╝ ╚═════╝ \n";
    std::cout << "Segmentation de maillage 3D - ISIMA ZZ3\n\n";
    
    Config config(argc, argv);
    config.print();
    
    Mesh mesh;
    if (!readOFF(config.meshPath, mesh)) {
        std::cerr << "ERREUR : Impossible de charger " << config.meshPath << "\n";
        return 1;
    }

    std::cout << "Maillage chargé : " << mesh.size_of_vertices() << " sommets, "
              << mesh.size_of_facets() << " faces\n\n";
    
    // ========== PARTIE I : PROPRIÉTÉS LOCALES ==========
    
    //partie_I0_Genre(mesh);
    
    //partie_I1_Aires(mesh, config.outputDir);
    
    //partie_I2_ExportColore(mesh, config.outputDir);//, "nefertiti_Area_Map.off");  
    
    //partie_I4_AutresMesures(mesh, config.outputDir);//, "nefertiti_Normal_Map.off");
    
    // ========== PARTIE II : SEGMENTATION ==========
    
    partie_II1_Segmentation(mesh, config.outputDir, config.robotPosition);
    
    partie_II2_ApplicationRobot(mesh, config.outputDir, config.robotPosition);
    
    // ========== RÉCAPITULATIF ==========
    
    afficherStatistiques(mesh, config.robotPosition);
    
    std::cout << "Tous les fichiers de sortie sont dans : " << config.outputDir << "\n";
    
    return 0;
}
