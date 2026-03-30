#pragma once
#include <igl/readOFF.h>
#include <Eigen/Dense>
#include <filesystem>
#include <iostream>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../../external/stb/stb_image_write.h"

namespace fs = std::filesystem;

// Définit un type d'image Row-Major (RGB entrelacés en mémoire)
typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, 3, Eigen::RowMajor> ImageRGB;

// Bresenham line
void draw_line(Eigen::Ref<ImageRGB> image, 
               int x0, int y0, int x1, int y1, int w, int h, unsigned char r, unsigned char g, unsigned char b) {
    int dx = std::abs(x1 - x0), dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    
    while (true) {
        if (x0 >= 0 && x0 < w && y0 >= 0 && y0 < h) {
            int idx = y0 * w + x0;
            if (idx >= 0 && idx < w * h) {
                image.row(idx) << r, g, b;
            }
        }
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

// Wireframe raster
ImageRGB raster_wireframe(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int w, int h) {
    ImageRGB image(w * h, 3);
    image.setConstant(150); // Fond gris uni
    
    double xmin = V.col(0).minCoeff(), xmax = V.col(0).maxCoeff();
    double ymin = V.col(1).minCoeff(), ymax = V.col(1).maxCoeff();
    
    // Scale UNIFORME à 80% de l'image
    double scale = 0.8 * std::min(w / (xmax - xmin + 1e-6), h / (ymax - ymin + 1e-6));
    
    // Centre géométrique absolu de l'écran
    double tx = (w / 2.0) - scale * ((xmax + xmin) / 2.0);
    double ty = (h / 2.0) - scale * ((ymax + ymin) / 2.0);
    
    for (int f = 0; f < F.rows(); ++f) {
        for (int e = 0; e < 3; ++e) {
            int i0 = F(f, e), i1 = F(f, (e + 1) % 3);
            int x0 = static_cast<int>(scale * V(i0, 0) + tx);
            int y0 = static_cast<int>(scale * V(i0, 1) + ty);
            int x1 = static_cast<int>(scale * V(i1, 0) + tx);
            int y1 = static_cast<int>(scale * V(i1, 1) + ty);
            draw_line(image, x0, y0, x1, y1, w, h, 0, 0, 0);  // jaune
        }
    }
    return image;
}

void batch_render_off(const std::string& input_dir, const std::string& output_dir, int w=1024, int h=1024) {
    fs::create_directories(output_dir);
    int count = 0;
    
    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (entry.path().extension() != ".off") continue;
        
        Eigen::MatrixXd V; Eigen::MatrixXi F;
        if (!igl::readOFF(entry.path().string(), V, F)) continue;
        
        auto image = raster_wireframe(V, F, w, h);
        
        std::string base_name = entry.path().stem().string();
        std::string png = output_dir + "/" + base_name + ".png";
        
        // stbi s'attend à du RowMajor : R0G0B0 R1G1B1...
        stbi_write_png(png.c_str(), w, h, 3, image.data(), w * 3);
        std::cout << "Rendu #" << ++count << ": " << png << std::endl;
    }
    std::cout << "Total rendus : " << count << std::endl;
}
