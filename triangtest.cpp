#include "rupperts.hpp"

#include <vector>
#include <fstream>

using namespace Rupperts;

int main() {
    std::vector<Vertex> poly;
    poly.resize(10);
    poly[0].coordinates = Eigen::Vector2f(0.0f, 0.0f);
    poly[1].coordinates = Eigen::Vector2f(-0.5f, -1.0f);
    poly[2].coordinates = Eigen::Vector2f(-1.5f, 0.0f);
    poly[3].coordinates = Eigen::Vector2f(-1.5f, -2.5f);
    poly[4].coordinates = Eigen::Vector2f(-0.5f, -2.0f);
    poly[5].coordinates = Eigen::Vector2f(0.0f, -3.0f);
    poly[6].coordinates = Eigen::Vector2f(1.0f, -2.0f);
    poly[7].coordinates = Eigen::Vector2f(1.0f, -1.0f);
    poly[8].coordinates = Eigen::Vector2f(3.0f, -2.5f);
    poly[9].coordinates = Eigen::Vector2f(0.0f, 0.5f);

    Delaunay2D delaunay;
    delaunay.vertices = poly;
    delaunay.DelaunayTriangulation();
    delaunay.ToFile("delaunay.out");
}