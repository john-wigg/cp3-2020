#include "rupperts.hpp"

#include <vector>
#include <fstream>

using namespace Rupperts;

int main() {
    std::vector<Vertex> poly;
    poly.resize(12);
    poly[0].coordinates = Eigen::Vector2f(0.0f, 0.0f);
    poly[1].coordinates = Eigen::Vector2f(-0.5f, -1.0f);
    poly[2].coordinates = Eigen::Vector2f(-1.5f, 0.0f);
    poly[3].coordinates = Eigen::Vector2f(-1.5f, -2.5f);
    poly[4].coordinates = Eigen::Vector2f(-0.5f, -2.0f);
    poly[5].coordinates = Eigen::Vector2f(0.0f, -3.0f);
    poly[6].coordinates = Eigen::Vector2f(1.0f, -2.0f);
    poly[7].coordinates = Eigen::Vector2f(2.0f, -1.0f);
    poly[8].coordinates = Eigen::Vector2f(2.0f, 2.0f);
    poly[9].coordinates = Eigen::Vector2f(-0.5f, 3.0f);
    poly[10].coordinates = Eigen::Vector2f(-1.5f, 3.0f);
    poly[11].coordinates = Eigen::Vector2f(0.0f, 2.0f);

    std::vector<Edge> segments;
    segments.push_back(Edge(0, 1));
    segments.push_back(Edge(1, 2));
    segments.push_back(Edge(2, 3));
    segments.push_back(Edge(3, 4));
    segments.push_back(Edge(4, 5));
    segments.push_back(Edge(5, 6));
    segments.push_back(Edge(6, 7));
    segments.push_back(Edge(7, 8));
    segments.push_back(Edge(8, 9));
    segments.push_back(Edge(9, 10));
    segments.push_back(Edge(10, 11));
    segments.push_back(Edge(11, 0));

    Delaunay2D delaunay;
    delaunay.vertices = poly;
    delaunay.segments = segments;
    //delaunay.DelaunayTriangulation();
    delaunay.RefineRupperts(19.9f);
    delaunay.ToFile("delaunay.out");
}