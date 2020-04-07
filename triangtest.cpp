#include "rupperts.hpp"

#include <vector>
#include <fstream>

using namespace Rupperts;

int main() {
    std::vector<Vertex *> poly;
    poly.resize(12);
    poly[0] = new Vertex(Eigen::Vector2f(0.0f, 0.0f));
    poly[1] = new Vertex(Eigen::Vector2f(-0.5f, -1.0f));
    poly[2] = new Vertex(Eigen::Vector2f(-1.5f, 0.0f));
    poly[3] = new Vertex(Eigen::Vector2f(-1.5f, -2.5f));
    poly[4] = new Vertex(Eigen::Vector2f(-0.5f, -2.0f));
    poly[5] = new Vertex(Eigen::Vector2f(0.0f, -3.0f));
    poly[6] = new Vertex(Eigen::Vector2f(1.0f, -2.0f));
    poly[7] = new Vertex(Eigen::Vector2f(2.0f, -1.0f));
    poly[8] = new Vertex(Eigen::Vector2f(2.0f, 2.0f));
    poly[9] = new Vertex(Eigen::Vector2f(-0.5f, 3.0f));
    poly[10] = new Vertex(Eigen::Vector2f(-1.5f, 3.0f));
    poly[11] = new Vertex(Eigen::Vector2f(0.0f, 2.0f));

    std::vector<Edge> segments;
    segments.push_back(Edge(*poly[0], *poly[1]));
    segments.push_back(Edge(*poly[1], *poly[2]));
    segments.push_back(Edge(*poly[2], *poly[3]));
    segments.push_back(Edge(*poly[3], *poly[4]));
    segments.push_back(Edge(*poly[4], *poly[5]));
    segments.push_back(Edge(*poly[5], *poly[6]));
    segments.push_back(Edge(*poly[6], *poly[7]));
    segments.push_back(Edge(*poly[7], *poly[8]));
    segments.push_back(Edge(*poly[8], *poly[9]));
    segments.push_back(Edge(*poly[9], *poly[10]));
    segments.push_back(Edge(*poly[10], *poly[11]));
    segments.push_back(Edge(*poly[11], *poly[0]));

    Delaunay2D delaunay;
    delaunay.vertices = poly;
    delaunay.segments = segments;
    //delaunay.DelaunayTriangulation();
    delaunay.RefineRupperts(19.9f);
    delaunay.ToFile("delaunay.out");
}