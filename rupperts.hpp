#include <Eigen/Dense>
#include <vector>
#include <array>

namespace Rupperts {
    struct Vertex;
    struct Segment;

    struct Vertex {
        Vertex() { };
        Vertex(Eigen::Vector2f init_coordinates) { coordinates = init_coordinates; }
        Eigen::Vector2f coordinates; // x and y coordinates of the Vertex
        //Segment &incident_segment; // Could be useful later?
    };

    struct Edge {
        Edge() { };
        Edge(int init_orig, int init_dest) {
            orig = init_orig;
            dest = init_dest;
        }
        int orig;
        int dest;
    };

    struct Segment {
        Segment() { };
        Segment(Vertex &init_orig, Vertex &init_dest) {
            orig = &init_orig;
            dest = &init_dest;
        }
        Vertex *orig;
        Vertex *dest;
    };

    struct Circle {
        Circle() { };
        Circle(Eigen::Vector2f init_center, float init_radius) {
            center = init_center;
            radius = init_radius;
        }

        Eigen::Vector2f center;
        float radius;
    };

    /* Stores vertices, connected triangles and circle of a single triangle */
    struct Triangle {
        std::array<int, 3> v{-1, -1, -1};
        std::array<int, 3> t{-1, -1, -1};
        Circle c;
        bool is_bad = false;
    };

    class Delaunay2D {
    private:
        Circle CalculateCircle(const Triangle &t);
        bool inCircle(const Circle &c, const Vertex &v);
    public:
        std::vector<Vertex> vertices;
        std::vector<Triangle> triangles;
        void DelaunayTriangulation();
        void ToFile(std::string oname);
    };

    /*
    struct PSLG {
        std::vector<Vertex> vertices;
        std::vector<Segment> segments;
    };
    */

    std::vector<Triangle> DelauneyTriangulation(std::vector<Vertex> &P);
    //void RuppertsAlgorithm(PSLG X, float alpha);
}