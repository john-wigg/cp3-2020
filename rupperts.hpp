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
        bool is_helper = false;
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
        std::array<Triangle *, 3> t{nullptr, nullptr, nullptr};
        Circle c;
        bool is_bad = false;
    };

    class Delaunay2D {
    private:
        void SplitSeg(int s);
        void SplitTri(Triangle &t, Vertex p);
        int GetFirstEncroached();
        bool IsTriangleLowQuality(const Triangle &t, float alpha);
        Circle CalculateCircle(const Triangle &t);
        bool inCircle(const Circle &c, const Vertex &v);
        void DelaunayAddPoint(int i);
    public:
        ~Delaunay2D();
        std::vector<Vertex> vertices;
        std::vector<std::reference_wrapper<Triangle>> triangles;
        std::vector<Edge> segments;
        void DelaunayTriangulation();
        void ToFile(std::string oname);
        void RefineRupperts(float alpha);
    };
}