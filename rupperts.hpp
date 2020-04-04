#include <Eigen/Dense>
#include <vector>

namespace Rupperts {
    struct Vertex;
    struct Segment;

    struct Vertex {
        Vertex();
        Vertex(Eigen::Vector2f init_coordinates) { coordinates = init_coordinates; }
        Eigen::Vector2f coordinates; // x and y coordinates of the Vertex
        //Segment &incident_segment; // Could be useful later?
    };

    struct Segment {
        Segment();
        Segment(Vertex &init_orig, Vertex &init_dest) {
            orig = &init_orig;
            dest = &init_dest;
        }
        Vertex *orig;
        Vertex *dest;
    };

    struct PSLG {
        std::vector<Vertex> vertices;
        std::vector<Segment> segments;
    };

    void DelauneyTriangulation(std::vector<Eigen::Vector2f> P);
    void RuppertsAlgorithm(PSLG X, float alpha);
}