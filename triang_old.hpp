#include <vector>
#include <Eigen/Dense>

namespace Triang
{
    enum VertexType {
        START,
        END,
        REGULAR,
        MERGE,
        SPLIT
    };
    // see p. 32
    struct HalfEdge;
    struct Vertex;
    struct Face;

    struct Vertex {
        Eigen::Vector2f coordinates; // x and y coordinates of the Vertex
        HalfEdge *incident_edge;     // Arbitrary half-edge that has the vertex as its origin
        VertexType type;             // Type of the vertex
    };

    struct Face {
        HalfEdge *outer_component; // Some half-edge on its outer boundary (nullptr if unbounded)
        std::vector<HalfEdge *> inner_components; // One half-edge on its inner boundary for each hole
    };

    struct HalfEdge {
        Vertex *origin;      // Origin vertex (no need to score destination since it's twin->origin)
        HalfEdge *twin;      // Twin half-edge (half-edge between same vertices but opposite direction)
        Face *incident_face; // Face that the edge bounds
        HalfEdge *next;      // Next edge on boundary of incident_face
        HalfEdge *prev;      // Previous edge on boundary of incident_face
        Vertex *helper;      // Used to store helper in MakeMonotone

        // Comparator for set
        bool operator<(const HalfEdge& rhs) const {
            return origin->coordinates(0) < rhs.origin->coordinates(0);
        }
    };

    // Doubly connected edge list
    struct DCEL {
        std::vector<Vertex> vertices;
        std::vector<Face> faces;
        std::vector<HalfEdge> half_edges;
    };

    void MakeMonotone(DCEL &P);
}