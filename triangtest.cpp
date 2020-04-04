#include "triang.hpp"

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace Triang;

void polygonToDCEL(DCEL &out, const std::vector<Eigen::Vector2f> &in) {
    int N = in.size();
    out.vertices.resize(N);
    out.half_edges.resize(2 * N);
    out.faces.resize(2); // f1 outer, f2 inner

    for (int i = 0; i < N; ++i) {
        out.vertices[i].coordinates = in[i];
        out.vertices[i].incident_edge = &(out.half_edges[i]);
    }

    for (int i = 0; i < N; ++i) { // revered edges have +1 offset
        out.half_edges[2*i+0].origin = &(out.vertices[i+0]);
        out.half_edges[2*i+1].origin = &(out.vertices[i+1]);
        out.half_edges[2*i+0].twin = &(out.half_edges[2*i+1]);
        out.half_edges[2*i+1].twin = &(out.half_edges[2*i+0]);
        if (i < N - 1) {
            out.half_edges[2*i+0].next = &(out.half_edges[2*(i+1)+0]);
            out.half_edges[2*i+1].prev = &(out.half_edges[2*(i+1)+1]);
        } else {
            out.half_edges[2*i+0].next = &(out.half_edges[0]);
            out.half_edges[2*i+1].prev = &(out.half_edges[1]);
        }
        if (i > 0) {
            out.half_edges[2*i+0].prev = &(out.half_edges[2*(i-1)+0]);
            out.half_edges[2*i+1].next = &(out.half_edges[2*(i-1)+1]);
        } else {
            out.half_edges[2*i+0].prev = &(out.half_edges[2*(N-1)+0]);
            out.half_edges[2*i+1].next = &(out.half_edges[2*(N-1)+1]);
        }
        out.half_edges[2*i+0].incident_face = &(out.faces[1]); // inner
        out.half_edges[2*i+1].incident_face = &(out.faces[0]); // outer
    }

    out.faces[0].outer_component = nullptr;
    out.faces[0].inner_components.push_back(&(out.half_edges[0]));
    out.faces[1].outer_component = &(out.half_edges[1]);

    // Set vertex types
    for (int i = 0; i < N; ++i) {
        float y = out.vertices[i].coordinates(1);
        float x = out.vertices[i].coordinates(0);
        float x1 = out.vertices[i].incident_edge->prev->origin->coordinates(0);
        float y1 = out.vertices[i].incident_edge->prev->origin->coordinates(1);
        float x2 = out.vertices[i].incident_edge->twin->origin->coordinates(0);
        float y2 = out.vertices[i].incident_edge->twin->origin->coordinates(1);
        // Vertex "above" inner face (we assume that the polygon is defined counter-clockwise)
        if (x1 > x2) {
            if (y1 < y && y2 < y) out.vertices[i].type = VertexType::START;
            else if (y1 > y && y2 > y) out.vertices[i].type = VertexType::MERGE;
            else out.vertices[i].type = VertexType::REGULAR;
        } else { // Vertex "below" inner face
            if (y1 < y && y2 < y) out.vertices[i].type = VertexType::SPLIT;
            else if (y1 > y && y2 > y) out.vertices[i].type = VertexType::END;
            else out.vertices[i].type = VertexType::REGULAR;
        }
    }
}

void summary(const DCEL &P) {
    std::cout << std::setw(20) << "Vertex" << std::setw(20) << "Coordinates" << std::setw(20) << "IncidentEdge" << std::setw(20) << "Type" << std::endl;
    for (int i = 0; i < P.vertices.size(); ++i) {
        std::cout << std::setw(20) << &(P.vertices[i])
                  << std::setw(14) << P.vertices[i].coordinates(0) << ", " << std::setw(4) << P.vertices[i].coordinates(1)
                  << std::setw(20) << P.vertices[i].incident_edge
                  << std::setw(20) << P.vertices[i].type
        << std::endl;
    }

    std::cout << "==========================" << std::endl;

    std::cout << std::setw(20) << "Half-edge"  << std::setw(20) << "Origin" << std::setw(20) << "Twin" << std::setw(20) << "Incident Face" << std::setw(20) << "Next" << std::setw(20) << "Prev" << std::endl;
    for (int i = 0; i < P.half_edges.size(); ++i) {
        std::cout << std::setw(20) << &(P.half_edges[i])
                  << std::setw(20) << P.half_edges[i].origin
                  << std::setw(20) << P.half_edges[i].twin
                  << std::setw(20) << P.half_edges[i].incident_face
                  << std::setw(20) << P.half_edges[i].next
                  << std::setw(20) << P.half_edges[i].prev
        << std::endl;
    }
}

int main() {
    DCEL P;
    std::vector<Eigen::Vector2f> poly;
    poly.resize(3);
    poly[0] = Eigen::Vector2f(0.0f, 0.0f);
    poly[1] = Eigen::Vector2f(-1.0f, 0.5f);
    poly[2] = Eigen::Vector2f(-2.0f, -0.5f);
    polygonToDCEL(P, poly);
    summary(P);
    MakeMonotone(P);
    return 0;
}