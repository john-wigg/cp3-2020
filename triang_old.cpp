#include "triang.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <Eigen/Dense>

#include <iostream>

namespace Triang {
    void HandleStartVertex(Vertex &v, std::set<HalfEdge> &T) {
        v.incident_edge->helper = &v;
        T.insert(*(v.incident_edge));
    }

    void HandleEndVertex(Vertex &v) {
        if (v.incident_edge->helper->type == VertexType::MERGE) {
            // TODO: Insert diagonal from v to helper
        } 
        // TODO
    }

    void HandleSplitVertex(Vertex &v) {
        // TODO
    }

    void HandleMergeVertex(Vertex &v) {
        // TODO
    }

    void HandleRegularVertex(Vertex &v) {
        // TODO
    }

    // See p. 53
    void MakeMonotone(DCEL &P) {
        // Construct priority queue of vertices using their y-coordinates
        // as priority.
        // If two points have the same y-coordinate, the one with smaller
        // x-coordinate has priority
        std::vector<Vertex *> Q;
        Q.resize(P.vertices.size());
        // First, fill the queue with references to all vertices of P
        for (int i = 0; i < P.vertices.size(); ++i) {
            Q[i] = &(P.vertices[i]);
        }

        // Now, sort the queue in order of ascending y coordinates
        std::sort(Q.begin(), Q.end(), [](const Vertex *a, const Vertex *b) { return a->coordinates(0) < b->coordinates(0); });
        // TODO: x-coordinates

        // TODO: Initialize an empty binary search tree T
        std::set<HalfEdge> T;

        while (!Q.empty()) {
            Vertex *v = Q.back();
            Q.pop_back();
            switch (v->type)
            {
                case VertexType::START:
                    HandleStartVertex(*v, T);
                    break;
                case VertexType::END:
                    HandleEndVertex(*v);
                    break;
                case VertexType::SPLIT:
                    HandleSplitVertex(*v);
                    break;
                case VertexType::MERGE:
                    HandleMergeVertex(*v);
                    break;
                case VertexType::REGULAR:
                    HandleRegularVertex(*v);
                    break;
            }
        }
    }
}