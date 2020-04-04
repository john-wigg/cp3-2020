#include "triang.hpp"

#include <vector>
#include <algorithm>
#include <Eigen/Dense>

#include <iostream>

namespace Triang {
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

        for (int i = 0; i < Q.size(); ++i) {
            std::cout << Q[i]->coordinates << std::endl;
        }
        // TODO: x-coordinates

        // TODO: Initialize an empty binary search tree T

        while (!Q.empty()) {
            Vertex *v = Q.back();
            Q.pop_back();
            // TODO: Handle Vertex
        }
    }
}