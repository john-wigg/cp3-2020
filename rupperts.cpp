#include "rupperts.hpp"

#include <algorithm>

namespace Rupperts {
    void RuppertsAlgorithm(PSLG X, float alpha) {
        // Add a bounding square B to X:
        // * Compute extremes of X
        float xmax = std::max_element(X.vertices.begin(), X.vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(0);
        float xmin = std::min_element(X.vertices.begin(), X.vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(0);
        float ymax = std::max_element(X.vertices.begin(), X.vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(1) < b.coordinates(1); })->coordinates(1);
        float ymin = std::min_element(X.vertices.begin(), X.vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(1) < b.coordinates(1); })->coordinates(1);
        float span = fmax(xmax-xmin, ymax-ymin);
        Eigen::Vector2f center((xmin+xmax)/2.0f, (ymin+ymax)/2.0f);
        // * Let B be the square of side 3 * span(X) centered on X
        // * Add the four boundary segments of B to X
        Vertex v1(Vertex(center + Eigen::Vector2f(1.5f * span, 1.5f * span)));
        Vertex v2(Vertex(center + Eigen::Vector2f(1.5f * span, 1.5f * span)));
        Vertex v3(Vertex(center + Eigen::Vector2f(1.5f * span, 1.5f * span)));
        Vertex v4(Vertex(center + Eigen::Vector2f(1.5f * span, 1.5f * span)));
        X.vertices.push_back(v1);
        X.vertices.push_back(v2);
        X.vertices.push_back(v3);
        X.vertices.push_back(v4);
        X.segments.push_back(Segment(v1, v2));
        X.segments.push_back(Segment(v2, v3));
        X.segments.push_back(Segment(v3, v4));
        X.segments.push_back(Segment(v4, v4));

        // TODO: Compute initial Delauney triangulation
    }
}