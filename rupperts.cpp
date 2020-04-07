#include "rupperts.hpp"

#include <algorithm>
#include <fstream>
#include <string>
#include <iostream>

namespace Rupperts {
    int positive_modulo(int i, int n) {
        return (i % n + n) % n;
    }

    // https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    Circle Delaunay2D::CalculateCircle(const Triangle &t) {
        Eigen::Vector2f a = vertices[t.v[0]].coordinates;
        Eigen::Vector2f b = vertices[t.v[1]].coordinates;
        Eigen::Vector2f c = vertices[t.v[2]].coordinates;
        float ba2 = (b - a).squaredNorm();
        float ca2 = (c - a).squaredNorm();
        float bc2 = (b - c).squaredNorm();
        Eigen::Matrix2f A1, A2, A3;
        A1 << b(0) - a(0), b(1) - a(1),
              c(0) - a(0), c(1) - a(1);
        A2 << b(1) - a(1), ba2,
              c(1) - a(1), ca2;
        A3 << b(0) - a(0), ba2,
              c(0) - a(0), ca2;
        float s = sqrt(ba2 * ca2 * bc2);
        float det1 = 2 * A1.determinant();
        float det2 = A2.determinant();
        float det3 = A3.determinant();

        Circle cc;
        cc.radius = s / det1;
        cc.center(0) = a(0) - det2 / det1;
        cc.center(1) = a(1) + det3 / det1;

        return cc;
    }

    bool Delaunay2D::inCircle(const Circle &c, const Vertex &v) {
        // Fast but not very robust, see http://www.cs.cmu.edu/~quake/robust.html
        return (v.coordinates - c.center).norm() <= c.radius;
    }

    // Add vertex at index i to triangulation while maintaining Delaunay condition
    // (Vertex must be INSIDE triangles of other vertices)
    void Delaunay2D::DelaunayAddPoint(int i) {
        // Add point 

        int last_bad_index = -1;
        int num_bad_triangles = 0;

        // First find alle the triangles that are no longer valid due to the insrtion
        for (int k = 0; k < triangles.size(); ++k) {
            // Add triangles that are no longer valid to bad triangles
            if (inCircle(triangles[k]->c, vertices[i])) {
                triangles[k]->is_bad = true;
                last_bad_index = k;
                num_bad_triangles++;
            }
        }

        // Get boundary edge loop of bad triangles
        std::vector<Edge> boundary;
        std::vector<Triangle *> boundary_triangles;

        // start with first triangle in bad triangles
        Triangle *T = triangles[last_bad_index];
        int edge = 0;
        //while (i == 5){};
        while (true) {
            Triangle *tri_op = T->t[edge];
            if (tri_op == nullptr || !tri_op->is_bad) {
                // bad_triangles contains tri_op
                boundary.push_back(Edge(T->v[(3 + (edge+1) % 3) % 3], T->v[(3 + (edge-1) % 3) % 3]));
                boundary_triangles.push_back(tri_op);
                edge = (edge + 1) % 3;

                // If loop is closed
                if (boundary.front().orig == boundary.back().dest) {
                    break;
                }
            } else {
                // Move to next CCW edge in opposite triangle
                int index;
                for (int m = 0; m < 3; ++m) { // Find edge that borders to original triangle
                    if (tri_op->t[m] == T) {
                        index = m;
                    }
                }
                edge = (index + 1) % 3;
                T = tri_op;
            }
        }
        // Remove bad triangles
        for (int k = 0; k < triangles.size(); ++k) {
            if (triangles[k]->is_bad) {
                delete triangles[k];
                triangles.erase(triangles.begin() + k);
                --k;
            }
        }

        // Re-triangulate the hole
        std::vector<Triangle *> new_triangles;
        for (int k = 0; k < boundary.size(); ++k) {
            Triangle *new_T = new Triangle();
            new_T->v = { i, boundary[k].orig, boundary[k].dest };
            new_T->c = CalculateCircle(*new_T);
            new_T->t[0] = boundary_triangles[k];

            // Link boundary triangle (if existing) to this triangle
            if (boundary_triangles[k] != nullptr) {
                for (int m = 0; m < 3; ++m) {
                    if (boundary_triangles[k]->v[m] == boundary[k].orig) {
                        boundary_triangles[k]->t[positive_modulo(m-2, 3)] = new_T;
                    }
                }
            }
            new_triangles.push_back(new_T);
        }

        // Link new triangles to each other and add to list
        for (int k = 0; k < new_triangles.size(); ++k) {
            new_triangles[k]->t[1] = new_triangles[positive_modulo((k + 1), new_triangles.size())];
            new_triangles[k]->t[2] = new_triangles[positive_modulo((k - 1), new_triangles.size())];
            triangles.push_back(new_triangles[k]);
        }
    }

    void Delaunay2D::DelaunayTriangulation() {
        // Add super-triangle
        float xmax = std::max_element(vertices.begin(), vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(0);
        float xmin = std::min_element(vertices.begin(), vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(0);
        float ymax = std::max_element(vertices.begin(), vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(1);
        float ymin = std::min_element(vertices.begin(), vertices.end(), [](const Vertex& a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(1);
        float span = fmax(xmax-xmin, ymax-ymin);
        Eigen::Vector2f center((xmin+xmax)/2.0f, (ymin+ymax)/2.0f);
        
        // Add points CCW
        Vertex v1(Vertex(center + Eigen::Vector2f(0.0f, 3.0f * span)));
        Vertex v2(Vertex(center + Eigen::Vector2f(-3.0f * span, -3.0f * span)));
        Vertex v3(Vertex(center + Eigen::Vector2f(3.0f * span, -3.0f * span)));
        v1.is_helper = true;
        v2.is_helper = true;
        v3.is_helper = true;

        vertices.push_back(v1);
        vertices.push_back(v2);
        vertices.push_back(v3);

        int N = vertices.size();

        // Add super-triangle to the triangulation
        Triangle *t = new Triangle();
        t->v = { N - 3, N - 2, N - 1 };
        t->c = CalculateCircle(*t);
        triangles.push_back(t);

        for (int i = 0; i < N-3; ++i) {
            DelaunayAddPoint(i);
        }
        
        // TODO: Cleanup of the super-triangle
        /*
        for (int i = 0; i < triangles.size(); ++i) {
            for (int k = 0; k < 3; ++k) {
                if (triangles[i].get().v[k] >= (N - 3)) {
                    triangles.erase(triangles.begin() + i);
                    --i;
                    break;
                }
            }
        }

        vertices.pop_back();
        vertices.pop_back();
        vertices.pop_back();
        */
    }

    void Delaunay2D::ToFile(std::string oname) {
        std::ofstream outfile(oname);
        for (int i = 0; i < vertices.size(); ++i) {
            outfile << vertices[i].coordinates(0) << " " << vertices[i].coordinates(1) << std::endl;
        }

        outfile << std::endl;

        for (int i = 0; i < triangles.size(); ++i) {
            outfile << triangles[i]->v[0] << " " << triangles[i]->v[1] << " " << triangles[i]->v[2] << std::endl;
        }
        outfile << std::endl;

        for (int i = 0; i < segments.size(); ++i) {
            outfile << segments[i].orig << " " << segments[i].dest << std::endl;
        }
        outfile << std::endl;
    }

    Delaunay2D::~Delaunay2D() {
        for (int i = 0; i < triangles.size(); ++i) {
            delete triangles[i];
        }
    }

    void Delaunay2D::SplitTri(Triangle &t, Vertex p) {
        vertices.push_back(p);
        DelaunayAddPoint(vertices.size() - 1);
    }

    void Delaunay2D::SplitSeg(int s_idx) {
        Edge s = segments[s_idx];
        Eigen::Vector2f m = ( vertices[s.orig].coordinates + vertices[s.dest].coordinates ) / 2.0f; // midpoint
        vertices.push_back(Vertex(m));
        DelaunayAddPoint(vertices.size() - 1);

        segments.erase(segments.begin() + s_idx);
        segments.push_back(Edge(s.orig, vertices.size() - 1));
        segments.push_back(Edge(vertices.size() - 1, s.dest));
    }

    int Delaunay2D::GetFirstEncroached() {
        for (int i = 0; i < segments.size(); ++i) {
            Circle c;
            c.center = ( vertices[segments[i].orig].coordinates + vertices[segments[i].dest].coordinates ) / 2.0f;
            c.radius = ( vertices[segments[i].orig].coordinates - vertices[segments[i].dest].coordinates ).norm() / 2.0f;
            for (int k = 0; k < vertices.size(); ++k) {
                if (inCircle(c, vertices[k]) && k != segments[i].orig && k != segments[i].dest) {
                    return i;
                }
            }
        }
        return -1;
    }

    Triangle *Delaunay2D::GetFirstBadTriangle(float alpha) {
        for (int i = 0; i < triangles.size(); ++i) {
            if (IsTriangleLowQuality(*triangles[i], alpha))
            {
                return triangles[i];
            }
        }
        return nullptr;
    }

    bool Delaunay2D::IsTriangleLowQuality(const Triangle &t, float alpha) {
        if (vertices[t.v[1]].is_helper) return false;
        if (vertices[t.v[0]].is_helper) return false;
        if (vertices[t.v[2]].is_helper) return false;
        for (int k = 0; k < 3; ++k) {
            Eigen::Vector2f a = vertices[t.v[(k+1)%3]].coordinates - vertices[t.v[k]].coordinates;
            Eigen::Vector2f b = vertices[t.v[(k+1)%3]].coordinates - vertices[t.v[(k+2)%3]].coordinates;
            float dot = a.x() * b.x() + a.y() + b.y();
            float det = a.x() * b.y() - a.y() * b.x();
            float angle = fabs(atan2(det, dot)) * 180.0 / 3.1415926535;
            if (angle < alpha) {
                return true;
            }
        }
        return false;
    }

    void Delaunay2D::RefineRupperts(float alpha) {
        float xmax = std::max_element(vertices.begin(), vertices.end(), [](const Vertex&a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(0);
        float xmin = std::min_element(vertices.begin(), vertices.end(), [](const Vertex&a, const Vertex& b){ return a.coordinates(0) < b.coordinates(0); })->coordinates(0);
        float ymax = std::max_element(vertices.begin(), vertices.end(), [](const Vertex&a, const Vertex& b){ return a.coordinates(1) < b.coordinates(1); })->coordinates(1);
        float ymin = std::min_element(vertices.begin(), vertices.end(), [](const Vertex&a, const Vertex& b){ return a.coordinates(1) < b.coordinates(1); })->coordinates(1);
        float span = fmax(xmax-xmin, ymax-ymin);
        Eigen::Vector2f center((xmin+xmax)/2.0f, (ymin+ymax)/2.0f);
        // * Let B be the square of side 3 * span(X) centered on X
        // * Add the four boundary segments of B to X
        Vertex v1(Vertex(center + Eigen::Vector2f(1.5f * span, 1.5f * span)));
        Vertex v2(Vertex(center + Eigen::Vector2f(-1.5f * span, 1.5f * span)));
        Vertex v3(Vertex(center + Eigen::Vector2f(-1.5f * span, -1.5f * span)));
        Vertex v4(Vertex(center + Eigen::Vector2f(1.5f * span, -1.5f * span)));
        vertices.push_back(v1);
        vertices.push_back(v2);
        vertices.push_back(v3);
        vertices.push_back(v4);
        int N = vertices.size();
        segments.push_back(Edge(N-4, N-3));
        segments.push_back(Edge(N-3, N-2));
        segments.push_back(Edge(N-2, N-1));
        segments.push_back(Edge(N-1, N-4));

        DelaunayTriangulation();

        int first_encroached = GetFirstEncroached();

        Triangle *first_bad_triangle = GetFirstBadTriangle(alpha);

        int it = 0;
        while (first_encroached >= 0 || first_bad_triangle != nullptr) {
            if (first_encroached >= 0) {
                SplitSeg(first_encroached);
            } else {
                // Get first bad triangle
                first_bad_triangle = GetFirstBadTriangle(alpha);

                if (first_bad_triangle != nullptr) {
                    Vertex p;
                    p.coordinates = CalculateCircle(*first_bad_triangle).center;
                    bool encroached = false;
                    for (int i = 0; i < segments.size(); ++i) {
                        Circle c;
                        c.center = ( vertices[segments[i].orig].coordinates + vertices[segments[i].dest].coordinates ) / 2.0f;
                        c.radius = ( vertices[segments[i].orig].coordinates - vertices[segments[i].dest].coordinates ).norm() / 2.0f;
                        if (inCircle(c,p)) { // Circumcenter of bad triangle encroaches a segment
                            SplitSeg(i);
                            encroached = true;
                        }
                    }

                    if (!encroached) {
                        SplitTri(*first_bad_triangle, p);
                    }
                }
            }

            // Get encroached segment
            first_encroached = GetFirstEncroached();

            ++it;
            if (it > 1000) {
                std::cout << "WARNING: Did not converge in time!" << std::endl;
                return;
            }
        }
    }  
}