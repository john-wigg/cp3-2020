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

    void Delaunay2D::DelaunayTriangulation() {
        std::cout << "Starting Delaunay Triangulation..." << std::endl;
        std::cout << "Add super triangle..." << std::endl;
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

        vertices.push_back(v1);
        vertices.push_back(v2);
        vertices.push_back(v3);

        int N = vertices.size();

        // Add super-triangle to the triangulation
        Triangle *t = new Triangle();
        t->v = { N - 3, N - 2, N - 1 };
        t->c = CalculateCircle(*t);
        triangles.push_back(*t);

        std::cout << "* ADDED super triangle with edges (" << t->v[0] << ", " << t->v[1] << ", " << t->v[2] << ")" << std::endl;

        std::cout << "Start adding points..." << std::endl;

        for (int i = 0; i < N-3; ++i) {
            std::cout << "=========================================================" << std::endl;
            std::cout << "Add Point " << i << " of " << N-3 << " at " << vertices[i].coordinates.x() << ", " << vertices[i].coordinates.y() << std::endl;
            std::cout << "Number of triangles: " << triangles.size() << std::endl;
            // Add point 

            std::cout << "Mark bad triangles..." << std::endl;
            int last_bad_index = -1;
            int num_bad_triangles = 0;

            // First find alle the triangles that are no longer valid due to the insrtion
            for (int k = 0; k < triangles.size(); ++k) {
                // Add triangles that are no longer valid to bad triangles
                if (inCircle(triangles[k].get().c, vertices[i])) {
                    std::cout << "* MARKING bad triangle: " << k << " (" << triangles[k].get().v[0] << ", " << triangles[k].get().v[1] << ", " << triangles[k].get().v[2] << ")" << std::endl;
                    triangles[k].get().is_bad = true;
                    last_bad_index = k;
                    num_bad_triangles++;
                }
            }

            std::cout << "Get boundary edge loop of bad triangles..." << std::endl;
            // Get boundary edge loop of bad triangles
            std::vector<Edge> boundary;
            std::vector<Triangle *> boundary_triangles;

            // start with first triangle in bad triangles
            Triangle *T = &triangles[last_bad_index].get();
            int edge = 0;
            //while (i == 5){};
            while (true) {
                Triangle *tri_op = T->t[edge];
                std::cout << "-------" << std::endl;
                std::cout << "* Triangle " << T << " has vertices " << T->v[0] << ", "  << T->v[1] << ", " << T->v[2] << std::endl;
                std::cout << "On Edge " << edge << "with opposite triangle " << tri_op << std::endl;
                if (tri_op == nullptr || !tri_op->is_bad) {
                    // bad_triangles contains tri_op
                    std::cout << "* ADDING edge to boundary..." << std::endl;
                    boundary.push_back(Edge(T->v[(3 + (edge+1) % 3) % 3], T->v[(3 + (edge-1) % 3) % 3]));
                    std::cout << "b: " << boundary.back().orig << " --> " << boundary.back().dest << std::endl;
                    boundary_triangles.push_back(tri_op);
                    edge = (edge + 1) % 3;

                    // If loop is closed
                    if (boundary.front().orig == boundary.back().dest) {
                        std::cout << "* Boundary CLOSED with vertex " << boundary.back().dest << std::endl;
                        break;
                    }
                } else {
                    std::cout << "* MOVE to next CCW edge" << std::endl;
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

            std::cout << "Remove bad triangles..." << std::endl;
            // Remove bad triangles
            // TODO: This has n^2 complexity
            for (int k = 0; k < triangles.size(); ++k) {
                if (triangles[k].get().is_bad) {
                    std::cout << "* REMOVING bad triangle: " << k << " (" << triangles[k].get().v[0] << ", " << triangles[k].get().v[1] << ", " << triangles[k].get().v[2] << ")" << std::endl;
                    delete &(triangles[k].get());
                    triangles.erase(triangles.begin() + k);
                    --k;
                    /*
                    for (int m = 0; m < triangles.size(); ++m) { // Update triangle refs of 
                        if (triangles[m].get().t[0] > k && triangles[m].get().t[0] != -1) --triangles[m].get().t[0];
                        if (triangles[m].get().t[1] > k && triangles[m].get().t[1] != -1) --triangles[m].get().t[1];
                        if (triangles[m].get().t[2] > k && triangles[m].get().t[2] != -1) --triangles[m].get().t[2];
                    }
                    

                    for (int m = 0; m < boundary_triangles.size(); ++m) {
                        if (boundary_triangles[m] > k && boundary_triangles[m] != -1) --boundary_triangles[m];
                    }
                    */
                }
            }

            std::cout << "Boundary edge loop has length " << boundary.size() << " and is ";
            for (int k = 0; k < boundary.size(); ++k) {
                std::cout << boundary[k].orig << " --> " << boundary[k].dest << " ";
            }
            std::cout << std::endl;
            std::cout << "With boundary triangles ";
            for (int k = 0; k < boundary_triangles.size(); ++k) {
                std::cout << boundary_triangles[k] << ", ";
            }
            std::cout << std::endl;

            std::cout << "Remaining triangles:" << std::endl;
            for (int k = 0; k < triangles.size(); ++k) {
                std::cout << triangles[k].get().v[0] << ", " << triangles[k].get().v[1] << ", " << triangles[k].get().v[2] << std::endl;
            }

            std::cout << "Retriangulate hole..." << std::endl;
            // Re-triangulate the hole
            std::vector<std::reference_wrapper<Triangle>> new_triangles;
            for (int k = 0; k < boundary.size(); ++k) {
                Triangle *new_T = new Triangle();
                new_T->v = { i, boundary[k].orig, boundary[k].dest };
                new_T->c = CalculateCircle(*new_T);
                new_T->t[0] = boundary_triangles[k];

                // Link boundary triangle (if existing) to this triangle
                if (boundary_triangles[k] != nullptr) {
                    std::cout << "Reconnecting boundary triangle " << boundary_triangles[k] << std::endl;
                    std::cout << "Currently " << boundary_triangles[k]->t[0] << ", " << boundary_triangles[k]->t[1] << ", " << boundary_triangles[k]->t[2] << std::endl;
                    for (int m = 0; m < 3; ++m) {
                        /*
                        int neigh = triangles[boundary_triangles[k]].t[m];
                        if (neigh != -1) {
                            std::cout << "Testing " << neigh << std::endl;
                            bool flag1 = false;
                            bool flag2 = false;
                            for (int n = 0; n < 3; n++) {
                                if (triangles[neigh].v[n] == boundary[k].orig) {
                                    flag1 = true;
                                }
                                if (triangles[neigh].v[n] == boundary[k].dest) {
                                    flag2 = true;
                                }
                            }

                            if (flag1 && flag2) {
                                triangles[boundary_triangles[k]].t[m] = triangles.size() + k; // This will be the index of the new triangle
                                std::cout << "* CONNECTING triangle " << triangles.size() + k  << " to " << " boundary triangle " << boundary_triangles[k] << std::endl;
                            }
                        }
                        */
                        if (boundary_triangles[k]->v[m] == boundary[k].orig) {
                            boundary_triangles[k]->t[positive_modulo(m-2, 3)] = new_T;
                            std::cout << "* CONNECTING triangle " << triangles.size() + k  << " to " << " boundary triangle " << boundary_triangles[k] << std::endl;
                        }
                    }
                }
                new_triangles.push_back(*new_T);
            }

            // Link new triangles to each other and add to list
            for (int k = 0; k < new_triangles.size(); ++k) {
                new_triangles[k].get().t[1] = &new_triangles[positive_modulo((k + 1), new_triangles.size())].get();
                new_triangles[k].get().t[2] = &new_triangles[positive_modulo((k - 1), new_triangles.size())].get();
                triangles.push_back(new_triangles[k]);
            }

            std::cout << "Final triangles:" << std::endl;
            for (int k = 0; k < triangles.size(); ++k) {
                std::cout << triangles[k].get().v[0] << ", " << triangles[k].get().v[1] << ", " << triangles[k].get().v[2] << std::endl;
            }
        }

        
        // TODO: Cleanup of the super-triangle
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
    }

    void Delaunay2D::ToFile(std::string oname) {
        std::ofstream outfile(oname);
        for (int i = 0; i < vertices.size(); ++i) {
            outfile << vertices[i].coordinates(0) << " " << vertices[i].coordinates(1) << std::endl;
        }

        outfile << std::endl;

        for (int i = 0; i < triangles.size(); ++i) {
            outfile << triangles[i].get().v[0] << " " << triangles[i].get().v[1] << " " << triangles[i].get().v[2] << std::endl;
        }
        outfile << std::endl;

        /*

        for (int i = 0; i < triangles.size(); ++i) {
            outfile << triangles[i].t[0] << " " << triangles[i].t[1] << " " << triangles[i].t[2] << std::endl;
        }

        outfile << std::endl;

        for (int i = 0; i < triangles.size(); ++i) {
            outfile << triangles[i].c.center(0) << " " << triangles[i].c.center(1) << " " << triangles[i].c.radius << std::endl;
        }

        outfile << std::endl;

        for (int i = 0; i < triangles.size(); ++i) {
            outfile << triangles[i].t[0] << " " << triangles[i].t[1] << " " << triangles[i].t[2] << std::endl;
        }
        */
    }

    Delaunay2D::~Delaunay2D() {
        for (int i = 0; i < triangles.size(); ++i) {
            delete &triangles[i].get();
        }
    }

/*
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
        DelauneyTriangulation(X.vertices);
    }
    */
}