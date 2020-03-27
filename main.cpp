#include <iostream>
#include <fstream>
#include <vector>

#include "fem.hpp"

int main(int argc, char *argv[])
{
    if ( argc != 3 )
    {
        std::cout<<"usage: "<< argv[0] <<" <input file> <output file>\n";
        return 1;
    }
		
	std::ifstream infile(argv[1]);
	std::ofstream outfile(argv[2]);
    float poisson_ratio;
    float young_modulus;
    
    std::cout << "Load poisson ration and Young modulus." << std::endl;

    infile >> poisson_ratio >> young_modulus;

    std::cout << "Load node positions." << std::endl;

    int nodes_count;
    std::vector<float> nodes_x;
    std::vector<float> nodes_y;

    infile >> nodes_count;
    nodes_x.resize(nodes_count);
    nodes_y.resize(nodes_count);

    for (int i = 0; i < nodes_count; ++i) {
        infile >> nodes_x[i] >> nodes_y[i];
    }

    std::cout << "Load elements." << std::endl;

    int element_count;
    std::vector<Element> elements;
    infile >> element_count;
    elements.resize(element_count);

    for (int i = 0; i < element_count; ++i) {
        Element element;
        infile >> element.node_ids_[0] >> element.node_ids_[1] >> element.node_ids_[2];
        elements[i] = element;
    }

    std::cout << "Create solver object." << std::endl;

    FEMSolver solver(elements, nodes_x, nodes_y, poisson_ratio, young_modulus);

    std::cout << "Load constraints." << std::endl;

    int constraint_count;
    std::vector<Constraint> constraints;
    infile >> constraint_count;
    constraints.resize(constraint_count);

    for (int i = 0; i < constraint_count; ++i)
    {
        Constraint constraint;
        int type;
        infile >> constraint.node >> type;
        constraint.type = static_cast<Constraint::Type>(type);
        constraints[i] = constraint;
    }


    std::cout << "Load loads." << std::endl;

    int load_count;
    Eigen::VectorXf loads;
    infile >> load_count;
    loads.resize(2 * nodes_count);
    loads.setZero();

    for (int i = 0; i < load_count; ++i) {
        int node;
        float x, y;
        infile >> node >> x >> y;
        loads[2 * node + 0] = x;
        loads[2 * node + 1] = y;
    }


    std::cout << "Calculating elasticity matrix." << std::endl;
    solver.CalculateElasticityMatrix();

    std::cout << "Calculate stiffness matrix." << std::endl;
    solver.CalculateGlobalStiffnessMatrix();

    std::cout << "Setting constraints." << std::endl;
    solver.SetConstraints(constraints);
    solver.ApplyConstraints();

    std::cout << "Set loads." << std::endl;
    solver.SetLoads(loads);

    std::cout << "Solve." << std::endl;
    Eigen::VectorXf displacements = solver.Solve();

    std::cout << "Writing to file." << std::endl;
    outfile << displacements << std::endl;

    return 0;
}