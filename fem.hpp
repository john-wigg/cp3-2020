#include <array>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

class Element {
 private:
  Eigen::Matrix<float, 3, 6> B_;
 public:
  std::array<int, 3> node_ids_; // TODO: This should be private
  void CalculateStiffnessMatrix(const Eigen::Matrix3f& D, const std::vector<float>& nodes_x,
                                const std::vector<float>& nodes_y, std::vector<Eigen::Triplet<float> >& triplets);
};

struct Constraint
{
	enum Type
	{
		UX = 1 << 0,
		UY = 1 << 1,
		UXY = UX | UY
	};
	int node;
	Type type;
};

class FEMSolver {
 private:
  std::vector<Element> elements_;
  std::vector<float> nodes_x_;
  std::vector<float> nodes_y_;

  Eigen::VectorXf loads_;

  int nodes_count_;

  float poisson_ratio_;
  float young_modulus_;
  Eigen::Matrix3f D_;

  Eigen::SparseMatrix<float> K_;

  std::vector<Constraint> constraints_;

  void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator& it, int index);
 public:
  // Constructor
  FEMSolver(std::vector<Element> elements, const std::vector<float>& nodes_x, const std::vector<float>& nodes_y, float poisson_ratio, float yound_modulus);

  // Calculates and sets the elasticity matrix D_
  void CalculateElasticityMatrix();
  // Calculate global stiffness matrix
  void CalculateGlobalStiffnessMatrix();
  void ApplyConstraints();
  void SetLoads(Eigen::VectorXf loads) { loads_ = loads; }
  void SetConstraints(std::vector<Constraint> constraints) { constraints_ = constraints; }
  Eigen::VectorXf Solve();
};