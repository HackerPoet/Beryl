#pragma once
#pragma warning(disable : 26451)
#pragma warning(disable : 26495)
#pragma warning(disable : 26812)
#pragma warning(push, 0)
#include <Eigen/Dense>
#pragma warning(pop)
#include <map>
#include <vector>

using Eigen::Vector2f;
using Eigen::Vector3f;
using Eigen::Vector3d;
using Eigen::VectorXf;
using Eigen::Matrix3f;
using Eigen::Matrix3d;

using Verts2D = std::vector<Vector2f>;
using Verts3D = std::vector<Vector3f>;
using Face = std::vector<int>;
using Faces = std::vector<Face>;
using Edge = std::pair<int, int>;
using Edges = std::vector<Edge>;

using FaceMap = std::map<int, std::vector<int>>;
using EdgeMap = std::map<Edge, std::vector<int>>;
