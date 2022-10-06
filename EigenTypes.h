#ifndef _EIGENTYPES_H_
#define _EIGENTYPES_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

typedef double ScalarType;

typedef Eigen::Matrix<ScalarType, 2, 2> Matrix2x2;
typedef Eigen::Matrix<ScalarType, 3, 3> Matrix3x3;
typedef Eigen::Matrix<ScalarType, 4, 4> Matrix4x4;
typedef Eigen::Matrix<ScalarType, 3, 4> Matrix3x4;
typedef Eigen::Matrix<ScalarType, 4, 3> Matrix4x3;

typedef Eigen::Matrix<ScalarType, 2, 1> Vector2;
typedef Eigen::Matrix<ScalarType, 3, 1> Vector3;
typedef Eigen::Matrix<ScalarType, 4, 1> Vector4;
typedef Eigen::Matrix<ScalarType, 4, 1> Vector5;
typedef Eigen::Matrix<ScalarType, 4, 1> Vector6;
typedef Eigen::Matrix<ScalarType, 4, 1> Vector7;
typedef Eigen::Matrix<ScalarType, 4, 1> Vector8;
typedef Eigen::Matrix<ScalarType, 4, 1> Vector9;
typedef Eigen::Matrix<ScalarType, 4, 1> Vector10;
template <size_t N>
using VectorN = Eigen::Matrix<ScalarType, N, 1>;

typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> Vector;

typedef Eigen::AngleAxis<ScalarType> AngleAxis;
typedef Eigen::Quaternion<ScalarType> Quaternion;

#endif
