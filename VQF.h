#ifndef _VQF_H_
#define _VQF_H_

#include "EigenTypes.h"

#include <assert.h>
#include <limits>
#define _USE_MATH_DEFINES
#include <math.h>

template <size_t N>
class LowPass
{
public:
  LowPass()
    : tau_(0)
    , dt_(0)
    , a_(Vector2::Zero())
    , b_(Vector3::Zero())
    , state0_(VectorN<N>::Zero())
    , state1_(VectorN<N>::Zero())
    , avg_(VectorN<N>::Zero())
    , nInit_(0)
    , init_(false)
  {}

  LowPass(const ScalarType& tau, const ScalarType& dt)
    : LowPass()
  {
    assert(tau > 0);
    assert(dt > 0);
    tau_ = tau;
    dt_ = dt;
    const ScalarType fc = (M_SQRT2 / (2.0 * M_PI)) / tau_;
    const ScalarType C = std::tan(M_PI * fc * dt_);
    const ScalarType D = C * C + std::sqrt(2) * C + 1;
    const ScalarType b0 = C * C / D;
    a_ << 2 * (C * C - 1) / D, (1 - std::sqrt(2) * C + C * C) / D;
    b_ << b0, 2 * b0, b0;
  }
  
  ~LowPass() = default;

  VectorN<N> operator() (const VectorN<N>& x)
  {
    VectorN<N> y;
    if (!init_)
    {
      avg_ += x;
      y = avg_ / ++nInit_;
      if (nInit_ * dt_ >= tau_)
      {
        state0_ = (1 - b_(0)) * y;
        state1_ = (b_(2) - a_(1)) * y;
        init_ = true;
      }
      return y;
    }
    y = b_(0) * x + state0_;
    state0_ = b_(1) * x - a_(0) * y + state1_;
    state1_ = b_(2) * x - a_(1) * y;
    return y;
  }

protected:
  ScalarType dt_;
  ScalarType tau_;
  Vector2 a_;
  Vector3 b_;
  VectorN<N> state0_;
  VectorN<N> state1_;
  VectorN<N> avg_;
  size_t nInit_;
  bool init_;
};

class BasicVQF
{
public:
  BasicVQF(const ScalarType& dt, const ScalarType& taua = 3, const ScalarType& taum = 9)
    : dt_(dt)
    , qw_(Quaternion(1, 0, 0, 0))
    , qa_(Quaternion(1, 0, 0, 0))
    , deltam_(0)
    , taua_(taua)
    , taum_(taum)
    , km_(tau2Gain(taum_))
    , kmInit_(1)
    , aILp_(LowPass<3>(taua_, dt_))
  {}

  ~BasicVQF() = default;

  Quaternion operator() (const Vector3& w, const Vector3& a)
  {
    updateGyr(w);
    updateAcc(a);
    return qa_ * qw_;
  }

  Quaternion operator() (const Vector3& w, const Vector3& a, const Vector3& m)
  {
    updateGyr(w);
    updateAcc(a);
    updateMag(m);
    return rotateByHeadingDelta(qa_ * qw_, deltam_);
  }

protected:
  void updateGyr(const Vector3& w)
  {
    const ScalarType wN = w.norm();
    if (wN <= std::numeric_limits<ScalarType>::epsilon())
      return;

    const ScalarType phi = wN * dt_;
    const ScalarType c = std::cos(0.5 * phi);
    const ScalarType s = std::sin(0.5 * phi) / wN;
    qw_ = (qw_ * Quaternion(c, s * w.x(), s * w.y(), s * w.z())).normalized();
  }

  void updateAcc(const Vector3& a)
  {
    if (a.norm() <= std::numeric_limits<ScalarType>::epsilon())
      return;

    const Vector3 aI = aILp_(rotateByQuat(a, qw_)); // Transform acc to I frame, c.f. algo 1, line 9 and then low-pass, c.f. algo 1, line 10
    const Vector3 aE = rotateByQuat(aI, qa_).normalized(); // Transform low-passed acc to E frame, c.f. algo 1, line 11 and normalize, c.f. algo 1, line 12
    const ScalarType sqaw = std::sqrt(0.5 * (aE.z() + 1)); // Inclination correction update, c.f. algo 1, line 14 and normalize 
    qa_ = ((sqaw > 1e-6 ? Quaternion(sqaw, (0.5 / sqaw) * aE.y() / sqaw, -(0.5 / sqaw) * aE.x(), 0) : Quaternion(0, 1, 0, 0)) * qa_).normalized();
  }

  void updateMag(const Vector3& m)
  {
    if (m.norm() <= std::numeric_limits<ScalarType>::epsilon())
      return;

    const Vector3 mE = rotateByQuat(m, qa_ * qw_); // Transform mag measurement to earth frame, c.f. algo 1, line 17
    ScalarType k = km_;
    if (kmInit_ != 0) // Ensure fast initial convergence
    {
      if (k < kmInit_)
        k = kmInit_;
      kmInit_ = kmInit_ / (kmInit_ + 1);
      if (kmInit_ * taum_ < dt_)
        kmInit_ = 0;
    }
    deltam_ += k * wrapToPi(std::atan2(mE.x(), mE.y()) - deltam_); // Calculate heading offset, c.f. algo 1, line 18
    deltam_ = wrapToPi(deltam_); // Ensure it to stay within [-Pi, Pi]
  }

  ScalarType tau2Gain(const ScalarType& tau)
  {
    if (tau < 0)
      return 0;
    if (tau == 0)
      return 1;
    return 1 - std::exp(-dt_ / tau);
  }

  ScalarType wrapToPi(const ScalarType& phi)
  {
    if (phi > M_PI)
      return phi - 2 * M_PI;
    if (phi < -M_PI)
      return phi + 2 * M_PI;
    return phi;
  }

  Vector3 rotateByQuat(const Vector3& v, const Quaternion& q)
  {
    const Vector3 u = q.vec();
    const ScalarType w = q.w();
    return 2 * u.dot(v) * u + (w * w - u.dot(u)) * v + 2 * w * u.cross(v);
  }

  Quaternion rotateByHeadingDelta(const Quaternion& q, const ScalarType& delta)
  {
    const ScalarType c = std::cos(0.5 * delta);
    const ScalarType s = std::sin(0.5 * delta);
    return Quaternion(c * q.w() - s * q.z(), c * q.x() - s * q.y(), c * q.y() + s * q.x(), c * q.z() + s * q.w());
  }

  ScalarType dt_;
  Quaternion qw_;
  Quaternion qa_;
  ScalarType deltam_;
  ScalarType taua_;
  ScalarType taum_;
  ScalarType km_;
  ScalarType kmInit_;
  LowPass<3> aILp_;
};

#endif
