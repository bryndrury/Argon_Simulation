#include <iostream>
#include <cmath>
#include "vector3d.hpp"

using std::cout, std::endl;

///////////////////////////
// Implement the vec class
///////////////////////////
vec::vec():
        x_(0),
        y_(0),
        z_(0)
{
}

vec::vec(const vec &v):
        x_(v.x_),
        y_(v.y_),
        z_(v.z_)
{
}

vec::vec(double x_init, double y_init, double z_init):
        x_(x_init),
        y_(y_init),
        z_(z_init)
{
}

void vec::set(double x_set, double y_set, double z_set)
{
    x_ = x_set;
    y_ = y_set;
    z_ = z_set;
}

vec &vec::operator=(const vec &v)
{
    x_ = v.x_;
    y_ = v.y_;
    z_ = v.z_;
    return *this;
}

vec &vec::operator+=(const vec &v)
{
    x_ += v.x_;
    y_ += v.y_;
    z_ += v.z_;
    return *this;
}

vec &vec::operator-=(const vec &v)
{
    x_ -= v.x_;
    y_ -= v.y_;
    z_ -= v.z_;
    return *this;
}

vec &vec::operator*=(double f)
{
    x_ *= f;
    y_ *= f;
    z_ *= f;
    return *this;
}

vec &vec::operator/=(double f)
{
    x_ /= f;
    y_ /= f;
    z_ /= f;
    return *this;
}

vec vec::operator+(const vec &v) const
{
    return vec(
            x_ + v.x_,
            y_ + v.y_,
            z_ + v.z_);
}

vec vec::operator-(const vec &v) const
{
    return vec(
            x_ - v.x_,
            y_ - v.y_,
            z_ - v.z_);
}

vec vec::operator-() const
{
    return vec(
            -x_,
            -y_,
            -z_);
}


vec vec::operator*(const double f) const
{
    return vec(
            x_ * f,
            y_ * f,
            z_ * f);
}

vec vec::operator/(const double f) const
{
    return vec(
            x_ / f,
            y_ / f,
            z_ / f);
}

void vec::output(std::ostream &os) const
{
    os << x_ << "," << y_ << "," << z_;
}

std::ostream &operator<<(std::ostream &os, const vec &v)
{
    v.output(os);
    return os;
}

vec operator*(double f, const vec &rhs)
{
    return rhs * f; // Return vec * f when f * vec is used
}

double vec::length() const
{
    return sqrt(lengthsq());
}

double vec::lengthsq() const
{
    return dot(*this); // Dot product with itself
}

vec vec::normalize() const
{
    double len = length();
    return *this/len;
}

vec vec::cross(const vec &v) const
{
    vec ret;
    // Compute components of u x v
    ret.x_ = y_ * v.z_ - v.y_ * z_;
    ret.y_ = z_ * v.x_ - v.z_ * x_;
    ret.z_ = x_ * v.y_ - v.x_ * y_;
    return ret;
}

double vec::dot(const vec &v) const
{
    return x_*v.x_ + y_*v.y_ + z_*v.z_;
}

double vec::angle(const vec &v) const
{
    double udotv = dot(v); // Compute u.v where u = this vector
    double u_len = length(); // Compute |u|
    double v_len = v.length(); // Compute |v|
    return acos(udotv/(u_len*v_len)); // Invert u.v = |u|*|v|*cos(theta)
}

double vec::x() const
{
    return x_;
}
double vec::y() const
{
    return y_;
}
double vec::z() const
{
    return z_;
}