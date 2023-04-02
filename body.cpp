#include <iostream>
#include <cmath>
#include "vector3d.hpp"
#include "body.hpp"

using std::cout, std::endl;

///////////////////////////
// Implement the body class
///////////////////////////
body::body():
        mass_(1), pos_(vec(0,0,0)), vel_(vec(0,0,0)), acc_(vec(0,0,0))
{
}

body::body(const body &b):
        name_(b.name_), mass_(b.mass_), pos_(b.pos_), vel_(b.vel_), acc_(vec(0,0,0))
{
}

body::body(std::string n, double m, vec pos, vec vel):
        name_(n), mass_(m), pos_(pos), vel_(vel), acc_(vec(0,0,0))
{
}

void body::set(std::string n, double m, vec pos, vec vel)
{
    name_ = n;
    mass_= m;
    pos_ = pos;
    vel_ = vel;
    acc_ = vec(0,0,0);
}

body &body::operator=(const body &b)
{
    mass_ = b.mass_;
    pos_ = b.pos_;
    vel_ = b.vel_;
    acc_ = b.acc_;
    return *this;
}

void body::output(std::ostream &os) const
{
    os << name_ << ": m = " << mass_ << ", pos = (" << pos_ << "), vel = (" << vel_ << ")";
}

// Implement get() and set() functions:

std::string body::get_name() const
{
    return name_;
}

double body::get_mass() const
{
    return mass_;
}

vec body::get_pos() const
{
    return pos_;
}

vec body::get_vel() const
{
    return vel_;
}

vec body::get_acc() const
{
    return acc_;
}

vec body::get_L() const
{
    return L_;
}

double body::get_ke() const
{
    return ke_;
}

double body::get_gpe() const
{
    return gpe_;
}

void body::set_name(const std::string n)
{
    name_ = n;
}

void body::set_mass(const double m)
{
    mass_ = m;
}

void body::set_pos(const vec &v)
{
    pos_ = v;
}

void body::set_vel(const vec &v)
{
    vel_ = v;
}

void body::set_acc(const vec &v)
{
    acc_ = v;
}

void body::set_L(const vec &v)
{
    L_ = v;
}

void body::set_ke(const double E)
{
    ke_ = E;
}

void body::set_gpe(const double E)
{
    gpe_ = E;
}

double body::distance(const body &b) const
{
    vec r = pos_ - b.pos_; // Vector connecting the two bodies
    return r.length(); // Return length of vector
}

vec body::direction(const body &b) const
{
    vec r = pos_ - b.pos_; // Vector connecting the two bodies
    return r.normalize(); // Return normalized vector
}

vec body::angular_momentum(const vec &point) const
{
    vec r = pos_ - point; // Vector connecting body to point
    return mass_*r.cross(vel_); // Return m r x v
}

std::ostream &operator<<(std::ostream &os, const body &b)
{
    b.output(os);
    return os;
}