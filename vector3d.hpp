#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>

// 3d vector class definition:
class vec
{
public:
    // Initialize to 0,0,0
    vec();
    // Initialize to the components of another vector
    vec(const vec &v);
    // Initialize to x, y, z
    vec(double x, double y, double z);
    // Set to x, y, z
    void set(double x, double y, double z);

    // Assign to the values of another vector
    vec &operator=(const vec &v);
    // Add to this vector, the values of another vector.
    vec &operator+=(const vec &v);
    // Subtract from this vector the values of another vector.
    vec &operator-=(const vec &v);
    // Multiply and store the values of this vector by a double.
    vec &operator*=(double f);
    // Divide and store the values of this vector by a double.
    vec &operator/=(double f);
    // Return the values of this vector plus another
    vec operator+(const vec &v) const;
    // Return the values of this vector minus another
    vec operator-(const vec &v) const;
    // Return a copy of this vector, negated
    vec operator-() const;
    // Return the values of this vector times a double
    vec operator*(const double f) const;
    // Return the values of this vector divided by a double
    vec operator/(const double f) const;

    // Formatted output function for displaying this vector's components
    void output(std::ostream &os) const;

    // Return the length of this vector
    double length() const;

    // Return the squared length of this vector
    double lengthsq() const;

    // Make this vector a unit-vector in the same direction
    vec normalize() const;

    // Cross product between this vector and another vector v
    vec cross(const vec &v) const;

    // Dot product between this vector and another vector v
    double dot(const vec &v) const;

    //  Angle between this vector another vector v
    double angle(const vec &v) const;

    double x() const;
    double y() const;
    double z() const;

private:
    // Components of the 3d vector
    double x_;
    double y_;
    double z_;
};

std::ostream &operator<<(std::ostream &os, const vec &v); // Overload insertion operator for displaying vectors
vec operator*(double f, const vec &rhs); // Overload * for doubles to allow for the more natural mutiplication ordering f*vec. Above we handled the case vec*f.

#endif  // VECTOR3D_H