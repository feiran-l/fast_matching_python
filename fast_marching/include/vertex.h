#ifndef VERTEX_H
#define VERTEX_H

#include "include.h"



/* The vertex class represents a 3D point and implements 3D vector operations. */
class vertex {
public:
    double x;
    double y;
    double z;

public:
    vertex(const double &x_ = 0, const double &y_ = 0, const double &z_ = 0);
    ~vertex() = default;
    double &operator[](const index_t &i);
    const double &operator[](const index_t &i) const;
    vertex unit() const;
    double operator*() const;                        // norm
    double operator,(const vertex &v) const;        // dot product
    vertex operator*(const vertex &v) const;        // cross product
    vertex operator/(const double &v) const;        // scalar division
    vertex operator+(const vertex &v) const;
    vertex operator-(const vertex &v) const;
    vertex operator-() const;

    void operator*=(const double &v);            // scalar produc
    void operator/=(const double &v);
    void operator+=(const vertex &v);
    void operator-=(const vertex &v);

    bool operator<(const vertex &v);
    bool operator==(const vertex &v);
};

vertex operator*(const double &a, const vertex &v);

std::ostream &operator<<(std::ostream &os, const vertex &v);

std::istream &operator>>(std::istream &is, vertex &v);



#endif // VERTEX_H

