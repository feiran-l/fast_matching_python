#include "vertex.h"



vertex::vertex(const double &x_, const double &y_, const double &z_) {
    x = x_;
    y = y_;
    z = z_;
}

double &vertex::operator[](const index_t &i) {
    return (&x)[i];
}

const double &vertex::operator[](const index_t &i) const {
    return (&x)[i];
}

vertex vertex::unit() const {
    return *this / **this;
}

double vertex::operator*() const {
    return sqrt(x * x + y * y + z * z);
}

double vertex::operator,(const vertex &v) const {
    return x * v.x + y * v.y + z * v.z;
}

vertex vertex::operator*(const vertex &v) const {
    return vertex(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

vertex vertex::operator/(const double &a) const {
    return (1. / a) * (*this);
}

vertex vertex::operator+(const vertex &v) const {
    return vertex(x + v.x, y + v.y, z + v.z);
}

vertex vertex::operator-(const vertex &v) const {
    return vertex(x - v.x, y - v.y, z - v.z);
}

vertex vertex::operator-() const {
    return vertex(-x, -y, -z);
}

void vertex::operator*=(const double &a) {
    x *= a;
    y *= a;
    z *= a;
}

void vertex::operator/=(const double &a) {
    (*this) *= (1. / a);
}

void vertex::operator+=(const vertex &v) {
    x += v.x;
    y += v.y;
    z += v.z;
}

void vertex::operator-=(const vertex &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
}

bool vertex::operator<(const vertex &v) {
    if (x != v.x) return x < v.x;
    if (y != v.y) return y < v.y;
    return z < v.z;
}

bool vertex::operator==(const vertex &v) {
    return x == v.x && y == v.y && z == v.z;
}

vertex operator*(const double &a, const vertex &v) {
    return vertex(a * v.x, a * v.y, a * v.z);
}

ostream &operator<<(ostream &os, const vertex &v) {
    os << v.x << " " << v.y << " " << v.z;
    return os;
}

istream &operator>>(istream &is, vertex &v) {
    is >> v.x >> v.y >> v.z;
    return is;
}

