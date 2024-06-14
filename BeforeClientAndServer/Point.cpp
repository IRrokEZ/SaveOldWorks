#include "Point.hpp"

Point& Point::operator= (const Point& other) { // remake coordinates
    x_ = other.x_;
    y_ = other.y_;
    label_ = other.label_;
    return *this;
}