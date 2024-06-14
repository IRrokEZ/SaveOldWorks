#pragma once

#include <fstream>
#include <string>
#include <vector>

class Point { // points
private:
    double x_, y_; // coordinates
    int label_; // label

public:
    Point () = default; //Default constructor
    Point (double new_x, double new_y) : x_(new_x), y_(new_y) {} // Initialize fields with coordinates
    Point (double new_x, double new_y, int new_label) // Initialize fields with new arguments
        : x_(new_x), y_(new_y), label_(new_label) {} 
    Point (const Point &other) 
        : x_(other.x_), y_(other.y_), label_(other.label_) {} // Copy constructor

    double GetX () const { // get X
        return x_;
    }
    double GetY () const { // get Y
        return y_;
    }
    int GetLabel () const { // get label_
        return label_;
    }

    void SetX (double new_x) { // set X coordinate
        x_ = new_x;
    }
    void SetY (double new_y) { // set Y coordinate
        y_ = new_y;
    }
    void SetLabel (int new_label) { // set label
        label_ = new_label;
    }

    void AddX (double x) { // move X coordinate
        x_ += x;
    }
    void AddY (double y) { // move Y coordinate
        y_ += y;
    }

    Point& operator= (const Point& other); // remake coordinates
    
    ~Point () = default; // default 
};