#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "globals.hpp"
#include "Point.hpp"

class Control {
private:
    std::vector<double> vec1_, vec2_; //vectors of coordinates
    std::vector<Point> group_; //vector of points
    size_t num_point1_, num_point2_, num_point_group_;
public:
    Control () = default; //default constructor

    int GetGroupSize () const {
        return num_point_group_;
    }

    int GetGroupLabel () const {
        return (group_.empty() ? 0 : group_[0].GetLabel());
    }
    
    std::vector<double> GetRavn () const { //ravn generation cpy
        return vec1_;
    }

    std::vector<double> GetNorm () const { //norm generation cpy
        return vec2_;
    }

    std::vector<Point> GetGroup () const { //group cpy
        return group_;
    }

    void MakeLabel (int new_label); //group label

    std::vector<double> CreateNorm (size_t number_of_points, double min_val, double max_val) const; //norm gisto
    
    void GenRnd (size_t number_of_points, double min_val, double max_val); //ravn gisto
    
    void GenNorm (size_t number_of_points, double min_val, double max_val); //norm generation
    
    void GenGroup (size_t number_of_points, double min_val_x, double max_val_x, double min_val_y, double max_val_y, int label); //group creating

    void GenGroup (const std::vector<Point>& new_values); //group creating copy
    
    void FileRavn () const; //printing ravn std::vector in file

    void FileNorm () const; //printing norm std::vector in file

    void FileGroup () const; //printing group in file
    
    void TurnNULL (double phi); //rotarion group relative (0;0)

    void TurnCenter (double phi); //rotarion group relative to group center

    void MoveX (double delta_x); //X axis moving
    
    void MoveY (double delta_y); //Y axis moving

    ~Control () = default; // default destructor
};