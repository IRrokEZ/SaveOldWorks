#pragma once

#include <memory>

#include "Control.hpp"

class Group{
private:
    std::unique_ptr<Control> group_; // array of groups
    int group_label_;
public:
    Group () = default; // Default constructor, required for std::vector::resize

    Group (int new_group_size, double new_min_x, double new_min_y,
           double new_max_x, double new_max_y, int new_group_label); // constructor

    Group (std::unique_ptr<Control> new_group); // constructor (means that label included into Point' vector)

    Group (std::unique_ptr <Control> new_group, int new_label); // constructor

    Group (const Group &new_group); // copy constructor

    Group (std::unique_ptr <Group> &new_group); // copy constructor

    int GetGroupLabel () const {
        return group_label_;
    }

    std::unique_ptr<Control> GetFieldGroup () const; // returning of a group
    
    void Regrupp (const std::unique_ptr <Control> &newgroup); // we modified this object by new Control object

    void Regrupp (const std::unique_ptr <Control> &newgroup, int new_label); // we modified this object by new Control object

    void ReMakeLabel (); // remake label_ of a group

    void ReMakeLabel (int LB); // remake label_ of a group

    ~Group () = default; // default destructor
};