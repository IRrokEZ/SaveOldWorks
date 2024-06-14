#include "Group.hpp"

Group::Group (int new_group_size, double new_min_x, double new_min_y,
              double new_max_x, double new_max_y, int new_group_label) { // constructor
    group_.reset(new Control());
    group_label_ = new_group_label;
    group_->GenGroup(new_group_size, new_min_x, new_max_x, new_min_y, new_max_y, new_group_label);
}

Group::Group (std::unique_ptr<Control> new_group) { // constructor (means that label included into Point' vector)
    group_label_ = new_group->GetGroupLabel();
    group_.reset(new Control());
    group_->GenGroup(new_group->GetGroup());
}

Group::Group (std::unique_ptr <Control> new_group, int new_label) { // constructor
    group_label_ = new_label;
    group_.reset(new Control());
    group_->GenGroup(new_group->GetGroup());
    group_->MakeLabel(new_label);
}

Group::Group (const Group &new_group) { // copy constructor
    std::unique_ptr<Control> GR = new_group.GetFieldGroup();
    std::vector<Point> arr = GR -> GetGroup();
    group_.reset(new Control());
    group_label_ = new_group.GetGroupLabel();
    group_->GenGroup(arr);
}

Group::Group (std::unique_ptr <Group> &new_group) { // copy constructor
    std::unique_ptr <Control> work_group = new_group->GetFieldGroup();
    std::vector<Point> work_points_array = work_group->GetGroup();
    group_.reset(new Control());
    group_label_ = new_group->GetGroupLabel();
    group_->GenGroup(work_points_array);
}

std::unique_ptr<Control> Group::GetFieldGroup () const { // returning of a group
    std::unique_ptr <Control> result_group(new Control());
    std::vector<Point> work_points_array = group_ -> GetGroup();
    result_group -> GenGroup(work_points_array);
    return result_group;
}

void Group::Regrupp (const std::unique_ptr <Control> &newgroup) { // we modified this object by new Control object
    std::vector<Point> work_points_array = newgroup->GetGroup();
    group_.reset(new Control());
    group_->GenGroup(work_points_array);
}

void Group::Regrupp (const std::unique_ptr <Control> &newgroup, int new_label) { // we modified this object by new Control object
    std::vector<Point> work_points_array = newgroup -> GetGroup();
    group_.reset(new Control());
    group_->GenGroup(work_points_array);
    group_label_ = new_label;
    group_->MakeLabel(new_label);
}

void Group::ReMakeLabel () {//remake label_ of a group
    group_->MakeLabel(group_label_);
}

void Group::ReMakeLabel (int LB) {//remake label_ of a group
    group_label_ = LB;
    group_->MakeLabel(group_label_);
}