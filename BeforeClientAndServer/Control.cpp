#include "Control.hpp"

void Control::MakeLabel (int new_label) { // set group label
    for (size_t i = 0; i < num_point_group_; ++ i) {
        group_[i].SetLabel(new_label);
    }
}

std::vector<double> Control::CreateNorm (size_t number_of_points, double min_val, double max_val) const { // norm gisto
    if (number_of_points < 1) {
        Error(1);
        return {}; // Empty, default-constructed std::vector
    }

    if (min_val < max_val) {
        std::swap(min_val, max_val);
    }

    std::vector<double> result(number_of_points); // Create a std::vector of predefined size
    double sum, temp;
    int min = static_cast<int>(min_val) * 10,
        max = (static_cast<int>(max_val) + 1) * 10,
        range;
    range = max - min + 1;

    for (int i = 0; i < static_cast<int>(number_of_points); ++ i) {
        sum = 0;
        for (int j = 0; j < 1000; ++ j) {
            temp = ((std::rand () % range) + min) / 10.0;
            if ((temp < min_val) || (temp > max_val)) {
                -- j;
            } else {
                sum += temp;
            }
        }
        sum /= 1000;
        if ((sum < min_val) || (sum > max_val)) {
            -- i;
        } else {
            result[i] = sum;
        }
    }
    return result;
}

void Control::GenRnd (size_t number_of_points, double min_val, double max_val) { // ravn gisto
    if (number_of_points < 1) {
        Error(1);
        return void();
    }

    num_point1_ = number_of_points;
    vec1_.resize(num_point1_);
    for (size_t i = 0; i < num_point1_; ++ i) {
        int64_t bottlom_line = static_cast<int64_t>(100.0 * min_val),
            top_line = static_cast<int64_t>(100.0 * max_val);
        int64_t random_result = rand() % (bottlom_line - top_line) + bottlom_line;
        
        vec1_[i] = static_cast<double>(random_result + static_cast<int64_t>(100.0 * min_val)) / 100.0;
    }
}

void Control::GenNorm (size_t number_of_points, double min_val, double max_val) { // norm generation
    vec2_ = CreateNorm(num_point2_, min_val, max_val);
}

void Control::GenGroup (size_t number_of_points, double min_val_x, double max_val_x,
                        double min_val_y, double max_val_y, int label) { // group creating
    if (number_of_points < 1) {
        Error(2);
        return;
    }

    num_point_group_ = number_of_points;
    group_.resize(num_point_group_);
    std::vector<double> x_coord_list = CreateNorm(num_point_group_, min_val_x, max_val_x);
    std::vector<double> y_coord_list = CreateNorm(num_point_group_, min_val_y, max_val_y);
    for(size_t i = 0; i < num_point_group_; ++ i){
        group_[i] = Point(x_coord_list[i], y_coord_list[i], label);
    }
}

void Control::GenGroup (const std::vector<Point>& new_values) { // group creating copy
    num_point_group_ = new_values.size();
    group_ = new_values;
}

void Control::FileRavn () const { // printing ravn std::vector in file
    std::fstream text, script;
    script.open("plotravn.plt", std::fstream::trunc); // writing script to the file
    script << "width=1" << std::endl << "bin(x, s) = s*int(x/s) + width/2"
            << std::endl << "set boxwidth width" << std::endl << "plot 'ravn.txt' u (bin($1,width)):(1.0) \\"
            << std::endl << "s f w boxes fs solid 0.5 title 'Ravn Gisto'" << std::endl;
    script.close();
    text.open("ravn.txt", std::fstream::trunc); // writing values to the file
    for (size_t i = 0; i < num_point1_; ++ i) {
        text << vec1_[i] << std::endl;
    }
    text.close();
}

void Control::FileNorm () const { // printing norm std::vector in file
    std::fstream text, script;
    script.open("plotnorm.plt", std::fstream::trunc);
    script << "width=0.1" << std::endl << "bin(x, s) = s*int(x/s) + width/2"
        << std::endl << "set boxwidth width" << std::endl << "plot 'norm.txt' u (bin($1,width)):(1.0) \\"
        << std::endl << "s f w boxes fs solid 0.5 title 'Norm Gisto'" << std::endl;
    script.close();
    text.open("norm.txt", std::fstream::trunc);
    for (size_t i = 0; i < num_point2_; ++ i) {
        text << vec2_[i] << std::endl;
    }
    text.close();
}

void Control::FileGroup () const { // printing group in file
    std::fstream text, script;
    script.open("plotgroup.plt", std::fstream::trunc);
    script << "plot 'group.txt'";
    script.close();
    text.open("group.txt", std::fstream::trunc); //base file
    for (size_t i = 0; i < num_point_group_; ++ i) {
        text << group_[i].GetX() << " " << group_[i].GetY() << std::endl;
    }
    text.close();
}

void Control::TurnNULL (double phi) { // rotarion group relative (0;0)
    double newx, newy;
    std::vector<Point> newvect;
    newvect.reserve(num_point_group_);
    for (size_t i = 0; i < num_point_group_; ++ i) {
        newx = group_[i].GetX();
        newy = group_[i].GetY();
        newvect.push_back(Point(
            newx * std::cos(phi) - newy * std::sin(phi),
            newx * std::sin(phi) + newy * std::cos(phi)));
    }
    group_.clear();
    GenGroup(newvect);
}

void Control::TurnCenter (double phi) { // rotarion group relative to group center
    double midx = 0, midy = 0, newx, newy;
    std::vector<Point> new_points_list(num_point_group_);
    for (size_t i = 0; i < num_point_group_; ++ i) {
        midx += group_[i].GetX();
        midy += group_[i].GetY();
    }
    midx /= static_cast<double>(num_point_group_);
    midy /= static_cast<double>(num_point_group_);
    for (size_t i = 0; i < num_point_group_; ++ i) {
        newx = group_[i].GetX() - midx;
        newy = group_[i].GetY() - midy;
        new_points_list[i] = Point(
            (newx * std::cos(phi) - newy * std::sin(phi)) + midx,
            (newx * std::sin(phi) + newy * std::cos(phi)) + midy);
    }
    group_.clear();
    GenGroup(new_points_list);
}

void Control::MoveX (double delta_x) { // X axis moving
    for (size_t i = 0; i < num_point_group_; ++ i) {
        group_[i].SetX(group_[i].GetX() + delta_x);
    }
}

void Control::MoveY (double delta_y) { // Y axis moving
    for (size_t i = 0; i < num_point_group_; ++ i) {
        group_[i].SetY(group_[i].GetY() + delta_y);
    }
}