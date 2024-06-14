#include "Finds.hpp"

void FindByWave::BeforeWork () { // set null object
    work_points_.clear();
    max_RO_ = 0;
    label_.clear();
    binary_table_.clear();
    RO_.clear();
}

void FindByWave::BeforeStart () { // saving
    for (size_t i = 0; i < work_points_.size(); ++ i) {
        work_points_[i].SetLabel(0);
    }
    max_RO_ = 0;
    label_.resize(work_points_.size());
    binary_table_.resize(work_points_.size());
    RO_.resize(work_points_.size());
    for (size_t i = 0; i < work_points_.size(); ++ i) {
        binary_table_[i].resize(work_points_.size());
        RO_[i].resize(work_points_.size());
    }
}

FindByWave::FindByWave (const std::vector<Point>& new_points) { // constructor
    BeforeWork();
    work_points_ = new_points;
    BeforeStart();
}

void FindByWave::CreateRO () { // creating of matrix of distances
    double x1, y1, x2, y2;
    for (size_t i = 0; i < work_points_.size(); ++ i) {
        label_[i] = -1;
        x1 = work_points_[i].GetX();
        y1 = work_points_[i].GetY();
        for (size_t j = 0; j < work_points_.size(); ++ j) {
            x2 = work_points_[j].GetX();
            y2 = work_points_[j].GetY();
            RO_[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
            if (RO_[i][j] > max_RO_) {
                max_RO_ = RO_[i][j];
            }
        }
    }
}

void FindByWave::GenBinary (double threshold) { // creating of binary matrix
    for (size_t i = 0; i < work_points_.size(); ++ i) {
        for (size_t j = 0; j < work_points_.size(); ++ j) {
            if ((RO_[i][j] > 0) && (RO_[i][j] < (threshold + eps))) {
                binary_table_[i][j] = 1;
            } else {
                binary_table_[i][j] = 0;
            }
        }
    }
}

void FindByWave::Wave () { // wave clasters finding
    int tag = -1, contRO_lsumm, add_elem;
    std::vector<int> positions;
    for (size_t q = 0; q < work_points_.size(); q ++) {
        positions.clear();
        contRO_lsumm = 0; // If we check all strings, where first elems was included into cluster
        add_elem = 0; // If we add new element to cluster
        for (size_t i = q; i < work_points_.size(); ++ i) {
            if (label_[i] == -1) {
                ++ tag;
                positions.push_back(i);
                ++ add_elem;
                label_[i] = tag;
                for (size_t j = 0; j < work_points_.size(); ++ j) {
                    if (binary_table_[i][j] == 1) {
                        positions.push_back(j);
                        add_elem ++;
                    }
                }
                break;
            }
        }
        while (contRO_lsumm != add_elem) {
            for (size_t i = 0; i < positions.size(); ++ i) {
                ++ contRO_lsumm;
                for (size_t j = 0; j < work_points_.size(); ++ j) {
                    if (binary_table_[positions[i]][j] == 1) {
                        if (!IsInVector(positions, j)) {
                            ++ add_elem;
                            positions.push_back(j);
                        }
                    }
                }
            }
        }
        for (const auto &strnumber : positions) {
            label_[strnumber] = tag;
        }
    }
    for (size_t i = 0; i < work_points_.size(); ++ i) {
        work_points_[i].SetLabel(label_[i]);
        if(label_[i] < 0){
            Error(6);
        }
    }
}

std::vector<Point> FindByWave::FindClusters (double threshold) {
    BeforeStart();
    CreateRO();
    GenBinary(threshold);
    Wave();
    return work_points_;
}



void FindByKM::BeforeWork () { // clear object
    work_points_.clear();
    field_koord_.clear();
    center_koords_.clear();
    start_points_.clear();
    k = 0;
    optimal_k_ = 0;
    best_weight_ = 0;
    fsize = 0;
}

void FindByKM::BeforeStart () {
    for (unsigned int i = 0; i < work_points_.size(); ++ i) {
        work_points_[i].SetLabel(0);
    }
    field_koord_.clear();
    for (unsigned int i = 0; i < work_points_.size(); ++ i) {
        field_koord_.push_back(Point(work_points_[i]));
    }
}

FindByKM::FindByKM (const std::vector<Point>& new_koords) { // constructor
    BeforeWork();
    work_points_ = new_koords;
    BeforeStart();
}

void FindByKM::KMeans (int userK) { // kmeans with fixed K
    int min_id, q_start, q_end, ptr = 0;
    double min_val;
    bool change, esc = false;
    std::vector<double> mu_function;
    if (userK < 0) {
        q_start = 1;
        q_end = fsize;
    } else {
        q_start = userK;
        q_end = userK + 1;
    }
    mu_function.resize(fsize - 1);
    for (int q = q_start; (q < q_end) && (!esc); q ++) {
        std::vector<Point> mid_point;
        std::vector<int> number_of_points_in_each_cluster(k);
        std::vector<std::vector<double>> RO__selected, RO__centers;
        mu_function[ptr] = 0;
        k = q;
        change = true;
        start_points_.clear();
        start_points_.reserve(k);
        RO__selected.resize(k);
        center_koords_.resize(k);

        for (int i = 0; i < k; ++ i) { // create first matrix of ro between all points by pairs
            start_points_.push_back(Point(field_koord_[i]));
            mid_point.push_back(Point(0, 0, -5)); // -5 for mid points
            number_of_points_in_each_cluster[i] = 0;
            RO__selected[i].resize(fsize);
            for (int j = 0; j < fsize; j++) {
                RO__selected[i][j] = Correct(sqrt(Sqr(start_points_[i].GetX() - field_koord_[j].GetX())
                                                    + Sqr(start_points_[i].GetY() - field_koord_[j].GetY())));
            }
        }
        for (int j = 0; j < fsize; ++ j) { // mark each point (which center was the nearest)
            min_val = RO__selected[0][j];
            min_id = 0;
            for (int i = 1; i < k; ++ i) { // find minimal ro from j point to k ros of each cluster
                if (RO__selected[i][j] < min_val) {
                    min_val = RO__selected[i][j];
                    min_id = i;
                }
            }
            field_koord_[j].SetLabel(min_id);
            number_of_points_in_each_cluster[min_id] ++;
        }
        for (int i = 0; i < fsize; ++ i) { // caclulate mid points
            int id = field_koord_[i].GetLabel();
            mid_point[id].AddX(Correct(field_koord_[i].GetX() / number_of_points_in_each_cluster[id]));
            mid_point[id].AddY(Correct(field_koord_[i].GetY() / number_of_points_in_each_cluster[id]));

        }
        for (int i = 0; i < k; ++ i) { // make new koord centers
            center_koords_[i].SetX(mid_point[i].GetX());
            center_koords_[i].SetY(mid_point[i].GetY());
        }
        while (change) {
            RO__centers.resize(k);
            change = false;
            mid_point.clear();
            number_of_points_in_each_cluster.clear();
            mid_point.resize(k);
            number_of_points_in_each_cluster.resize(k);
            for (int i = 0; i < k; ++ i) { // calc all rors from centers
                RO__centers[i].resize(fsize);
                mid_point[i] = Point(0,0,-5);
                number_of_points_in_each_cluster[i] = 0;
                for (int j = 0; j < fsize; ++ j) { // find new rors from new centers to each point
                    RO__centers[i][j] = Correct(sqrt(Sqr(center_koords_[i].GetX() - field_koord_[j].GetX())
                                                    + Sqr(center_koords_[i].GetY() - field_koord_[j].GetY())));
                }
            }
            for (int j = 0; j < fsize; ++ j) {
                min_val = RO__centers[0][j];
                min_id = 0;
                for (int i = 1; i < k; ++ i) {
                    if (RO__centers[i][j] < min_val) {
                        min_val = RO__centers[i][j];
                        min_id = i;
                    }
                }
                field_koord_[j].SetLabel(min_id);
                number_of_points_in_each_cluster[min_id] ++;
            }
            for (int i = 0; i < fsize; ++ i) { // caclulate mid points
                int id = field_koord_[i].GetLabel();
                mid_point[id].AddX(Correct(field_koord_[i].GetX() / number_of_points_in_each_cluster[id]));
                mid_point[id].AddY(Correct(field_koord_[i].GetY() / number_of_points_in_each_cluster[id]));
            }
            for (int i = 0; i < k; ++ i) { // check definitions between new and old center koords
                double xOld = center_koords_[i].GetX(),
                    yOld = center_koords_[i].GetY(),
                    xNew = mid_point[i].GetX(),
                    yNew = mid_point[i].GetY();
                if ((ABS(xOld - xNew) > eps) || (ABS(yOld - yNew) > eps)) {
                    change = true;
                    i = k;
                }
            }
            if (change) {
                center_koords_.resize(k);
                for(int i = 0; i < k; ++ i){
                    center_koords_[i].SetX(mid_point[i].GetX());
                    center_koords_[i].SetY(mid_point[i].GetY());
                }
            }
        }
        for (int i = 0; i < k; ++ i) { // here we find optimal k
            for (int j = 0; j < fsize; ++ j) {
                if (field_koord_[j].GetLabel() != i) {
                    continue;
                }
                for (int p = j + 1; p < fsize; ++ p) {
                    if (field_koord_[p].GetLabel() == i) {
                        mu_function[ptr] += Correct(sqrt(Sqr(field_koord_[j].GetX() - field_koord_[p].GetX()) +
                                                        Sqr(field_koord_[j].GetY() - field_koord_[p].GetY())));
                    }
                }
            }
        }
        for (int i = 0; i < q; ++ i) { // here we find optimal k
            for (int j = i + 1; j < q; ++ j) {
                mu_function[ptr] += Correct(sqrt(Sqr(center_koords_[i].GetX() - center_koords_[j].GetX()) +
                                                Sqr(center_koords_[i].GetY() - center_koords_[j].GetY())));
            }
        }
        if ((q != 1) && ((q_end - q_start) > 1)) {
            if(mu_function[ptr - 1] < mu_function[ptr]){
                esc = true;
                ptr -= 2;
            }
        }
        ++ ptr;
    }
    optimal_k_ = ptr + 1;
    best_weight_ = mu_function[ptr];
    work_points_.clear();
    for(unsigned int i = 0; i < field_koord_.size(); ++ i){
        work_points_.push_back(field_koord_[i]);
    }
}

std::vector<Point> FindByKM::FindClusters (int user_number_of_clusters) {
    BeforeStart();
    KMeans(user_number_of_clusters);
    if(user_number_of_clusters == -1){
        std::cout << std::endl << "Optimal number of clusters is " << optimal_k_ << std::endl;
        BeforeStart();
        KMeans(optimal_k_);
    }
    return work_points_;
}