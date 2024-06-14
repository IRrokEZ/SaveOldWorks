#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include "globals.hpp"
#include "Point.hpp"

class FindByWave{
private:
    std::vector<int> label_; // labels
    std::vector<Point> work_points_; // current points
    std::vector<std::vector<int>> binary_table_; // binary matrix
    std::vector<std::vector<double>> RO_; // distances between points
    double max_RO_;
public:
    void BeforeWork (); // set null object

    void BeforeStart (); // saving

    FindByWave (const std::vector<Point>& new_points); // constructor

    void CreateRO (); // creating of matrix of distances

    void GenBinary (double threshold); // creating of binary matrix
    
    void Wave (); // wave clasters finding
    
    std::vector<Point> FindClusters (double threshold);

    ~FindByWave () = default;
};

class FindByKM{ // Refactor later
private:
    int k, fsize, optimal_k_;
    double best_weight_;
    std::vector<Point> work_points_, field_koord_, center_koords_, start_points_;
public:
    void BeforeWork (); // clear object

    void BeforeStart ();

    FindByKM (const std::vector<Point>& new_koords); // constructor

    void KMeans (int userK); // kmeans with fixed K

    std::vector<Point> FindClusters (int user_number_of_clusters);

    ~FindByKM () = default;
};

class FindBySPTR{
    private:
        std::vector<Point> work_points_;
        int launchNumber;
        std::vector<std::vector<int>> binary_table_;//binary matrix
        std::vector<double> allLengths;
        std::vector<std::vector<double>> RO_;//distances between points
        double max_RO_;
        std::vector<int> label_;//label_s
    public:
        void BeforeWork(){//default constructor
            std::fstream f;
            launchNumber = -1;
            f.open("launch.log");
            f >> launchNumber;
            if(this -> launchNumber <= 1){
                f.close();
                f.open("launch.log", std::fstream::trunc);
                f << 1;
            }
            f.close();
            work_points_.clear();
            allLengths.clear();
            max_RO_ = 0;
            label_.clear();
            binary_table_.clear();
            RO_.clear();
        }
        void BeforeStart(){//saving
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                work_points_[i].SetLabel(0);
            }
            allLengths.clear();
            max_RO_ = 0;
            label_.resize(this -> work_points_.size());
            binary_table_.resize(this -> work_points_.size());
            RO_.resize(this -> work_points_.size());
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                binary_table_[i].resize(this -> work_points_.size());
                RO_[i].resize(this -> work_points_.size());
            }
        }
        FindBySPTR(const std::vector<Point>& koords){//constructor
            BeforeWork();
            work_points_ = koords;
            BeforeStart();
        }
        void CreateRO_(){//creating of matrix of distances
            double x1, y1, x2, y2;
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                label_[i] = -1;
                x1 = work_points_[i].GetX();
                y1 = work_points_[i].GetY();
                for(unsigned int j = 0; j < work_points_.size(); ++ j){
                    x2 = work_points_[j].GetX();
                    y2 = work_points_[j].GetY();
                    RO_[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                    if(this -> RO_[i][j] > max_RO_){
                        max_RO_ = RO_[i][j];
                    }
                }
            }
        }
        void GenBinary(double threshold){//creating of binary matrix
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                for(unsigned int j = 0; j < work_points_.size(); ++ j){
                    if((this -> RO_[i][j] > 0) && (this -> RO_[i][j] < (threshold + eps))){
                        binary_table_[i][j] = 1;
                    } else {
                        binary_table_[i][j] = 0;
                    }
                }
            }
        }
        void Wave(){//wave clasters finding
            int tag = -1, contRO_lsumm, add_elem;
            std::vector<int> positions;
            for(unsigned int q = 0; q < work_points_.size(); q ++){
                positions.clear();
                contRO_lsumm = 0; //If we check all std::strings, where first elems was included into cluster
                add_elem = 0; //If we add new element to cluster
                for(unsigned int i = q; i < work_points_.size(); ++ i){
                    if(this -> label_[i] == -1){
                        tag ++;
                        positions.push_back(i);
                        add_elem ++;
                        label_[i] = tag;
                        for(unsigned int j = 0; j < work_points_.size(); ++ j){
                            if(this -> binary_table_[i][j] == 1){
                                positions.push_back(j);
                                add_elem ++;
                            }
                        }
                        i = work_points_.size();//break;
                    }
                }
                while(contRO_lsumm != add_elem){
                    for(unsigned int i = 0; i < positions.size(); ++ i){
                        contRO_lsumm ++;
                        for(unsigned int j = 0; j < work_points_.size(); ++ j){
                            if(this -> binary_table_[positions[i]][j] == 1){
                                if(!(IsInVector(positions, j))){
                                    add_elem ++;
                                    positions.push_back(j);
                                }
                            }
                        }
                    }
                }
                for(const auto &strnumber : positions){
                    label_[strnumber] = tag;
                }
            }
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                work_points_[i].SetLabel(this -> label_[i]);
                if(this -> label_[i] < 0){
                    Error(6);
                }
            }
        }
        void CalculateRO_(){//calculate all RO_ between all points
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                for(unsigned int j = 0; j < work_points_.size(); ++ j){
                    if(i == j){
                        RO_[i][j] = -1;
                    } else if(j == 0){
                        RO_[i][j] = -1;
                    } else {
                        RO_[i][j] = sqrt(Sqr(this -> work_points_[i].GetX() - work_points_[j].GetX()) +
                                                Sqr(this -> work_points_[i].GetY() - work_points_[j].GetY()));
                    }
                }
            }
        }
        void FillLengths(){//making std::vector of distances for GNU gistogram
            std::vector<int> index;
            std::fstream tree, script, text, launchlog;
            std::vector<std::vector<double>> TempRO_ = RO_;
            double minRO_, mr;
            int minId = -1, mi = -1, saveIND = -1;
            launchlog.open("launch.log");
            launchlog >> launchNumber;
            launchNumber ++;
            launchlog << launchNumber;
            launchlog.close();
            std::string filename = "Tree" + std::to_string(this -> launchNumber) + ".txt", scriptname = "plotTree" + std::to_string(this -> launchNumber) + ".plt",
                        sstr = "TempGisto" + std::to_string(this -> launchNumber) + ".txt", ssstr = "plotTempGisto" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, std::fstream::trunc);
            script << "plot '" << filename << "' using 1:2 with lines lc rgb \"black\" lw 1 notitle";
            script.close();
            tree.open(filename, std::fstream::trunc);
            index.clear();
            index.push_back(0);
            while(this -> allLengths.size() < work_points_.size() - 1){
                mr = -1;
                for(const auto &id : index){//for each id fRO_m index array
                    minRO_ = TempRO_[id][0];
                    for(unsigned int i = 0; i < work_points_.size(); ++ i){
                        if((minRO_ < 0) && (TempRO_[id][i] > 0)){
                            minRO_ = TempRO_[id][i];
                            minId = i;
                        } else if ((TempRO_[id][i] < minRO_) && (TempRO_[id][i] > 0)){
                            minRO_ = TempRO_[id][i];
                            minId = i;
                        }
                    }
                    if((mr < 0) || (mr > minRO_)){
                        mr = minRO_;
                        mi = minId;
                        saveIND = id;
                    }
                }
                tree << work_points_[saveIND].GetX() << " " <<this -> work_points_[saveIND].GetY() << std::endl
                    << work_points_[mi].GetX() << " " <<this -> work_points_[mi].GetY() << std::endl << std::endl;
                allLengths.push_back(mr);
                index.push_back(mi);
                for(unsigned int i = 0; i < work_points_.size(); ++ i){
                    TempRO_[i][mi] = -1;
                }
            }
            tree.close();
            script.open(ssstr, std::fstream::trunc);
            script << "width=0.01" << std::endl << "bin(x, s) = s*int(x/s) + width/2"
                   << std::endl << "set boxwidth width" << std::endl << "plot '" << sstr << "' u (bin($1,width)):(1.0) \\"
                   << std::endl << "s f w boxes fs solid 0.5 title 'PoRO_g Gisto'" << std::endl;
            script.close();
            text.open(sstr, std::fstream::trunc); //writing values to the file
            for(unsigned int i = 0; i < (this -> work_points_.size() - 1); ++ i){
                text << allLengths[i] << std::endl;
            }
            text.close();
        }
        std::vector<Point> FindClusters(double threshold){//current points
            std::string scn1 = "TempGisto" + std::to_string(this -> launchNumber);
            BeforeStart();
            CalculateRO_();
            FillLengths();
            launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            BeforeStart();
            CreateRO_();
            GenBinary(threshold);
            Wave();
            return work_points_;
        }
};

class FindByHierarchy{
    private:
        int launchNumber;
        std::vector<Point> work_points_;
    public:
        void BeforeStart(){//saving
            std::fstream f;
            f.open("launch.log");
            f >> launchNumber;
            if(this -> launchNumber < 1){
                f.close();
                f.open("launch.log", std::fstream::trunc);
                f << 1;
            }
            f.close();
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                work_points_[i].SetLabel(0);
            }
        }
        FindByHierarchy(const std::vector<Point>& koords){//constructor
            work_points_ = koords;
        }
        void Hierarchy(int userK){
            double x1, y1, x2, y2, minRO_, midx, midy;
            std::fstream tree, script, launchlog;
            launchlog.open("launch.log");
            launchlog >> launchNumber;
            launchlog.close();
            std::string filename = "HieTree" + std::to_string(this -> launchNumber) + ".txt";
            std::string scriptname = "plotHieTree" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, std::fstream::trunc);
            script << "plot '" << filename << "' using 1:2 with lines lc rgb \"black\" lw 2 notitle";
            script.close();
            launchlog.open("launch.log", std::fstream::trunc);
            launchlog.close();
            tree.open(filename, std::fstream::trunc);
            std::vector<Point> newPoints = work_points_;
            int p1, p2, numpoints = work_points_.size(), numclusters = userK;
            std::vector<std::vector<double>> RO_, rRO_;
            while(numpoints != numclusters){
                RO_.resize(newPoints.size());
                minRO_ = -1;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    RO_[i].resize(newPoints.size());
                    x1 = newPoints[i].GetX();
                    y1 = newPoints[i].GetY();
                    for(unsigned int j = 0; j < newPoints.size(); ++ j){
                            if(i != j){
                            x2 = newPoints[j].GetX();
                            y2 = newPoints[j].GetY();
                            RO_[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                            if(((minRO_ < 0) || (minRO_ > RO_[i][j])) && (RO_[i][j] > 0)){
                                minRO_ = RO_[i][j];
                                p1 = i;
                                p2 = j;
                            }
                        }
                    }
                }
                x1 = newPoints[p1].GetX();
                x2 = newPoints[p2].GetX();
                y1 = newPoints[p1].GetY();
                y2 = newPoints[p2].GetY();
                tree << x1 << " " << y1 << std::endl << x2 << " " << y2 << std::endl << std::endl;
                midx = (x1 + x2) / 2.0;
                midy = (y1 + y2) / 2.0;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    RO_[i].erase(RO_[i].begin() + std::max(p1, p2));
                    RO_[i].erase(RO_[i].begin() + std::min(p1, p2));
                }
                RO_.erase(RO_.begin() + std::max(p1, p2));
                RO_.erase(RO_.begin() + std::min(p1, p2));
                newPoints.erase(newPoints.begin() + std::max(p1, p2));
                newPoints.erase(newPoints.begin() + std::min(p1, p2));
                newPoints.push_back(Point(midx, midy, 0));
                numpoints --;
            }
            tree.close();
            rRO_.resize(newPoints.size());
            for(unsigned int i = 0; i < newPoints.size(); ++ i){
                rRO_[i].resize(this -> work_points_.size());
                x1 = newPoints[i].GetX();
                y1 = newPoints[i].GetY();
                for(unsigned int j = 0; j < work_points_.size(); ++ j){
                    x2 = work_points_[j].GetX();
                    y2 = work_points_[j].GetY();
                    rRO_[i][j] = Correct(sqrt(Sqr(x1 - x2) + Sqr(y1 - y2)));
                }
            }
            for(unsigned int j = 0; j < work_points_.size(); ++ j){
                minRO_ = rRO_[0][j];
                p1 = 0;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    if(rRO_[i][j] < minRO_){
                        minRO_ = rRO_[i][j];
                        p1 = i;
                    }
                }
                work_points_[j].SetLabel(p1);
            }
        }
        std::vector<Point> FindClusters(int user_number_of_clusters){//current points
            BeforeStart();
            Hierarchy(user_number_of_clusters);
            return work_points_;
        }
};

class FindByForel{
    private:
        int launchNumber;
        std::vector<Point> work_points_;
    public:
        FindByForel(const std::vector<Point>& koords){//constructor
            work_points_ = koords;
        }
        void Forel(double threshold){
            std::fstream circle, script, text;//, launchlog;
            std::vector<double> RO_;//matrix of distances
            std::vector<int> binary, id;
            bool change;
            int startID, counterPoints, ptr, clustID = 0;
            double midx, midy;
            launchNumber = ReadLaunch();
            RewriteLaunch(this -> launchNumber + 1);
            std::string filename = "CircleF" + std::to_string(this -> launchNumber) + ".txt";
            std::string scriptname = "plotCircleF" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, std::fstream::trunc);
            script << "plot '" << filename << "' using 1:2:3 with circles lc rgb \"black\" lw 2 notitle";
            script.close();
            circle.open(filename, std::fstream::trunc);
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                work_points_[i].SetLabel(-i);
            }
            std::vector<Point> tempPoints = work_points_;
            std::unique_ptr<Point> center (new Point(tempPoints[0]));
            while(tempPoints.size() > 0){
                startID = rand () % tempPoints.size();
                RO_.resize(tempPoints.size());
                binary.resize(tempPoints.size());
                change = true;
                center.reset(new Point(tempPoints[startID]));
                while(change){
                    midx = 0;
                    midy = 0;
                    counterPoints = 0;
                    for(unsigned int i = 0; i < tempPoints.size(); ++ i){
                        RO_[i] = sqrt(Sqr(center -> GetX() - tempPoints[i].GetX()) + Sqr(center -> GetY() - tempPoints[i].GetY()));
                        if(RO_[i] < threshold){
                            binary[i] = 1;
                            counterPoints ++;
                            midx += tempPoints[i].GetX();
                            midy += tempPoints[i].GetY();
                        } else {
                            binary[i] = 0;
                        }
                    }
                    midx /= counterPoints;
                    midy /= counterPoints;
                    if((ABS(center -> GetX() - midx) < eps) && (ABS(center -> GetY() - midy) < eps)){//when center stops (in that iteration)
                        change = false;
                        circle << midx << " " << midy << " " << threshold << std::endl;
                        id.resize(counterPoints);
                        ptr = 0;
                        for(unsigned int i = 0; ((i < tempPoints.size()) && (ptr < counterPoints)); ++ i){
                            if(binary[i] == 1){
                                id[ptr] = tempPoints[i].GetLabel();
                                ptr ++;
                            }
                        }
                        ptr = counterPoints - 1;
                        for(unsigned int i = work_points_.size() - 1; /*((i >= 0) && (*/ptr >= 0/*))*/; i --){//making label_s
                            if(this -> work_points_[i].GetLabel() == id[ptr]){
                                work_points_[i].SetLabel(clustID);
                                ptr --;
                            }
                        }
                        ptr = counterPoints - 1;
                        for(unsigned int i = tempPoints.size() - 1; /*((i >= 0) && (*/ptr >= 0/*))*/; i --){//delete found cluster
                            if(tempPoints[i].GetLabel() == id[ptr]){
                                tempPoints.erase(tempPoints.begin() + i);
                                ptr --;
                            }
                        }
                        clustID ++;
                    } else {
                        center.reset(new Point(midx, midy, -5));
                    }
                }
            }
            circle.close();
        }
        std::vector<Point> FindClusters(double threshold){//current points
            Forel(threshold);
            std::string filename = "Forel" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            launchNumber ++;
            RewriteLaunch(this ->  launchNumber);
            return work_points_;
        }
};

class FindByDbscan{
    private:
        int launchNumber;
        std::vector<int> counters, marks;//label_s
        std::vector<Point> work_points_, field_koord_;
        std::vector<std::vector<double>> RO_;//distances between points
    public:
        void BeforeStart(){//saving
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                work_points_[i].SetLabel(0);
            }
            field_koord_ = work_points_;
            RO_.resize(this -> work_points_.size());
            counters.clear();
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                RO_[i].resize(this -> work_points_.size());
            }
        }
        FindByDbscan(const std::vector<Point>& koords){//constructor
            work_points_ = koords;
        }
        void CalcRO_(){//matrix of distances
            double x1, y1, x2, y2;
            RO_.resize(this -> work_points_.size());
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                RO_[i].resize(this -> work_points_.size());
                x1 = work_points_[i].GetX();
                y1 = work_points_[i].GetY();
                for(unsigned int j = 0; j < work_points_.size(); ++ j){
                    x2 = work_points_[j].GetX();
                    y2 = work_points_[j].GetY();
                    RO_[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                }
            }
        }
        std::vector<int> CountNearPoints(std::vector<std::vector<double>> RO_Matrix, double radius){
            std::vector<int> ctr;
            ctr.resize(this -> work_points_.size());
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                ctr[i] = 0;
                for(unsigned int j = 0; j < work_points_.size(); ++ j){
                    if(RO_Matrix[i][j] < (radius + eps)){
                        ctr[i] ++;
                    }
                }
            }
            return ctr;
        }
        std::vector<int> CreateMarks(std::vector<std::vector<double>> RO_, double radius, std::vector<int> cnt, int EntNum){//label_s for clusters
            std::vector<int> marker;
            marker.resize(this -> work_points_.size());
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                if(cnt[i] > EntNum){
                    marker[i] = 1;
                } else {
                    marker[i] = -1;
                }
            }
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                if(marker[i] == -1){
                    for(unsigned int j = 0; j < work_points_.size(); ++ j){
                        if((RO_[i][j] < (radius + eps)) && (marker[j] == 1)){
                            marker[i] = 0;
                        }
                    }
                }
            }
            return marker;
        }
        void DelPoints(double radius, int EntryNum){//delete rubbish
            field_koord_.clear();
            for(unsigned int i = 0; i < work_points_.size(); ++ i){
                field_koord_.push_back(Point(this -> work_points_[i]));
            }
            CalcRO_();
            counters = CountNearPoints(this -> RO_, radius);
            marks = CreateMarks(this -> RO_, radius, counters, EntryNum);
            for(int i = field_koord_.size() - 1; i >= 0; i --){
                if(this -> marks[i] == -1){
                    field_koord_.erase(this -> field_koord_.begin() + i);
                    counters.erase(this -> counters.begin() + i);
                    for(unsigned int j = 0; j < field_koord_.size(); ++ j){
                        RO_[j].erase(this -> RO_[j].begin() + i);
                    }
                    RO_.erase(this -> RO_.begin() + i);
                    marks.erase(this -> marks.begin() + i);
                    work_points_[i].SetLabel(-9999);
                }
            }
        }
        void DBSCANAlghorithm(double radius, int numpointsincircle){//wave on points without rubbish
            BeforeStart();
            DelPoints(radius, numpointsincircle);
            work_points_.clear();
            for(unsigned int i = 0; i < field_koord_.size(); ++ i){
                work_points_.push_back(Point(this -> field_koord_[i]));
            }
            std::unique_ptr <FindByWave> Waw(new FindByWave(this -> work_points_));
            work_points_ = Waw -> FindClusters(radius + eps);
        }
        std::vector<Point> FindClusters(double radius, int NumNearPoints){//current points
            std::string filename = "DBSCAN" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            DBSCANAlghorithm(radius, NumNearPoints);
            launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            return work_points_;
        }
};

class FindByEM{
    private:
        int k, launchNumber;//launch No.
        std::vector<Point> work_points_;
        std::vector<Point> CentersA, CentersB;
    public:
        void BeforeWork(){//default constructor
            work_points_.clear();
            CentersA.clear();
            CentersB.clear();
            k = 0;
        }
        FindByEM(const std::vector<Point>& koords){//constructor
            BeforeWork();
            work_points_ = koords;
        }
        double CentRO_(){
            double rz = 0;
            for(unsigned int i = 0; i < CentersB.size(); ++ i){
                rz += sqrt(Sqr(CentersB[i].GetX() - CentersA[i].GetX()) + Sqr(CentersB[i].GetY() - CentersA[i].GetY()));
            }
            if(rz < eps){
                rz =  -999;
            }
            return rz;
        }
        void GNUMultOut(int ctr, std::vector<std::vector<int>> clrs){
            std::fstream mult;
            mult.open("mult.txt", std::fstream::ate);
            for(int i = 0; i < k; ++ i){
                mult << work_points_[i].GetX() << " "
                     << work_points_[i].GetY() << " " << clrs[this -> work_points_[i].GetLabel()][0] << " "
                     << clrs[this -> work_points_[i].GetLabel()][1] << " " << clrs[this -> work_points_[i].GetLabel()][2] << std::endl;
            }
            mult << std::endl << std::endl;
            mult.close();
        }
        void EM(int userK){//EM with fixed K
            int itercount = 0;
            k = work_points_.size();
            std::vector<std::vector<int>> Mycolors;
            Mycolors.resize(userK);
            srand(time(NULL));
            for(int i = 0; i < userK; ++ i){
                Mycolors[i].resize(3);
                Mycolors[i][0] = rand() % 256;
                Mycolors[i][1] = rand() % 256;
                Mycolors[i][2] = rand() % 256;
            }
            std::vector<std::vector<std::vector<double>>> Sigma, SigmaM;
            std::vector<std::vector<double>> d, p;
           double mx;
           int index, sUk = userK;
            d.resize(this -> k);
            p.resize(this -> k);
            for(int i = 0; i < k; ++ i){
                d[i].resize(userK);
                p[i].resize(userK);
            }
            CentersB.clear();
            for(int i = 0; i < userK; ++ i){
                CentersB.push_back(Point(this -> work_points_[i]));
            }
            Sigma.resize(userK);
            SigmaM.resize(userK);
            for(int i = 0; i < userK; ++ i){
                Sigma[i].resize(2);
                SigmaM[i].resize(2);
                Sigma[i][0].resize(2);
                Sigma[i][1].resize(2);
                SigmaM[i][0].resize(2);
                SigmaM[i][1].resize(2);
                Sigma[i][0][0] = 1;
                Sigma[i][0][1] = 0;
                Sigma[i][1][0] = 0;
                Sigma[i][1][1] = 1;
                SigmaM[i][0][0] = 1;
                SigmaM[i][0][1] = 0;
                SigmaM[i][1][0] = 0;
                SigmaM[i][1][1] = 1;
            }
            double summ, sumpx, sump, sumpy;
            std::fstream mult;
            mult.open("mult.txt", std::fstream::trunc);
            mult.close();
            mult.open("mult.plt", std::fstream::trunc);
            mult << "set terminal gif animate delay 100" << std::endl;
            mult << "set output 'MULT.gif'" << std::endl;
            mult << "stats 'mult.txt' nooutput" << std::endl << std::endl;
            mult << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)" << std::endl << std::endl;
            mult << "do for[i=1:int(STATS_blocks)]{" << std::endl;
                mult << "plot 'mult.txt' index(i-1) using 1:2:(rgb($3,$4,$5)) with points lc rgb variable" << std::endl;
            mult << "}" << std::endl;
            mult.close();
            bool esc = true;
                while(esc){
                itercount ++;
                    for(int i = 0; i < k; ++ i){
                        summ = 0;
                        for(int j = 0; j < userK; ++ j){
                            d[i][j] = sqrt(SigmaM[j][0][0] * Sqr(this -> work_points_[i].GetX() - CentersB[j].GetX()) +
                                   (SigmaM[j][1][0] + SigmaM[j][0][1]) * (this -> work_points_[i].GetX() -
                                    CentersB[j].GetX()) * (this -> work_points_[i].GetY() - CentersB[j].GetY()) +
                                    SigmaM[j][1][1] * Sqr(this -> work_points_[i].GetY() - CentersB[j].GetY()));
                            summ += d[i][j];
                        }
                        for(int j = 0; j < userK; ++ j){
                            p[i][j] = 1.0 / (sUk - 1) - d[i][j] / (sUk - 1) / summ;
                        }
                    }
                    CentersA.clear();
                    for(int j = 0; j < userK; ++ j){
                        sumpx = 0;
                        sumpy = 0;
                        sump = 0;
                        for(int i = 0; i < k; ++ i){
                            sumpx += work_points_[i].GetX() * p[i][j];
                            sumpy += work_points_[i].GetY() * p[i][j];
                            sump += p[i][j];
                        }
                        CentersA.push_back(Point(sumpx / sump, sumpy / sump));
                    }
                    if(this -> CentRO_() > 0){
                        CentersB.clear();
                        for(int i = 0; i < userK; ++ i){
                            CentersB.push_back(CentersA[i]);
                        }
                        for(int j = 0; j < userK; ++ j){
                            Sigma[j][0][0] = 0;
                            Sigma[j][0][1] = 0;
                            Sigma[j][1][0] = 0;
                            Sigma[j][1][1] = 0;
                            for(int i = 0; i < k; ++ i){
                                Sigma[j][0][0] += p[i][j] * Sqr(this -> work_points_[i].GetX() - CentersB[j].GetX());
                                Sigma[j][0][1] += p[i][j] * (this -> work_points_[i].GetX() - CentersB[j].GetX()) *
                                                            (this -> work_points_[i].GetY() - CentersB[j].GetY());
                                Sigma[j][1][1] += p[i][j] * Sqr(this -> work_points_[i].GetY() - CentersB[j].GetY());
                            }
                            Sigma[j][1][0] = Sigma[j][0][1];

                            SigmaM[j][0][0] = Sigma[j][1][1] / ABS(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][0][1] = - Sigma[j][0][1] / ABS(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][1][0] = - Sigma[j][1][0] / ABS(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][1][1] = Sigma[j][0][0] / ABS(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                        }
                        if((itercount % 2) == 1){
                            for(int i = 0; i < k; i++){
                                mx = p[i][0];
                                index = 0;
                                for(int j = 0; j < userK; ++ j){
                                    if(p[i][j] > mx){
                                        mx = p[i][j];
                                        index = j;
                                    }
                                }
                                work_points_[i].SetLabel(index);
                            }
                            GNUMultOut(itercount, Mycolors);
                        }
                    } else {
                        esc = false;
                        for(int i = 0; i < k; i++){
                            mx = p[i][0];
                            index = 0;
                            for(int j = 0; j < userK; ++ j){
                                if(p[i][j] > mx){
                                    mx = p[i][j];
                                    index = j;
                                }
                            }
                            work_points_[i].SetLabel(index);
                        }
                        GNUMultOut(itercount, Mycolors);
                    }
                    CentersB.clear();
                    for(int i = 0; i < userK; ++ i){
                        CentersB.push_back(Point(this -> CentersA[i]));
                    }
                }
        }
        std::vector<Point> FindClusters(int user_number_of_clusters){//current points
            EM(user_number_of_clusters);
            launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            return work_points_;
        }
};
