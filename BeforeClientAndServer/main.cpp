#include <clocale>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

const double PI = 3.1415926535897932384626433832795;
const double eps = 0.0000001;

std::fstream Mylog;

void RewriteLaunch (int val) { //launcher
    std::fstream launcher;
    launcher.open("launch.log", std::fstream::trunc);
    launcher << val;
    launcher.close();
}

int ReadLaunch () { //launch reading
    std::fstream launcher;
    int val;
    launcher.open("launch.log");
    launcher >> val;
    launcher.close();
    return val;
}

double DegreesToRadians (double Grad) { //convertion of degrees to radians
    return (Grad * PI / 180);
}

void Error (int type) { //error warning
    switch (type) {
        case 1:
            std::cout << std::endl << "Size of std::vector less than 1" << std::endl;
            Mylog << "Error! Size of std::vector less than 1" << std::endl;
            break;
        case 2:
            std::cout << std::endl << "Size of group less than 1" << std::endl;
            Mylog << "Error! Size of group less than 1" << std::endl;
            break;
        case 3:
            std::cout << std::endl << "Couldn't open file" << std::endl;
            Mylog << "Error! Couldn't open file" << std::endl;
            break;
        case 4:
            std::cout << std::endl << "Worksize mustn't be less than 1" << std::endl;
            Mylog << "Error! Worksize mustn't be less than 1" << std::endl;
            break;
        case 5:
            std::cout << std::endl << "Incorrect input. Please try again" << std::endl;
            Mylog << "Error! Incorrect input. Please try again" << std::endl;
            break;
        case 6:
            std::cout << std::endl << "Incorrect value of point's label_" << std::endl;
            Mylog << "Error! Incorrect value of point's label_" << std::endl;
            break;
        case 7:
            std::cout << std::endl << "Incorrect file data" << std::endl;
            Mylog << "Error! Incorrect file data" << std::endl;
            break;
        default:
            std::cout << std::endl << "Unnamed Error" << std::endl;
            Mylog << "Error! Unnamed Error" << std::endl;
            break;
    }
}

bool IsInVector (const std::vector<int> &vect, int elem) { //checking of element in std::vector
    bool rez = false;
    for (size_t i = 0; i < vect.size(); ++ i){
        if (vect[i] == elem){
            rez = true;
            return true;
        }
    }
    return rez;
}

double Sqr (double val) {
    return val * val;
}

double Correct (double val) { //change too little value to 0
    if (val < eps) {
        return 0.0;
    }
    return val;
}

class Point { //points
private:
    double x_, y_; //coordinates
    int label_; //label

public:
    Point () = default; //Default constructor
    Point (double new_x, double new_y) : x_(new_x), y_(new_y) {} // Initialize fields with coordinates
    Point (double new_x, double new_y, int new_label) // Initialize fields with new arguments
        : x_(new_x), y_(new_y), label_(new_label) {} 
    Point (const Point &other) 
        : x_(other.x_), y_(other.y_), label_(other.label_) {} // Copy constructor

    double GetX () { //get X
        return x_;
    }
    double GetY () { //get Y
        return y_;
    }
    int GetLabel () { //get label_
        return label_;
    }

    void SetX (double new_x) { //set X coordinate
        x_ = new_x;
    }
    void SetY (double new_y) { //set Y coordinate
        y_ = new_y;
    }
    void SetLabel (int new_label) { //set label
        label_ = new_label;
    }

    void AddX (double x) { //move X coordinate
        x_ += x;
    }
    void AddY (double y) { //move Y coordinate
        y_ += y;
    }
    Point& operator = (const Point& other) { //remake coordinates
        x_ = other.x_;
        y_ = other.y_;
        label_ = other.label_;
        return *this;
    }
};

std::vector<Point> arrr (const std::string& filename) { //read Point std::vector from file
    std::fstream f;
    double x, y;
    int l;
    f.open(filename);
    std::vector<Point> rez;
    while (!f.eof()) {
        f >> x >> y >> l;
        rez.push_back(Point(x, y, l));
    }
    return rez;
}

class Control {
private:
    std::vector<double> vec1_, vec2_; //std::vectors of coordinates
    std::vector<Point> group_; //std::vector of groups
    size_t num_point1_, num_point2_, num_point_group_;
public:
    Control () = default; //default constructor

    void MakeLabel (int lb) { //group label_
        for (size_t i = 0; i < num_point_group_; ++ i) {
            group_[i].SetLabel(lb);
        }
    }
    std::vector<double> CreateNorm (size_t number_of_points, double min_val, double max_val) { //norm gisto
        if (number_of_points < 1) {
            Error(1);
            return {}; // Empty, default-constructed std::vector
        }

        if (min_val < max_val) {
            std::swap(min_val, max_val);
        }

        std::vector<double> arr (number_of_points); // Create a std::vector of predefined size
        double sum, temp;
        int mmn = static_cast<int>(min_val) * 10,
            mmx = (static_cast<int>(max_val) + 1) * 10,
            range;
        range = mmx - mmn + 1;

        for (int i = 0; i < static_cast<int>(number_of_points); ++ i) {
            sum = 0;
            for (int j = 0; j < 1000; ++ j) {
                temp = ((rand () % range) + mmn) / 10.0;
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
                arr[i] = sum;
            }
        }
        return arr;
    }
    void GenRnd (size_t number_of_points, double min_val, double max_val) { //ravn gisto
        if (number_of_points < 1) {
            Error(1);
            return void();
        }

        vec1_.clear();
        num_point1_ = number_of_points;
        for (size_t i = 0; i < num_point1_; ++ i) {
            vec1_.push_back(((rand() % (static_cast<int>(100 * max_val)
                                        - static_cast<int>(100 * min_val)))
                            + static_cast<int>(100 * min_val)) / 100.0);
        }
    }
    void GenNorm (size_t number_of_points, double min_val, double max_val) { //norm generation
        if (number_of_points < 1) {
            Error(1);
            return void();
        }

        vec2_.clear();
        num_point2_ = number_of_points;
        if (max_val < min_val) {
            std::swap(min_val, max_val);
        }
        std::vector<double> arr = CreateNorm(num_point2_, min_val, max_val);
        for (size_t i = 0; i < num_point2_; ++ i) {
            vec2_.push_back(arr[i]);
        }
    }
    void GenGroup (size_t number_of_points, double min_val_x, double max_val_x, double min_val_y, double max_val_y, int label) { //group creating
        if (number_of_points < 1) {
            Error(2);
            return void();
        }

        num_point_group_ = number_of_points;
        group_.clear();
        std::vector<double> arrx = CreateNorm(num_point_group_, min_val_x, max_val_x);
        std::vector<double> arry = CreateNorm(num_point_group_, min_val_y, max_val_y);
        for(size_t i = 0; i < num_point_group_; ++ i){
            group_.push_back(Point(arrx[i], arry[i], label));
        }
    }
    void GenGroup (const std::vector<Point>& values) { //group creating copy
        num_point_group_ = values.size();
        group_.clear();
        for (size_t i = 0; i < num_point_group_; ++ i) {
            group_.push_back(Point(values[i]));
        }
    }
    void FileRavn () { //printing ravn std::vector in file
        std::fstream text, script;
        script.open("plotravn.plt", std::fstream::trunc); //writing script to the file
        script << "width=1" << std::endl << "bin(x, s) = s*int(x/s) + width/2"
               << std::endl << "set boxwidth width" << std::endl << "plot 'ravn.txt' u (bin($1,width)):(1.0) \\"
               << std::endl << "s f w boxes fs solid 0.5 title 'Ravn Gisto'" << std::endl;
        script.close();
        text.open("ravn.txt", std::fstream::trunc); //writing values to the file
        for (size_t i = 0; i < num_point1_; ++ i) {
            text << vec1_[i] << std::endl;
        }
        text.close();
    }
    void FileNorm () { //printing norm std::vector in file
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
    void FileGroup () { //printing group in file
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
    std::vector<double> RetRavn () { //ravn generation cpy
        std::vector<double> arr(num_point1_);
        for (size_t i = 0; i < num_point1_; ++ i) {
            arr[i] = vec1_[i];
        }
        return arr;
    }
    std::vector<double> RetNorm () { //norm generation cpy
        std::vector<double> arr(num_point2_);
        for (size_t i = 0; i < num_point2_; ++ i) {
            arr[i] = vec2_[i];
        }
        return arr;
    }
    std::vector<Point> RetGroup () { //group cpy
        std::vector<Point> arr; // No default ctor for Point, so no resizes
        arr.reserve(num_point_group_); // At least we'll hint about our dataset size
        for (size_t i = 0; i < num_point_group_; ++ i) {
            arr.push_back(group_[i]);
        }
        return arr;
    }
    void TurnNULL (double phi) { //rotarion group relative (0;0)
        double newx, newy;
        std::vector<Point> newvect;
        newvect.reserve(num_point_group_);
        for (size_t i = 0; i < num_point_group_; ++ i) {
            newx = group_[i].GetX();
            newy = group_[i].GetY();
            newvect.push_back(Point(
                newx * cos(phi) - newy * sin(phi),
                newx * sin(phi) + newy * cos(phi)));
        }
        group_.clear();
        GenGroup(newvect);
    }
    void turnCenter (double phi) { //rotarion group relative to group center
        double midx = 0, midy = 0, newx, newy;
        std::vector<Point> newvect;
        newvect.reserve(num_point_group_);
        for (size_t i = 0; i < num_point_group_; ++ i) {
            midx += group_[i].GetX();
            midy += group_[i].GetY();
        }
        midx /= num_point_group_;
        midy /= num_point_group_;
        for (size_t i = 0; i < num_point_group_; ++ i) {
            newx = group_[i].GetX() - midx;
            newy = group_[i].GetY() - midy;
            newvect.push_back(Point(
                (newx * cos(phi) - newy * sin(phi)) + midx,
                (newx * sin(phi) + newy * cos(phi)) + midy));
        }
        group_.clear();
        GenGroup(newvect);
    }
    void MoveX (double deltax) { //X axis moving
        std::vector<Point> newvect = RetGroup();
        for (size_t i = 0; i < num_point_group_; ++ i) {
            newvect[i].SetX(newvect[i].GetX() + deltax);
        }
        group_.clear();
        GenGroup(newvect);
    }
    void MoveY (double deltay) { //Y axis moving
        std::vector<Point> newvect = RetGroup();
        for (size_t i = 0; i < num_point_group_; ++ i) {
            newvect[i].SetY(newvect[i].GetY() + deltay);
        }
        group_.clear();
        GenGroup(newvect);
    }
    int RetGroupSize () {
        return num_point_group_;
    }
};

class Group{
    private:
        std::unique_ptr<Control> group;//array of groups
        int Grouplabel_;
    public:
        Group() = default; // Default constructor, required for std::vector::resize

        Group(int groupsize, double minx, double miny, double maxx, double maxy, int GL){//constructor
            this -> group.reset(new Control());
            this -> Grouplabel_ = GL;
            this -> group -> GenGroup(groupsize, minx, maxx, miny, maxy, GL);
        }
        Group(std::unique_ptr<Control> newgroup){//constructor //it means that label_ included into Point' std::vector
            std::vector<Point> arr = newgroup -> RetGroup();
            this -> Grouplabel_ = arr[0].GetLabel();
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
        }
        Group(std::unique_ptr <Control> newgroup, int LB){//constructor
            std::vector<Point> arr = newgroup -> RetGroup();
            this -> Grouplabel_ = LB;
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
            this -> group -> MakeLabel(LB);
        }
        int RetGrouplabel_(){
            return this -> Grouplabel_;
        }
        std::unique_ptr<Control> Ret_Field_Group() const {//returning of a group
            std::unique_ptr <Control> rgroup(new Control());
            std::vector<Point> arr = this -> group -> RetGroup();
            rgroup -> GenGroup(arr);
            return rgroup;
        }
        Group(const Group &newgroup) {//cpy constructor
            std::unique_ptr<Control> GR = newgroup.Ret_Field_Group();
            std::vector<Point> arr = GR -> RetGroup();
            this -> group.reset(new Control());
            this -> Grouplabel_ = arr[0].GetLabel();
            this -> group -> GenGroup(arr);
        }
        Group(std::unique_ptr <Group> &newgroup) {//cpy constructor
            std::unique_ptr <Control> GR = newgroup->Ret_Field_Group();
            std::vector<Point> arr = GR -> RetGroup();
            this -> group.reset(new Control());
            this -> Grouplabel_ = arr[0].GetLabel();
            this -> group -> GenGroup(arr);
        }
        void Regrupp(const std::unique_ptr <Control> &newgroup){ //we modified this object by new Control object
            std::vector<Point> arr = newgroup -> RetGroup();
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
        }
        void Regrupp(const std::unique_ptr <Control> &newgroup, int LB){ //we modified this object by new Control object
            std::vector<Point> arr = newgroup -> RetGroup();
            this -> group.reset(new Control());
            this -> group -> GenGroup(arr);
            this -> Grouplabel_ = LB;
            this -> group -> MakeLabel(LB);
        }
        void ReMakeLabel(){//remake label_ of a group
            this -> group -> MakeLabel(this -> Grouplabel_);
        }
        void ReMakeLabel(int LB){//remake label_ of a group
            this -> Grouplabel_ = LB;
            this -> group -> MakeLabel(this -> Grouplabel_);
        }
};

class FindByWave{
    private:
        std::vector<int> label_;//label_s
        std::vector<Point> workPoints;//current points
        std::vector<std::vector<int>> binary_table;//binary matrix
        std::vector<std::vector<double>> RO;//distances between points
        double maxRO;
    public:
        void BeforeWork(){//default constructor
            this -> workPoints.resize(0);
            this -> maxRO = 0;
            this -> label_.resize(0);
            this -> binary_table.resize(0);
            this -> RO.resize(0);
        }
        void BeforeStart(){//saving
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(0);
            }
            this -> maxRO = 0;
            this -> label_.resize(this -> workPoints.size());
            this -> binary_table.resize(this -> workPoints.size());
            this -> RO.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> binary_table[i].resize(this -> workPoints.size());
                this -> RO[i].resize(this -> workPoints.size());
            }
        }
        FindByWave(const std::vector<Point>& koords){//constructor
            this -> BeforeWork();
            this -> workPoints = koords;
            this -> BeforeStart();
        }
        void CreateRo(){//creating of matrix of distances
            double x1, y1, x2, y2;
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> label_[i] = -1;
                x1 = this -> workPoints[i].GetX();
                y1 = this -> workPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    this -> RO[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                    if(this -> RO[i][j] > this -> maxRO){
                        this -> maxRO = this -> RO[i][j];
                    }
                }
            }
        }
        void GenBinary(double threshold){//creating of binary matrix
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    if((this -> RO[i][j] > 0) && (this -> RO[i][j] < (threshold + eps))){
                        this -> binary_table[i][j] = 1;
                    } else {
                        this -> binary_table[i][j] = 0;
                    }
                }
            }
        }
        void Wave(){//wave clasters finding
            int tag = -1, controlsumm, add_elem;
            std::vector<int> positions;
            for(unsigned int q = 0; q < this -> workPoints.size(); q ++){
                positions.clear();
                controlsumm = 0; //If we check all std::strings, where first elems was included into cluster
                add_elem = 0; //If we add new element to cluster
                for(unsigned int i = q; i < this -> workPoints.size(); ++ i){
                    if(this -> label_[i] == -1){
                        tag ++;
                        positions.push_back(i);
                        add_elem ++;
                        this -> label_[i] = tag;
                        for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                            if(this -> binary_table[i][j] == 1){
                                positions.push_back(j);
                                add_elem ++;
                            }
                        }
                        i = this -> workPoints.size();//break;
                    }
                }
                while(controlsumm != add_elem){
                    for(unsigned int i = 0; i < positions.size(); ++ i){
                        controlsumm ++;
                        for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                            if(this -> binary_table[positions[i]][j] == 1){
                                if(!(IsInVector(positions, j))){
                                    add_elem ++;
                                    positions.push_back(j);
                                }
                            }
                        }
                    }
                }
                for(const auto &strnumber : positions){
                    this -> label_[strnumber] = tag;
                }
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(this -> label_[i]);
                if(this -> label_[i] < 0){
                    Error(6);
                }
            }
        }
        std::vector<Point> FindClusters(double threshold){
            this -> BeforeStart();
            this -> CreateRo();
            this -> GenBinary(threshold);
            this -> Wave();
            return this -> workPoints;
        }
};

class FindByKM{
    private:
        int k, fsize, optimalK;
        double bestWeight;
        std::vector<Point> workPoints, fieldkoord, centerKoords, startPoints;
    public:
        void BeforeWork(){//default constructor
            this -> workPoints.resize(0);
            this -> fieldkoord.resize(0);
            this -> centerKoords.resize(0);
            this -> startPoints.resize(0);
            this -> k = 0;
            this -> optimalK = 0;
            this -> bestWeight = 0;
            this -> fsize = 0;
        }
        void BeforeStart(){
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(0);
            }
            this -> fieldkoord.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> fieldkoord.push_back(Point(this -> workPoints[i]));
            }
        }
        FindByKM(const std::vector<Point>& koords){//constructor
            this -> BeforeWork();
            this -> workPoints = koords;
            this -> BeforeStart();
        }
        void KMeans(int userK){//kmeans with fixed K
            int minid, qStart, qEnd, ptr = 0;
            double minVal;
            bool change, esc = false;
            std::vector<double> MuFunction;
            if(userK < 0){
                qStart = 1;
                qEnd = this -> fsize;
            } else {
                qStart = userK;
                qEnd = userK + 1;
            }
            MuFunction.resize(this -> fsize - 1);
            for(int q = qStart; ((q < qEnd) && (!esc)); q ++){
                MuFunction[ptr] = 0;
                this -> k = q;
                std::vector<std::vector<double>> RO_selected, RO_centers;
                change = true;
                this -> startPoints.clear();
                this -> startPoints.reserve(this -> k);
                RO_selected.resize(this -> k);
                this -> centerKoords.resize(this -> k);
                std::vector<Point> midPoint;
                std::vector<int> numberOfPointsInEachCluster(this -> k);
                for(int i = 0; i < this -> k; ++ i){ //create first matrix of RO between all points by pairs
                    this -> startPoints.push_back(Point(this -> fieldkoord[i]));
                    midPoint.push_back(Point(0, 0, -5));//-5 for mid points
                    numberOfPointsInEachCluster[i] = 0;
                    RO_selected[i].resize(this -> fsize);
                    for(int j = 0; j < this -> fsize; j++){
                        RO_selected[i][j] = Correct(sqrt(Sqr(this -> startPoints[i].GetX() - this -> fieldkoord[j].GetX())
                                                         + Sqr(this -> startPoints[i].GetY() - this -> fieldkoord[j].GetY())));
                    }
                }
                for(int j = 0; j < this -> fsize; ++ j){//mark each point (which center was the nearest)
                    minVal = RO_selected[0][j];
                    minid = 0;
                    for(int i = 1; i < this -> k; ++ i){//find minimal ro from j point to k ros of each cluster
                        if(RO_selected[i][j] < minVal){
                            minVal = RO_selected[i][j];
                            minid = i;
                        }
                    }
                    this -> fieldkoord[j].SetLabel(minid);
                    numberOfPointsInEachCluster[minid] ++;
                }
                for(int i = 0; i < this -> fsize; ++ i){ //caclulate mid points
                    int id = this -> fieldkoord[i].GetLabel();
                    midPoint[id].AddX(Correct(this -> fieldkoord[i].GetX() / numberOfPointsInEachCluster[id]));
                    midPoint[id].AddY(Correct(this -> fieldkoord[i].GetY() / numberOfPointsInEachCluster[id]));

                }
                for(int i = 0; i < this -> k; ++ i){ //make new koord centers
                    this -> centerKoords[i].SetX(midPoint[i].GetX());
                    this -> centerKoords[i].SetY(midPoint[i].GetY());
                }
                while(change){
                    RO_centers.resize(this -> k);
                    change = false;
                    midPoint.clear();
                    numberOfPointsInEachCluster.clear();
                    midPoint.resize(this -> k);
                    numberOfPointsInEachCluster.resize(this -> k);
                    for(int i = 0; i < this -> k; ++ i){ //calc all rors from centers
                        RO_centers[i].resize(this -> fsize);
                        midPoint[i] = Point(0,0,-5);
                        numberOfPointsInEachCluster[i] = 0;
                        for(int j = 0; j < this -> fsize; ++ j){ //find new ros from new centers to each point
                            RO_centers[i][j] = Correct(sqrt(Sqr(this -> centerKoords[i].GetX() - this -> fieldkoord[j].GetX())
                                                            + Sqr(this -> centerKoords[i].GetY() - this -> fieldkoord[j].GetY())));
                        }
                    }
                    for(int j = 0; j < this -> fsize; ++ j){
                        minVal = RO_centers[0][j];
                        minid = 0;
                        for(int i = 1; i < this -> k; ++ i){
                            if(RO_centers[i][j] < minVal){
                                minVal = RO_centers[i][j];
                                minid = i;
                            }
                        }
                        this -> fieldkoord[j].SetLabel(minid);
                        numberOfPointsInEachCluster[minid] ++;
                    }
                    for(int i = 0; i < this -> fsize; ++ i){ //caclulate mid points
                        int id = this -> fieldkoord[i].GetLabel();
                        midPoint[id].AddX(Correct(this -> fieldkoord[i].GetX() / numberOfPointsInEachCluster[id]));
                        midPoint[id].AddY(Correct(this -> fieldkoord[i].GetY() / numberOfPointsInEachCluster[id]));
                    }
                    for(int i = 0; i < this -> k; ++ i){//check definitions between new and old center koords
                        double xOld, yOld, xNew, yNew;
                        xOld = this -> centerKoords[i].GetX();
                        yOld = this -> centerKoords[i].GetY();
                        xNew = midPoint[i].GetX();
                        yNew = midPoint[i].GetY();
                        if((abs(xOld - xNew) > eps) || (abs(yOld - yNew) > eps)){
                            change = true;
                            i = this -> k;
                        }
                    }
                    if(change){
                        this -> centerKoords.clear();
                        this -> centerKoords.resize(this -> k);
                        for(int i = 0; i < this -> k; ++ i){
                            this -> centerKoords[i].SetX(midPoint[i].GetX());
                            this -> centerKoords[i].SetY(midPoint[i].GetY());
                        }
                    }
                }
                for(int i = 0; i < this -> k; ++ i){// here we find optimal k
                    for(int j = 0; j < this -> fsize; ++ j){
                        if(this -> fieldkoord[j].GetLabel() == i){
                            for(int p = j + 1; p < this -> fsize; p ++){
                                if(this -> fieldkoord[p].GetLabel() == i){
                                   MuFunction[ptr] += Correct(sqrt(Sqr(this -> fieldkoord[j].GetX() - this -> fieldkoord[p].GetX()) +
                                                                    Sqr(this -> fieldkoord[j].GetY() - this -> fieldkoord[p].GetY())));
                                }
                            }
                        }
                    }
                }
                for (int i = 0; i < q; ++ i){//here we find optimal k
                    for(int j = i + 1; j < q; ++ j){
                        MuFunction[ptr] += Correct(sqrt(Sqr(this -> centerKoords[i].GetX() - this -> centerKoords[j].GetX()) +
                                                      Sqr(this -> centerKoords[i].GetY() - this -> centerKoords[j].GetY())));
                    }
                }
                if((q != 1) && ((qEnd - qStart) > 1)){
                    if(MuFunction[ptr - 1] < MuFunction[ptr]){
                        esc = true;
                        ptr -= 2;
                    }
                }
                ptr ++;
            }
            this -> optimalK = ptr + 1;
            this -> bestWeight = MuFunction[ptr];
            this -> workPoints.resize(0);
            for(unsigned int i = 0; i < this -> fieldkoord.size(); ++ i){
                this -> workPoints.push_back(this -> fieldkoord[i]);
            }
        }
        std::vector<Point> FindClusters(int userNumberOfClusters){
            this -> BeforeStart();
            this -> KMeans(userNumberOfClusters);
            std::cout << "3" << std::endl;
            if(userNumberOfClusters == -1){
                std::cout << std::endl << "Optimal number of clusters is " << this -> optimalK << std::endl;
                this -> BeforeStart();
                this -> KMeans(this -> optimalK);
            }
            return this -> workPoints;
        }
};

class FindBySPTR{
    private:
        std::vector<Point> workPoints;
        int launchNumber;
        std::vector<std::vector<int>> binary_table;//binary matrix
        std::vector<double> allLengths;
        std::vector<std::vector<double>> RO;//distances between points
        double maxRO;
        std::vector<int> label_;//label_s
    public:
        void BeforeWork(){//default constructor
            std::fstream f;
            this -> launchNumber = -1;
            f.open("launch.log");
            f >> this -> launchNumber;
            if(this -> launchNumber <= 1){
                f.close();
                f.open("launch.log", std::fstream::trunc);
                f << 1;
            }
            f.close();
            this -> workPoints.resize(0);
            this -> allLengths.resize(0);
            this -> maxRO = 0;
            this -> label_.resize(0);
            this -> binary_table.resize(0);
            this -> RO.resize(0);
        }
        void BeforeStart(){//saving
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(0);
            }
            this -> allLengths.resize(0);
            this -> maxRO = 0;
            this -> label_.resize(this -> workPoints.size());
            this -> binary_table.resize(this -> workPoints.size());
            this -> RO.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> binary_table[i].resize(this -> workPoints.size());
                this -> RO[i].resize(this -> workPoints.size());
            }
        }
        FindBySPTR(const std::vector<Point>& koords){//constructor
            this -> BeforeWork();
            this -> workPoints = koords;
            this -> BeforeStart();
        }
        void CreateRo(){//creating of matrix of distances
            double x1, y1, x2, y2;
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> label_[i] = -1;
                x1 = this -> workPoints[i].GetX();
                y1 = this -> workPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    this -> RO[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                    if(this -> RO[i][j] > this -> maxRO){
                        this -> maxRO = this -> RO[i][j];
                    }
                }
            }
        }
        void GenBinary(double threshold){//creating of binary matrix
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    if((this -> RO[i][j] > 0) && (this -> RO[i][j] < (threshold + eps))){
                        this -> binary_table[i][j] = 1;
                    } else {
                        this -> binary_table[i][j] = 0;
                    }
                }
            }
        }
        void Wave(){//wave clasters finding
            int tag = -1, controlsumm, add_elem;
            std::vector<int> positions;
            for(unsigned int q = 0; q < this -> workPoints.size(); q ++){
                positions.clear();
                controlsumm = 0; //If we check all std::strings, where first elems was included into cluster
                add_elem = 0; //If we add new element to cluster
                for(unsigned int i = q; i < this -> workPoints.size(); ++ i){
                    if(this -> label_[i] == -1){
                        tag ++;
                        positions.push_back(i);
                        add_elem ++;
                        this -> label_[i] = tag;
                        for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                            if(this -> binary_table[i][j] == 1){
                                positions.push_back(j);
                                add_elem ++;
                            }
                        }
                        i = this -> workPoints.size();//break;
                    }
                }
                while(controlsumm != add_elem){
                    for(unsigned int i = 0; i < positions.size(); ++ i){
                        controlsumm ++;
                        for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                            if(this -> binary_table[positions[i]][j] == 1){
                                if(!(IsInVector(positions, j))){
                                    add_elem ++;
                                    positions.push_back(j);
                                }
                            }
                        }
                    }
                }
                for(const auto &strnumber : positions){
                    this -> label_[strnumber] = tag;
                }
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(this -> label_[i]);
                if(this -> label_[i] < 0){
                    Error(6);
                }
            }
        }
        void CalculateRo(){//calculate all Ro between all points
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    if(i == j){
                        this -> RO[i][j] = -1;
                    } else if(j == 0){
                        this -> RO[i][j] = -1;
                    } else {
                        this -> RO[i][j] = sqrt(Sqr(this -> workPoints[i].GetX() - this -> workPoints[j].GetX()) +
                                                Sqr(this -> workPoints[i].GetY() - this -> workPoints[j].GetY()));
                    }
                }
            }
        }
        void FillLengths(){//making std::vector of distances for GNU gistogram
            std::vector<int> index;
            std::fstream tree, script, text, launchlog;
            std::vector<std::vector<double>> TempRo = this -> RO;
            double minRo, mr;
            int minId = -1, mi = -1, saveIND = -1;
            launchlog.open("launch.log");
            launchlog >> this -> launchNumber;
            this -> launchNumber ++;
            launchlog << this -> launchNumber;
            launchlog.close();
            std::string filename = "Tree" + std::to_string(this -> launchNumber) + ".txt", scriptname = "plotTree" + std::to_string(this -> launchNumber) + ".plt",
                        sstr = "TempGisto" + std::to_string(this -> launchNumber) + ".txt", ssstr = "plotTempGisto" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, std::fstream::trunc);
            script << "plot '" << filename << "' using 1:2 with lines lc rgb \"black\" lw 1 notitle";
            script.close();
            tree.open(filename, std::fstream::trunc);
            index.resize(0);
            index.push_back(0);
            while(this -> allLengths.size() < this -> workPoints.size() - 1){
                mr = -1;
                for(const auto &id : index){//for each id from index array
                    minRo = TempRo[id][0];
                    for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                        if((minRo < 0) && (TempRo[id][i] > 0)){
                            minRo = TempRo[id][i];
                            minId = i;
                        } else if ((TempRo[id][i] < minRo) && (TempRo[id][i] > 0)){
                            minRo = TempRo[id][i];
                            minId = i;
                        }
                    }
                    if((mr < 0) || (mr > minRo)){
                        mr = minRo;
                        mi = minId;
                        saveIND = id;
                    }
                }
                tree << this -> workPoints[saveIND].GetX() << " " <<this -> workPoints[saveIND].GetY() << std::endl
                    << this -> workPoints[mi].GetX() << " " <<this -> workPoints[mi].GetY() << std::endl << std::endl;
                this -> allLengths.push_back(mr);
                index.push_back(mi);
                for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                    TempRo[i][mi] = -1;
                }
            }
            tree.close();
            script.open(ssstr, std::fstream::trunc);
            script << "width=0.01" << std::endl << "bin(x, s) = s*int(x/s) + width/2"
                   << std::endl << "set boxwidth width" << std::endl << "plot '" << sstr << "' u (bin($1,width)):(1.0) \\"
                   << std::endl << "s f w boxes fs solid 0.5 title 'Porog Gisto'" << std::endl;
            script.close();
            text.open(sstr, std::fstream::trunc); //writing values to the file
            for(unsigned int i = 0; i < (this -> workPoints.size() - 1); ++ i){
                text << this -> allLengths[i] << std::endl;
            }
            text.close();
        }
        std::vector<Point> FindClusters(double threshold){//current points
            std::string scn1 = "TempGisto" + std::to_string(this -> launchNumber);
            this -> BeforeStart();
            this -> CalculateRo();
            this -> FillLengths();
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            this -> BeforeStart();
            this -> CreateRo();
            this -> GenBinary(threshold);
            this -> Wave();
            return this -> workPoints;
        }
};

class FindByHierarchy{
    private:
        int launchNumber;
        std::vector<Point> workPoints;
    public:
        void BeforeStart(){//saving
            std::fstream f;
            f.open("launch.log");
            f >> this -> launchNumber;
            if(this -> launchNumber < 1){
                f.close();
                f.open("launch.log", std::fstream::trunc);
                f << 1;
            }
            f.close();
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(0);
            }
        }
        FindByHierarchy(const std::vector<Point>& koords){//constructor
            this -> workPoints = koords;
        }
        void Hierarchy(int userK){
            double x1, y1, x2, y2, minro, midx, midy;
            std::fstream tree, script, launchlog;
            launchlog.open("launch.log");
            launchlog >> this -> launchNumber;
            launchlog.close();
            std::string filename = "HieTree" + std::to_string(this -> launchNumber) + ".txt";
            std::string scriptname = "plotHieTree" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, std::fstream::trunc);
            script << "plot '" << filename << "' using 1:2 with lines lc rgb \"black\" lw 2 notitle";
            script.close();
            launchlog.open("launch.log", std::fstream::trunc);
            launchlog.close();
            tree.open(filename, std::fstream::trunc);
            std::vector<Point> newPoints = this -> workPoints;
            int p1, p2, numpoints = this -> workPoints.size(), numclusters = userK;
            std::vector<std::vector<double>> ro, rro;
            while(numpoints != numclusters){
                ro.resize(newPoints.size());
                minro = -1;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    ro[i].resize(newPoints.size());
                    x1 = newPoints[i].GetX();
                    y1 = newPoints[i].GetY();
                    for(unsigned int j = 0; j < newPoints.size(); ++ j){
                            if(i != j){
                            x2 = newPoints[j].GetX();
                            y2 = newPoints[j].GetY();
                            ro[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                            if(((minro < 0) || (minro > ro[i][j])) && (ro[i][j] > 0)){
                                minro = ro[i][j];
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
                    ro[i].erase(ro[i].begin() + std::max(p1, p2));
                    ro[i].erase(ro[i].begin() + std::min(p1, p2));
                }
                ro.erase(ro.begin() + std::max(p1, p2));
                ro.erase(ro.begin() + std::min(p1, p2));
                newPoints.erase(newPoints.begin() + std::max(p1, p2));
                newPoints.erase(newPoints.begin() + std::min(p1, p2));
                newPoints.push_back(Point(midx, midy, 0));
                numpoints --;
            }
            tree.close();
            rro.resize(newPoints.size());
            for(unsigned int i = 0; i < newPoints.size(); ++ i){
                rro[i].resize(this -> workPoints.size());
                x1 = newPoints[i].GetX();
                y1 = newPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    rro[i][j] = Correct(sqrt(Sqr(x1 - x2) + Sqr(y1 - y2)));
                }
            }
            for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                minro = rro[0][j];
                p1 = 0;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    if(rro[i][j] < minro){
                        minro = rro[i][j];
                        p1 = i;
                    }
                }
                this -> workPoints[j].SetLabel(p1);
            }
        }
        std::vector<Point> FindClusters(int userNumberOfClusters){//current points
            this -> BeforeStart();
            this -> Hierarchy(userNumberOfClusters);
            return this -> workPoints;
        }
};

class FindByForel{
    private:
        int launchNumber;
        std::vector<Point> workPoints;
    public:
        FindByForel(const std::vector<Point>& koords){//constructor
            this -> workPoints = koords;
        }
        void Forel(double threshold){
            std::fstream circle, script, text;//, launchlog;
            std::vector<double> ro;//matrix of distances
            std::vector<int> binary, id;
            bool change;
            int startID, counterPoints, ptr, clustID = 0;
            double midx, midy;
            this -> launchNumber = ReadLaunch();
            RewriteLaunch(this -> launchNumber + 1);
            std::string filename = "CircleF" + std::to_string(this -> launchNumber) + ".txt";
            std::string scriptname = "plotCircleF" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, std::fstream::trunc);
            script << "plot '" << filename << "' using 1:2:3 with circles lc rgb \"black\" lw 2 notitle";
            script.close();
            circle.open(filename, std::fstream::trunc);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(-i);
            }
            std::vector<Point> tempPoints = this -> workPoints;
            std::unique_ptr <Point> center (new Point(tempPoints[0]));
            while(tempPoints.size() > 0){
                startID = rand () % tempPoints.size();
                ro.resize(tempPoints.size());
                binary.resize(tempPoints.size());
                change = true;
                center.reset(new Point(tempPoints[startID]));
                while(change){
                    midx = 0;
                    midy = 0;
                    counterPoints = 0;
                    for(unsigned int i = 0; i < tempPoints.size(); ++ i){
                        ro[i] = sqrt(Sqr(center -> GetX() - tempPoints[i].GetX()) + Sqr(center -> GetY() - tempPoints[i].GetY()));
                        if(ro[i] < threshold){
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
                    if((abs(center -> GetX() - midx) < eps) && (abs(center -> GetY() - midy) < eps)){//when center stops (in that iteration)
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
                        for(unsigned int i = this -> workPoints.size() - 1; /*((i >= 0) && (*/ptr >= 0/*))*/; i --){//making label_s
                            if(this -> workPoints[i].GetLabel() == id[ptr]){
                                this -> workPoints[i].SetLabel(clustID);
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
            this -> Forel(threshold);
            std::string filename = "Forel" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            this -> launchNumber ++;
            RewriteLaunch(this ->  launchNumber);
            return this -> workPoints;
        }
};

class FindByDbscan{
    private:
        int launchNumber;
        std::vector<int> counters, marks;//label_s
        std::vector<Point> workPoints, fieldkoord;
        std::vector<std::vector<double>> RO;//distances between points
    public:
        void BeforeStart(){//saving
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(0);
            }
            this -> fieldkoord = this -> workPoints;
            this -> RO.resize(this -> workPoints.size());
            this -> counters.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> RO[i].resize(this -> workPoints.size());
            }
        }
        FindByDbscan(const std::vector<Point>& koords){//constructor
            this -> workPoints = koords;
        }
        void CalcRo(){//matrix of distances
            double x1, y1, x2, y2;
            this -> RO.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> RO[i].resize(this -> workPoints.size());
                x1 = this -> workPoints[i].GetX();
                y1 = this -> workPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    this -> RO[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                }
            }
        }
        std::vector<int> CountNearPoints(std::vector<std::vector<double>> RoMatrix, double radius){
            std::vector<int> ctr;
            ctr.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                ctr[i] = 0;
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    if(RoMatrix[i][j] < (radius + eps)){
                        ctr[i] ++;
                    }
                }
            }
            return ctr;
        }
        std::vector<int> CreateMarks(std::vector<std::vector<double>> ro, double radius, std::vector<int> cnt, int EntNum){//label_s for clusters
            std::vector<int> marker;
            marker.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                if(cnt[i] > EntNum){
                    marker[i] = 1;
                } else {
                    marker[i] = -1;
                }
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                if(marker[i] == -1){
                    for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                        if((ro[i][j] < (radius + eps)) && (marker[j] == 1)){
                            marker[i] = 0;
                        }
                    }
                }
            }
            return marker;
        }
        void DelPoints(double radius, int EntryNum){//delete rubbish
            this -> fieldkoord.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> fieldkoord.push_back(Point(this -> workPoints[i]));
            }
            this -> CalcRo();
            this -> counters = CountNearPoints(this -> RO, radius);
            this -> marks = CreateMarks(this -> RO, radius, this -> counters, EntryNum);
            for(int i = this -> fieldkoord.size() - 1; i >= 0; i --){
                if(this -> marks[i] == -1){
                    this -> fieldkoord.erase(this -> fieldkoord.begin() + i);
                    this -> counters.erase(this -> counters.begin() + i);
                    for(unsigned int j = 0; j < this -> fieldkoord.size(); ++ j){
                        this -> RO[j].erase(this -> RO[j].begin() + i);
                    }
                    this -> RO.erase(this -> RO.begin() + i);
                    this -> marks.erase(this -> marks.begin() + i);
                    this -> workPoints[i].SetLabel(-9999);
                }
            }
        }
        void DBSCANAlghorithm(double radius, int numpointsincircle){//wave on points without rubbish
            this -> BeforeStart();
            this -> DelPoints(radius, numpointsincircle);
            this -> workPoints.resize(0);
            for(unsigned int i = 0; i < this -> fieldkoord.size(); ++ i){
                this -> workPoints.push_back(Point(this -> fieldkoord[i]));
            }
            std::unique_ptr <FindByWave> Waw(new FindByWave(this -> workPoints));
            this -> workPoints = Waw -> FindClusters(radius + eps);
        }
        std::vector<Point> FindClusters(double radius, int NumNearPoints){//current points
            std::string filename = "DBSCAN" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            this -> DBSCANAlghorithm(radius, NumNearPoints);
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            return this -> workPoints;
        }
};

class FindByEM{
    private:
        int k, launchNumber;//launch No.
        std::vector<Point> workPoints;
        std::vector<Point> CentersA, CentersB;
    public:
        void BeforeWork(){//default constructor
            this -> workPoints.resize(0);
            this -> CentersA.resize(0);
            this -> CentersB.resize(0);
            this -> k = 0;
        }
        FindByEM(const std::vector<Point>& koords){//constructor
            this -> BeforeWork();
            this -> workPoints = koords;
        }
        double CentRo(){
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
            for(int i = 0; i < this -> k; ++ i){
                mult << this -> workPoints[i].GetX() << " "
                     << this -> workPoints[i].GetY() << " " << clrs[this -> workPoints[i].GetLabel()][0] << " "
                     << clrs[this -> workPoints[i].GetLabel()][1] << " " << clrs[this -> workPoints[i].GetLabel()][2] << std::endl;
            }
            mult << std::endl << std::endl;
            mult.close();
        }
        void EM(int userK){//EM with fixed K
            int itercount = 0;
            this -> k = this -> workPoints.size();
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
            for(int i = 0; i < this -> k; ++ i){
                d[i].resize(userK);
                p[i].resize(userK);
            }
            this -> CentersB.resize(0);
            for(int i = 0; i < userK; ++ i){
                this -> CentersB.push_back(Point(this -> workPoints[i]));
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
                    for(int i = 0; i < this -> k; ++ i){
                        summ = 0;
                        for(int j = 0; j < userK; ++ j){
                            d[i][j] = sqrt(SigmaM[j][0][0] * Sqr(this -> workPoints[i].GetX() - this -> CentersB[j].GetX()) +
                                   (SigmaM[j][1][0] + SigmaM[j][0][1]) * (this -> workPoints[i].GetX() -
                                    this -> CentersB[j].GetX()) * (this -> workPoints[i].GetY() - this -> CentersB[j].GetY()) +
                                    SigmaM[j][1][1] * Sqr(this -> workPoints[i].GetY() - this -> CentersB[j].GetY()));
                            summ += d[i][j];
                        }
                        for(int j = 0; j < userK; ++ j){
                            p[i][j] = 1.0 / (sUk - 1) - d[i][j] / (sUk - 1) / summ;
                        }
                    }
                    this -> CentersA.resize(0);
                    for(int j = 0; j < userK; ++ j){
                        sumpx = 0;
                        sumpy = 0;
                        sump = 0;
                        for(int i = 0; i < this -> k; ++ i){
                            sumpx += this -> workPoints[i].GetX() * p[i][j];
                            sumpy += this -> workPoints[i].GetY() * p[i][j];
                            sump += p[i][j];
                        }
                        this -> CentersA.push_back(Point(sumpx / sump, sumpy / sump));
                    }
                    if(this -> CentRo() > 0){
                        this -> CentersB.resize(0);
                        for(int i = 0; i < userK; ++ i){
                            this -> CentersB.push_back(CentersA[i]);
                        }
                        for(int j = 0; j < userK; ++ j){
                            Sigma[j][0][0] = 0;
                            Sigma[j][0][1] = 0;
                            Sigma[j][1][0] = 0;
                            Sigma[j][1][1] = 0;
                            for(int i = 0; i < this -> k; ++ i){
                                Sigma[j][0][0] += p[i][j] * Sqr(this -> workPoints[i].GetX() - this -> CentersB[j].GetX());
                                Sigma[j][0][1] += p[i][j] * (this -> workPoints[i].GetX() - this -> CentersB[j].GetX()) *
                                                            (this -> workPoints[i].GetY() - this -> CentersB[j].GetY());
                                Sigma[j][1][1] += p[i][j] * Sqr(this -> workPoints[i].GetY() - this -> CentersB[j].GetY());
                            }
                            Sigma[j][1][0] = Sigma[j][0][1];

                            SigmaM[j][0][0] = Sigma[j][1][1] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][0][1] = - Sigma[j][0][1] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][1][0] = - Sigma[j][1][0] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][1][1] = Sigma[j][0][0] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                        }
                        if((itercount % 2) == 1){
                            for(int i = 0; i < this -> k; i++){
                                mx = p[i][0];
                                index = 0;
                                for(int j = 0; j < userK; ++ j){
                                    if(p[i][j] > mx){
                                        mx = p[i][j];
                                        index = j;
                                    }
                                }
                                this -> workPoints[i].SetLabel(index);
                            }
                            this -> GNUMultOut(itercount, Mycolors);
                        }
                    } else {
                        esc = false;
                        for(int i = 0; i < this -> k; i++){
                            mx = p[i][0];
                            index = 0;
                            for(int j = 0; j < userK; ++ j){
                                if(p[i][j] > mx){
                                    mx = p[i][j];
                                    index = j;
                                }
                            }
                            this -> workPoints[i].SetLabel(index);
                        }
                        this -> GNUMultOut(itercount, Mycolors);
                    }
                    this -> CentersB.resize(0);
                    for(int i = 0; i < userK; ++ i){
                        this -> CentersB.push_back(Point(this -> CentersA[i]));
                    }
                }
        }
        std::vector<Point> FindClusters(int userNumberOfClusters){//current points
            this -> EM(userNumberOfClusters);
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            return this -> workPoints;
        }
};

class Field{
    private:
        std::vector<Group> field;//std::vector field
        int fieldsize/*number of groups*/, pointer, totalsize,/*number of points*/k, fsize, optimalK, launchNumber;//launch No.
        std::vector<int> arrsz, label_, counters, marks;//label_s
        std::vector<Point> allkoord, centerKoords,/*centers of groups*/ startPoints,/*start coordinate std::vector*/ save;/*for save*/
        std::vector<Point> workPoints,/*std::vector of current coordinates*/ fieldkoord,/*field coordinates*/ forClusters, saveData;/*for save*/
        std::vector<Point> CentersA, CentersB; //need foe em alg
        std::vector<double> allLengths;
        std::vector<std::vector<int>> binary_table;//binary matrix
        std::vector<std::vector<double>> RO;//distances between points
        double bestWeight, maxRO;
    public:
        void BeforeWork(){//default constructor
            this -> workPoints.resize(0);
            this -> fieldkoord.resize(0);
            this -> allLengths.resize(0); //we will set them after
            this -> centerKoords.resize(0);
            this -> startPoints.resize(0);
            this -> save.resize(0);
            this -> CentersA.resize(0);
            this -> CentersB.resize(0);
            this -> k = 0;
            this -> maxRO = 0;
            this -> optimalK = 0;
            this -> bestWeight = 0;
            this -> fsize = 0;
            this -> label_.resize(0);
            this -> binary_table.resize(0);
            this -> RO.resize(0);
            this -> counters.resize(0);
            this -> forClusters.resize(0);
            this -> saveData.resize(0);
        }
        void BeforeStart(){//saving
            if(this -> saveData.size() > 0){
                this -> workPoints.resize(0);
                for(unsigned int i = 0; i < this -> saveData.size(); ++ i){
                    this -> workPoints.push_back(Point(this -> saveData[i]));
                }
                this -> saveData.resize(0);
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(0);
            }
            this -> fieldkoord = this -> workPoints;
            this -> allLengths.resize(0); //we will set them after
            this -> centerKoords.resize(0);
            this -> startPoints.resize(0);
            this -> save.resize(0);
            this -> CentersA.resize(0);
            this -> CentersB.resize(0);
            this -> k = 0;
            this -> maxRO = 0;
            this -> optimalK = 0;
            this -> bestWeight = 0;
            this -> fsize = this -> workPoints.size();
            this -> label_.resize(this -> workPoints.size());
            this -> binary_table.resize(this -> workPoints.size());
            this -> RO.resize(this -> workPoints.size());
            this -> counters.resize(0);
            this -> forClusters.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> binary_table[i].resize(this -> workPoints.size());
                this -> RO[i].resize(this -> workPoints.size());
            }
        }
        Field(int numofgroups){//constructor
            this -> BeforeWork();
            if(numofgroups > 0){
                this -> fieldsize = numofgroups;
            } else {
                this -> fieldsize = 0;
            }
            this -> field.resize(0); //Reserve memory for our object
            this -> pointer = 0;
            this -> totalsize = 0;
            this -> CentersA.resize(0);
            this -> CentersB.resize(0);
            this -> arrsz.resize(this -> fieldsize); // int default-constructs to 0
        }
        Field(const std::vector<Point>& koords){//constructor
            this -> BeforeWork();
            this -> totalsize = koords.size();
            this -> allkoord = koords;
            this -> field.resize(this -> totalsize); //Reserve memory for our object
            this -> arrsz.resize(this -> totalsize); // int default-constructs to 0
            this -> fieldsize = 0;
            this -> pointer = 0;
            this -> CentersA.resize(0);
            this -> CentersB.resize(0);
            this -> workPoints = koords;
        }
        Field(const std::string& filename){//constructor from file
            this -> BeforeWork();
            this -> allkoord = arrr(filename);
            this -> totalsize = allkoord.size();
            this -> field.resize(this -> totalsize); //Reserve memory for our object
            this -> arrsz.resize(this -> totalsize); // int default-constructs to 0
            this -> fieldsize = 0;
            this -> CentersA.resize(0);
            this -> CentersB.resize(0);
            this -> pointer = 0;
        }
        int retTOTSIZE(){//size of field
            return this -> totalsize;
        }
        std::vector<Point> retALLKOORD(){//current points
            if(this -> workPoints.size() < 1){
                return this -> allkoord;
            }
            return this -> workPoints;
        }
        void AddGroup(const std::unique_ptr <Group> &newgroup){//adding of a group to the field
            std::unique_ptr <Control> tc = newgroup -> Ret_Field_Group();
            int gsz = tc -> RetGroupSize();
            this -> arrsz.push_back(gsz);
            std::unique_ptr <Group> TempGroup(new Group());
            TempGroup -> Regrupp(tc, this -> pointer);
            this -> field.push_back(Group(TempGroup));
            this -> pointer ++;
        }
        void MakeAllKoord(){//list of field coordinates
            int totsize = 0;
            this -> allkoord.clear();
            std::vector<Point> tmp;
            for(int i = 0; i < this -> fieldsize; ++ i){
                tmp.clear();
                std::unique_ptr <Control> tempcontrol = this -> field[i].Ret_Field_Group();
                tmp = tempcontrol -> RetGroup();
                totsize += tmp.size();
                for(unsigned int j = 0; j < tmp.size(); ++ j){
                    this -> allkoord.push_back(Point(tmp[j]));
                }
            }
            this -> totalsize = totsize;
            this -> workPoints.resize(0);
            for(unsigned int i = 0; i < this -> allkoord.size(); ++ i){
                this -> workPoints.push_back(Point(this -> allkoord[i]));
            }
        }
        void ToTxt(std::string answ){//printing of the field to file
            std::fstream text;
            bool addcolor = false;
            int r, g, b;
            if((answ != "yes") && (answ != "no")){
                Error(5);
            }
            if(answ == "yes"){
                addcolor = true;
                std::fstream script;
                script.open("plotField.plt", std::fstream::trunc);
                script << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b) " << std::endl;
                script << "plot 'Field.txt' using 1:2:(rgb($3,$4,$5)) with points lc rgb variable";
                script.close();
            }
            text.open("Field.txt", std::fstream::trunc);
            r = rand() % 256;
            g = rand() % 256;
            b = rand() % 256;
            for(int i = 0; i < this -> totalsize; ++ i){
                text << this -> allkoord[i].GetX() << " " << this -> allkoord[i].GetY();
                if(addcolor){
                    if((i == 0) || (this -> allkoord[i].GetLabel() != this -> allkoord[i - 1].GetLabel())){
                        r = rand() % 256;
                        g = rand() % 256;
                        b = rand() % 256;
                    }
                    text << " " << r << " " << g << " " << b;
                }
                text << std::endl;
            }
            text.close();
            text.open("FieldBackup.txt", std::fstream::trunc);
            text << "FIELD" << std::endl << "FIELD" << std::endl;
            for(int i = 0; i < this -> totalsize; ++ i){
                text << this -> allkoord[i].GetX() << " " << this -> allkoord[i].GetY() << " "
                    << this -> allkoord[i].GetLabel();
                if(i != this -> totalsize - 1){
                    text << std::endl;
                }
            }
            text.close();
        }
/*
        //WAVE
        void CreateRo(){//creating of matrix of distances
            double x1, y1, x2, y2;
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> label_[i] = -1;
                x1 = this -> workPoints[i].GetX();
                y1 = this -> workPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    this -> RO[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                    if(this -> RO[i][j] > this -> maxRO){
                        this -> maxRO = this -> RO[i][j];
                    }
                }
            }
        }
        void GenBinary(double threshold){//creating of binary matrix
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    if((this -> RO[i][j] > 0) && (this -> RO[i][j] < (threshold + eps))){
                        this -> binary_table[i][j] = 1;
                    } else {
                        this -> binary_table[i][j] = 0;
                    }
                }
            }
        }
        void Wave(){//wave clasters finding
            int tag = -1, controlsumm, add_elem;
            std::vector<int> positions;
            for(unsigned int q = 0; q < this -> workPoints.size(); q ++){
                positions.clear();
                controlsumm = 0; //If we check all std::strings, where first elems was included into cluster
                add_elem = 0; //If we add new element to cluster
                for(unsigned int i = q; i < this -> workPoints.size(); ++ i){
                    if(this -> label_[i] == -1){
                        tag ++;
                        positions.push_back(i);
                        add_elem ++;
                        this -> label_[i] = tag;
                        for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                            if(this -> binary_table[i][j] == 1){
                                positions.push_back(j);
                                add_elem ++;
                            }
                        }
                        i = this -> workPoints.size();//break;
                    }
                }
                while(controlsumm != add_elem){
                    for(unsigned int i = 0; i < positions.size(); ++ i){
                        controlsumm ++;
                        for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                            if(this -> binary_table[positions[i]][j] == 1){
                                if(!(Instd::vector(positions, j))){
                                    add_elem ++;
                                    positions.push_back(j);
                                }
                            }
                        }
                    }
                }
                for(const auto &strnumber : positions){
                    this -> label_[strnumber] = tag;
                }
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(this -> label_[i]);
                if(this -> label_[i] < 0){
                    Error(6);
                }
            }
        }
*/
        //TXT
        void ToTxtCluster(std::string name){//printing clasters in file. Files set by user
            std::vector<std::vector<int>> rgb;
            std::vector<Point> bpoint;
            int kkk = 0, l, counterPoints;
            std::fstream script, myFile, bup;
            std::string filename = "", scriptname = "plot", backupname;
            double midx, midy, leftx, rightx, topy, bottomy;
            bool start;
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                if(this -> workPoints[i].GetLabel() > kkk){
                    kkk = workPoints[i].GetLabel();
                }
            }
            rgb.resize(kkk + 1);
            for(int i = 0; i < (kkk + 1); ++ i){
                rgb[i].resize(3);
                rgb[i][0] = rand() % 256;
                rgb[i][1] = rand() % 256;
                rgb[i][2] = rand() % 256;
            }
            for(unsigned i = 0; i < name.size(); ++ i){
                if(name[i] == '.'){
                    i = name.size();
                } else {
                    filename += name[i];
                }
            }
            scriptname += filename + ".plt";
            backupname = filename + "Backup.txt";
            filename = name;
            script.open(scriptname, std::fstream::trunc);
            script << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b) " << std::endl
                << "plot '" << filename << "' using 1:2:(rgb($3,$4,$5)) with points lc rgb variable";
            script.close();
            myFile.open(filename, std::fstream::trunc);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                l = this -> workPoints[i].GetLabel();
                if(l >= 0){
                    myFile << this -> workPoints[i].GetX() << " " << this -> workPoints[i].GetY() << " "
                        << rgb[l][0] << " " << rgb[l][1] << " " << rgb[l][2] << std::endl;
                }
            }
            myFile.close();
            bup.open(backupname, std::fstream::trunc);
            if(name != "Field.txt"){
                bup << "CLUSTER" << std::endl << "Total clusters: " << kkk + 1 << std::endl
                    << "Total points: " << this -> workPoints.size() << std::endl;
                for(int i = 0; i <= kkk; ++ i){
                    bpoint.resize(0);
                    counterPoints = 0;
                    midx = 0;
                    midy = 0;
                    start = true;
                    for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                        if(this -> workPoints[j].GetLabel() == i){
                            bpoint.push_back(Point(this -> workPoints[i]));
                            counterPoints ++;
                            midx += this -> workPoints[j].GetX();
                            midy += this -> workPoints[j].GetY();
                            if(start){
                                leftx = this -> workPoints[j].GetX();
                                rightx = this -> workPoints[j].GetX();
                                topy = this -> workPoints[j].GetY();
                                bottomy = this -> workPoints[j].GetY();
                                start = false;
                            } else {
                                if(leftx > this -> workPoints[j].GetX()){
                                    leftx = this -> workPoints[j].GetX();
                                }
                                if(rightx < this -> workPoints[j].GetX()){
                                    rightx = this -> workPoints[j].GetX();
                                }
                                if(topy < this -> workPoints[j].GetY()){
                                    topy = this -> workPoints[j].GetY();
                                }
                                if(bottomy > this -> workPoints[j].GetY()){
                                    bottomy = this -> workPoints[j].GetY();
                                }
                            }
                        }
                    }
                    bup << std::endl << "id: " << i << std::endl << "left border: " << leftx << std::endl << "right border: "
                        << rightx << std::endl << "top border: " << topy << std::endl << "bottom border: " << bottomy << std::endl
                        << "center point: " << midx / counterPoints << "; " << midy / counterPoints << std::endl
                        << "total points: " << counterPoints << std::endl << std::endl;
                }
                bup << "CLUSTER" << std::endl;
            } else {
                bup << "FIELD" << std::endl << "FIELD" << std::endl;
            }
            for(int i = 0; i < this -> fsize; ++ i){
                bup << this -> workPoints[i].GetX() << " " << this -> workPoints[i].GetY() << " "
                    << this -> workPoints[i].GetLabel() << std::endl;
            }
            bup.close();
        }
/*
        //KMEANS
        void KMeans(int userK){//kmeans with fixed K
            int minid, qStart, qEnd, ptr = 0;
            double minVal;
            bool change, esc = false;
            std::vector<double> MuFunction;
            if(userK < 0){
                qStart = 1;
                qEnd = this -> fsize;
            } else {
                qStart = userK;
                qEnd = userK + 1;
            }
            MuFunction.resize(this -> fsize - 1);
            for(int q = qStart; ((q < qEnd) && (!esc)); q ++){
                MuFunction[ptr] = 0;
                this -> k = q;
                std::vector<std::vector<double>> RO_selected, RO_centers;
                change = true;
                this -> startPoints.clear();
                this -> startPoints.reserve(this -> k);
                RO_selected.resize(this -> k);
                this -> centerKoords.resize(this -> k);
                std::vector<Point> midPoint;
                std::vector<int> numberOfPointsInEachCluster(this -> k);
                for(int i = 0; i < this -> k; ++ i){ //create first matrix of RO between all points by pairs
                    this -> startPoints.push_back(Point(this -> fieldkoord[i]));
                    midPoint.push_back(Point(0, 0, -5));//-5 for mid points
                    numberOfPointsInEachCluster[i] = 0;
                    RO_selected[i].resize(this -> fsize);
                    for(int j = 0; j < this -> fsize; j++){
                        RO_selected[i][j] = Correct(sqrt(Sqr(this -> startPoints[i].GetX() - this -> fieldkoord[j].GetX())
                                                         + Sqr(this -> startPoints[i].GetY() - this -> fieldkoord[j].GetY())));
                    }
                }
                for(int j = 0; j < this -> fsize; ++ j){//mark each point (which center was the nearest)
                    minVal = RO_selected[0][j];
                    minid = 0;
                    for(int i = 1; i < this -> k; ++ i){//find minimal ro from j point to k ros of each cluster
                        if(RO_selected[i][j] < minVal){
                            minVal = RO_selected[i][j];
                            minid = i;
                        }
                    }
                    this -> fieldkoord[j].SetLabel(minid);
                    numberOfPointsInEachCluster[minid] ++;
                }
                for(int i = 0; i < this -> fsize; ++ i){ //caclulate mid points
                    int id = this -> fieldkoord[i].GetLabel();
                    midPoint[id].AddX(Correct(this -> fieldkoord[i].GetX() / numberOfPointsInEachCluster[id]));
                    midPoint[id].AddY(Correct(this -> fieldkoord[i].GetY() / numberOfPointsInEachCluster[id]));

                }
                for(int i = 0; i < this -> k; ++ i){ //make new koord centers
                    this -> centerKoords[i].SetX(midPoint[i].GetX());
                    this -> centerKoords[i].SetY(midPoint[i].GetY());
                }
                while(change){
                    RO_centers.resize(this -> k);
                    change = false;
                    midPoint.clear();
                    numberOfPointsInEachCluster.clear();
                    midPoint.resize(this -> k);
                    numberOfPointsInEachCluster.resize(this -> k);
                    for(int i = 0; i < this -> k; ++ i){ //calc all rors from centers
                        RO_centers[i].resize(this -> fsize);
                        midPoint[i] = Point(0,0,-5);
                        numberOfPointsInEachCluster[i] = 0;
                        for(int j = 0; j < this -> fsize; ++ j){ //find new ros from new centers to each point
                            RO_centers[i][j] = Correct(sqrt(Sqr(this -> centerKoords[i].GetX() - this -> fieldkoord[j].GetX())
                                                            + Sqr(this -> centerKoords[i].GetY() - this -> fieldkoord[j].GetY())));
                        }
                    }
                    for(int j = 0; j < this -> fsize; ++ j){
                        minVal = RO_centers[0][j];
                        minid = 0;
                        for(int i = 1; i < this -> k; ++ i){
                            if(RO_centers[i][j] < minVal){
                                minVal = RO_centers[i][j];
                                minid = i;
                            }
                        }
                        this -> fieldkoord[j].SetLabel(minid);
                        numberOfPointsInEachCluster[minid] ++;
                    }
                    for(int i = 0; i < this -> fsize; ++ i){ //caclulate mid points
                        int id = this -> fieldkoord[i].GetLabel();
                        midPoint[id].AddX(Correct(this -> fieldkoord[i].GetX() / numberOfPointsInEachCluster[id]));
                        midPoint[id].AddY(Correct(this -> fieldkoord[i].GetY() / numberOfPointsInEachCluster[id]));
                    }
                    for(int i = 0; i < this -> k; ++ i){//check definitions between new and old center koords
                        double xOld, yOld, xNew, yNew;
                        xOld = this -> centerKoords[i].GetX();
                        yOld = this -> centerKoords[i].GetY();
                        xNew = midPoint[i].GetX();
                        yNew = midPoint[i].GetY();
                        if((abs(xOld - xNew) > eps) || (abs(yOld - yNew) > eps)){
                            change = true;
                            i = this -> k;
                        }
                    }
                    if(change){
                        this -> centerKoords.clear();
                        this -> centerKoords.resize(this -> k);
                        for(int i = 0; i < this -> k; ++ i){
                            this -> centerKoords[i].SetX(midPoint[i].GetX());
                            this -> centerKoords[i].SetY(midPoint[i].GetY());
                        }
                    }
                }
                for(int i = 0; i < this -> k; ++ i){// here we find optimal k
                    for(int j = 0; j < this -> fsize; ++ j){
                        if(this -> fieldkoord[j].GetLabel() == i){
                            for(int p = j + 1; p < this -> fsize; p ++){
                                if(this -> fieldkoord[p].GetLabel() == i){
                                   MuFunction[ptr] += Correct(sqrt(Sqr(this -> fieldkoord[j].GetX() - this -> fieldkoord[p].GetX()) +
                                                                    Sqr(this -> fieldkoord[j].GetY() - this -> fieldkoord[p].GetY())));
                                }
                            }
                        }
                    }
                }
                for (int i = 0; i < q; ++ i){//here we find optimal k
                    for(int j = i + 1; j < q; ++ j){
                        MuFunction[ptr] += Correct(sqrt(Sqr(this -> centerKoords[i].GetX() - this -> centerKoords[j].GetX()) +
                                                      Sqr(this -> centerKoords[i].GetY() - this -> centerKoords[j].GetY())));
                    }
                }
                if((q != 1) && ((qEnd - qStart) > 1)){
                    if(MuFunction[ptr - 1] < MuFunction[ptr]){
                        esc = true;
                        ptr -= 2;
                    }
                }
                ptr ++;
            }
            this -> optimalK = ptr + 1;
            this -> bestWeight = MuFunction[ptr];
            this -> workPoints.clear();
            for(unsigned int i = 0; i < this -> fieldkoord.size(); ++ i){
                this -> workPoints.push_back(this -> fieldkoord[i]);
            }
        }
        int RetBestK(){//optimal K for kmeans
            return this -> optimalK;
        }
        //SPTR
        void CalculateRo(){//calculate all Ro between all points
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    if(i == j){
                        this -> RO[i][j] = -1;
                    } else if(j == 0){
                        this -> RO[i][j] = -1;
                    } else {
                        this -> RO[i][j] = sqrt(Sqr(this -> workPoints[i].GetX() - this -> workPoints[j].GetX()) +
                                                Sqr(this -> workPoints[i].GetY() - this -> workPoints[j].GetY()));
                    }
                }
            }
        }
        void FillLengths(){//making std::vector of distances for GNU gistogram
            std::vector<int> index;
            std::fstream tree, script, text, launchlog;
            std::vector<std::vector<double>> TempRo = this -> RO;
            double minRo, mr;
            int minId = -1, mi = -1, saveIND = -1;
            launchlog.open("launch.log", ios::in | ios::out);
            launchlog >> this -> launchNumber;
            this -> launchNumber ++;
            launchlog << this -> launchNumber;
            launchlog.close();
            std::string filename = "Tree" + std::to_string(this -> launchNumber) + ".txt", scriptname = "plotTree" + std::to_string(this -> launchNumber) + ".plt",
                sstr = "TempGisto" + std::to_string(this -> launchNumber) + ".txt", ssstr = "plotTempGisto" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, ios::out | ios::trunc);
            script << "plot '" << filename << "' using 1:2 with lines lc rgb \"black\" lw 1 notitle";
            script.close();
            tree.open(filename, ios::out | ios::trunc);
            index.resize(0);
            index.push_back(0);
            while(this -> allLengths.size() < this -> workPoints.size() - 1){
                mr = -1;
                for(const auto &id : index){//for each id from index array
                    minRo = TempRo[id][0];
                    for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                        if((minRo < 0) && (TempRo[id][i] > 0)){
                            minRo = TempRo[id][i];
                            minId = i;
                        } else if ((TempRo[id][i] < minRo) && (TempRo[id][i] > 0)){
                            minRo = TempRo[id][i];
                            minId = i;
                        }
                    }
                    if((mr < 0) || (mr > minRo)){
                        mr = minRo;
                        mi = minId;
                        saveIND = id;
                    }
                }
                tree << this -> workPoints[saveIND].GetX() << " " <<this -> workPoints[saveIND].GetY() << endl
                    << this -> workPoints[mi].GetX() << " " <<this -> workPoints[mi].GetY() << std::endl << endl;
                this -> allLengths.push_back(mr);
                index.push_back(mi);
                for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                    TempRo[i][mi] = -1;
                }
            }
            tree.close();
            script.open(ssstr, ios::out | ios::trunc);
            script << "width=0.01" << std::endl << "bin(x, s) = s*int(x/s) + width/2"
                << std::endl << "set boxwidth width" << std::endl << "plot '" << sstr << "' u (bin($1,width)):(1.0) \\"
                << std::endl << "s f w boxes fs solid 0.5 title 'Porog Gisto'" << endl;
            script.close();
            text.open(sstr, ios::out | ios::trunc); //writing values to the file
            for(unsigned int i = 0; i < (this -> workPoints.size() - 1); ++ i){
                text << this -> allLengths[i] << endl;
            }
            text.close();
        }
        //HIERARCHY
        void Hierarchy(int userK){
            double x1, y1, x2, y2, minro, midx, midy;
            std::fstream tree, script, launchlog;
            launchlog.open("launch.log", ios::in | ios::out);
            launchlog >> this -> launchNumber;
            launchlog.close();
            std::string filename = "HieTree" + std::to_string(this -> launchNumber) + ".txt";
            std::string scriptname = "plotHieTree" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, ios::out | ios::trunc);
            script << "plot '" << filename << "' using 1:2 with lines lc rgb \"black\" lw 2 notitle";
            script.close();
            launchlog.open("launch.log", ios::trunc);
            launchlog.close();
            tree.open(filename, ios::out | ios::trunc);
            std::vector<Point> newPoints = this -> workPoints;
            int p1, p2id of point 1 and point 2, numpoints = this -> workPoints.size(), numclusters = userK;
            std::vector<std::vector<double>> ro, rro;
            while(numpoints != numclusters){
                ro.resize(newPoints.size());
                minro = -1;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    ro[i].resize(newPoints.size());
                    x1 = newPoints[i].GetX();
                    y1 = newPoints[i].GetY();
                    for(unsigned int j = 0; j < newPoints.size(); ++ j){
                            if(i != j){
                            x2 = newPoints[j].GetX();
                            y2 = newPoints[j].GetY();
                            ro[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                            if(((minro < 0) || (minro > ro[i][j])) && (ro[i][j] > 0)){
                                minro = ro[i][j];
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
                tree << x1 << " " << y1 << std::endl << x2 << " " << y2 << std::endl << endl;
                midx = (x1 + x2) / 2.0;
                midy = (y1 + y2) / 2.0;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    ro[i].erase(ro[i].begin() + std::max(p1, p2));
                    ro[i].erase(ro[i].begin() + std::min(p1, p2));
                }
                ro.erase(ro.begin() + std::max(p1, p2));
                ro.erase(ro.begin() + std::min(p1, p2));
                newPoints.erase(newPoints.begin() + std::max(p1, p2));
                newPoints.erase(newPoints.begin() + std::min(p1, p2));
                newPoints.push_back(Point(midx, midy, 0));
                numpoints --;
            }
            tree.close();
            rro.resize(newPoints.size());
            for(unsigned int i = 0; i < newPoints.size(); ++ i){
                rro[i].resize(this -> workPoints.size());
                x1 = newPoints[i].GetX();
                y1 = newPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    rro[i][j] = Correct(sqrt(Sqr(x1 - x2) + Sqr(y1 - y2)));
                }
            }
            for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                minro = rro[0][j];
                p1 = 0;
                for(unsigned int i = 0; i < newPoints.size(); ++ i){
                    if(rro[i][j] < minro){
                        minro = rro[i][j];
                        p1 = i;
                    }
                }
                this -> workPoints[j].SetLabel(p1);
            }
        }
        //FOREL
        void Forel(double threshold){
            std::fstream circle, script, text;//, launchlog;
            std::vector<double> ro;//matrix of distances
            std::vector<int> binary, id;
            bool change;
            int startID, counterPoints, ptr, clustID = 0;
            double midx, midy;
            this -> launchNumber = ReadLaunch();
            RewriteLaunch(this -> launchNumber + 1);
            std::string filename = "CircleF" + std::to_string(this -> launchNumber) + ".txt";
            std::string scriptname = "plotCircleF" + std::to_string(this -> launchNumber) + ".plt";
            script.open(scriptname, ios::out | ios::trunc);
            script << "plot '" << filename << "' using 1:2:3 with circles lc rgb \"black\" lw 2 notitle";
            script.close();
            circle.open(filename, ios::out | ios::trunc);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> workPoints[i].SetLabel(-i);
            }
            std::vector<Point> tempPoints = this -> workPoints;
            std::unique_ptr <Point> center (new Point(tempPoints[0]));
            while(tempPoints.size() > 0){
                startID = rand () % tempPoints.size();
                ro.resize(tempPoints.size());
                binary.resize(tempPoints.size());
                change = true;
                center.reset(new Point(tempPoints[startID]));
                while(change){
                    midx = 0;
                    midy = 0;
                    counterPoints = 0;
                    for(unsigned int i = 0; i < tempPoints.size(); ++ i){
                        ro[i] = sqrt(Sqr(center -> GetX() - tempPoints[i].GetX()) + Sqr(center -> GetY() - tempPoints[i].GetY()));
                        if(ro[i] < threshold){
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
                    if((abs(center -> GetX() - midx) < eps) && (abs(center -> GetY() - midy) < eps)){//when center stops (in that iteration)
                        change = false;
                        circle << midx << " " << midy << " " << threshold << endl;
                        id.resize(counterPoints);
                        ptr = 0;
                        for(unsigned int i = 0; ((i < tempPoints.size()) && (ptr < counterPoints)); ++ i){
                            if(binary[i] == 1){
                                id[ptr] = tempPoints[i].GetLabel();
                                ptr ++;
                            }
                        }
                        ptr = counterPoints - 1;
                        for(unsigned int i = this -> workPoints.size() - 1; ((i >= 0) && (ptr >= 0)); i --){//making label_s
                            if(this -> workPoints[i].GetLabel() == id[ptr]){
                                this -> workPoints[i].SetLabel(clustID);
                                ptr --;
                            }
                        }
                        ptr = counterPoints - 1;
                        for(unsigned int i = tempPoints.size() - 1; ((i >= 0) && (ptr >= 0)); i --){//delete found cluster
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
            }   circle.close();
        }
        */
        //RECOVERY
        void GetFromFile(){//getting clusters from file
            std::string filename, specialCommand, tempstr;
            std::cout << std::endl << "Please enter backup filename" << std::endl;
            std::cin >> filename;
            std::fstream bup;
            bool read;
            double x, y;
            int l;
            std::vector<Point> tempvec;
            tempvec.resize(0);
            bup.open(filename);
            if(!bup.is_open()){
                std::cout << std::endl << "File open error" << std::endl;
            } else {
                read = true;
                getline(bup, specialCommand);
                if((specialCommand != "CLUSTER") && (specialCommand != "FIELD")){
                    Error(7);
                } else {
                    while(!bup.eof()){
                        if(read){
                            getline(bup, tempstr);
                            if(tempstr == specialCommand){
                                read = false;
                            } else {
                                std::cout << tempstr << std::endl;
                            }
                        } else {
                            bup >> x >> y >> l;
                            tempvec.push_back(Point(x, y, l));
                        }
                    }
                }
            }
            bup.close();
            this -> workPoints = tempvec;
            this -> BeforeStart();
        }
        /*
        //DBSCAN
        void CalcRo(){//matrix of distances
            double x1, y1, x2, y2;
            this -> RO.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> RO[i].resize(this -> workPoints.size());
                x1 = this -> workPoints[i].GetX();
                y1 = this -> workPoints[i].GetY();
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    x2 = this -> workPoints[j].GetX();
                    y2 = this -> workPoints[j].GetY();
                    this -> RO[i][j] = sqrt(Sqr(x1 - x2) + Sqr(y1 - y2));
                }
            }
        }
        std::vector<int> CountNearPoints(std::vector<std::vector<double>> RoMatrix, double radius){
            std::vector<int> ctr;
            ctr.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                ctr[i] = 0;
                for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                    if(RoMatrix[i][j] < (radius + eps)){
                        ctr[i] ++;
                    }
                }
            }
            return ctr;
        }
        std::vector<int> CreateMarks(std::vector<std::vector<double>> ro, double radius, std::vector<int> cnt, int EntNum){//label_s for clusters
            std::vector<int> marker;
            marker.resize(this -> workPoints.size());
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                if(cnt[i] > EntNum){
                    marker[i] = 1;
                } else {
                    marker[i] = -1;
                }
            }
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                if(marker[i] == -1){
                    for(unsigned int j = 0; j < this -> workPoints.size(); ++ j){
                        if((ro[i][j] < (radius + eps)) && (marker[j] == 1)){
                            marker[i] = 0;
                        }
                    }
                }
            }
            return marker;
        }
        void DelPoints(double radius, int EntryNum){//delete rubbish
            this -> fieldkoord.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                this -> fieldkoord.push_back(Point(this -> workPoints[i]));
            }
            this -> CalcRo();
            this -> counters = CountNearPoints(this -> RO, radius);
            this -> marks = CreateMarks(this -> RO, radius, this -> counters, EntryNum);
            for(int i = this -> fieldkoord.size() - 1; i >= 0; i --){
                if(this -> marks[i] == -1){
                    this -> fieldkoord.erase(this -> fieldkoord.begin() + i);
                    this -> counters.erase(this -> counters.begin() + i);
                    for(unsigned int j = 0; j < this -> fieldkoord.size(); ++ j){
                        this -> RO[j].erase(this -> RO[j].begin() + i);
                    }
                    this -> RO.erase(this -> RO.begin() + i);
                    this -> marks.erase(this -> marks.begin() + i);
                    this -> workPoints[i].SetLabel(-9999);
                }
            }
        }
        void DBSCANAlghorithm(double radius, int numpointsincircle){//wave on points without rubbish
            this -> BeforeStart();
            saveData.resize(0);
            for(unsigned int i = 0; i < this -> workPoints.size(); ++ i){
                saveData.push_back(Point(this -> workPoints[i]));
            }
            this -> DelPoints(radius, numpointsincircle);
            this -> workPoints.resize(0);
            for(unsigned int i = 0; i < this -> fieldkoord.size(); ++ i){
                this -> workPoints.push_back(Point(this -> fieldkoord[i]));
            }
            this -> BeforeStart();
            this -> CreateRo();
            this -> GenBinary(radius + eps);
            this -> Wave();
        }
        //EM Algorithm
        double CentRo(){
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
            mult.open("mult.txt", ios::ate | ios::in | ios::out);
            for(int i = 0; i < this -> k; ++ i){
                mult << this -> workPoints[i].GetX() << " " <<
                this -> workPoints[i].GetY() << " " << clrs[this -> workPoints[i].GetLabel()][0] << " " <<
                clrs[this -> workPoints[i].GetLabel()][1] << " " << clrs[this -> workPoints[i].GetLabel()][2] << endl;
            }
            mult << std::endl << endl;
            mult.close();
        }
        void EM(int userK){//EM with fixed K
            int itercount = 0;
            this -> k = this -> workPoints.size();
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
            for(int i = 0; i < this -> k; ++ i){
                d[i].resize(userK);
                p[i].resize(userK);
            }
            this -> CentersB.resize(0);
            for(int i = 0; i < userK; ++ i){
                this -> CentersB.push_back(Point(this -> workPoints[i]));
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
            mult.open("mult.txt", ios::trunc | ios::in | ios::out);
            mult.close();
            mult.open("mult.plt", ios::trunc | ios::in | ios::out);
            mult << "set terminal gif animate delay 100" << endl;
            mult << "set output 'MULT.gif'" << endl;
            mult << "stats 'mult.txt' nooutput" << std::endl << endl;
            mult << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)" << std::endl << endl;
            mult << "do for[i=1:int(STATS_blocks)]{" << endl;
                mult << "plot 'mult.txt' index(i-1) using 1:2:(rgb($3,$4,$5)) with points lc rgb variable" << endl;
            mult << "}" << endl;
            mult.close();
            bool esc = true;
                while(esc){
                itercount ++;
                    for(int i = 0; i < this -> k; ++ i){
                        summ = 0;
                        for(int j = 0; j < userK; ++ j){
                            d[i][j] = sqrt(SigmaM[j][0][0] * Sqr(this -> workPoints[i].GetX() - this -> CentersB[j].GetX()) +
                                   (SigmaM[j][1][0] + SigmaM[j][0][1]) * (this -> workPoints[i].GetX() -
                                    this -> CentersB[j].GetX()) * (this -> workPoints[i].GetY() - this -> CentersB[j].GetY()) +
                                    SigmaM[j][1][1] * Sqr(this -> workPoints[i].GetY() - this -> CentersB[j].GetY()));
                            summ += d[i][j];
                        }
                        for(int j = 0; j < userK; ++ j){
                            p[i][j] = 1.0 / (sUk - 1) - d[i][j] / (sUk - 1) / summ;
                        }
                    }
                    this -> CentersA.resize(0);
                    for(int j = 0; j < userK; ++ j){
                        sumpx = 0;
                        sumpy = 0;
                        sump = 0;
                        for(int i = 0; i < this -> k; ++ i){
                            sumpx += this -> workPoints[i].GetX() * p[i][j];
                            sumpy += this -> workPoints[i].GetY() * p[i][j];
                            sump += p[i][j];
                        }
                        this -> CentersA.push_back(Point(sumpx / sump, sumpy / sump));
                    }
                    if(this -> CentRo() > 0){
                        this -> CentersB.resize(0);
                        for(int i = 0; i < userK; ++ i){
                            this -> CentersB.push_back(CentersA[i]);
                        }
                        for(int j = 0; j < userK; ++ j){
                            Sigma[j][0][0] = 0;
                            Sigma[j][0][1] = 0;
                            Sigma[j][1][0] = 0;
                            Sigma[j][1][1] = 0;
                            for(int i = 0; i < this -> k; ++ i){
                                Sigma[j][0][0] += p[i][j] * Sqr(this -> workPoints[i].GetX() - this -> CentersB[j].GetX());
                                Sigma[j][0][1] += p[i][j] * (this -> workPoints[i].GetX() - this -> CentersB[j].GetX()) *
                                                            (this -> workPoints[i].GetY() - this -> CentersB[j].GetY());
                                Sigma[j][1][1] += p[i][j] * Sqr(this -> workPoints[i].GetY() - this -> CentersB[j].GetY());
                            }
                            Sigma[j][1][0] = Sigma[j][0][1];

                            SigmaM[j][0][0] = Sigma[j][1][1] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][0][1] = - Sigma[j][0][1] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][1][0] = - Sigma[j][1][0] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                            SigmaM[j][1][1] = Sigma[j][0][0] / abs(Sigma[j][0][0] * Sigma[j][1][1] -
                                        Sigma[j][1][0] * Sigma[j][0][1]);
                        }
                        if((itercount % 2) == 1){
                            for(int i = 0; i < this -> k; i++){
                                mx = p[i][0];
                                index = 0;
                                for(int j = 0; j < userK; ++ j){
                                    if(p[i][j] > mx){
                                        mx = p[i][j];
                                        index = j;
                                    }
                                }
                                this -> workPoints[i].SetLabel(index);
                            }
                            this -> GNUMultOut(itercount, Mycolors);
                        }
                    } else {
                        esc = false;
                        for(int i = 0; i < this -> k; i++){
                            mx = p[i][0];
                            index = 0;
                            for(int j = 0; j < userK; ++ j){
                                if(p[i][j] > mx){
                                    mx = p[i][j];
                                    index = j;
                                }
                            }
                            this -> workPoints[i].SetLabel(index);
                        }
                        this -> GNUMultOut(itercount, Mycolors);
                    }
                    this -> CentersB.resize(0);
                    for(int i = 0; i < userK; ++ i){
                        this -> CentersB.push_back(Point(this -> CentersA[i]));
                    }
                }
        }
        */
};

class Find{
    private:
        std::string algname, filename;//names of algorithms and files
        std::vector<Point> findrez;//std::vector of findings
        unsigned int numclust, numrubb;
        std::vector<std::vector<double>> centers;//centers of clusters
        std::vector<double> leftborders, rightborders, topborders, butttomborders;//borders of clusters
        std::vector<int>  numpointsincluster;//number of points in clusters
        std::vector<std::vector<double>> Amatrix, Lambda, R, Z;

    public:
        Find(){//default constroctor
            this -> algname = "";
            this -> filename = "";
            this -> numclust = 0;
            this -> numrubb = 0;
            this -> findrez.resize(0);
            this -> centers.resize(0);
            this -> leftborders.resize(0);
            this -> rightborders.resize(0);
            this -> topborders.resize(0);
            this -> butttomborders.resize(0);
            this -> numpointsincluster.resize(0);
            this -> Amatrix.resize(0);
            this -> Lambda.resize(0);
            this -> R.resize(0);
            this -> Z.resize(0);
        }
        Find(std::string alg, const std::vector<Point>& values){//cluster finding with factors
            this -> findrez.resize(0);
            for(unsigned int i = 0; i < values.size(); i++){
                this -> findrez.push_back(Point(values[i]));
            }
            this -> algname = alg;
            int kkk = -1;
            for(unsigned int i = 0; i < this -> findrez.size(); ++ i){
                if(this -> findrez[i].GetLabel() > kkk){
                    kkk = findrez[i].GetLabel();
                }
            }
            kkk ++;
            this -> numclust = kkk;
            this -> Amatrix.resize(kkk);
            this -> Lambda.resize(kkk);
            this -> R.resize(kkk);
            this -> Z.resize(kkk);
            this -> centers.resize(kkk);
            this -> leftborders.resize(0);
            this -> rightborders.resize(0);
            this -> topborders.resize(0);
            this -> butttomborders.resize(0);
            this -> numpointsincluster.resize(0);
            double midx, midy, leftx, rightx, topy, bottomy, sqX, xy, sqY;
            int counterPoints;
            bool start;
            for(int i = 0; i < kkk; ++ i){
                counterPoints = 0;
                midx = 0;
                midy = 0;
                start = true;
                this -> centers[i].resize(2);
                for(unsigned int j = 0; j < this -> findrez.size(); ++ j){
                    if(this -> findrez[j].GetLabel() == i){
                        counterPoints ++;
                        midx += this -> findrez[j].GetX();
                        midy += this -> findrez[j].GetY();
                        if(start){
                            leftx = this -> findrez[j].GetX();
                            rightx = this -> findrez[j].GetX();
                            topy = this -> findrez[j].GetY();
                            bottomy = this -> findrez[j].GetY();
                            start = false;
                        } else {
                            if(leftx > this -> findrez[j].GetX()){
                                leftx = this -> findrez[j].GetX();
                            }
                            if(rightx < this -> findrez[j].GetX()){
                                rightx = this -> findrez[j].GetX();
                            }
                            if(topy < this -> findrez[j].GetY()){
                                topy = this -> findrez[j].GetY();
                            }
                            if(bottomy > this -> findrez[j].GetY()){
                                bottomy = this -> findrez[j].GetY();
                            }
                        }
                    }
                }
                sqX = 0;
                sqY = 0;
                xy = 0;
                this -> centers[i][0] = midx / counterPoints;
                this -> centers[i][1] = midy / counterPoints;
                for(unsigned int j = 0; j < this -> findrez.size(); ++ j){
                    if(this -> findrez[j].GetLabel() == i){
                        sqX += Sqr(this -> findrez[j].GetX() - this -> centers[i][0]);
                        sqY += Sqr(this -> findrez[j].GetY() - this -> centers[i][1]);
                        xy += (this -> findrez[j].GetX() - this -> centers[i][0]) *
                            (this -> findrez[j].GetY() - this -> centers[i][1]);
                    }
                }
                this -> Amatrix[i].resize(3);
                sqX /= counterPoints;
                xy /= counterPoints;
                sqY /= counterPoints;
                this -> Amatrix[i][0] = sqX;
                this -> Amatrix[i][1] = xy;
                this -> Amatrix[i][2] = sqY;
                this -> Lambda[i].resize(2);
                this -> Lambda[i][0] = (sqX + sqY + sqrt(Sqr(sqX - sqY) + 4 * Sqr(xy))) / 2;
                this -> Lambda[i][1] = (sqX + sqY - sqrt(Sqr(sqX - sqY) + 4 * Sqr(xy))) / 2;
                this -> R[i].resize(2);
                this -> R[i][0] = 2 * xy * sqrt(this -> Lambda[i][0]) / sqrt(Sqr(sqX - this -> Lambda[i][0]) + Sqr(xy));
                this -> R[i][1] = 2 * xy * sqrt(this -> Lambda[i][1]) / sqrt(Sqr(sqX - this -> Lambda[i][1]) + Sqr(xy));
                this -> Z[i].resize(2);
                this -> Z[i][0] = - 2 * (sqX  - this -> Lambda[i][0]) * sqrt(this -> Lambda[i][0])
                    / sqrt(Sqr(sqX - this -> Lambda[i][0]) + Sqr(xy));
                this -> Z[i][1] = - 2 *  (sqX  - this -> Lambda[i][1]) * sqrt(this -> Lambda[i][1])
                    / sqrt(Sqr(sqX - this -> Lambda[i][1]) + Sqr(xy));
                this -> leftborders.push_back(leftx);
                this -> rightborders.push_back(rightx);
                this -> topborders.push_back(topy);
                this -> butttomborders.push_back(bottomy);
                this -> numpointsincluster.push_back(counterPoints);
            }
            this -> algname = alg;
            this -> filename = "";
        }
        void totalBackup(){//backup of all findings
            std::fstream totbup;
            totbup.open("totalBackup.log", std::fstream::app);
            totbup << this -> algname << std::endl << this -> findrez.size() << std::endl;
            for(unsigned int i = 0; i < this -> findrez.size(); ++ i){
                totbup << this -> findrez[i].GetX() << " " << this -> findrez[i].GetY() << " "
                    << this -> findrez[i].GetLabel() << std::endl;
            }
            totbup.close();
        }
        void ShowAlgorithm(){
            std::cout << this -> algname << std::endl;
        }
        void ClearTotBackup(){//clearing
            std::fstream totbup;
            totbup.open("totalBackup.log", std::fstream::trunc);
            totbup.close();
        }
        void ToTxtCluster(std::string name){//printing clasters in file. Files set by user
            std::vector<std::vector<int>> rgb;
            int kkk = this -> numclust, l;
            std::fstream script, myFile, bup;
            std::string filename = "", scriptname = "plot", backupname;
            rgb.resize(kkk);
            for(int i = 0; i < kkk; ++ i){
                rgb[i].resize(3);
                rgb[i][0] = rand() % 256;
                rgb[i][1] = rand() % 256;
                rgb[i][2] = rand() % 256;
            }
            for(unsigned i = 0; i < name.size(); ++ i){
                if(name[i] == '.'){
                    i = name.size();
                } else {
                    filename += name[i];
                }
            }
            scriptname += filename + ".plt";
            backupname = filename + "Backup.txt";
            filename = name;
            script.open(scriptname, std::fstream::trunc);
            script << "rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b) " << std::endl
                << "plot '" << filename << "' using 1:2:(rgb($3,$4,$5)) with points lc rgb variable";
            script.close();
            myFile.open(filename, std::fstream::trunc);
            for(unsigned int i = 0; i < this -> findrez.size(); ++ i){
                l = this -> findrez[i].GetLabel();
                if(l >= 0){
                    myFile << this -> findrez[i].GetX() << " " << this -> findrez[i].GetY() << " "
                        << rgb[l][0] << " " << rgb[l][1] << " " << rgb[l][2] << std::endl;
                }
            }
            myFile.close();
            bup.open(backupname, std::fstream::trunc);
            if((name != "Field.txt") && (name != "Field")){
                bup << "CLUSTER" << std::endl << "Total clusters: " << kkk << std::endl
                    << "Total points: " << this -> findrez.size() << std::endl;
                for(unsigned int i = 0; i < this -> numclust; ++ i){
                    bup << std::endl << "id: " << i << std::endl << "left border: " << this -> leftborders[i] << std::endl
                        << "right border: " << this -> rightborders[i] << std::endl
                        << "top border: " << this -> topborders[i] << std::endl
                        << "bottom border: " << this -> butttomborders[i] << std::endl
                        << "center point: " << this -> centers[i][0] << "; " << this -> centers[i][0] << std::endl
                        << "total points: " << this ->numpointsincluster[i] << std::endl << std::endl;
                }
                bup << "CLUSTER" << std::endl;
            } else {
                bup << "FIELD" << std::endl << "FIELD" << std::endl;
            }
            for(unsigned int i = 0; i < this -> findrez.size(); ++ i){
                bup << this -> findrez[i].GetX() << " " << this -> findrez[i].GetY() << " "
                    << this -> findrez[i].GetLabel() << std::endl;
            }
            bup.close();
        }
        std::vector<Point> retPoints(){
            return this -> findrez;
        }
        std::string retName(){
            return this -> algname;
        }
        void Factors(int ID){//making file for factors
            std::fstream script, data;
            std::string sfilename = "Factors" + std::to_string(ID) + ".plt", filename = "Factors" + std::to_string(ID) + ".txt";
            script.open(sfilename, std::fstream::trunc);
            script << "plot '" + filename + "' using 1:2:3:4 with std::vectors filled head lw 3";
            script.close();
            data.open(filename, std::fstream::trunc);
            for(unsigned int i = 0; i < this -> Lambda.size(); ++ i){
                data << this -> centers[i][0] << " " << this -> centers[i][1] << " " << this -> R[i][0]
                    << " " << this -> Z[i][0] << std::endl << this -> centers[i][0] << " " << this -> centers[i][1]
                    << " " << this -> R[i][1] << " " << this -> Z[i][1] << std::endl << std::endl;
            }
            data.close();
        }
};

class Cluster{
    private:
        //std::unique_ptr <Field> MyField;
        std::vector<Point> resultFind, workPoints;
        int launchNumber;
    public:
        Cluster(){//default constructor
            std::fstream launcher;
            launcher.open("launch.log",  std::fstream::trunc);
            launcher << 1;
            launcher.close();
            //this -> MyField.reset(new Field(0));
            this -> resultFind.resize(0);
            this -> workPoints.resize(0);
        }
        Cluster(Field *mf){//constructor
            std::fstream launcher;
            launcher.open("launch.log", std::fstream::out);
            launcher >> this -> launchNumber;
            launcher.close();
            this -> resultFind.resize(0);
            this -> workPoints = mf -> retALLKOORD();
        }
        std::vector<Point> retResult(){
            return this -> resultFind;
        }
        void FindByWaveAlgorithm(){//RUN WAVE
            /*this -> MyField -> BeforeStart();
            this -> MyField -> CreateRo();
            double threshold = -1;
            std::string filename = "Wave" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            while(threshold < 0){
                cout << std::endl << "Please enter your threshold value" << endl;
                cin >> threshold;
                if(threshold < 0){
                    Error(5);
                }
            }
            this -> MyField -> GenBinary(threshold);
            this -> MyField -> Wave();
            this -> resultFind = this -> MyField -> retALLKOORD();*/
            double threshold = -1;
            std::string filename = "Wave" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            while(threshold < 0){
                std::cout << std::endl << "Please enter your threshold value" << std::endl;
                std::cin >> threshold;
                if(threshold < 0){
                    Error(5);
                }
            }
            std::unique_ptr <FindByWave> findcl(new FindByWave(this -> workPoints));
            this -> resultFind = findcl -> FindClusters(threshold);
        }
        void FindByKMeansAlgorithm(){//RUN KMEANS
            /*this -> MyField -> BeforeStart();*/
            int userNumberOfClusters = -100;
            while((userNumberOfClusters != -1) && (userNumberOfClusters < 0)){
                std::cout << std::endl << "Please enter your number of clusters or '-1' to find optimal number" << std::endl;
                std::cin >> userNumberOfClusters;
                if(userNumberOfClusters > int(this -> workPoints.size())){
                    userNumberOfClusters = -1000;
                    std::cout << std::endl << "Too big number" << std::endl;
                }
                if((userNumberOfClusters != -1) && (userNumberOfClusters < 0)){
                    Error(5);
                }
            }
            /*this -> MyField -> KMeans(userNumberOfClusters);
            if(userNumberOfClusters == -1){
                cout << std::endl << "Optimal number of clusters is " << this -> MyField -> RetBestK() << endl
                    << "You could rerun k-means algorithm with that value" << endl;
            } else {
                std::string filename = "K-Means" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
                this -> launchNumber ++;
                RewriteLaunch(this -> launchNumber);
                this -> resultFind = this -> MyField -> retALLKOORD();
            }*/
            std::unique_ptr <FindByKM> findcl(new FindByKM(this -> workPoints));
            this -> resultFind = findcl -> FindClusters(userNumberOfClusters);
        }
        void RunSpainningTreeAlgorithm(){//RUN SPTR
            /*std::string scn1 = "TempGisto" + std::to_string(this -> launchNumber);
            this -> MyField -> BeforeStart();
            this -> MyField -> CalculateRo();
            this -> MyField -> FillLengths();
            cout << std::endl << "You could plot '" << scn1 << "' and choose correct threshold value" << endl;
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            this -> FindByWaveAlgorithm();*/
            double threshold = -1;
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            while(threshold < 0){
                std::cout << std::endl << "Please enter your threshold value" << std::endl;
                std::cin >> threshold;
                if(threshold < 0){
                    Error(5);
                }
            }
            std::unique_ptr <FindBySPTR> findcl(new FindBySPTR(this -> workPoints));
            this -> resultFind = findcl -> FindClusters(threshold);
        }
        void FindByHierarchyAlgorithm(){
            /*this -> MyField -> BeforeStart();*/
            int userNumberOfClusters = -100;
            while((userNumberOfClusters != -1) && (userNumberOfClusters < 0)){
                std::cout << std::endl << "Please enter your number of clusters" << std::endl;
                std::cin >> userNumberOfClusters;
                if(userNumberOfClusters > int(this -> workPoints.size())){
                    userNumberOfClusters = -1000;
                }
                if(userNumberOfClusters < 0){
                    Error(5);
                }
            }
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            std::unique_ptr <FindByHierarchy> findcl(new FindByHierarchy(this -> workPoints));
            this -> resultFind = findcl -> FindClusters(userNumberOfClusters);
        }
        void FindByForelAlghorithm(){//RUN FOREL
           /* this -> MyField -> BeforeStart();*/
            double threshold = -1;
            while(threshold < 0){
                std::cout << std::endl << "Please enter your threshold" << std::endl;
                std::cin >> threshold;
                if(threshold < 0){
                    Error(5);
                }
            }
            /*this -> MyField -> Forel(threshold);
            std::string filename = "Forel" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            this -> launchNumber ++;
            RewriteLaunch(this ->  launchNumber);
            std::unique_ptr <Field> TempField;
            TempField.reset(new Field(this -> MyField -> retALLKOORD()));
            this -> resultFind = TempField -> retALLKOORD();*/
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            std::unique_ptr <FindByForel> findcl(new FindByForel(this -> workPoints));
            this -> resultFind = findcl -> FindClusters(threshold);
        }
        void FindByDBSCANAlghorithm(){//RUN DBSCAN
            double radius = -1;
            std::string filename = "DBSCAN" + std::to_string(this -> launchNumber) + ".txt", accept = "zzz";
            int NumNearPoints = -1;
            while(radius < 0){
                std::cout << std::endl << "Please enter your radius" << std::endl;
                std::cin >> radius;
                if(radius < 0){
                    Error(5);
                }
            }
            while(NumNearPoints < 0){
                std::cout << std::endl << "Please enter number of points in circle" << std::endl;
                std::cin >> NumNearPoints;
                if(NumNearPoints < 0){
                    Error(5);
                }
            }
            /*this -> MyField -> DBSCANAlghorithm(radius, NumNearPoints);
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            this -> resultFind = this -> MyField -> retALLKOORD();*/
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            std::unique_ptr <FindByDbscan> findcl(new FindByDbscan(this -> workPoints));
            this -> resultFind = findcl -> FindClusters(radius, NumNearPoints);
        }
        void FindByEMAlghorithm(){ //RUN EM
            /*this -> MyField -> BeforeStart();*/
            int userNumberOfClusters = -100;
            while(userNumberOfClusters < 0){
                std::cout << std::endl << "Please enter your number of clusters" << std::endl;
                std::cin >> userNumberOfClusters;
                if(userNumberOfClusters > int(this -> workPoints.size())){
                    userNumberOfClusters = -1000;
                    std::cout << std::endl << "Too big number" << std::endl;
                }
                if(userNumberOfClusters < 0){
                    Error(5);
                }
            }
           /* this -> MyField -> EM(userNumberOfClusters);
            std::string filename = "EM" + std::to_string(this -> launchNumber) + ".txt";
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            this -> resultFind = this -> MyField -> retALLKOORD();*/
            this -> launchNumber ++;
            RewriteLaunch(this -> launchNumber);
            std::unique_ptr <FindByEM> findcl(new FindByEM(this -> workPoints));
            this -> resultFind = findcl -> FindClusters(userNumberOfClusters);
        }
};

class Interface{//hand interface
    private:
        int sz;//size of field
        std::vector<Find> FindList;//std::vector of findings
    public:
        std::string command;

        Interface(){//constructor
            this -> sz = 1000;
            this -> FindList.resize(0);
        }
        void PrintHelp(){//help printing in console
            std::cout << std::endl << " Help:" << std::endl << " 'genrnd' for DEMO ravn generation of the std::vector" << std::endl
                      << " 'gennorm' for DEMO norm generation of the std::vector" << std::endl << " 'genfield' to make field of groups"
                      << std::endl << " 'gengroup' to make one DEMO group" << std::endl << " 'moveX' to move X DEMO" << std::endl
                      << " 'moveY' to move Y DEMO" << std::endl << " 'rollN' to turn group(0;0) DEMO" << std::endl
                      << " 'rollC' to turn group (center) DEMO" << std::endl << " 'exit' to exit" << std::endl << std::endl
                      << " 'setsize' to set new basic size" << std::endl << std::endl << " 'clear' to clear total backup file" << std::endl
                      << std::endl << " 'clog' to clear log file" << std::endl << std::endl;
        }
        void run(){//method of using of the commands by the user
            std::string command = "start";
            int trg;
            this -> PrintHelp();
            while(command != "exit"){
                trg = -1;
                std::cout << std::endl << " Please enter the command" << std::endl;
                std::cin >> command;
                Mylog << "command: '" << command << "'" << std::endl;
                if((command == "setsize")){//generation of ravn std::vector (DEMO)
                    int newsz = 0;
                    while(newsz < 1){
                        std::cout << std::endl << "Please enter size (number of points)" << std::endl;
                        std::cin >> newsz;
                        Mylog << newsz << " points in one group" << std::endl;
                        if(newsz < 1){
                            Error(4);
                        }
                    }
                    this -> sz = newsz;
                    trg = 0;
                }
                if((command == "genrnd")){//generation of ravn std::vector (DEMO)
                    double mn, mx;
                    std::cout << std::endl << " Please enter min and max value of random" << std::endl;
                    std::cin >> mn >> mx;
                    Mylog << mn << " " << mx << std::endl;
                    std::unique_ptr <Control> RND(new Control());
                    RND -> GenRnd(this -> sz, mn, mx);
                    RND -> FileRavn();
                    trg = 0;
                }
                else if((command == "clear")){
                    std::unique_ptr <Find> MyFind(new Find());
                    MyFind -> ClearTotBackup();
                    trg = 0;
                }
                else if((command == "clog")){
                    std::fstream f;
                    f.open("log.log", std::fstream::trunc);
                    f.close();
                    trg = 0;
                }
                else if(command == "gennorm"){//generation of norm std::vector (DEMO)
                    double mn, mx;
                    std::cout << std::endl << " Please enter min and max value of random" << std::endl;
                    std::cin >> mn >> mx;
                    Mylog << mn << " " << mx << std::endl;
                    std::unique_ptr <Control> NORM(new Control());
                    NORM -> GenNorm(this -> sz, mn, mx);
                    NORM -> FileNorm();
                    trg = 0;
                }
                else if(command == "gengroup"){//generation of a group of points (DEMO)
                    double mnx, mxx, mny, mxy;
                    std::cout << std::endl << " Please enter min and max value of random (x)" << std::endl;
                    std::cin >> mnx >> mxx;
                    Mylog << mnx << " " << mxx << std::endl;
                    std::cout << std::endl << " Please enter min and max value of random (y)" << std::endl;
                    std::cin >> mny >> mxy;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> GROUP(new Control());
                    GROUP -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0);
                    GROUP -> FileGroup();
                    trg = 0;
                }
                else if(command == "rollN"){//turning of the group (0;0) (DEMO)
                    double phi, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    std::cout << std::endl << "Please create new group to show this function" << std::endl
                              << std::endl << " Please enter min and max value of random (x)" << std::endl;
                    std::cin >> mnx >> mxx;
                    Mylog << mnx << " " << mxx << std::endl;
                    std::cout << std::endl << " Please enter min and max value of random (y)" << std::endl;
                    std::cin >> mny >> mxy;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream rotated;
                    rotated.open("rotateN.txt", std::fstream::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    std::cout << std::endl << " Please enter angle" << std::endl;
                    std::cin >> phi;
                    Mylog << phi << std::endl;
                    phi = DegreesToRadians(phi);
                    work -> TurnNULL(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        rotated << gr[i].GetX() << " " << gr[i].GetY() << std::endl;
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "rollC"){//turning of the group (Center) (DEMO)
                    double phi, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    std::cout << std::endl << "Please create new group to show this function" << std::endl
                              << std::endl << " Please enter min and max value of random (x)" << std::endl;
                    std::cin >> mnx >> mxx;
                    Mylog << mnx << " " << mxx << std::endl;
                    std::cout << std::endl << " Please enter min and max value of random (y)" << std::endl;
                    std::cin >> mny >> mxy;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream rotated;
                    rotated.open("rotS.txt", std::fstream::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    std::cout << std::endl << " Please enter angle" << std::endl;
                    std::cin >> phi;
                    Mylog << phi << std::endl;
                    phi = DegreesToRadians(phi);
                    work -> turnCenter(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        rotated << gr[i].GetX() << " " << gr[i].GetY() << std::endl;
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "moveX"){//moving X (DEMO)
                    double delx, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    std::cout << std::endl << " Please create new group to show this function" << std::endl
                              << std::endl << " Please enter min and max value of random (x)" << std::endl;
                    std::cin >> mnx >> mxx;
                    Mylog << mnx << " " << mxx << std::endl;
                    std::cout << std::endl << " Please enter min and max value of random (y)" << std::endl;
                    std::cin >> mny >> mxy;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream moved;
                    moved.open("moveX.txt", std::fstream::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    std::cout << std::endl << " Please enter delta X" << std::endl;
                    std::cin >> delx;
                    Mylog << delx << std::endl;
                    work -> MoveX(delx);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        moved << gr[i].GetX() << " " << gr[i].GetY() << std::endl;
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "moveY"){//moving Y (DEMO)
                    double dely, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    std::cout << std::endl << " Please create new group to show this function" << std::endl
                              << std::endl << " Please enter min and max value of random (x)" << std::endl;
                    std::cin >> mnx >> mxx;
                    Mylog << mnx << " " << mxx << std::endl;
                    std::cout << std::endl << " Please enter min and max value of random (y)" << std::endl;
                    std::cin >> mny >> mxy;
                    Mylog << mny << " " << mxx << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream moved;
                    moved.open("moveY.txt", std::fstream::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    std::cout << std::endl << " Please enter delta Y" << std::endl;
                    std::cin >> dely;
                    Mylog << dely << std::endl;
                    work -> MoveY(dely);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        moved << gr[i].GetX() << " " << gr[i].GetY() << std::endl;
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "genfield"){//field creating
                    int numofgroups = 0, tr;
                    double mnx, mny, mxx, mxy;
                    std::string com = "rrr";
                    std::unique_ptr <Field> MyField;
                    std::unique_ptr <Group> temp;
                    while((com != "get") && (com!="work") && (com != "rec")){
                        std::cout << std::endl << "Please enter 'get' to read field from file" << std::endl
                                  << std::endl << "Please enter 'work' to work in basik mode" << std::endl
                                  << std::endl << "Please enter 'rec' to run total recovery" << std::endl;
                        std::cin >> com;
                        if((com != "get") && (com != "work") && (com != "rec")){
                            Mylog << "command: " << com << std::endl;
                            Error(5);//incorrect value
                        }
                    }
                    if(com == "rec"){
                        Mylog << "command: " << com << std::endl;
                        std::fstream recfile;
                        this -> FindList.resize(0);
                        recfile.open("totalBackup.log");
                        while(!recfile.eof()){
                            std::string algname, numpts;
                            int npts, l;
                            double x, y;
                            std::vector<Point> tempvec;
                            getline(recfile, algname);
                            if(algname.size() > 1){
                                recfile >> npts;
                                for(int i = 0; i < npts; ++ i){
                                    recfile >> x >> y >> l;
                                    tempvec.push_back(Point(x, y, l));
                                }
                                this -> FindList.push_back(Find(algname, tempvec));
                            }
                        }
                        std::vector<Point> newKoord = this -> FindList[this -> FindList.size() - 1].retPoints();
                        MyField.reset(new Field(newKoord));
                    }
                    else if(com == "work"){
                        Mylog << "command: 'work'" << std::endl;
                        while(numofgroups < 1){
                            std::cout << std::endl << " Please enter number of groups in field (>=1)" << std::endl;
                            std::cin >> numofgroups;
                            Mylog <<"  "<< numofgroups << " groups generated: " << std::endl;
                            if(numofgroups < 1){
                                std::cout << std::endl << " Error! Please enter correctly value" << std::endl;
                                Mylog << "Error! Incorrect value" << std::endl;
                            }
                        }
                        MyField.reset(new Field(numofgroups));
                        for(int i = 0; i < numofgroups; ++ i){
                            com = "start";
                            std::cout << std::endl << " Please enter min and max value of random (x)" << std::endl;
                            std::cin >> mnx >> mxx;
                            Mylog << "  X adges: " << mnx << " " << mxx << std::endl;
                            std::cout << std::endl << " Please enter min and max value of random (y)" << std::endl;
                            std::cin >> mny >> mxy;
                            Mylog << "  Y adges: " << mny << " " << mxy << std::endl;
                            temp.reset(new Group(this -> sz, mnx, mny, mxx, mxy, i));
                            while(com != "esc"){
                                tr = -1;
                                std::cout << std::endl << " Enter 'RN' to turn group(0;0)" << std::endl << std::endl
                                          << " Enter 'RC' to turn group (center)" << std::endl << std::endl << " Enter 'MX' to move X"
                                          << std::endl << std::endl << " Enter 'MY' to move Y" << std::endl << std::endl
                                          << " Enter 'esc' to finish field creating" << std::endl;
                                std::cin >> com;
                                if(com == "RN"){//turning of the group (0;0)
                                    Mylog << "  RN - turning (0,0)" << std::endl;
                                    double alf;
                                    std::cout << std::endl << " Please enter angle" << std::endl;
                                    std::cin >> alf;
                                    Mylog << "  angle: "<< alf << std::endl;
                                    alf = DegreesToRadians(alf);
                                    std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    work -> TurnNULL(alf);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                else if(com == "RC"){//turning of the group (Center)
                                    Mylog << "  RC - turning (Center)" << std::endl;
                                    double alf;
                                    std::cout << std::endl << " Please enter angle" << std::endl;
                                    std::cin >> alf;
                                    Mylog << "  angle: "<< alf << std::endl;
                                    alf = DegreesToRadians(alf);
                                    std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    work -> turnCenter(alf);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                else if(com == "MX"){//moving X
                                    Mylog << "  MX - moving X" << std::endl;
                                    double dx;
                                    std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    std::cout << std::endl << " Enter delta x" << std::endl;
                                    std::cin >> dx;
                                    Mylog <<"  delta x: "<<   dx << std::endl;
                                    work -> MoveX(dx);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                else if(com == "MY"){//moving Y
                                    Mylog << "  MY - moving Y" << std::endl;
                                    double dy;
                                    std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                    std::cout << std::endl << " Enter delta y" << std::endl;
                                    std::cin >> dy;
                                    Mylog << "  delta y: "<< dy << std::endl;
                                    work -> MoveY(dy);
                                    temp -> Regrupp(work);
                                    tr = 0;
                                }
                                if(tr == -1){
                                    Mylog << "command: 'esc'" << std::endl;
                                    if(com != "esc"){
                                        Mylog << "command: " << com << std::endl;
                                        Error(5); //incorrect value
                                        com = "start";
                                    }
                                }
                            }
                            MyField -> AddGroup(temp);
                        }
                    } else if (com == "get"){
                        Mylog << "command: 'get'" << std::endl;
                        MyField -> GetFromFile();
                        MyField -> ToTxtCluster("Field.txt");
                    }
                    if((com != "get") && (com != "rec")){
                        MyField -> MakeAllKoord();
                        std::string answ = "wait";
                        while((answ != "yes") && (answ != "no")){
                            std::cout << std::endl << "Do you want to add color? ('yes'/'no') " << std::endl;
                            std::cin >> answ;
                            if((answ != "yes") && (answ != "no")){
                                Error(5);
                            }
                        }
                        if(answ == "yes"){
                            MyField -> ToTxt("yes");
                        } else {
                            MyField -> ToTxt("no");
                        }
                    }
                    std::string q = "NO";
                    while((q != "yes") && (q != "no")){
                        std::cout << std::endl << "Do you want to find clusters? ('yes'/'no')" << std::endl;
                        std::cin >> q;
                        if((q != "yes") && (q != "no")){
                            Mylog << "command: " << q << std::endl;
                            Error(5);
                        }
                    }
                    if(q == "no"){
                        Mylog << "No clusters found" << std::endl;
                        q = "esc";
                    }
                    std::unique_ptr <Cluster> MyCluster(new Cluster(MyField.get()));
                    while(q != "esc"){
                        std::cout << std::endl << "To Find by Wave enter 'wave'" << std::endl << "To Find by K-means enter 'km'" << std::endl
                                  << "To Spainning Tree enter 'sptr'" << std::endl << "To Hierarchy enter 'ie'" << std::endl
                                  << "To Forel alghorithm enter 'fish'" << std::endl << "To DBSCAN enter 'dbs'" << std::endl
                                  << "To EM enter 'em'" << std::endl << "To Exit enter 'esc'" << std::endl;
                        std::cin >> q;
                        if(q == "wave"){
                            Mylog << "WAVE algorithm" << std::endl;
                            MyCluster -> FindByWaveAlgorithm();
                            this -> FindList.push_back(Find("WAVE algorithm", MyCluster -> retResult()));
                            this -> FindList[0].totalBackup();
                        } else if (q == "km"){
                            Mylog << "K-MEANS algorithm" << std::endl;
                            MyCluster -> FindByKMeansAlgorithm();
                            std::cout << "3.0" << std::endl;
                            this -> FindList.push_back(Find("K-MEANS algorithm", MyCluster -> retResult()));
                            std::cout << "4.0" << std::endl;
                            this -> FindList[0].totalBackup();
                            std::cout << "5.0" << std::endl;
                        } else if (q == "sptr"){
                            Mylog << "SPANNING TREE algorithm" << std::endl;
                            MyCluster -> RunSpainningTreeAlgorithm();
                            this -> FindList.push_back(Find("SPANNING TREE algorithm", MyCluster -> retResult()));
                            this -> FindList[0].totalBackup();
                        } else if (q == "ie"){
                            Mylog << "HIERARCHY algorithm" << std::endl;
                            MyCluster -> FindByHierarchyAlgorithm();
                            this -> FindList.push_back(Find("HIERARCHY algorithm", MyCluster -> retResult()));
                            this -> FindList[0].totalBackup();
                        } else if (q == "fish"){
                            Mylog << "FOREL algorithm" << std::endl;
                            MyCluster -> FindByForelAlghorithm();
                            this -> FindList.push_back(Find("FOREL algorithm", MyCluster -> retResult()));
                            this -> FindList[0].totalBackup();
                        } else if (q == "dbs"){
                            Mylog << "DBSCAN algorithm" << std::endl;
                            MyCluster -> FindByDBSCANAlghorithm();
                            this -> FindList.push_back(Find("DBSCAN algorithm", MyCluster -> retResult()));
                            this -> FindList[0].totalBackup();
                        } else if (q == "em"){
                            Mylog << "EM algorithm" << std::endl;
                            MyCluster -> FindByEMAlghorithm();
                            this -> FindList.push_back(Find("EM algorithm", MyCluster -> retResult()));
                            this -> FindList[0].totalBackup();
                        } else if (q == "esc"){
                            Mylog << "command: esc" << std::endl;
                            std::string question = "zzz";
                            while(question != "esc"){
                                while((question != "no") && (question != "yes") && (question != "esc")){
                                    std::cout << std::endl << "Do you want to write some results into file? ('yes'/'no')" << std::endl;
                                    std::cin >> question;
                                    if((question!= "yes") && (question != "no")){
                                            Mylog << "command: " << question << std::endl;
                                            Error(5);
                                    }
                                }
                                if(question == "yes"){
                                    std::cout << std::endl << "  Which algorithm write into file?" << std::endl;
                                    for(unsigned int i = 0; i < this -> FindList.size(); ++ i){
                                        std::cout << std::endl << i + 1 << " ";
                                        this -> FindList[i].ShowAlgorithm();
                                    }
                                    int ID = -1;
                                    while((ID < 1) || (unsigned(ID) > (this -> FindList.size() + 1))){
                                        std::cout << std::endl << "  Please enter algorithm ID" << std::endl;
                                        std::cin >> ID;
                                        Mylog << "To file: alg " << ID << std::endl;
                                    }
                                    std::string accept = "zzz", filename;
                                    while((accept != "yes") && (accept != "no")){
                                        std::cout << std::endl << "Do you want to set your file name? ('yes'/'no')" << std::endl;
                                        std::cin >> accept;
                                        if((accept != "yes") && (accept != "no")){
                                            Mylog << "command: " << accept << std::endl;
                                            Error(5);
                                        }
                                    }
                                    if(accept == "yes"){
                                        std::cout << std::endl << "Please enter your filename" << std::endl;
                                        std::cin >> filename;
                                        Mylog << "filename is: " << filename << std::endl;
                                    } else {
                                        int vers = ReadLaunch();
                                        filename = this -> FindList[ID].retName() + "_" + std::to_string(vers) + ".txt";
                                        vers ++;
                                        RewriteLaunch(vers);
                                    }
                                    this -> FindList[ID - 1].ToTxtCluster(filename);
                                    question = "zzz";
                                    std::string fact = "zzz";
                                    while((fact != "yes") && (fact != "no")){
                                        std::cout << std::endl << "Do you want to add Factors to that algorithm? ('yes'/'no')" << std::endl;
                                        std::cin >> fact;
                                        if((fact != "yes") && (fact != "no")){
                                            Mylog << "command: " << fact << std::endl;
                                            Error(5);
                                        }
                                    }
                                    if(fact == "yes"){
                                        Mylog << " with factors" << std::endl;
                                        this -> FindList[ID - 1].Factors(ID);
                                    }
                                }
                                if(question == "no"){
                                    question = "esc";
                                }
                            }
                        } else if (q != "esc"){
                            Mylog << "command: " << q << std::endl;
                            Error(5);
                        }
                    }
                    trg = 0;
                }
                if(trg == -1){
                    if(command == "exit"){
                        break;
                    } else{
                        Error(5);
                    }
                }
                else if(trg == 0){
                    std::cout << std::endl << " Please enter the command" << std::endl;
                    PrintHelp();
                    getline(std::cin, command);
                    trg = -1;
                }
            }
        }
};

class InterfaceSTR{//file interface
    private:
        int sz;
        std::vector<std::string> comlist;
        std::vector<Find> FindList;
    public:
        std::string command;

        InterfaceSTR(){//constructor
            this -> sz = 1000;
            this -> FindList.resize(0);
            std::string fname, tempstr;
            this -> comlist.resize(0);
            std::cout << std::endl << "Please enter master filename" << std::endl;
            std::cin >> fname;
            Mylog << "Master filename is: " << fname << std::endl;
            std::fstream master;
            master.open(fname);
            if(!master.is_open()){
                std::cout << std::endl << "File open error" << std::endl;
                Mylog << "File open error" << std::endl;
            } else {
                while(!master.eof()){
                    master >> fname;
                    if(fname != ""){
                        comlist.push_back(std::string(fname));
                    }
                }
            }
        }
        void run(){//method of using of the commands by the user
            std::string command = "start";
            int trg, pos = 0;
            while(command != "exit"){
                trg = -1;
                command = this -> comlist[pos];
                pos ++;
                if((command == "setsize")){//generation of ravn std::vector (DEMO)
                    Mylog << "command: 'setsize'"<< std::endl;
                    int newsz = 0;
                    while(newsz < 1){
                        newsz = stoi(this -> comlist[pos]);
                        pos ++;
                        Mylog << "  "<< newsz << " points in one group" << std::endl;
                        if(newsz < 1){
                            Error(4);
                        }
                    }
                    this -> sz = newsz;
                    trg = 0;
                }
                if((command == "genrnd")){//generation of ravn std::vector (DEMO)
                    Mylog << "command: 'genrnd'" << std::endl;
                    double mn, mx;
                    mn = stod(this -> comlist[pos]);
                    pos ++;
                    mx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mn << " "<< mx << std::endl;
                    std::unique_ptr <Control> RND(new Control());
                    RND -> GenRnd(this -> sz, mn, mx);
                    RND -> FileRavn();
                    trg = 0;
                }
                else if((command == "clear")){
                    std::unique_ptr <Find> MyFind(new Find());
                    MyFind -> ClearTotBackup();
                    trg = 0;
                }
                else if(command == "gennorm"){//generation of norm std::vector (DEMO)
                    Mylog << "command: 'gennorm'"<< std::endl;
                    double mn, mx;
                    mn = stod(this -> comlist[pos]);
                    pos ++;
                    mx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mn << " "<< mx << std::endl;
                    std::unique_ptr <Control> NORM(new Control());
                    NORM -> GenNorm(this -> sz, mn, mx);
                    NORM -> FileNorm();
                    trg = 0;
                }
                else if(command == "gengroup"){//generation of a group of points (DEMO)
                    Mylog << "command: 'gengroup'" << std::endl;
                    double mnx, mxx, mny, mxy;
                    mnx = stod(this -> comlist[pos]);
                    pos ++;
                    mxx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mnx << " " << mxx << std::endl;
                    mny = stod(this -> comlist[pos]);
                    pos ++;
                    mxy = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> GROUP(new Control());
                    GROUP -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0);
                    GROUP -> FileGroup();
                    trg = 0;
                }
                else if(command == "rollN"){//turning of the group (0;0) (DEMO)
                    Mylog << "command: 'rollN'" << std::endl;
                    double phi, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    mnx = stod(this -> comlist[pos]);
                    pos ++;
                    mxx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mnx << " " << mxx << std::endl;
                    mny = stod(this -> comlist[pos]);
                    pos ++;
                    mxy = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream rotated;
                    rotated.open("rotateN.txt", std::fstream::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    phi = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << "  angle: "<< phi << std::endl;
                    phi = DegreesToRadians(phi);
                    work -> TurnNULL(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        rotated << gr[i].GetX() << " ";
                        rotated << gr[i].GetY() << "\n";
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "rollC"){//turning of the group (Center) (DEMO)
                    Mylog << "command: 'rollC'"<< std::endl;
                    double phi, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    mnx = stod(this -> comlist[pos]);
                    pos ++;
                    mxx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mnx << " " << mxx << std::endl;
                    mny = stod(this -> comlist[pos]);
                    pos ++;
                    mxy = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream rotated;
                    rotated.open("rotS.txt", std::fstream::trunc);
                    if(!rotated.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    phi = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << "  angle: "<< phi << std::endl;
                    phi = DegreesToRadians(phi);
                    work -> turnCenter(phi);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        rotated << gr[i].GetX() << " ";
                        rotated << gr[i].GetY() << "\n";
                    }
                    rotated.close();
                    trg = 0;
                }
                else if(command == "moveX"){//moving X (DEMO)
                    Mylog << "command: 'moveX'"<< std::endl;
                    double delx, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    mnx = stod(this -> comlist[pos]);
                    pos ++;
                    mxx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mnx << " " << mxx << std::endl;
                    mny = stod(this -> comlist[pos]);
                    pos ++;
                    mxy = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream moved;
                    moved.open("moveX.txt", std::fstream::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    delx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << delx << std::endl;
                    work -> MoveX(delx);
                    gr = work -> RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        moved << gr[i].GetX() << " ";
                        moved << gr[i].GetY() << "\n";
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "moveY"){//moving Y (DEMO)
                    Mylog << "command: 'moveY'"<< std::endl;
                    double dely, mnx, mny, mxx, mxy;
                    std::vector<Point> gr;
                    mnx = stod(this -> comlist[pos]);
                    pos ++;
                    mxx = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mnx << " " << mxx << std::endl;
                    mny = stod(this -> comlist[pos]);
                    pos ++;
                    mxy = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << mny << " " << mxy << std::endl;
                    std::unique_ptr <Control> work (new Control());
                    work -> GenGroup(this -> sz, mnx, mxx, mny, mxy, 0); //base group
                    std::fstream moved;
                    moved.open("moveY.txt", std::fstream::trunc);
                    if(!moved.is_open()){
                        Error(3);
                    }
                    work -> FileGroup();
                    dely = stod(this -> comlist[pos]);
                    pos ++;
                    Mylog << dely << std::endl;
                    work -> MoveY(dely);
                    gr = work ->RetGroup();
                    for(int i = 0; i < this -> sz; ++ i){
                        moved << gr[i].GetX() << " ";
                        moved << gr[i].GetY() << "\n";
                    }
                    moved.close();
                    trg = 0;
                }
                else if(command == "genfield"){//field creating
                    Mylog << "command: 'genfield'"<< std::endl;
                    int numofgroups = 0, tr;
                    double mnx, mny, mxx, mxy;
                    std::string com = this -> comlist[pos];
                    pos ++;
                    std::unique_ptr <Field> MyField;
                    std::unique_ptr <Cluster> MyCluster(new Cluster());
                    std::unique_ptr <Group> temp;
                    if((com != "rec") && (com != "get") && (com!="work")){
                        Mylog << "command: " << com << std::endl;
                        Error(5);//incorrect value
                    } else {
                        if(com == "rec"){
                            Mylog << "command: 'rec'" << std::endl;
                            std::fstream recfile;
                            this -> FindList.resize(0);
                            recfile.open("totalBackup.txt");
                            while(!recfile.eof()){
                                std::string algname;
                                int npts, l;
                                double x, y;
                                std::vector<Point> tempvec;
                                tempvec.resize(0);
                                recfile >> algname;
                                recfile >> npts;
                                for(int i = 0; i < npts; ++ i){
                                    recfile >> x >> y >> l;
                                    tempvec.push_back(Point(x, y, l));
                                }
                                this -> FindList.push_back(Find(algname, tempvec));
                            }
                        } else if(com == "work"){
                            Mylog << "command: 'work'" << std::endl;
                            numofgroups = stoi(this -> comlist[pos]);
                            pos ++;
                            Mylog << "  " << numofgroups << " groups will be created:" << std::endl;
                            if(numofgroups < 1){
                                Error(5);
                            }
                            MyField.reset(new Field(numofgroups));
                            for(int i = 0; i < numofgroups; ++ i){
                                mnx = stod(this -> comlist[pos]);
                                pos ++;
                                mxx = stod(this -> comlist[pos]);
                                pos ++;
                                Mylog << "  X adges: " << mnx << " " << mxx << std::endl;
                                mny = stod(this -> comlist[pos]);
                                pos ++;
                                mxy = stod(this -> comlist[pos]);
                                pos ++;
                                Mylog << "  Y adges: " << mny << " " << mxy << std::endl;
                                temp.reset(new Group(this -> sz, mnx, mny, mxx, mxy, i));
                                com = "start";
                                while(com != "esc"){
                                    tr = -1;
                                    com = this -> comlist[pos];
                                    pos ++;
                                    if(com == "RN"){//turning of the group (0;0)
                                        Mylog << "  'RN' - turning (0,0)" << std::endl;
                                        double alf = DegreesToRadians(stod(this -> comlist[pos]));
                                        pos ++;
                                        Mylog << "  angle: " << alf;
                                        Mylog << ";  angle in radians: " << alf << std::endl;
                                        std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        work -> TurnNULL(alf);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    else if(com == "RC"){//turning of the group (Center)
                                        Mylog << "  'RC' - turning (Center)" << std::endl;
                                        double alf = DegreesToRadians(stod(this -> comlist[pos]));
                                        pos ++;
                                        Mylog << "  angle: " << alf;
                                        Mylog << ";  angle in radians: " << alf << std::endl;
                                        std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        work -> turnCenter(alf);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    else if(com == "MX"){//moving X
                                        Mylog << "  'MX' - moving x" << std::endl;
                                        double dx = stod(this -> comlist[pos]);
                                        std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        pos ++;
                                        Mylog <<"  delta x: "<< dx << std::endl;
                                        work -> MoveX(dx);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    else if(com == "MY"){//moving Y
                                        Mylog << "  'MY' - moving y" << std::endl;
                                        double dy = stod(this -> comlist[pos]);
                                        std::unique_ptr <Control> work = temp -> Ret_Field_Group();
                                        pos ++;
                                        Mylog <<"  delta y: "<< dy << std::endl;
                                        work -> MoveY(dy);
                                        temp -> Regrupp(work);
                                        tr = 0;
                                    }
                                    if(tr == -1){
                                        Mylog << "command: 'esc'" << std::endl;
                                        if(com != "esc"){
                                            Mylog << "command: " << com << std::endl;
                                            Error(5); //incorrect value
                                        }
                                    }
                                }
                                MyField -> AddGroup(temp);
                            }
                        }
                        else if (com == "get"){
                            Mylog << "command: 'get'" << std::endl;
                            MyField -> GetFromFile();
                            MyField -> ToTxtCluster("Field.txt");
                            MyCluster.reset(new Cluster(MyField.get()));
                        }
                        if(com != "get"){
                            MyField -> MakeAllKoord();
                            MyField -> ToTxt(this -> comlist[pos]);
                            pos ++;
                            MyCluster.reset(new Cluster(MyField.get()));
                        }
                        std::string q = this -> comlist[pos];
                        pos ++;
                        if((q != "yes") && (q != "no")){
                            Mylog << "command: " << q << std::endl;
                            Error(5);
                        }
                        if(q == "no"){
                            Mylog << "No clusters found" << std::endl;
                            q = "esc";
                        }
                        while(q != "esc"){
                            std::cout << std::endl << "To Find by Wave enter 'wave'" << std::endl << "To Find by K-means enter 'km'" << std::endl
                                      << "To Spainning Tree enter 'sptr'" << std::endl << "To Hierarchy enter 'ie'" << std::endl
                                      << "To Forel alghorithm enter 'fish'" << std::endl << "To DBSCAN enter 'dbs'" << std::endl
                                      << "To Exit enter 'esc'" << std::endl;
                            std::cin >> q;
                            if(q == "wave"){
                                Mylog << "WAVE algorithm" << std::endl;
                                MyCluster -> FindByWaveAlgorithm();
                                this -> FindList.push_back(Find("WAVE algorithm", MyCluster -> retResult()));
                                this -> FindList[0].totalBackup();
                            } else if (q == "km"){
                                Mylog << "K-MEANS algorithm" << std::endl;
                                MyCluster -> FindByKMeansAlgorithm();
                                this -> FindList.push_back(Find("K-MEANS algorithm", MyCluster -> retResult()));
                                this -> FindList[0].totalBackup();
                            } else if (q == "sptr"){
                                Mylog << "SPANNING TREE algorithm" << std::endl;
                                MyCluster -> RunSpainningTreeAlgorithm();
                                this -> FindList.push_back(Find("SPANNING TREE algorithm", MyCluster -> retResult()));
                                this -> FindList[0].totalBackup();
                            } else if (q == "ie"){
                                Mylog << "HIERARCHY algorithm" << std::endl;
                                MyCluster -> FindByHierarchyAlgorithm();
                                this -> FindList.push_back(Find("HIERARCHY algorithm", MyCluster -> retResult()));
                                this -> FindList[0].totalBackup();
                            } else if (q == "fish"){
                                Mylog << "FOREL algorithm" << std::endl;
                                MyCluster -> FindByForelAlghorithm();
                                this -> FindList.push_back(Find("FOREL algorithm", MyCluster -> retResult()));
                                this -> FindList[0].totalBackup();
                            } else if (q == "dbs"){
                                Mylog << "DBSCAN algorithm" << std::endl;
                                MyCluster -> FindByDBSCANAlghorithm();
                                this -> FindList.push_back(Find("DBSCAN algorithm", MyCluster -> retResult()));
                                this -> FindList[0].totalBackup();
                            } else if (q == "em"){
                                Mylog << "EM algorithm" << std::endl;
                                MyCluster -> FindByEMAlghorithm();
                                this -> FindList.push_back(Find("EM algorithm", MyCluster -> retResult()));
                                this -> FindList[0].totalBackup();
                            } else if (q == "esc"){
                                Mylog << "esc" << std::endl;
                                std::string question = "zzz";
                                while(question != "esc"){
                                    while((question != "no") && (question != "yes") && (question != "esc")){
                                        std::cout << std::endl << "Do you want to write some results into file? 'yes'/'no'" << std::endl;
                                        std::cin >> question;
                                        if((question != "yes") && (question != "no")){
                                                Mylog << "command: " << question << std::endl;
                                                Error(5);
                                        }
                                    }
                                    if(question == "yes"){
                                        std::cout << std::endl << "  Which algorithm write into file?" << std::endl;
                                        for(unsigned int i = 0; i < this -> FindList.size(); ++ i){
                                            std::cout << std::endl << i + 1 << " ";
                                            this -> FindList[i].ShowAlgorithm();
                                        }
                                        int ID = -1;
                                        while((ID < 1) || (unsigned(ID) > (this -> FindList.size() + 1))){
                                            std::cout << std::endl << "  Please enter algorithm ID" << std::endl;
                                            std::cin >> ID;
                                            Mylog << "To file: alg " << ID << std::endl;
                                        }
                                        std::string accept = "zzz", filename;
                                        while((accept != "yes") && (accept != "no")){
                                            std::cout << std::endl << "Do you want to set your file name? ('yes'/'no')" << std::endl;
                                            std::cin >> accept;
                                            if((accept != "yes") && (accept != "no")){
                                                Mylog << "command: " << accept << std::endl;
                                                Error(5);
                                            }
                                        }
                                        if(accept == "yes"){
                                            std::cout << std::endl << "Please enter your filename" << std::endl;
                                            std::cin >> filename;
                                            Mylog << "filename is: " << filename << std::endl;
                                        } else {
                                            int vers = ReadLaunch();
                                            filename = this -> FindList[ID].retName() + "_" + std::to_string(vers) + ".txt";
                                            Mylog << "filename is: " << filename << std::endl;
                                            std::cout << "filename is: " << filename << std::endl;
                                            vers ++;
                                            RewriteLaunch(vers);
                                        }
                                        this -> FindList[ID - 1].ToTxtCluster(filename);
                                        question = "zzz";
                                        std::string fact = "zzz";
                                        while((fact != "yes") && (fact != "no")){
                                            std::cout << std::endl << "Do you want to add Factors to that algorithm? ('yes'/'no')" << std::endl;
                                            std::cin >> fact;
                                            if((fact != "yes") && (fact != "no")){
                                                Mylog << "command: " << fact << std::endl;
                                                Error(5);
                                            }
                                        }
                                        if(fact == "yes"){
                                            Mylog << " with factors" << std::endl;
                                            this -> FindList[ID - 1].Factors(ID);
                                        }
                                    }
                                    if(question == "no"){
                                        question = "esc";
                                    }
                                }
                            } else if (q != "esc"){
                                Mylog << "command: " << q << std::endl;
                                Error(5);
                            }
                        }
                        command = this -> comlist[pos];
                        trg = 0;
                    }
                    if(trg == -1){
                        if(command == "exit"){
                            Mylog << "command: 'exit'" << std::endl;
                            break;
                        } else{
                            Mylog << "command: " << command << std::endl;
                            Error(5);
                        }
                    }
                    else if(trg == 0){
                        command = this -> comlist[pos];
                        Mylog << "command: " << command << std::endl;
                        pos ++;
                        trg = -1;
                    }
                }
            }
        }
};

int main(){
    srand(time(NULL));
    setlocale(LC_ALL, "");
    Mylog.open("log.log", std::fstream::app);
    std::string var = "zzz";
    while((var != "file") && (var != "hand")){
        std::cout << std::endl << "Please choose mode ('file'/'hand')" << std::endl;
        getline(std::cin, var);
        if((var != "file") && (var != "hand")){
            Mylog << "command: " << var << std::endl;
            Error(5);
        }
    }
    if(var == "hand"){
        Mylog << "command: 'hand'" << std::endl;
        std::unique_ptr <Interface> MyProject (new Interface());
        MyProject -> run();
    } else {
        Mylog << "command: 'file'" << std::endl;
        std::unique_ptr <InterfaceSTR> MyProject (new InterfaceSTR());
        MyProject -> run();
    }
    Mylog.close();
    return 0;
}
