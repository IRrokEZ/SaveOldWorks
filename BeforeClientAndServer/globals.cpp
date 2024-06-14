#include "globals.hpp"

const double PI = 3.1415926535897932384626433832795;
const double eps = 0.0000001;

std::fstream Mylog;

void RewriteLaunch (int val) { // launcher
    std::fstream launcher;
    launcher.open("launch.log", std::fstream::trunc);
    launcher << val;
    launcher.close();
}

int ReadLaunch () { // launch reading
    std::fstream launcher;
    int val;
    launcher.open("launch.log");
    launcher >> val;
    launcher.close();
    return val;
}

double DegreesToRadians (double degree) { // convertion of degrees to radians
    return (degree * PI / 180);
}

void Error (int type) { // error warning // REWRITE to enum
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

bool IsInVector (const std::vector<int> &vect, int search_elem) { // checking of element in std::vector
   return std::find(vect.begin(), vect.end(), search_elem) != std::end(vect);
}

double Sqr (double val) {
    return val * val;
}

double Correct (double val) { // change too little value to 0
    if (val < eps) {
        return 0.0;
    }
    return val;
}


double ABS (double val) { // calculate abs of double value
    if (val < 0) {
        return val * -1.0;
    }
    return val;
}