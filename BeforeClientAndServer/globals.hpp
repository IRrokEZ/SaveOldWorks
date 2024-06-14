#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

const double PI = 3.1415926535897932384626433832795;
const double eps = 0.0000001;

std::fstream Mylog;

void RewriteLaunch (int val); // launcher

int ReadLaunch (); // launch reading

double DegreesToRadians (double degree); // convertion of degrees to radians

void Error (int type); // error warning

bool IsInVector (const std::vector<int> &vect, int search_elem); // checking of element in vector

double Sqr (double val);

double Correct (double val); // change too little value to 0

double ABS (double val); // calculate abs of double value