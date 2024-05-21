#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>

using namespace std;

class Person {
public:
    void ChangeFirstName(const int& inp_year, const string& inp_first_name) {
        first_name[inp_year] = inp_first_name;
    }

    void ChangeLastName(const int& inp_year, const string& inp_last_name) {
        last_name[inp_year] = inp_last_name;
    }

    string GetFullName(const int& inp_year) {
        string f_name = GetNameByYear(inp_year, first_name);
        string l_name = GetNameByYear(inp_year, last_name);
        string result;
        if ((l_name == "") && (f_name == "")) {
            result = "Incognito";
            return result;
        }
        if (l_name == "") {
            result = f_name + " with unknown last name";
            return result;
        }
        if (f_name == "") {
            result = l_name + " with unknown first name";
            return result;
        }
        result = f_name + " " + l_name;
        return result;
    }
private:
    map <int, string> first_name, last_name;

    string GetNameByYear (const int& inp_year, const map <int, string>& name_list) {
        string name = "";
        for (const auto &pair : name_list) {
            if (pair.first <= inp_year) {
                name = pair.second;
            } else {
                break;
            }
        }
        return name;
    }
};

int main() {
    return 0;
}
