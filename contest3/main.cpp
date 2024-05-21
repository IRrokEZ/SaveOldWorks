#include <iostream>
#include <map>
#include <set>
#include <algorithm>

struct Person {
    bool inTeam;
    long timeIn;
    long previousExp;
};

std::map<std::string, Person> persons;
std::set<std::pair<long, std::string>> expSet;

std::string getMax(std::string& newName, long newExp, std::string& oldName, long time) {
    auto oldExpPerson = persons.at(oldName);
    long oldExp = oldExpPerson.previousExp + (time - oldExpPerson.timeIn);
    if (newExp > oldExp) {
        return newName;
    }
    if (newExp == oldExp) {
        return std::min(newName, oldName);
    }
    return oldName;
}

int main() {
    std::string mostExp = "";
    long expSum = 0;
    long peopleInTeam = 0;
    long prevTime = 0;

    long count;
    std::cin >> count;
    for(long i = 0; i < count; ++i) {
        std::string name;
        long time;
        std::cin >> name >> time;
        expSum += peopleInTeam * (time - prevTime);

        if (persons.find(name) != persons.end()) {
            auto person = persons.at(name);
            if (person.inTeam) {
                peopleInTeam--;
                expSet.erase(std::make_pair(-1 * (person.previousExp - person.timeIn), name));
                person.inTeam = false;
                person.previousExp += time - person.timeIn;
                if (expSet.empty()) {
                    mostExp = "";
                } else {
                    mostExp = expSet.begin()->second;
                    //mostExp = std::min_element(expSet.begin(), expSet.end())->second;
                }
                expSum -= person.previousExp;
                persons[name] = person;

            } else {
                peopleInTeam ++;
                person.inTeam = true;
                person.timeIn = time;
                mostExp = getMax(name, person.previousExp, mostExp, time);
                expSet.insert(std::make_pair(-1 * (person.previousExp - time), name));
                persons[name] = person;
                expSum += person.previousExp;
            }
        } else {
            persons[name] = {true, time, 0};
            peopleInTeam++;
            if (mostExp.empty()) {
                mostExp = name;
            } else {
                mostExp = getMax(name, 0, mostExp, time);
            }
            expSet.insert(std::make_pair(time, name));
        }

        auto mostExpPerson = persons.at(mostExp);
        long totalExp = mostExpPerson.previousExp + (time - mostExpPerson.timeIn);
        std::cout << mostExp << ' ' << expSum - 2 * totalExp << std::endl;
        prevTime = time;
    }
}
