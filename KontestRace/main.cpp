#include <iostream>
#include <vector>

using namespace std;

size_t CountOvertakes(vector<double>& speeds, double mintimer, size_t n, double t, double s) {
    size_t overtakes = 0;
    vector<vector<double>> dist(2, vector<double>(n + 1));
    for(size_t i = 1; i <= n; ++ i)
    {
        dist[1][i] = speeds[i] * mintimer;
        if(!(dist[1][i] < s)){
            dist[1][i] -= s;
        }
    }
    for(double timer = mintimer; timer < t; timer += mintimer){
        dist[0][1] = dist[1][1];
        dist[1][1] += speeds[1] * mintimer;
        for(size_t i = 2; i <= n; ++ i){
            dist[0][i] = dist[1][i];
            dist[1][i] += speeds[i] * mintimer;
            if(!(dist[0][1] > dist[0][i]) && (dist[1][1] > dist[1][i])){
                ++ overtakes;
            }
            if(!(dist[1][i] < s)){
                dist[1][i] -= s;
            }
        }

        if(!(dist[1][1] < s)){
                dist[1][1] -= s;
            }
    }
    return overtakes;
}

int main() {
    double t, s, mintimer;
    size_t n, num = 1;
    cin >> n >> t >> s;
    vector<double> speeds(n + 1);
    vector<double> newspeeds(n + 1);
    cin >> speeds[1];
    newspeeds[1] = speeds[1];
    mintimer = s / speeds[1];
    for (size_t i = 2; i <= n; ++ i) {
        cin >> speeds[i];
        if(speeds[1] > speeds[i]){
            ++ num;
            newspeeds[num] = speeds[i];
        }
    }
    if (num < n)
    {
        newspeeds.resize(num);
    }
    size_t overtakes = CountOvertakes(newspeeds, mintimer, num, t, s);
    cout << overtakes;
    return 0;
}
