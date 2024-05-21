#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

// Функция для подсчета инверсий в перестановке
int countInversions(const std::vector<int>& arr) {
    int n = arr.size();
    int inversions = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (arr[i] > arr[j]) {
                inversions += 1;
            }
        }
    }
    return inversions;
}

// Функция для перестановки элементов по заданной паре индексов
void permutePair(std::vector<int>& arr, const std::pair<int, int>& pair) {
    std::swap(arr[pair.first - 1], arr[pair.second - 1]);
}

int main() {
    int n;
    std::cout << "Enter the size of the permutation (n): ";
    std::cin >> n;

    std::vector<int> permutation(n);
    std::cout << "Enter the permutation: ";
    for (int i = 0; i < n; ++i) {
        std::cin >> permutation[i];
    }

    std::vector<int> original_permutation = permutation;

    // Генерация случайной пары индексов
    int i = rand() % n + 1;
    int j = rand() % n + 1;
    while (i == j) {
        j = rand() % n + 1;
    }
    std::pair<int, int> randomPair = std::make_pair(std::min(i, j), std::max(i, j));

    // Перестановка элементов по случайной паре индексов
    permutePair(permutation, randomPair);

    // Подсчет инверсий
    int inversions_after_permutation = countInversions(permutation);
    int original_inversions = countInversions(original_permutation);

    // Вывод среднего количества инверсий в формате обыкновенной дроби
    std::cout << "Inversions after permuting a random pair: " << inversions_after_permutation << "/" << original_inversions << std::endl;

    return 0;
}
