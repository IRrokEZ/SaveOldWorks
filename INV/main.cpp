#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

// ������� ��� �������� �������� � ������������
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

// ������� ��� ������������ ��������� �� �������� ���� ��������
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

    // ��������� ��������� ���� ��������
    int i = rand() % n + 1;
    int j = rand() % n + 1;
    while (i == j) {
        j = rand() % n + 1;
    }
    std::pair<int, int> randomPair = std::make_pair(std::min(i, j), std::max(i, j));

    // ������������ ��������� �� ��������� ���� ��������
    permutePair(permutation, randomPair);

    // ������� ��������
    int inversions_after_permutation = countInversions(permutation);
    int original_inversions = countInversions(original_permutation);

    // ����� �������� ���������� �������� � ������� ������������ �����
    std::cout << "Inversions after permuting a random pair: " << inversions_after_permutation << "/" << original_inversions << std::endl;

    return 0;
}
