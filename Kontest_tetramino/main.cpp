#include <iostream>
#include <vector>

using namespace std;

size_t N = 8; // boardsize

void CreateBoard(vector<vector<char>>& board){
    for(size_t i = 0; i < N; ++ i){
        for(size_t j = 0; j < N; ++ j){
            cin >> board[i][j];
        }
    }
}

bool CheckLTetramino(vector<vector<char>>& board, size_t startrow, size_t startcol){
    if(board[startrow][startcol] == '*')
    {
        return false;
    }
    if(startrow + 1 >= N)
    {
        return false;
    }
    if(startrow + 2 >= N){
        return false;
    }
    if(startcol + 1 >= N){
        return false;
    }
    if(board[startrow + 1][startcol] == '*')
    {
        return false;
    }
    if(board[startrow + 2][startcol] == '*')
    {
        return false;
    }
    if(board[startrow + 1][startcol + 1] == '*')
    {
        return false;
    }
    return true;
}

bool CheckRTetramino(vector<vector<char>>& board, size_t startrow, size_t startcol){
    if(board[startrow][startcol] == '*')
    {
        return false;
    }
    if(startcol == 0){
        return false;
    }
    if(startrow + 1 >= N)
    {
        return false;
    }
    if(startrow + 2 >= N){
        return false;
    }
    if(board[startrow + 1][startcol] == '*')
    {
        return false;
    }
    if(board[startrow + 2][startcol] == '*')
    {
        return false;
    }
    if(board[startrow + 1][startcol - 1] == '*')
    {
        return false;
    }
    return true;
}

bool CheckTTetramino(vector<vector<char>>& board, size_t startrow, size_t startcol){
    if(board[startrow][startcol] == '*')
    {
        return false;
    }
    if(startrow == 0){
        return false;
    }
    if(startcol + 1 >= N)
    {
        return false;
    }
    if(startcol + 2 >= N){
        return false;
    }
    if(board[startrow][startcol + 1] == '*')
    {
        return false;
    }
    if(board[startrow][startcol + 2] == '*')
    {
        return false;
    }
    if(board[startrow - 1][startcol + 1] == '*')
    {
        return false;
    }
    return true;
}

bool CheckBTetramino(vector<vector<char>>& board, size_t startrow, size_t startcol){
    if(board[startrow][startcol] == '*')
    {
        return false;
    }
    if(startrow + 1 >= N){
        return false;
    }
    if(startcol + 1 >= N)
    {
        return false;
    }
    if(startcol + 2 >= N){
        return false;
    }
    if(board[startrow][startcol + 1] == '*')
    {
        return false;
    }
    if(board[startrow][startcol + 2] == '*')
    {
        return false;
    }
    if(board[startrow + 1][startcol + 1] == '*')
    {
        return false;
    }
    return true;
}

int CheckBoard(vector<vector<char>>& board){
    int counter = 0;
    for(size_t i = 0; i < N; ++ i){
        for(size_t j = 0; j < N; ++ j){
            if(CheckLTetramino(board, i, j)){
                ++ counter;
            }
            if(CheckRTetramino(board, i, j)){
                ++ counter;
            }
            if(CheckTTetramino(board, i, j)){
                ++ counter;
            }
            if(CheckBTetramino(board, i, j)){
                ++ counter;
            }
        }
    }
    return counter;
}

int main()
{
    /*
    vector<vector<char>> myboard(N, vector<char>(N));
    CreateBoard(myboard);
    int result = CheckBoard(myboard);
    cout << result;
    */
    cout << sizeof(size_t);
    return 0;
}
