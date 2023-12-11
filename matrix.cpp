//
// Created by zhest on 11.12.2023.
//
#include "vector"
#include <iostream>
#include "cmath"

std::ostream& operator << (std::ostream &os, std::vector<std::vector<double>> matr)
{
    for (int i = 0; i < matr.size(); i++){
        for (int j = 0; j < matr[i].size(); j++)
            std::cout << matr[i][j] << " ";
        std::cout << std::endl;
    }
    return os;
}
void TransponMtx(std::vector<std::vector<double>> &matr, std::vector<std::vector<double>> &tMatr, int n){//
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            tMatr[j][i] = matr[i][j];
}
void Get_matr(std::vector<std::vector<double>> &matr, int n, std::vector<std::vector<double>> &temp_matr, int indRow, int indCol)
{
    int ki = 0;
    for (int i = 0; i < n; i++){
        if(i != indRow){
            for (int j = 0, kj = 0; j < n; j++){
                if (j != indCol){
                    temp_matr[ki][kj] = matr[i][j];
                    kj++;
                }
            }
            ki++;
        }
    }
}
double Det(std::vector<std::vector<double>> &matr, int n)
{
    setlocale(LC_ALL, "Rus");
    double temp = 0;   //временная переменная для хранения определителя
    int k = 1;      //степень
    if(n < 1){
        std::cout<<"Не верный размер матрицы!!!" << std::endl;
        return 0;
    }
    else if (n == 1)
        temp = matr[0][0];
    else if (n == 2)
        temp = matr[0][0] * matr[1][1] - matr[1][0] * matr[0][1];
    else{
        for(int i = 0; i < n; i++){
            int m = n - 1;
            std::vector<std::vector<double>> temp_matr(m, std::vector<double> (n, 0));
            Get_matr(matr, n, temp_matr, 0, i);
            temp = temp + k * matr[0][i] * Det(temp_matr, m);
            k = -k;
        }
    }
    return temp;
}
std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> &matr){
    int n = matr.size();
    double det = Det(matr, n);
    std::vector<std::vector<double>> obr_matr(matr[0].size(), std::vector<double> (n, 0));
    std::vector<std::vector<double>> tobr_matr(matr[0].size(), std::vector<double> (n, 0));
    if(det){
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                int m = n - 1;
                std::vector<std::vector<double>> temp_matr(m, std::vector<double> (n, 0));
                Get_matr(matr, n, temp_matr, i, j);
                obr_matr[i][j] = pow(-1.0, i + j + 2) * Det(temp_matr, m) / det;
            }
        }
    }
    else std::cout << "Т.к. определитель матрицы = 0,\nто матрица вырожденная и обратной не имеет!!!" << std::endl;
    TransponMtx(obr_matr, tobr_matr, n);
    return tobr_matr;
}