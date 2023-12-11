//
// Created by zhest on 11.12.2023.
//
#include "vector"
#include <iostream>

#ifndef OPTIMISATION_MATRIX_H
#define OPTIMISATION_MATRIX_H
double Det(std::vector<std::vector<double>> &matr, int n);
void TransponMtx(std::vector<std::vector<double>> &matr, std::vector<std::vector<double>> &tMatr, int n);
void Get_matr(std::vector<std::vector<double>> &matr, int n, std::vector<std::vector<double>> &temp_matr, int indRow, int indCol);
std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> &matr);
std::ostream& operator << (std::ostream &os, std::vector<std::vector<double>> matr);


#endif //OPTIMISATION_MATRIX_H
