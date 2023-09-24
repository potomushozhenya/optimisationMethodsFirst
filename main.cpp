#include <iostream>
#include "cmath"

double f(double &x){
    return std::exp(x*0.6)+std::exp(x*0.1)+std::exp(x*0.12)-2.8*std::sin(x);
}
double polyMean(double *x, double *poly){
    return poly[0]*pow(x[0],2)+poly[1]*pow(x[1],2)+poly[2]*pow(x[2],2)+poly[3]*x[0]*x[1]+poly[4]*x[0]*x[2]+poly[5]*x[1]*x[2]+poly[6]*x[0]+poly[7]*x[1]+poly[8]*x[2]+poly[9];
}
double** derivativeQuadraticFunc(double *coef){
    auto **resultCoef = new double*[3];
    for (int i = 0; i < 3; ++i) {
        resultCoef[i] = new double[4];
    }
    resultCoef[0][0]=2*coef[0];
    resultCoef[0][1]=coef[3];
    resultCoef[0][2]=coef[4];
    resultCoef[0][3]=coef[6];
    resultCoef[1][0]=coef[3];
    resultCoef[1][1]=2*coef[1];
    resultCoef[1][2]=coef[5];
    resultCoef[1][3]=coef[7];
    resultCoef[2][0]=coef[4];
    resultCoef[2][1]=coef[5];
    resultCoef[2][2]=2*coef[2];
    resultCoef[2][3]=coef[8];
    return resultCoef;
}
double goldenRatio(double &leftBorder, double &rightBorder, double &eps){
    double a=leftBorder;
    double b=rightBorder;
    double c=((3- pow(5,0.5))/2)*(b-a)+a;
    double d=((pow(5,0.5)-1)/2)*(b-a)+a;
    do {
        if (f(c)<=f(d)){
            b=d;
            d=c;
            c=((3- pow(5,0.5))/2)*(b-a)+a;
        } else{
            a=c;
            c=d;
            d=((pow(5,0.5)-1)/2)*(b-a)+a;
        }
    } while (((b-a)/2)>=eps);
    return (a+b)/2;
}
double dichotomy(double &leftBorder, double &rightBorder, double &eps){
    double a=leftBorder;
    double b=rightBorder;
    double c,d;
    do {
        c=(a+b)/2-(eps/2);
        d=(a+b)/2+(eps/2);
        if (f(c)<=f(d)){
            b=d;
        } else{
            a=c;
        }
    } while (((b-a)/2)>=eps);
    return (a+b)/2;
}
int main() {
    //Начальные границы для метода дихотомии a=0 и b=2
    //Золотое сечение, метод касательных, Ньютона-Рафсона
    //Для квадратичной функции используется следующий порядок ячейка-моном
    //0-x1^2    1-x2^2    2-x3^2    3-x1x2    4-x1x3    5-x2x3    6-x1    7-x2    8-x3    9-c
    double a=0;
    double b=2;
    double eps=0.001;
    double myPoly[10] = {3,4,5,2,-1,-2,1,0,-3,0};
    std::cout << goldenRatio(a,b,eps) << std::endl;
    return 0;
}
