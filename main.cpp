#include <iostream>
#include "cmath"

double f(double &x){
    return std::exp(x*0.6)+std::exp(x*0.1)+std::exp(x*0.12)-2.8*std::sin(x);
}
double passiveSearch(double &x){

}
double dichotomy(double &leftBorder, double &rightBorder, double &eps){
    double a=leftBorder;
    double b=rightBorder;
    do {
        double c=(a+b)/2-(eps/2);
        double d=(a+b)/2+(eps/2);
        if (f(c)<=f(d)){
            b=d;
        } else{
            a=c;
        }
    } while (((b-a)/2)>=eps);
    return (a+b)/2;
}
int main() {
    //Начальные границы для метода дихотомии 0 и 2
    double a=0;
    double b=2;
    double eps=0.001;
    std::cout << dichotomy(a,b,eps) << std::endl;
    return 0;
}
