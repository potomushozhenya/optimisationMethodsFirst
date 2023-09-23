#include <iostream>
#include "cmath"

double f(double &x){
    return std::exp(x*0.6)+std::exp(x*0.1)+std::exp(x*0.12)-2.8*std::sin(x);
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
    double a=0;
    double b=2;
    double eps=0.001;
    std::cout << goldenRatio(a,b,eps) << std::endl;
    return 0;
}
