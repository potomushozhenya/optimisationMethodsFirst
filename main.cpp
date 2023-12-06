#include <iostream>
#include "cmath"
#include "vector"
#include "iomanip"
/* check-list
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */
std::vector<double> operator+(std::vector<double> a, std::vector<double> b) {
    for (int i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}
std::vector<double> operator-(std::vector<double> a, std::vector<double> b) {
    for (int i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
    return a;
}
std::vector<double> operator*(double a, std::vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] *= a;
    }
    return b;
}
double f(double &x){
    return std::exp(x*0.6)+std::exp(x*0.1)+std::exp(x*0.12)-2.8*std::sin(x);
}
double fDer(double &x){
    return 0.6*std::exp(x*0.6)+0.1*std::exp(x*0.1)+0.12*std::exp(x*0.12)-2.8*std::cos(x);
}
double f2Der(double &x){
    return 0.36*std::exp(x*0.6)+0.01*std::exp(x*0.1)+0.0144*std::exp(x*0.12)+2.8*std::sin(x);
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
double passiveSearch(double (*f)(double&), double &leftBorder, double &rightBorder, double k){
    double res = INFINITY;
    double resPoint;
    for (int i = 0; i < k; ++i) {
        double point = leftBorder + i*((rightBorder-leftBorder)/k);
        double fMean = f(point);
        if (fMean < res){
            res = fMean;
            resPoint = point;
        }
    }
    return resPoint;
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
double fibonacci(double (*f)(double&), double &leftBorder, double &rightBorder, double eps){
    std::vector<int> fibSeq = {1, 1, 2};
    int n = 2;
    while (fibSeq[n] < ((rightBorder - leftBorder)/eps)){
        fibSeq.push_back(fibSeq[n]+fibSeq[n-1]);
        ++n;
    }
    fibSeq.push_back(fibSeq[n]+fibSeq[n-1]);
    ++n;
    fibSeq.push_back(fibSeq[n]+fibSeq[n-1]);
    --(--n);
    double a = leftBorder, b = rightBorder;
    for (int i = 0; i < n; ++i) {
        double c = a + (b - a)*((double)fibSeq[n+1-i]/(double)fibSeq[n+3-i]);
        double d = a + (b - a)*((double)fibSeq[n+2-i]/(double)fibSeq[n+3-i]);
        if (f(c) < f(d)) {
            b = d;
        } else {
            a = c;
        }
    }
    return (a+b)/2;
}
double tangentMethod(double (*f)(double&), double (*fDer)(double&), double& a, double& b, double eps){
    if (fDer(a) >= 0){
        return a;
    }
    if (fDer(b) <= 0){
        return b;
    }
    double cCurr;
    double k1 = fDer(a);
    double c1 = f(a) - k1*a;
    double k2 = fDer(b);
    double c2 = f(b) - k2*b;
    double kCurr;
    do {
        cCurr = (c2-c1)/(k1-k2);
        kCurr = fDer(cCurr);
        if (kCurr>0){
            k2 = kCurr;
            c2 = f(cCurr) - k2*cCurr;
        } else{
            k1 = kCurr;
            c1 = f(cCurr) -k1*cCurr;
        }
    } while (std::abs(fDer(cCurr))>= eps);
    return cCurr;
}
double newtonRaphson(double (*fDer)(double&), double (*f2Der)(double&), double& x0, double eps){
    double x = x0;
    double kCurr = fDer(x0);
    double k2Curr = f2Der(x0);
    while (std::abs(kCurr)>= eps){
        x = x - kCurr/k2Curr;
        kCurr = fDer(x);
        k2Curr = f2Der(x);
    }
    return x;
}
double der2Apr(double (*fDer)(double&), double &x, double xPrev){
    double res = (fDer(x)-fDer(xPrev))/(x-xPrev);
    return res;
}
double secantMethod(double (*fDer)(double&), double& x0, double eps){
    double x = x0;
    double xPrev = 0;
    double kCurr = fDer(x0);
    double k2Curr = der2Apr(fDer, x, xPrev);
    while (std::abs(kCurr)>= eps){
        xPrev = x;
        x = x - kCurr/k2Curr;
        kCurr = fDer(x);
        k2Curr = f2Der(x);
    }
    return x;
}
double coordinateDescent(double (*f)(double&), std::vector<double> &x0, double step){
    bool flag;
    std::vector<double> x = x0;
    while (true){
        flag = true;
        for (int i = 0; i < x0.size(); ++i) {


        }
        if (flag){
            break;
        }
    }

    return 0;
}
int main() {
    //Начальные границы для метода дихотомии a=0 и b=2
    //Золотое сечение, метод касательных, Ньютона-Рафсона
    //Для квадратичной функции используется следующий порядок ячейка-моном
    //0-x1^2    1-x2^2    2-x3^2    3-x1x2    4-x1x3    5-x2x3    6-x1    7-x2    8-x3    9-c
    double a=1;
    double b=2;
    double x0 = 0.5;
    double eps=0.0000001;
    std::cout << std::fixed;
    std::cout << std::setprecision(16);
    std::vector<double> myPoly = {3,4,5,2,-1,-2,1,0,-3,0};
    //std::cout << goldenRatio(a,b,eps) << std::endl;
    //std::cout << fibonacci(f, a, b, 0.00001)<< std::endl;
    //std::cout << passiveSearch(f, a, b, 100000) << std::endl;
    //std::cout << newtonRaphson(fDer, f2Der,a , 0.0000001) << std::endl;
    //std::cout << secantMethod(fDer, x0, 0.001) << std::endl;

    return 0;
}
