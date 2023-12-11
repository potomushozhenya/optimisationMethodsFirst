#include <iostream>
#include "cmath"
#include "vector"
#include "iomanip"
#include "matrix.h"

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
double operator*(std::vector<double> a, std::vector<double> b) {
    double res = 0;
    unsigned int vecSize = a.size();
    if (vecSize != b.size()){
        throw "Invalid vector arguments";
    }
    for (int i = 0; i < vecSize; i++) {
        res += a[i]*b[i];
    }
    return res;
}
std::vector<double> operator*(std::vector<std::vector<double>> a, std::vector<double> b) {
    unsigned int vecSize = b.size();
    std::vector<double> res(b.size(), 0);
    if (vecSize != a.size()){
        throw "Invalid vector arguments";
    }
    for (int i = 0; i < vecSize; i++) {
        double curr = 0;
        for (int j = 0; j < a[0].size(); ++j) {
            curr += a[i][j]*b[j];
        }
        res[i] = curr;
    }
    return res;
}
std::vector<std::vector<double>> operator+(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b){
    unsigned int vecSize = b.size();
    std::vector<std::vector<double>> res(vecSize, std::vector<double> (vecSize, 0));
    for (int i = 0; i < vecSize; ++i) {
        for (int j = 0; j < vecSize; ++j) {
            res[i][j] = a[i][j]+b[i][j];
        }
    }
    return res;
};
std::vector<std::vector<double>> operator*(double a, std::vector<std::vector<double>> b){
    unsigned int vecSize = b.size();
    std::vector<std::vector<double>> res(vecSize, std::vector<double> (vecSize, 0));
    for (int i = 0; i < vecSize; ++i) {
        for (int j = 0; j < vecSize; ++j) {
            res[i][j] = a*b[i][j];
        }
    }
    return res;
}
double norm(std::vector<double> x){
    double res = 0;
    for (int i = 0; i < x.size(); ++i) {
        res += pow(x[i],2);
    }
    return sqrt(res);
}
std::vector<double> unitVec(std::vector<double> &x){
    std::vector<double> res;
    unsigned int vecSize = x.size();
    res.resize(vecSize);
    for (int i = 0; i < vecSize; ++i) {
        res[i] = 1;
    }
    return res;
}
std::vector<double> baseVec(int &i, double &h, std::vector<double> &x){
    std::vector<double> res;
    unsigned int vecSize = x.size();
    res.resize(vecSize);
    for (int j = 0; j < vecSize; ++j) {
        if (j!=i){
            res[j] = 0;
        } else{
            res[i] = h;
        }
    }
    return res;
}
std::vector<std::vector<double>> E(int n){
    std::vector<std::vector<double>> res(n, std::vector<double> (n, 0));
    for (int i = 0; i < n; ++i) {
        res[i][i] = 1;
    }
    return res;
}
double f(double &x){
    return std::exp(x*0.6)+std::exp(x*0.1)+std::exp(x*0.12)-2.8*std::sin(x);
}
double f1761(double &x){
    return -pow(x, 3)+ 3*(1+x)*(log(1+x)-1);
}
double fMul(std::vector<double> &x){
    double res;
    res = pow(x[0]-1, 2);
    res += pow(x[1]+1, 2);
    return res;
}
double f17138(std::vector<double> &x){
    double res = 0;
    res += 3*x[0]*x[0];
    res += 4*x[1]*x[1];
    res += 5*x[2]*x[2];
    res += 2*x[0]*x[1];
    res += -x[0]*x[2];
    res += -2*x[1]*x[2];
    res += x[0];
    res += -3*x[2];
    return res;
}
std::vector<double> fDer17138(std::vector<double> &x){
    std::vector<double> res(3,0);
    res[0] += 6*x[0]+2*x[1]-x[2]+1;
    res[1] += 8*x[1]+2*x[0]-2*x[2];
    res[2] += 10*x[2]-x[0]-2*x[1]-3;
    return res;
}
std::vector<double> fDerMul(std::vector<double> &x){
    std::vector<double> res;
    res.resize(2);
    res[0] = 2*(x[0]-1);
    res[1] = 2*(x[1]+1);
    return res;
}
std::vector<std::vector<double>> f2DerMul(std::vector<double> &x){
    std::vector<std::vector<double>> res(2, std::vector<double> (2, 0));
    res[0][0] = 2;
    res[0][1] = 0;
    res[1][0] = 0;
    res[1][1] = 2;
    return res;
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
std::ostream& operator << (std::ostream &os, std::vector<double> x)
{
    for (int i = 0; i < x.size(); ++i) {
        std::cout<<x[i]<<" ";
    }
    std::cout<<std::endl;
    return os;
}
/*
double** derivativeQuadraticFunc(double *x){
    auto **res = new double*[3];
    for (int i = 0; i < 3; ++i) {
        res[i] = new double[4];
    }
    res[0][0]=2*x[0];
    res[0][1]=x[3];
    res[0][2]=x[4];
    res[0][3]=x[6];
    res[1][0]=x[3];
    res[1][1]=2*x[1];
    res[1][2]=x[5];
    res[1][3]=x[7];
    res[2][0]=x[4];
    res[2][1]=x[5];
    res[2][2]=2*x[2];
    res[2][3]=x[8];
    return res;
}
*/
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
double goldenRatio(double (*f)(double&), double &leftBorder, double &rightBorder, double &eps){
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
double goldenRatioGradientDuplicate(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> &xk, double &eps){
    double a=0.001;
    double b=10;
    double c=((3- pow(5,0.5))/2)*(b-a)+a;
    double d=((pow(5,0.5)-1)/2)*(b-a)+a;
    do {
        std::vector<double> point1 = xk - c*fDer(xk);
        double fMean1 = f(point1);
        std::vector<double> point2 = xk - d*fDer(xk);
        double fMean2 = f(point2);
        if (fMean1<=fMean2){
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
double goldenRatioGradientYDuplicate(double (*f)(std::vector<double>&), std::vector<double> &yk, std::vector<double> &xk, double &eps){
    double a=0.0001;
    double b=10;
    double c=((3- pow(5,0.5))/2)*(b-a)+a;
    double d=((pow(5,0.5)-1)/2)*(b-a)+a;
    do {
        std::vector<double> point1 = xk - c*(yk - xk);
        double fMean1 = f(point1);
        std::vector<double> point2 = xk - d*(yk - xk);
        double fMean2 = f(point2);
        if (fMean1<=fMean2){
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
double goldenRatioFletcherDuplicate(double (*f)(std::vector<double>&), std::vector<double> &xk, std::vector<double> dk, double &eps){
    double a=0.001;
    double b=10;
    double c=((3- pow(5,0.5))/2)*(b-a)+a;
    double d=((pow(5,0.5)-1)/2)*(b-a)+a;
    do {
        std::vector<double> point1 = xk + c*dk;
        double fMean1 = f(point1);
        std::vector<double> point2 = xk + d*dk;
        double fMean2 = f(point2);
        if (fMean1<=fMean2){
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
double goldenRatioNewtonDuplicate(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<std::vector<double>> &matrix, std::vector<double> &xk, double &eps){
    double a=-10;
    double b=10;
    double c=((3- pow(5,0.5))/2)*(b-a)+a;
    double d=((pow(5,0.5)-1)/2)*(b-a)+a;
    std::vector<double> vec = inverse(matrix)*fDer(xk);
    do {
        std::vector<double> point1 = xk + c*(vec);
        double fMean1 = f(point1);
        std::vector<double> point2 = xk + d*(vec);
        double fMean2 = f(point2);
        if (fMean1<=fMean2){
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
double goldenRatioQuasiNewtonDuplicate(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<std::vector<double>> &matrix, std::vector<double> &xk, double &eps){
    double a=-10;
    double b=10;
    double c=((3- pow(5,0.5))/2)*(b-a)+a;
    double d=((pow(5,0.5)-1)/2)*(b-a)+a;
    std::vector<double> vec = matrix*fDer(xk);
    do {
        std::vector<double> point1 = xk - c*(vec);
        double fMean1 = f(point1);
        std::vector<double> point2 = xk - d*(vec);
        double fMean2 = f(point2);
        if (fMean1<=fMean2){
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
double dichotomy(double (*f)(double&), double &leftBorder, double &rightBorder, double &eps){
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
std::vector<double> coordinateDescent(double (*f)(std::vector<double>&), std::vector<double> &x0, double step, double localTol){
    bool flag;
    bool flagLocal;
    std::vector<double> x = x0;
    while (true){
        flag = true;
        for (int i = 0; i < x.size(); ++i) {
            flagLocal = true;
            while (flagLocal) {
                std::vector<double> curr = baseVec(i, step, x);
                std::vector<double> aPoint = x + curr;
                double a = f(aPoint);
                std::vector<double> bPoint = x - curr;
                double b = f(bPoint);
                if (std::abs(a-b) <= localTol){
                    flagLocal = false;
                } else{
                    flag = false;
                    if (a < b){
                        x = aPoint;
                    } else{
                        x = bPoint;
                    }
                }
            }
        }
        if (flag){
            break;
        }
    }
    return x;
}
std::vector<double> gradientStepSplitting(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> x0, double initialStep, double multiplier, double eps){
    std::vector<double> xk = x0;
    double step;
    while (norm(fDer(xk)) >= eps){
        step = initialStep;
        std::vector<double> xPoint = xk - step*fDer(xk);
        double fMean = f(xPoint);
        while ((fMean-f(xk)) > (multiplier*step* pow(norm(fDer(xk)),2))){
            step *= multiplier;
            xPoint = xk - step*fDer(xk);
            fMean = f(xPoint);
        }
        xk = xPoint;
    }
    return xk;
}
std::vector<double> gradientConstStep(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> x0, double initialStep, double eps){
    std::vector<double> xk = x0;
    while (norm(fDer(xk)) >= eps){
        xk = xk - initialStep*fDer(xk);
    }
    return xk;
}
std::vector<double> gradientFixedStep(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> x0, double initialStep, double eps){
    std::vector<double> xk = x0;
    double i = 0;
    while (norm(fDer(xk)) >= eps){
        i += 1;
        xk = xk - (1/i)*fDer(xk);
    }
    return xk;
}
double minimizationGradientArgument(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> xk, double &eps){
    return goldenRatioGradientDuplicate(f, fDer, xk, eps);
}
std::vector<double> gradientSteepestDescent(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> &x0,  double eps){
    std::vector<double> xk = x0;
    double step;
    while (norm(fDer(xk)) >= eps){
        xk = xk - minimizationGradientArgument(f, fDer, xk, eps)*fDer(xk);
    }
    return xk;
}
double minimizationGradientY(double (*f)(std::vector<double>&), std::vector<double> &yk, std::vector<double> &xk, double &eps){
    return goldenRatioGradientYDuplicate(f, yk, xk, eps);
}
std::vector<double> gradientPthOrder(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> &x0,  double eps){
    std::vector<double> x = x0;
    while (norm(fDer(x)) >= eps){
        std::vector<double> xk = x;
        for (int i = 0; i < xk.size(); ++i) {
            x = x - minimizationGradientArgument(f, fDer, x, eps)*fDer(x);
        }
        x = x + minimizationGradientY(f, x, xk, eps)*(x - xk);
    }
    return x;
}
std::vector<double> gullyMethod(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> &x0,  double eps, double diff){
    //diff - параметр, который задает разницу между значениям x для овражного метода
    std::vector<double> x = x0;
    std::vector<double> _x = x0 + diff*unitVec(x);
    while (norm(fDer(x)) >= eps){
        for (int i = 0; i < x.size(); ++i) {
            x = x - minimizationGradientArgument(f, fDer, x, eps)*fDer(x);
            _x = _x - minimizationGradientArgument(f, fDer, _x, eps)*fDer(_x);
        }
        x = x + minimizationGradientY(f, _x, x, eps)*(_x - x);
    }
    return x;
}
std::vector<double> conjugateDirectionsMethod(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> &x0,  double eps, bool isFletcherReeves){
    std::vector<double> d = -1*fDer(x0);
    std::vector<double> x = x0;
    std::vector<double> grad = fDer(x);
    unsigned int n = x0.size();
    unsigned int k = 0;
    while (norm(grad) >= eps){
        double a = goldenRatioFletcherDuplicate(f, x, d, eps);
        std::vector<double> xPrev = x;
        std::vector<double> gradPrev = grad;
        x = x + a*d;
        grad = fDer(x);
        if ((k+1) % n == 0){
            d = -1*grad;
        } else{
            double b;
            if (isFletcherReeves) {
                b = (pow(norm(grad), 2)) / (pow(norm(gradPrev), 2));
            } else{
                b = (grad*(grad-gradPrev)) / (pow(norm(gradPrev), 2));
            }
            d = -1*grad+b*d;
        }
        ++k;
    }
    return x;
}
std::vector<double> newtonMethod(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<std::vector<double>>(*f2Der)(std::vector<double>&), std::vector<double> &x0,  double eps){
    std::vector<double> xk = x0;
    while (norm(fDer(xk)) >= eps){
        std::vector<std::vector<double>> currMatrix = f2Der(xk);
        xk = xk + goldenRatioNewtonDuplicate(f, fDer, currMatrix, xk, eps)*(inverse(currMatrix)*fDer(xk));
    }
    return xk;
}
std::vector<double> quasiNewtonianMethod(double (*f)(std::vector<double>&), std::vector<double> (*fDer)(std::vector<double>&), std::vector<double> &x0,  double eps){
    std::vector<double> xk = x0;
    unsigned n = x0.size();
    std::vector<std::vector<double>> H = E(n);
    std::vector<double> grad = fDer(xk);
    unsigned int k = 0;
    while (norm(grad) >= eps){
        double a = goldenRatioQuasiNewtonDuplicate(f, fDer, H, xk, eps);
        std::vector<double> xPrev = xk;
        std::vector<double> gradPrev = grad;
        xk = xk - a*(H*fDer(xk));
        grad = fDer(xk);
        if ((k+1) % n == 0){
            H = E(n);
        } else{
            std::vector<double> delta = xk- xPrev;
            std::vector<double> gamma = grad - gradPrev;
            H = H + (1/((delta - H*gamma)*gamma))*((delta - H*gamma)*(delta - H*gamma))*E(n);
        }
        ++k;
    }
    return xk;
}
int main() {
    //Начальные границы для метода дихотомии a=0 и b=2
    //Золотое сечение, метод касательных, Ньютона-Рафсона
    //Для квадратичной функции используется следующий порядок ячейка-моном
    //0-x1^2    1-x2^2    2-x3^2    3-x1x2    4-x1x3    5-x2x3    6-x1    7-x2    8-x3    9-c
    double a=1;
    double b=2;
    double a1 = -0.5;
    double b1 = 0.5;
    double x0 = 0.5;
    std::vector<double> xVec = {-13, 100};
    std::vector<double> xVec1 = {0,0,0};
    double eps=0.000001;
    std::cout << std::fixed;
    std::cout << std::setprecision(16);
    //std::vector<double> myPoly = {3,4,5,2,-1,-2,1,0,-3,0};
    //std::cout << goldenRatio(f,a,b,eps) << std::endl;
    //std::cout << goldenRatio(f1761, a1, b1, eps)<< std::endl;
    //std::cout << fibonacci(f, a, b, 0.00001)<< std::endl;
    //std::cout << passiveSearch(f, a, b, 100000) << std::endl;
    //std::cout << newtonRaphson(fDer, f2Der,a , 0.0000001) << std::endl;
    //std::cout << secantMethod(fDer, x0, 0.001) << std::endl;
    //std::cout << coordinateDescent(fMul, xVec, 0.05, 0.00000001);
    //std::cout << gradientStepSplitting(fMul, fDerMul, xVec, 0.5, 0.5, 0.00001);
    std::cout << gradientConstStep(fMul, fDerMul, xVec, 0.5, 0.00001);
    //std::cout << gradientFixedStep(fMul, fDerMul, xVec, 0.5, 0.00001);
    //std::cout << gradientSteepestDescent(fMul, fDerMul, xVec, 0.0000001);
    std::cout << gradientSteepestDescent(f17138, fDer17138, xVec1, 0.000001);
    //std::cout << gradientPthOrder(fMul, fDerMul, xVec, 0.000001);
    //std::cout << gullyMethod(fMul, fDerMul, xVec, 0.00001, 0.1);
    //std::cout << conjugateDirectionsMethod(fMul, fDerMul, xVec, 0.000001, true);
    //std::cout << conjugateDirectionsMethod(fMul, fDerMul, xVec, 0.000001, false);
    //std::cout << newtonMethod(fMul, fDerMul, f2DerMul, xVec, 0.00001);
    //std::cout << quasiNewtonianMethod(fMul, fDerMul, xVec, 0.00001);
    //std::cout << quasiNewtonianMethod(f17138, fDer17138, xVec1, 0.000001);


    return 0;
}
