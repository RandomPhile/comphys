#include <iostream>
#include <cmath>

double log1(double x);
float log21(float x);

using namespace std;

int main() {
    float x0=1.f;
    double x1=1.l, err_rel, err_rel1;
    
    for(int i=0; i<50; i++){
        err_rel=fabs((log(1.l+x1)-log(1.f+x0))/log(1.l+x1));
        err_rel1=fabs((log1(1.l+x1)-log21(1.f+x0))/log1(1.l+x1));
        cout<<x0<<"\t"<<err_rel<<"\t"<<err_rel1<<endl;
        x0=x0/2.f;
        x1=x1/2.l;
    }
    
    return 0;
}

double log1(double x){
    if(1.l+x==1.l){
        return x;
    }
    else{
        double y;
        y=x*(log(1.l+x)/((1.l+x)-1.l));
        return y;
    }
}
float log21(float x){
    if(1.f+x==1.f){
        return x;
    }
    else{
        double y;
        y=x*(log(1.f+x)/((1.f+x)-1.f));
        return y;
    }
}
