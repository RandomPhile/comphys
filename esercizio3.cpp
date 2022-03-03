#include <iostream>
#include <cmath>

using namespace std;

int main() {
    float a=1, b, c=1, x1, x2;
    double ad, bd, cd, x1d, x2d, err_rel1, err_rel2;
    ad=a;
    b=0;
    bd=b;
    cd=c;
    
    
    
    for(int i=3; i<100000000;i=i*1.4){
        b=i;
        bd=i;
        x1=-b+sqrt(b*b-4*a*c)/(2*a);
        x2=2*c/(-b-sqrt(b*b-4*a*c));
        
        x1d=-bd+sqrt(pow(bd,2)-4*ad*cd);
        x2d=2*cd/(-bd-sqrt(pow(bd,2)-4*ad*cd));
        
        
        err_rel1=fabs((x1d-x1)/x1d)/(2*ad);
        err_rel2=fabs((x2d-x2)/x2d)/(2*ad);
        
        cout<<i<<"\t"<<err_rel1<<"\t"<<err_rel2<<endl;
    }
    return 0;
}
