#include <iostream>
#include <cmath>

using namespace std;

int main() {
    float a, b, c, x1, x2, x1c;
    double ad=1, bd, cd=1, x1d, x1dc, x2d, err_rel1, err_rel2;
    a=ad;
    bd=1;//per avere una inizializzazione, poi inutile
    b=bd;
    c=cd;
    
    
    
    for(double i=3; i<10000000;i=i*1.4l){
        b=i;
        bd=i;
        printf("%f\n", b);
        
        x1=(-b+sqrt(b*b-4*a*c))/(2*a);//guardo la prima soluzione perche instabile e valuto l'errore per b grande
        x1c=2.f*c/(-b-sqrt(b*b-4.f*a*c));
        x2=(-b-sqrt(b*b-4.f*a*c))/(2.f*a);
        
        x1d=-bd+sqrt(pow(bd,2)-4*ad*cd);//controllo la soluzione 1 con la versione stabile e controllo l'errore
        x1dc=2*cd/(-bd-sqrt(pow(bd,2)-4*ad*cd));
        x2d=(-bd-sqrt(pow(bd,2)-4*ad*cd))/(2*ad);
        
        
        err_rel1=fabs((x1d-x1)/x1d)/(2.f*a);
        err_rel2=fabs((x1dc-x1c)/x1dc)/(2*ad);
    
        cout<<i<<"\t"<<err_rel1<<"\t"<<err_rel2<<endl;
    }
    return 0;
}
