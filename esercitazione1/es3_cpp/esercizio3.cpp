#include <iostream>
#include <cmath>

using namespace std;

int main() {
    float a=1, b, c=1, x1, x2, x1c;
    double ad, bd, cd, x1d, x1dc, x2d, err_rel1, err_rel2;
    ad=a;//stai assegnando il valore di un float ad un double, sarebbe meglio fare viceversa se proprio
    b=0;//questo valore di b non viene usato
    bd=b;
    cd=c;
    
    
    
    for(int i=3; i<10000000;i=i*1.4){//i è un intero e gli assegni il valore di i*1.4
        b=i;//valore di un intero assegnato ad un float
        bd=i;//valore di un intero ad un double
        printf("%f\n", b);
        x1=-b+sqrt(b*b-4*a*c)/(2*a);//mancano parentesi prima di -b
        x1c=2*c/(-b-sqrt(b*b-4*a*c));
        x2=-b-sqrt(b*b-4*a*c)/(2*a);//la soluzione x2 non viene mai usata... invece è quella interessante perchè è stabile
        
        x1d=-bd+sqrt(pow(bd,2)-4*ad*cd);
        x1dc=2*cd/(-bd-sqrt(pow(bd,2)-4*ad*cd));
        x2=-bd-sqrt(pow(bd,2)-4*ad*cd);//qua immagino fosse x2d
        
        
        err_rel1=fabs((x1d-x1)/x1d)/(2*ad);
        err_rel2=fabs((x1dc-x1c)/x1dc)/(2*ad);
        //manca appunto la soluzione x2...
        cout<<i<<"\t"<<err_rel1<<"\t"<<err_rel2<<endl;
    }
    return 0;
}
