#include <iostream>
#include <cmath>

using namespace std;

int main() {
    float a=1, b, c=1, x1, x2, x1c;//definisco arbitrariamente a e c come uguali a 1 cosi da poter vedere eventualmente un grafico in funzione di b
    double ad, bd, cd, x1d, x1dc, x2d, err_rel1, err_rel2;
    ad=a;
    b=0;
    bd=b;
    cd=c;
    
    
    
    for(int i=3; i<100000000;i=i*1.4){//trovo i vari valori in funzione di b e li stampo a video
        b=i;
        bd=i;
        x1=-b+sqrt(b*b-4*a*c)/(2*a);
        x1c=2*c/(-b-sqrt(b*b-4*a*c)); //correzione per stabilità
        x2=-b-sqrt(b*b-4*a*c)/(2*a);
        
        x1d=-bd+sqrt(pow(bd,2)-4*ad*cd);
        x1dc=2*cd/(-bd-sqrt(pow(bd,2)-4*ad*cd)); //correzione per stabilità in double
        x2=-bd-sqrt(pow(bd,2)-4*ad*cd);
        
        
        err_rel1=fabs((x1d-x1)/x1d)/(2*ad);
        err_rel2=fabs((x1dc-x1c)/x1dc)/(2*ad);
        
        cout<<i<<"\t"<<err_rel1<<"\t"<<err_rel2<<endl;
    }
    return 0;
}
