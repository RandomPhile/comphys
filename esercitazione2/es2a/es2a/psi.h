//
//  psi.h
//  es2a
//
//  Created by Mattia Lupi on 20/03/22.
//

#ifndef psi_h
#define psi_h
#define _USE_MATH_DEFINES

using namespace std;

extern double k;
extern double q;
extern double R;

double psi2(double r, double A){//definisco il modulo quadro di psi che essendo reale corrisponde al quadrato
    if(r<=R){
        double y=sin(k*r)/(r*sqrt(4*M_PI));
        return y*y;
    }
    else{
        double y=sin(k*R)*exp(q*(R-r))/(r*sqrt(4*M_PI));
        return y*y;
    }
}
double rpsi2(double r, double A){//definisco il modulo quadro di psi che essendo reale corrisponde al quadrato
    if(r<=R){
        double y=sin(k*r)/(r*sqrt(4*M_PI));
        return r*y*y;
    }
    else{
        double y=sin(k*R)*exp(q*(R-r))/(r*sqrt(4*M_PI));
        return r*y*y;
    }
}
double r2psi2(double r, double A){//come sopra ma moltiplico per r^2 per avere l'integrale di mio interesse
    if(r<=R){
        double y=r*sin(k*r)/(r*sqrt(4*M_PI));
        return y*y;
    }
    else{
        double y=r*sin(k*R)*exp(q*(R-r))/(r*sqrt(4*M_PI));
        return y*y;
    }
}

#endif /* psi_h */
