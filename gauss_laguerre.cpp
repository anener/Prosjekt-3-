#include "gauss_laguerre.h"
#include "lib.h"

//funfor aa regne ut Gauss-Laguerre
void gaussLaguerre(int n) {
    double pi = M_PI;

    //finner r verdiene og de vektede verdiene
    double *r = new double [n+1];
    double *wr = new double [n+1];
    double alpha = 2.0;
    gaulag(r, wr, n, alpha); //Laguerre

    //finner theta, pi og de vektede verdiene
    double *t = new double [n];
    double *wt = new double [n];
    gauleg(0.0, pi, t, wt, n); //Legendre

    double *p = new double [n];
    double *wp = new double [n];
    gauleg(0.0, 2.0*pi, p, wp, n); //Legendre

    double gaussint = 0.0;

    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++) {
            for(int k=0; k<n; k++) {
                for(int l=0; l<n; l++) {
                    for(int x=0; x<n; x++) {
                        for(int y=0; y<n; y++) {
                            gaussint += wr[i]*wr[j]*wt[k]*wt[l]*wp[x]*wp[y]*intFunction_laguerre(r[i], r[j], t[k], t[l], p[x], p[y]); //for-loop for aa summere opp alle vekt verdiene og funksjonen med de forskjellige koordinatene
                        }
                    }
                }
            }
        }
    }
    cout << "integrasjon = " << gaussint << "\n";
    delete [] r;
    delete [] wr;
    delete [] t;
    delete [] wt;
    delete [] p;
    delete [] wp;
}


//Laguerre funksjon for aa finne x og vektede verdier med Laguerre polynome
void gaulag(double *x, double *w, int n, double alpha) {
    double a;
    //int t;
    double p1, p2, p3, p, z, z1;
    p2 = 0.0; z = 0.0; p = 0.0;

    for(int i=1; i <= n; i++) {
        if(i == 1) {
            z = (1.0 + alpha)*(3.0 + 0.92*alpha)/(1.0 + 2.4*n + 1.8*alpha);
        }
        else if(i == 2) {
            z += (15.0 + 6.25*alpha)/(1.0 + 0.9*alpha + 2.5*n);
        }
        else {
            a = i-2;
            z += ((1.0 + 2.55*a)/(1.9*a) + 1.26*a*alpha/(1.0 + 3.5*a))*(z - x[i-2])/(1.0 + 0.3*alpha);
        }

        for(int t=1; t<= MAXIT; t++) {
            p1 = 1.0;
            p2 = 0.0;
            for (int j=1; j<= n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ( (2*j - 1 + alpha - z)*p2 - (j - 1 + alpha)*p3 )/j;
            }
            p = (n*p1 - (n+alpha)*p2)/z;
            z1 = z;
            z = z1 - p1/p;
            if(fabs(z - z1) <= EPS) break;
        }
        //if(t > MAXIT) cout << "too many iterations in gaulag" << endl;
        double N = static_cast<double>(n);
        x[i] = z;
        w[i] = -exp(gammln(alpha+n) - gammln(N))/(p*n*p2);
    }
}

//funksjonen for aa finne gamma funksjonen
double gammln(double x) {
    double v, w, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
    w = v = x;
    tmp = v + 5.5;
    tmp -= (v + 0.5)*log(tmp);
    ser = 1.000000000190015;

    for(int j=0; j <= 5; j++) {
        ser += cof[j]/++w;
    }

    return -tmp + log(2.5066282746310005*ser/v);
}

//funksjonene som vi onsker aa integrere
double intFunction_laguerre(double r1, double r2, double t1, double t2, double p1, double p2) {
    //double alpha = 2.0;
    double cos_b = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
    double r12 = (r1*r1) + (r2*r2) - (2.0*r1*r2*cos_b);
    //double dd = (r1*r1)*(r2*r2)*sin(t1)*sin(t2);
    double f;
    if(r12 > 1e-10) {
        //f = (exp(-2.0*alpha*(r1+r2))*dd)/sqrt(r12);
        f = (sin(t1)*sin(t2))/(1024.0*sqrt(r12));
    }
    else {
        f = 0.0;
    }
    return f;
}

#undef EPS
#undef MAXIT
