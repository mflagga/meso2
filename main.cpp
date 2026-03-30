#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;
typedef complex<double> cmp;
const cmp iu(0,1);
double H0(double x){return 1.0;}
double H1(double x){return 2.0*x;}
double H2(double x){return 4.0*x*x-2.0;}
double H3(double x){return 8.0*x*x*x-12.0*x;}

void psi_init(cmp **psi, int nx, double A, double sigma, double *x, double xc, double p0, int nt){
    for (int i=0;i<=nx;i++){
        for (int n=0;n<=nt;n++){
            psi[i][n]=0.0;
        }
    }
    for (int i=0;i<=nx;i++){
        psi[i][0] = A*exp(-sigma*pow(x[i]-xc,2))*exp(iu*p0*x[i]);
    }
    psi[0][0]=0.0;
    psi[nx][0]=0.0;
}

void psi_eig(cmp **psi, int nx, int nt, double *x){
    for (int i=0;i<=nx;i++){
        for (int n=0;n<=nt;n++){
            psi[i][n]=0.0;
        }
    }
    int n=0;
    for (int i=0;i<=nx;i++){
        psi[i][0] = (1.0/(sqrt(pow(2,n)*tgamma(n+1))))*pow(M_PI,-0.25)*exp(-0.5*x[i]*x[i])*H0(x[i]);
    }
    psi[0][0]=0.0;
    psi[nx][0]=0.0;
}

void init_barrier(cmp *V, double *x, int nx, double Vmax){
    for (int i=0;i<=nx;i++){
        if (x[i]<= 0.2 && x[i]>=0.0) V[i]=Vmax;
        else V[i]=0.0;
    }
}

void update_V(cmp *V, int nx, double *x, double t, double E0, double omega_ext) {
    double k = 40.0;
    for (int i = 0; i <= nx; i++) {
        V[i] = 0.5 * k * x[i] * x[i] - E0 * x[i] * sin(omega_ext * t);
    }
}

void thomas(cmp *nw, cmp *prawa, int n, cmp a, cmp *diag){
    cmp *cprim = new cmp[n-1];
    cmp *dprim = new cmp[n];

    cprim[0]=a/diag[0];
    for (int i=1;i<n-1;i++){
        cprim[i] = a/(diag[i]-a*cprim[i-1]);
    }
    dprim[0]=prawa[0]/diag[0];
    for (int i=1;i<n;i++){
        dprim[i]=(prawa[i]-a*dprim[i-1])/(diag[i]-a*cprim[i-1]);
    }
    nw[n-1]=dprim[n-1];
    for (int i=n-2;i>=0;i--){
        nw[i]=dprim[i]-cprim[i]*nw[i+1];
    }

    delete [] cprim;
    delete [] dprim;
}

void fillPrawa(cmp *prawa, int N, double dt, double theta, cmp **psi, int n, double dx, cmp *V){
    for (int i=0;i<N;i++){
        prawa[i] = iu*dt*(1.0-theta)*psi[i][n]+(2.0*dx*dx-2.0*iu*dt*(1-theta)-2.0*iu*dt*dx*dx*V[i+1]*(1.0-theta))*psi[i+1][n]+iu*dt*(1.0-theta)*psi[i+2][n];
    }
}

int main(){
    const double theta=0.5;
    const double xmin=-5.0;
    const double xmax=5.0;
    const double tmax=2.0;
    const int nx=350;
    const int nt=2000;
    const double dx=(xmax-xmin)/nx;
    const double dt=tmax/nt;
    const double xc=0.0;
    const double p0=0.0;
    const double sigma=15.0;
    const double A=sqrt(sigma/M_PI);
    const double Vmax=50.0;
    const bool bar = false;
    const int fps=20;
    const int co_ktora=8;
    const double E0 = 10.0;
    const double omega_ext = sqrt(40.0);
    // alokacja
    cmp **psi = new cmp*[nx+1];
    for (int i=0;i<=nx;i++){
        psi[i] = new cmp[nt+1];
    }
    cmp *V = new cmp[nx+1];
    double *x = new double[nx+1];
    double *t = new double[nt+1];
    // inicjalizacja
    for (int i=0;i<=nx;i++) x[i] = xmin+i*dx;
    if (bar) init_barrier(V,x,nx,Vmax);
    else {for (int i=0;i<=nx;i++) V[i] = 0.0;}
    for (int n=0;n<=nt;n++) t[n] = n*dt;
    update_V(V,nx,x,t[0],E0,omega_ext);
    psi_init(psi,nx,A,sigma,x,xc,p0,nt);
    // psi_eig(psi,nx,nt,x);
    // zmienne do pętli
    cmp *prawa = new cmp[nx-1];
    cmp *nw = new cmp[nx-1];
    fillPrawa(prawa,nx-1,dt,theta,psi,0,dx,V);
    cmp *diag = new cmp[nx-1];
    for (int i=0;i<nx-1;i++){
        diag[i] = 2.0*dx*dx+2.0*iu*dt*theta+2.0*iu*dt*dx*dx*V[i+1]*theta;
    }
    double expect_x=0.0;
    double expect_xx=0.0;
    ofstream Exx("Exx.dat");
    // pętla 
    for (int n=1;n<=nt;n++){
        update_V(V,nx,x,t[n],E0,omega_ext);
        for (int i=0;i<nx-1;i++) diag[i] = 2.0*dx*dx + 2.0*iu*dt*theta + 2.0*iu*dt*dx*dx*V[i+1]*theta;
        thomas(nw,prawa,nx-1,-iu*dt*theta,diag);
        for (int i=1;i<nx;i++) psi[i][n] = nw[i-1];
        for (int i=0;i<=nx;i++){
            expect_x += norm(psi[i][n])*x[i]*dx;
            expect_xx += norm(psi[i][n])*x[i]*x[i]*dx;
        }
        Exx<<t[n]<<'\t'<<expect_x<<'\t'<<sqrt(abs(expect_xx-expect_x*expect_x))<<'\n';
        expect_x=0.0;
        expect_xx=0.0;
        fillPrawa(prawa,nx-1,dt,theta,psi,n,dx,V);
    }
    // zapis 
    ofstream file("psi.dat");
    for (int n=0;n<=nt;n++){
        if (n%co_ktora==0) {for (int i=0;i<=nx;i++) file<<t[n]<<'\t'<<x[i]<<'\t'<<real(psi[i][n])<<'\t'<<imag(psi[i][n])<<'\t'<<norm(psi[i][n])<<'\n';}
    }
    file.close();
    // pliki pod przekaz danych
    ofstream fpsfile("fps.dat");
    fpsfile<<fps;
    fpsfile.close();
    ofstream misc("misc.dat");
    misc<<
    xmin<<'\n'<<
    xmax<<'\n'<<
    Vmax<<'\n'<<
    bar<<'\n';
    misc.close();
    ofstream Vfile("V.dat");
    for (int n=0;n<=nt;n++){
        if (n%co_ktora==0){
            update_V(V,nx,x,t[n],E0,omega_ext);
            for (int i=0;i<=nx;i++){
                Vfile<<t[n]<<'\t'<<x[i]<<'\t'<<real(V[i])<<'\n';
            }
        }
    }
    Vfile.close();
    // czystki
    for (int i=0;i<=nx;i++){
        delete [] psi[i];
    }
    delete [] psi;
    delete [] V;
    delete [] x;
    delete [] prawa;
    delete [] nw;
    delete [] diag;
    delete [] t;
    Exx.close();
    // return zero
    return 0;
}
