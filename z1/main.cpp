#include <iostream>
#include <cmath>
#include <complex>

using namespace std;
typedef complex<double> cmp;
const cmp iu(0,1);

void psi_init(cmp **psi, int nx, double A, double sigma, double *x, double xc, double p0, int nt){
    for (int i=0;i<=nx;i++){
        for (int n=0;n<=nt;n++){
            psi[i][n]=0.0;
        }
    }
    for (int i=0;i<=nx;i++){
        psi[i][0] = A*exp(-sigma*pow(x[i]-xc,2))*exp(iu*p0);
    }
    psi[0][0]=0.0;
    psi[nx][0]=0.0;
}

void thomas(cmp *nw, cmp *prawa, int n, cmp a, cmp b){
    cmp *cprim = new cmp[n-1];
    cmp *dprim = new cmp[n];

    cprim[0]=a/b;
    for (int i=1;i<n-1;i++){
        cprim[i] = a/(b-a*cprim[i-1]);
    }
    dprim[0]=prawa[0]/b;
    for (int i=1;i<n;i++){
        dprim[i]=(prawa[i]-a*dprim[i-1])/(b-a*cprim[i-1]);
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
    const double xmin=-1.89;
    const double xmax=1.89;
    const double tmax=1.0;
    const int nx=100;
    const int nt=100;
    const double dx=(xmax-xmin)/nx;
    const double dt=tmax/nt;
    const double xc=0.0;
    const double p0=10.0;
    const double A=1.758;
    const double sigma=15.0;
    // alokacja
    cmp **psi = new cmp*[nx+1];
    for (int i=0;i<=nx;i++){
        psi[i] = new cmp[nt+1];
    }
    cmp *V = new cmp[nx+1];
    double *x = new double[nx+1];
    // inicjalizacja
    for (int i=0;i<=nx;i++){
        x[i] = xmin+i*dx;
    }
    for (int i=0;i<=nx;i++){
        V[i]=0.0;
    }
    psi_init(psi,nx,A,sigma,x,xc,p0,nt);
    // zmienne do pętli
    cmp *prawa = new cmp[nx-1];
    cmp *nw = new cmp[nx-1];
    fillPrawa(prawa,nx-1,dt,theta,psi,0,dx,V);
    // pętla 
    for (int n=1;n<=nt;n++){
        thomas(nw,prawa,nx-1,-iu*dt*theta,2.0*dx*dx+2.0*iu*dt*theta+2.0*iu*dt*dx*dx*V*theta) // diagonala będzie sie zmieniala razem z V więc trzeba zmodyfikować algorytm thomasa
    }
    // czystki
    for (int i=0;i<=nx;i++){
        delete [] psi[i];
    }
    delete [] psi;
    delete [] V;
    delete [] x;
    delete [] prawa;
    delete [] nw;
    // return zero
    return 0;
}