#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
const cmplx ii = cmplx(0.0,1.0); 
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);
void LINstep(cmplx* const u1,  cmplx* const u0,  const double dt,
          const double dx, cmplx D,const int N);
void NLINste(cmplx* u0,const double dt,const int N);
//-----------------------------------
int main(){
	cmplx D = -ii;
	

	const int N = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (N - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);
	double t=0;
	const double eta = 0.2;
	
	cmplx* u0 = new cmplx[N];
	cmplx* u1 = new cmplx[N];
	cmplx* h;
  
	stringstream strm;	


	init(u0, eta, dx, dt,N);

	writeToFile(u0,"psi_0", dx,N,xmin);
	
	
	LINstep(u1,u0,dt/2,dx,D,N);
      
      
  


	for (int i = 1; i <= Na; i++) {
	
		for (int j = 1; j <= Nk-1; j++) {
		  LINstep(u1,u0,dt,dx,D,N);
		  NLINste(u0,dt,N);
		  h = u0;
		  u0 = u1;
		  u1 = h;
		  t +=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(u0,strm.str(), dx,N,xmin);
	}

	return 0;
}

//-----------------------------------
void LINstep(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx,
          cmplx D, const int N)
{

  cmplx* d=new cmplx[N];
  cmplx* u=new cmplx[N];
  cmplx* l=new cmplx[N];

  for(int i=0;i<N;i++) d[i] = 1.0 + 2.0*D*dt/(dx*dx);
  for(int i=0;i<N;i++) u[i] = - D*dt/(dx*dx);
  for(int i=0;i<N;i++) l[i] = - D*dt/(dx*dx);
    
    for (int i = 0; i< N-1; i++){
        d[i+1] -= (l[i+1]/d[i]) * u[i]; // hiermit elliminiere ich die Unterdiagonales arrays l
        f0[i+1] -=  (l[i+1]/d[i]) * f0[i];} // f0 wächst genau so wie d
        
    f1[N-1] = f0[N-1]/d[N-1]; // somit habe ich den letzten Eintrag des f1 Arrays. Nun kann ich die anderen davon ausgehend füllen
  for (int i = N-2; i > 0; i--){
      f1[i] = (f0[i] - u[i]*f1[i+1])/d[i];}    //warum gilt diese Formel?
        
        

  delete[] d;
  delete[] u;
  delete[] l;
}
void NLINste(cmplx* u0,const double dt,const int N)
{
  for (int j = 0; j < N; j++)
    u0[j] = u0[j]*exp(-ii*abs(u0[j])*abs(u0[j])*dt);
}
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
