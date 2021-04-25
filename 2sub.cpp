#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <boost/format.hpp>
#include <fmt/core.h>
#include "mfsw.hpp"

using namespace std;


struct param {
double J;
double K;
double hx;
double hy;
double hz;
double h;
double gamma;
};

void xbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J + 2 * p.K; // 相互作用行列をセット
  Js(1, 2) = p.gamma;
  Js(1, 1) = p.J;
  Js(2, 1) = p.gamma;
  Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void ybond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J; // 相互作用行列をセット
  Js(0, 2) = p.gamma;
  Js(1, 1) = p.J + 2 * p.K;
  Js(2, 0) = p.gamma;
  Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void zbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J; // 相互作用行列をセット
  Js(0, 1) = p.gamma;
  Js(1, 0) = p.gamma;
  Js(1, 1) = p.J;
  Js(2, 2) = p.J + 2 * p.K;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}



int main(int argc, char *argv[])
{
  if(argc!=9) exit(1);


  complex<double> im(0, 1);

  int Ns = 2; // 局所空間の次元
  CPPL::zhematrix Sxop(Ns);
  Sxop.zero();
  CPPL::zhematrix Syop(Ns);
  Syop.zero();
  CPPL::zhematrix Szop(Ns);
  Szop.zero();

  Sxop(0, 1) = 0.5;
  Syop(0, 1) = -0.5 * im;
  Szop(0, 0) = 0.5;
  Szop(1, 1) = -0.5;

  int N = 2*3;  // 平均場の数
  int Nc = 3; // サイトあたりの平均場の数
  int M = 2;

  int Nb = 2*3/2; //ボンドの種類
  int Ni = atoi(argv[1]);

  mfsw ms(N, Nc, Nb, Ni, 1e-15);

  for (size_t i = 0; i < N; i+=Nc)
  {
    ms.set_mat(i) = Sxop;
    ms.set_mat(i+1) = Syop;
    ms.set_mat(i+2) = Szop;
  }

  param p;

  double a = atof(argv[2]);  //K,Jの変数
  p.h = atof(argv[3]); //印加磁場の絶対値
  
  p.hx = atof(argv[4]);
  p.hy = atof(argv[5]);
  p.hz = atof(argv[6]);

  double th = atof(argv[7]);
  double phi = atof(argv[8]);

  p.J = cos(a);
  p.K = sin(a);
  p.gamma = 0.;

  
  ofstream ofs("spec2sub.txt");
  ofstream ofs1("ch2sub.txt");

  // 初期値のinputファイル(行数分だけ回る。列数はNと合致させる)を用いて平均場近似を実行    
    //double phi = 2.0;
  int n_=0;
  xbond(n_, 0, 1, ms, p); n_++;
  ybond(n_, 0, 1, ms, p); n_++;
  zbond(n_, 0, 1, ms, p); n_++;
  
  for(size_t i = 0; i < N; i += 3)
  {
    ms.set_H(i) = p.hx;
    ms.set_H(i+1) = p.hy;
    ms.set_H(i+2) = p.hz;
  }
  ms.exec_mf("cand2sub.txt");
  //cout << "# ";
  cout << ms.mf_out() << endl;   


  double gm = 2. / sqrt(3.0);
  double gk = 4. / 3.;

  for (double x = 0.; x < gk; x += 0.01)
  {
    double kx = x;
    double ky = 0;
    //cout << scientific << x << " " << kx << " " << ky << "   ";
    kx *= M_PI;
    ky *= M_PI;
    complex<double> gx = exp((-0.5 * kx - 0.5 / sqrt(3) * ky) * im);
    complex<double> gy = exp((0.5 * kx - 0.5 / sqrt(3) * ky) * im);
    complex<double> gz = exp((1. / sqrt(3) * ky) * im);
    complex<double> dkx_gx = -0.5 * im * gx;
    complex<double> dkx_gy = 0.5 * im * gy;
    complex<double> dkx_gz = 0.;
    complex<double> dky_gx = -0.5/sqrt(3) * im * gx;
    complex<double> dky_gy = -0.5/sqrt(3) * im * gy;
    complex<double> dky_gz = 1./sqrt(3) * im * gz;

    vector<complex<double>> gamma = {gx, gy, gz};
    vector<complex<double>> dkx_gamma = {dkx_gx, dkx_gy, dkx_gz};
    vector<complex<double>> dky_gamma = {dky_gx, dky_gy, dky_gz};
    ofs << x << " " << ms.exec_sw_out(gamma, dkx_gamma, dky_gamma) << endl;
    ofs1 << x << " " << ms.exec_bc_out(gamma,dkx_gamma,dky_gamma) << endl; // x n1, n2, n3, n4 , ..., n8 までのベリー曲率
  }

  for (double x = gk; x < 5./3; x += 0.01)
  {
    double kx = gk-(x-gk);
    double ky = -sqrt(3.0)*kx+4./sqrt(3.0);
    //cout << scientific << x + gk << " " << kx << " " << ky << "   ";
    kx *= M_PI;
    ky *= M_PI;
    complex<double> gx = exp((-0.5 * kx - 0.5 / sqrt(3) * ky) * im);
    complex<double> gy = exp((0.5 * kx - 0.5 / sqrt(3) * ky) * im);
    complex<double> gz = exp((1. / sqrt(3) * ky) * im);
    complex<double> dkx_gx = -0.5 * im * gx;
    complex<double> dkx_gy = 0.5 * im * gy;
    complex<double> dkx_gz = 0.;
    complex<double> dky_gx = -0.5/sqrt(3) * im * gx;
    complex<double> dky_gy = -0.5/sqrt(3) * im * gy;
    complex<double> dky_gz = 1./sqrt(3) * im * gz;

    vector<complex<double>> gamma = {gx, gy, gz};
    vector<complex<double>> dkx_gamma = {dkx_gx, dkx_gy, dkx_gz};
    vector<complex<double>> dky_gamma = {dky_gx, dky_gy, dky_gz};
    ofs << x << " " << ms.exec_sw_out(gamma, dkx_gamma, dky_gamma) << endl;
    ofs1 << x << " " << ms.exec_bc_out(gamma,dkx_gamma,dky_gamma) << endl; // x n1, n2, n3, n4 , ..., n8 までのベリー曲率
  }
  for (double x = 5./3.; x < 8./3.; x += 0.01)
  {
    double kx = 1. - (x - 5./3.);
    double ky = kx*1./sqrt(3.0);
   // cout << scientific << x + gk + 2 * gm << " " << kx << " " << ky << "   ";
    kx *= M_PI;
    ky *= M_PI;
    complex<double> gx = exp((-0.5 * kx - 0.5 / sqrt(3) * ky) * im);
    complex<double> gy = exp((0.5 * kx - 0.5 / sqrt(3) * ky) * im);
    complex<double> gz = exp((1. / sqrt(3) * ky) * im);
    complex<double> dkx_gx = -0.5 * im * gx;
    complex<double> dkx_gy = 0.5 * im * gy;
    complex<double> dkx_gz = 0.;
    complex<double> dky_gx = -0.5/sqrt(3) * im * gx;
    complex<double> dky_gy = -0.5/sqrt(3) * im * gy;
    complex<double> dky_gz = 1./sqrt(3) * im * gz;

    vector<complex<double>> gamma = {gx, gy, gz};
    vector<complex<double>> dkx_gamma = {dkx_gx, dkx_gy, dkx_gz};
    vector<complex<double>> dky_gamma = {dky_gx, dky_gy, dky_gz};
    ofs << x << " " << ms.exec_sw_out(gamma, dkx_gamma, dky_gamma) << endl;
    ofs1 << x << " " << ms.exec_bc_out(gamma,dkx_gamma,dky_gamma) << endl; // x n1, n2, n3, n4 , ..., n8 までのベリー曲率
  }

  
}
