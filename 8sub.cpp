#include <iostream>
#include <fstream>
#include <complex>
#include <boost/format.hpp>
#include <fmt/core.h>

#include "mfsw.hpp"

using namespace std;
using namespace fmt;


struct param {
double J;
double K;
double h;
double hx;
double hy;
double hz;
double Tmax;
double Tmin;
};

void xbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J + 2 * p.K; // 相互作用行列をセット
  Js(1, 1) = p.J; 
  Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void ybond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J; // 相互作用行列をセット
  Js(1, 1) = p.J + 2 * p.K;
  Js(2, 2) = p.J;
  ms.set_J(n) = std::tuple<int, int, CPPL::dsymatrix>{i, j, Js}; // bond index 0 として、サイト0 と サイト1 の相互作用をセット
}
void zbond(int n, int i, int j, mfsw& ms, const param& p)
{
  int Nc = 3; // サイトあたりの平均場の数
  CPPL::dsymatrix Js(Nc);
  Js.zero();
  Js(0, 0) = p.J; // 相互作用行列をセット
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

  int N = 8*3;  // 平均場の数
  int Nc = 3; // サイトあたりの平均場の数
  int M = 8; //副格子 

  int Nb = 8*3/2; //ボンドの種類

  int Ni = atoi(argv[1]);

  mfsw ms(N, Nc, Nb, Ni, 1e-8);

  //局所状態の演算子をセット SxA:0, SyA:1, SzA:2, SxB:3, SyB:4, SzB:5

  for (size_t i = 0; i < N; i+=Nc)
  {
  ms.set_mat(i) = Sxop;
  ms.set_mat(i+1) = Syop;
  ms.set_mat(i+2) = Szop;
  }

  param p;


  double a = atof(argv[2]);  //K,Jの変数
  p.h = atof(argv[3]); //印加磁場の絶対値
  
  //pythonから実行する時
  p.hx = atof(argv[4]);
  p.hy = atof(argv[5]);
  p.hz = atof(argv[6]);

  double th = atof(argv[7]);
  double phi = atof(argv[8]);

  p.J = cos(a);
  p.K = sin(a);


  int n_ =0;
  xbond(n_, 0, 1, ms, p); n_++;
  xbond(n_, 2, 7, ms, p); n_++;
  xbond(n_, 4, 3, ms, p); n_++;
  xbond(n_, 6, 5, ms, p); n_++;
  ybond(n_, 0, 5, ms, p); n_++;
  ybond(n_, 6, 1, ms, p); n_++;
  ybond(n_, 2, 3, ms, p); n_++;
  ybond(n_, 4, 7, ms, p); n_++;
  zbond(n_, 0, 7, ms, p); n_++;
  zbond(n_, 2, 1, ms, p); n_++;
  zbond(n_, 6, 3, ms, p); n_++;
  zbond(n_, 4, 5, ms, p); n_++;
  
  

  // 局所磁場のセット。 左辺の引数は、SxA:0, SyA:1, SzA:2, SxB:3, SyB:4, SzB:5
  
  for(size_t i = 0; i < N; i += 3)
  {
    ms.set_H(i) = p.hx;
    ms.set_H(i+1) = p.hy;
    ms.set_H(i+2) = p.hz;
  }
  //std::string t = fmt::format("aa");


  // 初期値のinputファイル(行数分だけ回る。列数はNと合致させる)を用いて平均場近似を実行
  // std::string s,S;
  // std::cout << "使用する読み込みファイル..." << endl;
  // std::cin >> s;
  // S = fmt::format(s + ".txt");
  ms.exec_mf("afstar_con.txt");


  std::string fn = fmt::format("qh8sub_α={}_h={}_θ={:.4f}_φ={:.4f}.txt",argv[2],argv[3],th,phi);

  ofstream ofs("spec8sub.txt");
  ofstream ofs1("ch8sub.txt");
  //ofstream ofs2("spin8sub.txt");
  //ofstream ofs2(fn);

  //cout << "# ";
  cout << ms.mf_out() << endl;
  //cout << ms.spin() << endl;

  double gm = 2. / sqrt(3.0);
  double gk = 4. / 3.;


  //Γ→K→M→Γ用
  for (double x = 0.; x < gk; x += 0.01)
    {
      double kx = x;
      double ky = 0;
      //cout << scientific << x << " " << kx << " " << ky << "   ";
      kx *= M_PI;
      ky *= M_PI;
      //vector<complex<double>> gamma = {exp((-0.5 * kx - 0.5 / sqrt(3) * ky) * im), exp((0.5 * kx - 0.5 / sqrt(3) * ky) * im), exp((1. / sqrt(3) * ky) * im)};
      complex<double> gx = exp((0.5 * kx + 0.5 / sqrt(3) * ky) * im);
      complex<double> gy = exp((-0.5 * kx + 0.5 / sqrt(3) * ky) * im);
      complex<double> gz = exp((-1. / sqrt(3) * ky) * im);
      complex<double> dkx_gx = 0.5 * im * gx;
      complex<double> dkx_gy = -0.5 * im * gy;
      complex<double> dkx_gz = 0.;
      complex<double> dky_gx = 0.5/sqrt(3) * im * gx;
      complex<double> dky_gy = 0.5/sqrt(3) * im * gy;
      complex<double> dky_gz = -1./sqrt(3) * im * gz;
      vector<complex<double>> gamma = {gx,gx,gx,gx, gy,gy,gy,gy, gz,gz,gz,gz};
      vector<complex<double>> dkx_gamma = {dkx_gx,dkx_gx,dkx_gx,dkx_gx, dkx_gy,dkx_gy,dkx_gy,dkx_gy, dkx_gz,dkx_gz,dkx_gz,dkx_gz};
      vector<complex<double>> dky_gamma = {dky_gx,dky_gx,dky_gx,dky_gx, dky_gy,dky_gy,dky_gy,dky_gy, dky_gz,dky_gz,dky_gz,dky_gz};
      ofs << x << " " << ms.exec_sw_out(gamma,dkx_gamma,dky_gamma) << endl;
      ofs1 << x << " " << ms.exec_bc_out(gamma,dkx_gamma,dky_gamma) << endl; // x n1, n2, n3, n4 , ..., n8 までのベリー曲率
      //cout << x << " " << ms.exec_bc_out(gamma,dkx_gamma,dky_gamma) << endl;
      //cout << ms.exec_sw_out(gamma) << endl; // energy, xx, yy, zz, Re xy, Re xz, Re yz, Im xy, Im xz, Im yz, それぞれM個ずつ出力
    }
    for (double x = gk; x < 5./3; x += 0.01)
    {
      double kx = gk-(x-gk);
      double ky = -sqrt(3.0)*kx+4./sqrt(3.0);
      //cout << scientific << x + gk << " " << kx << " " << ky << "   ";
      kx *= M_PI;
      ky *= M_PI;
      //vector<complex<double>> gamma = {exp((-0.5 * kx - 0.5 / sqrt(3) * ky) * im), exp((0.5 * kx - 0.5 / sqrt(3) * ky) * im), exp((1. / sqrt(3) * ky) * im)};
      complex<double> gx = exp((0.5 * kx + 0.5 / sqrt(3) * ky) * im);
      complex<double> gy = exp((-0.5 * kx + 0.5 / sqrt(3) * ky) * im);
      complex<double> gz = exp((-1. / sqrt(3) * ky) * im);
      complex<double> dkx_gx = 0.5 * im * gx;
      complex<double> dkx_gy = -0.5 * im * gy;
      complex<double> dkx_gz = 0.;
      complex<double> dky_gx = 0.5/sqrt(3) * im * gx;
      complex<double> dky_gy = 0.5/sqrt(3) * im * gy;
      complex<double> dky_gz = -1./sqrt(3) * im * gz;
      vector<complex<double>> gamma = {gx,gx,gx,gx, gy,gy,gy,gy, gz,gz,gz,gz};
      vector<complex<double>> dkx_gamma = {dkx_gx,dkx_gx,dkx_gx,dkx_gx, dkx_gy,dkx_gy,dkx_gy,dkx_gy, dkx_gz,dkx_gz,dkx_gz,dkx_gz};
      vector<complex<double>> dky_gamma = {dky_gx,dky_gx,dky_gx,dky_gx, dky_gy,dky_gy,dky_gy,dky_gy, dky_gz,dky_gz,dky_gz,dky_gz};
      ofs << x << " " << ms.exec_sw_out(gamma,dkx_gamma,dky_gamma) << endl;
      ofs1 << x << " " << ms.exec_bc_out(gamma,dkx_gamma,dky_gamma) << endl;
      //cout << ms.exec_sw_out(gamma) << endl;
    }
    for (double x = 5./3.; x < 8./3.; x += 0.01)
    {
      double kx = 1. - (x - 5./3.);
      double ky = kx*1./sqrt(3.0);
    // cout << scientific << x + gk + 2 * gm << " " << kx << " " << ky << "   ";
      kx *= M_PI;
      ky *= M_PI;
      //vector<complex<double>> gamma = {exp((-0.5 * kx - 0.5 / sqrt(3) * ky) * im), exp((0.5 * kx - 0.5 / sqrt(3) * ky) * im), exp((1. / sqrt(3) * ky) * im)};
      complex<double> gx = exp((0.5 * kx + 0.5 / sqrt(3) * ky) * im);
      complex<double> gy = exp((-0.5 * kx + 0.5 / sqrt(3) * ky) * im);
      complex<double> gz = exp((-1. / sqrt(3) * ky) * im);
      complex<double> dkx_gx = 0.5 * im * gx;
      complex<double> dkx_gy = -0.5 * im * gy;
      complex<double> dkx_gz = 0.;
      complex<double> dky_gx = 0.5/sqrt(3) * im * gx;
      complex<double> dky_gy = 0.5/sqrt(3) * im * gy;
      complex<double> dky_gz = -1./sqrt(3) * im * gz;
      vector<complex<double>> gamma = {gx,gx,gx,gx, gy,gy,gy,gy, gz,gz,gz,gz};
      vector<complex<double>> dkx_gamma = {dkx_gx,dkx_gx,dkx_gx,dkx_gx, dkx_gy,dkx_gy,dkx_gy,dkx_gy, dkx_gz,dkx_gz,dkx_gz,dkx_gz};
      vector<complex<double>> dky_gamma = {dky_gx,dky_gx,dky_gx,dky_gx, dky_gy,dky_gy,dky_gy,dky_gy, dky_gz,dky_gz,dky_gz,dky_gz};
      ofs << x << " " << ms.exec_sw_out(gamma,dkx_gamma,dky_gamma) << endl;
      ofs1 << x << " " << ms.exec_bc_out(gamma,dkx_gamma,dky_gamma) << endl;
      //cout << ms.exec_sw_out(gamma) << endl;
    }


  //ブリルアンゾーン全体(比熱用)
  
  vector <double> b1(2);
  vector <double> b2(2);
  vector <double> k(2);
  vector<complex<double>> ch;

  b1[0] = 2. * M_PI;
  b1[1] = -2. * M_PI / sqrt(3);
  b2[0] = 0.;
  b2[1] = 4. * M_PI / sqrt(3);
  p.Tmax = 4.;
  p.Tmin = 0.001;
  int cntx = 0;
  int cnty = 0;
  double T = 0;
  complex<double> ch_ = 0.;

  int n = 6; //固有エネルギーのラベル

  ch.push_back(0.);

  

  //for(int n = 0;n<8;n++)
  //{
    // for (double m1 = 0.; m1 <= 1.; m1 += 0.02)
    // {
    //   cntx += 1;
    //   cnty = 0;
    //   for (double m2 = 0.; m2 <= 1.; m2 += 0.02)
    //   {

    //     cnty += 1;
    //     double kx = m1 * b1[0] + m2 * b2[0];
    //     double ky = m1 * b1[1] + m2 * b2[1];

    //     complex<double> gx = exp((0.5 * kx + 0.5 / sqrt(3.) * ky) * im);
    //     complex<double> gy = exp((-0.5 * kx + 0.5 / sqrt(3.) * ky) * im);
    //     complex<double> gz = exp((-1. / sqrt(3.) * ky) * im);
    //     complex<double> dkx_gx = 0.5 * im * gx;
    //     complex<double> dkx_gy = -0.5 * im * gy;
    //     complex<double> dkx_gz = 0.;
    //     complex<double> dky_gx = 0.5/sqrt(3) * im * gx;
    //     complex<double> dky_gy = 0.5/sqrt(3) * im * gy;
    //     complex<double> dky_gz = -1./sqrt(3) * im * gz;
    //     vector<complex<double>> gamma = {gx,gx,gx,gx, gy,gy,gy,gy, gz,gz,gz,gz};
    //     vector<complex<double>> dkx_gamma = {dkx_gx,dkx_gx,dkx_gx,dkx_gx, dkx_gy,dkx_gy,dkx_gy,dkx_gy, dkx_gz,dkx_gz,dkx_gz,dkx_gz};
    //     vector<complex<double>> dky_gamma = {dky_gx,dky_gx,dky_gx,dky_gx, dky_gy,dky_gy,dky_gy,dky_gy, dky_gz,dky_gz,dky_gz,dky_gz};
    //     //ch_ += ms.bc(gamma, dkx_gamma, dky_gamma, n);
    //     //cout << ms.qholl(T,gamma,dkx_gamma,dky_gamma) << endl;
    //     ofs2 << ms.exec_qh_out(gamma,dkx_gamma,dky_gamma) << endl;

    //   }
    // }
  //}
  //cout << " charn number =  " << ch_ / ( double(cntx) * double(cnty) ) * 2. * M_PI / ( sqrt(3.)/2. )  << endl;


  // //K点
  // double kx = 0;
  // double ky = 0.;
  // kx *= M_PI;
  // complex<double> gx = exp((0.5 * kx + 0.5 / sqrt(3) * ky) * im);
  // complex<double> gy = exp((-0.5 * kx + 0.5 / sqrt(3) * ky) * im);
  // complex<double> gz = exp((-1. / sqrt(3) * ky) * im);
  // complex<double> dkx_gx = 0.5 * im * gx;
  // complex<double> dkx_gy = -0.5 * im * gy;
  // complex<double> dkx_gz = 0.;
  // complex<double> dky_gx = 0.5/sqrt(3) * im * gx;
  // complex<double> dky_gy = 0.5/sqrt(3) * im * gy;
  // complex<double> dky_gz = -1./sqrt(3) * im * gz;
  // vector<complex<double>> gamma = {gx,gx,gx,gx, gy,gy,gy,gy, gz,gz,gz,gz};
  // vector<complex<double>> dkx_gamma = {dkx_gx,dkx_gx,dkx_gx,dkx_gx, dkx_gy,dkx_gy,dkx_gy,dkx_gy, dkx_gz,dkx_gz,dkx_gz,dkx_gz};
  // vector<complex<double>> dky_gamma = {dky_gx,dky_gx,dky_gx,dky_gx, dky_gy,dky_gy,dky_gy,dky_gy, dky_gz,dky_gz,dky_gz,dky_gz};
  // std::complex<double> vc = ms.qholl(1,gamma,dkx_gamma,dky_gamma);
  // //cout << ms.qholl(T,gamma,dkx_gamma,dky_gamma) << endl;
  // cout << vc << endl;
}
