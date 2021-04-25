#ifndef MFSW_HPP
#define MFSW_HPP

//////////////////////////////
// v 1.3
// 2020/5/17
// modify calculations of eigenvalues

// v 1.2
// 2020/5/15
// modify SW.eg() part

// v 1.1
// 2020/5/13
// A bug for calculation of energy is fixed.

// v 1.0
// 2020/5/10
//////////////////////////////

#include <iostream>
#include <complex>
#include <boost/format.hpp>
#include "cpplapack/cpplapack.h"
#include "bogo.hpp"
#include <iterator>

using namespace std;

class mf
{
private:
  double val;
  double MIX;
  CPPL::zhematrix mat;

public:
  CPPL::zhematrix &assign(double MIX_)
  {
    MIX = MIX_;
    return mat;
  }
  void set_init(const double init_)
  {
    val = init_;
  }
  CPPL::zhematrix op()
  {
    return mat;
  }
  double calc(const CPPL::zcovector &v)
  {
    CPPL::zcovector vt = mat * v;
    complex<double> sum(0, 0);
    for (int i = 0; i < v.l; i++)
      sum += conj(v(i)) * vt(i);
    double val_old = val;
    val = MIX * std::real(sum) + (1 - MIX) * val_old;
    return fabs(val - val_old);
  }
  double ave() { return val; };
};

class mfsw
{

private:
  const int Nitr;
  const double MIX;
  const double EPS;

  int N;  // 平均場の数
  int M;  // 副格子の数
  int Nc; // 副格子あたりの平均場の数
  int Nb; //ボンドの種類

  vector<tuple<int, int, CPPL::dsymatrix>> bond; // i,jサイトを繋ぐ相互作用行列
  vector<double> hs;                             // iサイトに働く外場
  vector<mf> MF;

  vector<double> ave;                // 平均場の期待値
  vector<vector<double>> w;          // 平均場ハミルトニアンの固有値
  vector<vector<CPPL::zcovector>> v; //固有ベクトル

  CPPL::zhematrix dxA;
  CPPL::zhematrix A;

  CPPL::zgematrix Vx; //速さ行列
  CPPL::zgematrix Vy;
  CPPL::zhematrix Hsw;
  CPPL::zhematrix dMx;  //H微分
  CPPL::zhematrix dMy;  //H微分
  CPPL::zgematrix omega;

  complex<double> qh;


  vector<double> ene;
  double ene_min;

  template <typename T>
  void read_file(std::string filename, std::vector<std::vector<T>> &x)
  {
    x.clear();
    std::string line;
    std::ifstream ifs(filename.c_str());
    if (!ifs)
    {
      std::cerr << "Can't open the file : " << filename << std::endl;
      std::exit(1);
    }
    int i = 0;
    while (!ifs.eof())
    {
      std::vector<T> temp;
      getline(ifs, line);
      std::istringstream is(line);
      std::copy(std::istream_iterator<T>(is), std::istream_iterator<T>(), std::back_inserter(temp));
      if (temp.size())
      {
        x.resize(i + 1);
        x[i] = temp;
        i++;
      }
    }
  }

  void calc_element(const vector<CPPL::zcovector> &vec, const CPPL::zhematrix &mat, vector<complex<double>> &op)
  {
    op.resize(vec.size() - 1);
    {
      CPPL::zcovector tv = mat * vec[0];
      for (int j = 1; j < vec.size(); j++)
      {
        complex<double> sum(0, 0);
        for (int i = 0; i < vec[j].l; i++)
        {
          sum += conj(vec[j](i)) * tv(i);
        }
        op[j - 1] = sum;
      }
    }
  }

  vector<vector<complex<double>>> OP;

  void eval_op()
  {
    OP.resize(N);
    for (int i = 0; i < N; ++i)
    {
      calc_element(v[i / Nc], MF[i].op(), OP[i]);
    }
  }

public:
  mfsw(int N_, int Nc_, int Nb_, int Nitr_ = 300000, double EPS_ = 1e-20, double MIX_ = 0.05) : N(N_), Nc(Nc_), M(N_ / Nc_), Nb(Nb_),
                                                                                            bond(Nb_), hs(N_, 0), MF(N_), ave(N_), OP(0),
                                                                                            Nitr(Nitr_), MIX(MIX_), EPS(EPS_), ene_min(1e10) {}

  tuple<int, int, CPPL::dsymatrix> &set_J(const int n)
  {
    // CPPL::dsymatrix>> Js is Nc x Nc matrix
    return bond[n];
  }
  double &set_H(const int i) { return hs[i]; }

  CPPL::zhematrix &set_mat(const int n)
  {
    return MF[n].assign(MIX);
  }

  double mf_val(const int n)
  {
    return ave[n];
  }

  double mf_ene()
  {
    return ene_min;
  }

  double mf_diff()//(1,1,1)方向の相図の計算の時のみに使う
  {
    double diff = 0;
    for(int i = 0; i<N; i++)
    {
      diff += abs(ave[i] - 0.288675);   
    }
    return diff;
  }
  double mf_ave(int i)  //aveを他のファイルで使う
  {
    return ave[i];
  }

  string mf_out()
  {

    std::stringstream ss;
    ss << boost::format(" %12.5e   ") % mf_ene() << " " ;
    for (int i = 0; i < N; i++)
    {
      //ss << ", " <<  mf_val(i);
      ss << mf_val(i) << " ";
    }
    return ss.str();
  }

  string spin()
  {
      std::stringstream ss;
    for (int i = 0; i < N; i++)
    {
      //ss << ", " <<  mf_val(i);
      ss << scientific << boost::format(" %12.5e   ") % mf_val(i) << " ";
    }
    return ss.str();

  }

  string diff_out()
  {
    std::stringstream difference;
    difference << mf_diff();

    return difference.str();
  }

  void exec_mf(const string &fn)
  {

    vector<vector<double>> cand;
    read_file(fn, cand);

    int Ncand = cand.size();
    for (int i = 0; i < Ncand; ++i)
      exec_mf(cand[i]);
  }

  void exec_mf(vector<double> &init)
  {
  

    for (int j = 0; j < N; j++)
    {
      MF[j].set_init(init[j]);
    }

    int Ns = MF[0].op().m;

    bool conv = 0;
    double eps;

    vector<vector<double>> wtemp(M);          // 平均場ハミルトニアンの固有値
    vector<vector<CPPL::zcovector>> vtemp(M); //固有ベクトル

    for (int i = 0; i < Nitr; i++)
    {

      vector<CPPL::zhematrix> H(M);

      for (int j = 0; j < M; ++j)
      {
         wtemp[j].clear();
         vtemp[j].clear();
         H[j].resize(Ns);
         H[j].zero();
      }

      for (int j = 0; j < M; ++j)
      {

        for (int k = 0; k < Nc; ++k)
        {
          H[j] = H[j] - hs[k + Nc * j] * MF[k + Nc * j].op();
        }
      }

      for (int l = 0; l < bond.size(); ++l)
      {
        int j1 = get<0>(bond[l]);
        int j2 = get<1>(bond[l]);
        const CPPL::dsymatrix &Js = get<2>(bond[l]);
        for (int k1 = 0; k1 < Nc; ++k1)
        {
          for (int k2 = 0; k2 < Nc; ++k2)
          {
            H[j1] = H[j1] + Js(k1, k2) * MF[k1 + Nc * j2].ave() * MF[k2 + Nc * j1].op();
            H[j2] = H[j2] + Js(k1, k2) * MF[k1 + Nc * j1].ave() * MF[k2 + Nc * j2].op();
          }
        }
      }

      for (int j = 0; j < M; ++j)
      {
        H[j].zheev(wtemp[j], vtemp[j]);
      }

      eps = 0;
      for (int j = 0; j < N; ++j)
      {
        eps += MF[j].calc(vtemp[j / Nc][0]);
      }

      if (eps < EPS)
      {
        conv = 1;
        break;
      }
    }

    double total_ene = 0;
    for (int j = 0; j < M; ++j)
      total_ene += wtemp[j][0];

    for (int l = 0; l < bond.size(); ++l)
    {
      int j1 = get<0>(bond[l]);
      int j2 = get<1>(bond[l]);
      const CPPL::dsymatrix &Js = get<2>(bond[l]);
      for (int k1 = 0; k1 < Nc; ++k1)
      {
        for (int k2 = 0; k2 < Nc; ++k2)
        {
          total_ene += -Js(k1, k2) * MF[k1 + Nc * j2].ave() * MF[k2 + Nc * j1].ave();
        }
      }
    }

    if (conv)
    {
      double temp_ene = total_ene / double(M);
      if (temp_ene < ene_min)
      {
        ene_min = temp_ene;
        for (int i = 0; i < ave.size(); ++i)
        {
          ave[i] = MF[i].ave();
        }
        w = wtemp;
        v = vtemp;
      }
    }
  }

  vector<CPPL::zhematrix> exec_sw(const vector<complex<double>> &gamma, const vector<complex<double>> &dkx_gamma, 
                                  const vector<complex<double>> &dky_gamma)
  {

    if (OP.size() == 0)
    {
      eval_op();
    }

    int Ne = OP[0].size(); // 局所状態の数-1

    int Nm = Ne * M;

    ene.resize(Nm);

    //CPPL::zhematrix A(Nm);
    A.resize(Nm);
    CPPL::zgematrix B(Nm, Nm);
    CPPL::zhematrix C(Nm);

    //CPPL::zhematrix dxA(Nm);  //kx微分
    dxA.resize(Nm);
    CPPL::zgematrix dxB(Nm, Nm);
    CPPL::zhematrix dxC(Nm);

    CPPL::zhematrix dyA(Nm);  //ky微分
    CPPL::zgematrix dyB(Nm, Nm);
    CPPL::zhematrix dyC(Nm);

    A.zero();
    B.zero();
    C.zero();

    dxA.zero();
    dxB.zero();
    dxC.zero();

    dyA.zero();
    dyB.zero();
    dyC.zero();

    for (int k = 0; k < M; ++k)
    {
      for (int n = 0; n < Ne; ++n)
      {
        int j = k * Ne + n;
        A(j, j) = w[k][n + 1] - w[k][0];
        C(j, j) = w[k][n + 1] - w[k][0];
        // dxA(j, j) = 0.;
        // dxC(j, j) = 0.;
        // dyA(j, j) = 0.;
        // dyC(j, j) = 0.;
      }
    }

    for (int l = 0; l < Nb; ++l)
    {
      int j1 = get<0>(bond[l]);
      int j2 = get<1>(bond[l]);
      const CPPL::dsymatrix &Js = get<2>(bond[l]);
      for (int n1 = 0; n1 < Ne; ++n1)
      {
        for (int n2 = 0; n2 < Ne; ++n2)
        {
          complex<double> Jtt(0), JttmT(0), Jto(0), Jot(0);
          complex<double> dxJtt(0), dxJttmT(0), dxJto(0), dxJot(0);
          complex<double> dyJtt(0), dyJttmT(0), dyJto(0), dyJot(0);
          
          for (int k1 = 0; k1 < Nc; ++k1)
          {
            for (int k2 = 0; k2 < Nc; ++k2)
            {
              Jtt += Js(k1, k2) * gamma[l] * OP[k1 + Nc * j1][n1] * OP[k2 + Nc * j2][n2];
              JttmT += Js(k1, k2) * conj(gamma[l]) * OP[k1 + Nc * j1][n2] * OP[k2 + Nc * j2][n1];
              Jto += Js(k1, k2) * gamma[l] * OP[k1 + Nc * j1][n1] * conj(OP[k2 + Nc * j2][n2]);
              Jot += Js(k1, k2) * gamma[l] * conj(OP[k1 + Nc * j1][n1]) * OP[k2 + Nc * j2][n2];

              dxJtt += Js(k1, k2) * dkx_gamma[l] * OP[k1 + Nc * j1][n1] * OP[k2 + Nc * j2][n2];
              dxJttmT += Js(k1, k2) * conj(dkx_gamma[l]) * OP[k1 + Nc * j1][n2] * OP[k2 + Nc * j2][n1];
              dxJto += Js(k1, k2) * dkx_gamma[l] * OP[k1 + Nc * j1][n1] * conj(OP[k2 + Nc * j2][n2]);
              dxJot += Js(k1, k2) * dkx_gamma[l] * conj(OP[k1 + Nc * j1][n1]) * OP[k2 + Nc * j2][n2];

              dyJtt += Js(k1, k2) * dky_gamma[l] * OP[k1 + Nc * j1][n1] * OP[k2 + Nc * j2][n2];
              dyJttmT += Js(k1, k2) * conj(dky_gamma[l]) * OP[k1 + Nc * j1][n2] * OP[k2 + Nc * j2][n1];
              dyJto += Js(k1, k2) * dky_gamma[l] * OP[k1 + Nc * j1][n1] * conj(OP[k2 + Nc * j2][n2]);
              dyJot += Js(k1, k2) * dky_gamma[l] * conj(OP[k1 + Nc * j1][n1]) * OP[k2 + Nc * j2][n2];
            }
          }
          A(j1 * Ne + n1, j2 * Ne + n2) += Jto;
          C(j1 * Ne + n1, j2 * Ne + n2) += Jot;
          B(j1 * Ne + n1, j2 * Ne + n2) += Jtt;
          B(j2 * Ne + n1, j1 * Ne + n2) += JttmT;

          dxA(j1 * Ne + n1, j2 * Ne + n2) += dxJto;
          dxC(j1 * Ne + n1, j2 * Ne + n2) += dxJot;
          dxB(j1 * Ne + n1, j2 * Ne + n2) += dxJtt;
          dxB(j2 * Ne + n1, j1 * Ne + n2) += dxJttmT;

          dyA(j1 * Ne + n1, j2 * Ne + n2) += dyJto;
          dyC(j1 * Ne + n1, j2 * Ne + n2) += dyJot;
          dyB(j1 * Ne + n1, j2 * Ne + n2) += dyJtt;
          dyB(j2 * Ne + n1, j1 * Ne + n2) += dyJttmT;
        }
      }
    }

    bogo<complex<double>> SW;
    SW.calc(A, B, C);
    Hsw.resize(2*Nm);
    Hsw = SW.mcalc(A,B,C);
    


    

    CPPL::zgematrix Vx_(2*Nm, 2*Nm);  //dMx J^{-1}
    CPPL::zgematrix Vy_(2*Nm, 2*Nm);  //dMy J^{-1}
    CPPL::zgematrix omega_(2*Nm,2*Nm);
    omega.resize(2*Nm,2*Nm);  //対角化できているか確認用
    
    dMx.resize(2*Nm);
    dMy.resize(2*Nm);
    dMx = SW.mcalc(dxA, dxB, dxC);
    dMy = SW.mcalc(dyA, dyB, dyC);
    Vx.resize(2*Nm,2*Nm);
    Vy.resize(2*Nm,2*Nm);


    Vx_.zero();
    Vy_.zero();
    Vx.zero();
    Vy.zero();

    // for(size_t i=0;i<Nm;i++){
    //   for(size_t j=i;j<Nm;j++){
    //     dMx(i,j) = dxA(i,j);
    //     dMx(i+Nm,j+Nm) = dxC(i,j); 
    //     dMy(i,j) = dyA(i,j);
    //     dMy(i+Nm,j+Nm) = dyC(i,j);
    //   }
    //   for(size_t j=0;j<Nm;j++){
    //     dMx(i, j+Nm) = dxB(i, j);
    //     dMy(i, j+Nm) = dyB(i, j);
    //   }
    // }

    for(size_t i=0;i<2*Nm;i++){
      for(size_t j=0;j<2*Nm;j++){
        for(size_t k=0;k<2*Nm;k++){
          Vx_(i,j) += dMx(i,k) * SW.Jinv(k,j);
          Vy_(i,j) += dMy(i,k) * SW.Jinv(k,j);
          omega_(i,j) += Hsw(i,k) * SW.Jinv(k,j);
        }
        //cout << SW.Jinv(i,j) << " " ;
      }
      //cout << endl;
    }
    //cout << endl;

    // for(int i=0;i<Nm;i++)
    // {
    //   for(int j=0;j<Nm;j++)
    //   {
    //     cout << dxB(i,j) << " ";
    //   }
    //   cout << endl;
    // }
    // cout << endl;
    for(size_t i=0;i<2*Nm;i++){
      for(size_t j=0;j<2*Nm;j++){
        for(size_t k=0;k<2*Nm;k++){
          Vx(i,j) += conj(SW.Jinv(k,i)) * Vx_(k,j);
          Vy(i,j) += conj(SW.Jinv(k,i)) * Vy_(k,j);
          omega(i,j) += conj(SW.Jinv(k,i)) * omega_(k,j);
        }
      }
    }


    for (int i = 0; i < Nm; i++)
      ene[i] = SW.eg(i);

    

    CPPL::zgematrix spec(Nc, Nm);
    spec.zero();
    for (int i = 0; i < Nc; i++)
    {
      for (int j = 0; j < Nm; j++)
      {
        for (int n = 0; n < Ne; ++n)
        {
          for (int k = 0; k < M; ++k)
          {
            spec(i, j) += OP[i + Nc * k][n] * SW.Vel(k * Ne + n, j) + conj(OP[i + Nc * k][n]) * SW.Uel(k * Ne + n, j);
          }
        }
      }
    }

    vector<CPPL::zhematrix> W;

    W.assign(Nm, CPPL::zhematrix(Nc));
    for (int k = 0; k < Nm; k++)
    {
      for (int i = 0; i < Nc; i++)
      {
        for (int j = 0; j < i; j++)
        {
          W[k](i, j) = spec(i, k) * conj(spec(j, k)) / double(M);
        }
        W[k](i, i) = spec(i, k) * conj(spec(i, k)) / double(M);
      }
    }

    return W;
  }

  complex<double> Omega(size_t i, size_t j) //対角化ハミルトニアン(確認用)
  {
    return omega(i,j);
  }

  complex<double> dxa(int i,int j)
  {
    return dxA(i,j);
  }

  complex<double> a(int i, int j)
  {
    return A(i,j);
  }

  CPPL::zhecomplex ha(size_t i, size_t j)
  {
    return Hsw(i,j);
  }

  double sw_ene(int i)  //スピン波エネルギーを他のファイルで使う
  {
    return ene[i];   
  }

  CPPL::zhecomplex vx(size_t i, size_t j)
  {
    return dMx(i,j);
  }
  

  std::complex<double> bc(const vector<complex<double>> &gamma, const vector<complex<double>> &dkx_gamma, 
                        const vector<complex<double>> &dky_gamma, int n)
  {
    exec_sw_out(gamma,dkx_gamma,dky_gamma);
    int Ne = OP[0].size(); // 局所状態の数-1
    int Nm = Ne * M;
    //int cnt=0;
    double delta = 1e-6;
    double f_be;
    complex<double> d;
    complex<double> im(0.,1.);
    qh = 0.;

    CPPL::zgematrix sgm3(2*Nm,2*Nm);
    CPPL::zgematrix sgm3_(2*Nm,2*Nm);
    CPPL::zgematrix Jd(2*Nm,2*Nm);

    for(int i=0;i<Nm;i++)
    {
      sgm3(i,i) = 1.;
    }
    for(int i=Nm;i<2*Nm;i++)
    {
      sgm3(i,i) = -1.;
    }


    // for(size_t n=0;n<Nm;n++)
    // {
      //f_be = 1./( exp(ene[n]/T) - 1. );
      f_be = 1.;
      for(size_t m = 0;m < 2*Nm; m++)
      {
        if(m != n)  //
        {      
          d = sgm3(m,m)*omega(m,m) - sgm3(n,n)*omega(n,n);
          qh += im * f_be * sgm3(n,n) * sgm3(m,m) * ( Vx(n,m) * Vy(m,n) - Vy(n,m) * Vx(m,n) ) / ( d*d + delta*delta ); 
        }
      }
    //}



    return qh;
  }

  string exec_qh_out(const vector<complex<double>> &gamma, const vector<complex<double>> &dkx_gamma,
                     const vector<complex<double>> &dky_gamma)
  {
    vector<CPPL::zhematrix> W = exec_sw(gamma, dkx_gamma, dky_gamma);

    std::stringstream ss;

    int Ne = OP[0].size(); // 局所状態の数-1
    int Nm = Ne * M;
    vector<double> qh(Nm);
    vector<double> qh_n(Nm);

    exec_sw(gamma, dkx_gamma, dky_gamma);

    for (int i = 0; i < Nm; i++)
    {
      ss << boost::format(" %12.5e ") % ene[i];
    }
    ss << "  ";

    for(int n = 0; n < Nm; n++)
    {
      qh[n] = real( bc(gamma, dkx_gamma, dky_gamma, n) );
      qh_n[n] = qh[n];
    }

    

    // for(int n = 0; n<Nm; n++)
    // {
    //   for(int m=0; m<Nm; m++)
    //   {
    //     if(m != n)
    //     {
    //       if( abs( ene[n]-ene[m] ) < 1e-3 )
    //       {
    //         qh_n[n] += qh[m];
    //       }
    //     }
    //   }
    // }

    for(int n = 0; n< Nm; n++)
    {
      ss << qh_n[n] << " ";
    }

    return ss.str();
  }


///////////////



  string exec_bc_out(const vector<complex<double>> &gamma, const vector<complex<double>> &dkx_gamma, 
                    const vector<complex<double>> &dky_gamma)
  {
    vector<CPPL::zhematrix> W = exec_sw(gamma, dkx_gamma, dky_gamma);
    std::stringstream ss;

    int Ne = OP[0].size(); // 局所状態の数-1
    int Nm = Ne * M;

    

    vector<double> qh(Nm);
    vector<double> qh_n(Nm);

    for(int n = 0; n < Nm; n++)
    {
      qh[n] = real( bc(gamma, dkx_gamma, dky_gamma, n) );
      qh_n[n] = qh[n];
    }

    

    for(int n = 0; n<Nm; n++)
    {
      for(int m=0; m<Nm; m++)
      {
        if(m != n)
        {
          if( abs( ene[n]-ene[m] ) < 1e-3 )
          {
            qh_n[n] += qh[m];
          }
        }
      }
    }

    for(int n = 0; n< Nm; n++)
    {
      ss << qh_n[n] << " ";
    }

    return ss.str();
  }

  string exec_sw_out(const vector<complex<double>> &gamma, const vector<complex<double>> &dkx_gamma, 
                     const vector<complex<double>> &dky_gamma)
  {
    vector<CPPL::zhematrix> W = exec_sw(gamma, dkx_gamma, dky_gamma);

    std::stringstream ss;

    for (int i = 0; i < ene.size(); i++)
    {
      ss << boost::format(" %12.5e  ") % ene[i];
    }
    ss << "    ";

    for (int i = 0; i < W[0].m; ++i)
    {
      for (int k = 0; k < W.size(); k++)
      {
        ss << boost::format(" %12.5e  ") % real(W[k](i, i));
      }
      ss << "  ";
    }

    for (int i = 0; i < W[0].m; ++i)
    {
      for (int j = i + 1; j < W[0].m; ++j)
      {
        for (int k = 0; k < W.size(); k++)
        {
          ss << boost::format(" %12.5e  ") % real(W[k](i, j));
        }
        ss << "  ";
      }
    }

    for (int i = 0; i < W[0].m; ++i)
    {
      for (int j = i + 1; j < W[0].m; ++j)
      {
        for (int k = 0; k < W.size(); k++)
        {
          ss << boost::format(" %12.5e  ") % imag(W[k](i, j));
        }
        ss << "  ";
      }
    }

    return ss.str();
  }
};

#endif
