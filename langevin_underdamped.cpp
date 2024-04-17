#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <random>
#include <string>
#include <cstring>
#include <unistd.h>
/*#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>*/


using namespace std;


double Apot(double,double,double); //potential

double Fx(double,double,double); //力


long SEED; //乱数のシード

//ループ変数達
long n,step;


double dt,TRANSIENT_TIME; //時間幅、観測時間
long TRANSIENT_STEP; //時間ステップ数

long N; //サンプル数

double Gamma,kB,Tbath,coef1,coef; //環境パラメータ



double H; //Hamiltonian
double W,Q,dW,dQ; //仕事、熱

double mass; //質量

double rx1,rx2,vx1,vx2; //位置、速度


double fx1,fx2; //力
double fx1noise,fx2noise; //ノイズ

double X,k,dX,dk; //ポテンシャルパラメータ



int main(int argc, char *argv[]){
    
    
    SEED = 1; //乱数のシード
    
    N = 1; //サンプル数

    dt = 0.05; //時間ステップ幅
    TRANSIENT_TIME = 100.0; //1サンプルの全時間
    TRANSIENT_STEP = (long)TRANSIENT_TIME/dt; //時間ステップ数
    
    kB = 1.0; //Boltzmann定数
    Gamma = 1.0; //抵抗係数
    
    Tbath = 1.0; //熱浴温度
    
    mass = 1.0; //粒子の質量
    
    X = 0.0; //ポテンシャル中心
    k = 5.0; //ばね定数
    
    
    dX = -0.5;
    dk = 2.0;
    
    
    
    
    
    //ノイズの乱数(メルセンヌツイスターとそれを用いたガウスノイズ)
    mt19937_64 mt(SEED);
    normal_distribution<double> gasdev(0, 1);
    
    
START:
    
    
//サンプル数ループ
    
    for(n=1; n<=N; n++){
        
        
        W = 0.0;
        
        rx1 = sqrt(Tbath/k)*gasdev(mt);
        vx1 = sqrt(Tbath/mass)*gasdev(mt);
        fx1 = Fx(rx1,X,k);
        fx1noise = sqrt(2.0*Gamma*kB*Tbath/dt)*gasdev(mt);
        fx2 = 0.0;
        fx2noise = 0.0;
        
        
        //時間ステップループ
        for(step=0; step<TRANSIENT_STEP; step++){
            
                    
        //Verlet法による位置、速度発展
            rx2 = rx1 + dt*vx1 + (0.5/mass)*dt*dt*(fx1 + fx1noise - Gamma*vx1);
                    
            fx2 = Fx(rx2,X,k);
            fx2noise = sqrt(2.0*Gamma*kB*Tbath/dt)*gasdev(mt);
                            
                            
            vx2 = vx1 + (0.5/mass)*dt*(fx1 + fx2 + fx1noise + fx2noise - Gamma*vx1 - Gamma*vx1) - (0.25*Gamma/(mass*mass))*dt*dt*(fx1 + fx2 + fx1noise + fx2noise - Gamma*vx1 - Gamma*vx1);

                    
            rx1 = rx2;
            vx1 = vx2;
            fx1 = fx2;
            fx1noise = fx2noise;
            
            } //時間ステップループ終了
        
    } //サンプル数ループ終了

    
 END:

  return 0;

}


//ポテンシャル関数
double Apot(double rx,double X,double k){

  double apot;

    apot = 0.5*k*(rx-X)*(rx-X);
    
  return apot;

}

//ポテンシャル力
double Fx(double rx,double X,double k){
    
    double fx;
    
    fx = -k*(rx-X);
    
    return fx;
    
    
}
