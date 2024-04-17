#langevin_underdamped.cpp　のjulia版です。

using Random #初回起動前にパッケージの追加は事前に行ってください。
#=
パッケージの追加は
juliaを起動した後
import Pkg
Pkg.add("Random")
とするか
パッケージモードで
add Random
とすれば良いです。
=#

#julia langevin_underdamped.jl で使えます。


function main_sim()
  #行の終わりに;は無くても良いです。
  SEED = 1 #乱数のシード
  
    
  N = 1; #サンプル数

  dt = 0.05; #時間ステップ幅
  TRANSIENT_TIME = 100.0; #1サンプルの全時間
  TRANSIENT_STEP = TRANSIENT_TIME/dt; #時間ステップ数
  
  kB = 1.0; #Boltzmann定数
  Gamma = 1.0; #抵抗係数
  
  Tbath = 1.0; #熱浴温度
  
  mass = 1.0; #粒子の質量
  
  X = 0.0; #ポテンシャル中心
  k = 5.0; #ばね定数
  
  
  dX = -0.5;
  dk = 2.0;
  #乱数の生成
  rng = MersenneTwister(SEED)

  #= for i in a:b
      処理
     end
     で aからbまでループします。
     例:
     for i in 1:10
      println(i) 変数iを表示
     end
  =#
  for n in 1:N
    W = 0.0;
        
        rx1 = sqrt(Tbath/k)*randn(rng);
        vx1 = sqrt(Tbath/mass)*randn(rng);
        fx1 = Fx(rx1,X,k);
        fx1noise = sqrt(2.0*Gamma*kB*Tbath/dt)*randn(rng);
        fx2 = 0.0;
        fx2noise = 0.0;
        for step in 0:TRANSIENT_STEP-1
          rx2 = rx1 + dt*vx1 + (0.5/mass)*dt*dt*(fx1 + fx1noise - Gamma*vx1);
                    
          fx2 = Fx(rx2,X,k);
          fx2noise = sqrt(2.0*Gamma*kB*Tbath/dt)*randn(rng);
                          
                          
          vx2 = vx1 + (0.5/mass)*dt*(fx1 + fx2 + fx1noise + fx2noise - Gamma*vx1 - Gamma*vx1) - (0.25*Gamma/(mass*mass))*dt*dt*(fx1 + fx2 + fx1noise + fx2noise - Gamma*vx1 - Gamma*vx1);

                  
          rx1 = rx2;
          vx1 = vx2;
          fx1 = fx2;
          fx1noise = fx2noise;
          #=if step==0
            println("#"*string(N))
          end=#
          #println(rx1)
        end
  end





end

#関数の書き方が二つあります。
#数学的な表現
Fx(rx::Float64,X::Float64,k::Float64) = -k*(rx-X);

#プログラム的な表現
function Apot(rx::Float64,X::Float64,k::Float64) 
  #return をしないと最終行がこの関数の戻り値になります
   0.5*k*(rx-X)*(rx-X);
end

main_sim()#この行をコメントアウトして　include("langevin_underdamped.jl")　とすれば対話形式でもmain_sim()を使えます。
