import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# import pickle

# 313 Eular angle 

def func_satellite(x, t, J1, J2, J3):
    # 変数の取り出し
    omega1 = x[0]
    omega2 = x[1]
    omega3 = x[2]
    phi1 = x[3]
    phi2 = x[4]
    phi3 = x[5]

    # 運動方程式
    domega1 = (J2-J3) * omega2 * omega3/J1
    domega2 = (J3-J1) * omega3 * omega1/J2
    domega3 = (J1-J2) * omega1 * omega2/J3

    s1 = np.sin(phi1)
    s2 = np.sin(phi2)
    s3 = np.sin(phi3)
    c1 = np.cos(phi1)
    c2 = np.cos(phi2)
    c3 = np.cos(phi3)

    # 3-2-1 オイラー角を用いる場合のキネマティクス方程式
    dphi1 = omega2 * s3/c2 + omega3 * c3/c2
    dphi2 = omega2 * c3 - omega3 * s3
    dphi3 = omega1 + omega2 * s2 * s3/c2 + omega3 * s2*c3/c2
    
    # 3-1-3 オイラー角を用いる場合のキネマティクス方程式
    # dphi1 = omega1 * s3/s2 + omega2 * c3/s2
    # dphi2 = omega1 * c3 - omega2 * s3
    # dphi3 = -omega1 * s3*c2/s2 - omega2 * c3*c2/s2 + omega3


    return [domega1, domega2, domega3,   dphi1, dphi2, dphi3]



if (__name__ == '__main__'):
    t_list = np.linspace(0.0, 100.0, 30000)

    # 慣性モーメントの定義 (いろいろ変えてみよう)
    J1 = 30.0
    J2 = 40.0
    J3 = 50.0 
    # J1 = 30.0
    # J2 = 50.0
    # J3 = 40.0 
    
    # 状態変数 x の初期値 [ω1, ω2, ω3, φ1, φ2, φ3] の順で格納している
    # ω1 には微小擾乱を与えている
    # ω3 に最大の初期角速度を与えることで、Z軸(第3軸)をスピン軸にしている
    x_0 = [0.1, 0.0, 1,      0.0, 0.0, 0.0]  
    
    #微分方程式を解く
    x_result = odeint(func_satellite, x_0, t_list, args=(J1, J2, J3))

    # グラフ描画の準備
    fig = plt.figure()
    r=3
    c=2
    ax1 = fig.add_subplot(c, r, 1)
    ax2 = fig.add_subplot(c, r, 2)
    ax3 = fig.add_subplot(c, r, 3)
    ax4 = fig.add_subplot(c, r, 4)
    ax6 = fig.add_subplot(c, r, (5,6), projection='3d')
    c1,c2,c3,c4 = "blue","green","red","black"
    
    omega1 = x_result[:, 0]
    omega2 = x_result[:, 1]
    omega3 = x_result[:, 2]
    phi1 = x_result[:, 3]
    phi2 = x_result[:, 4]
    phi3 = x_result[:, 5]

    # φ1, φ2, φ3の時系列変化
    ax1.plot(t_list, phi1, color=c1)
    ax1.set_xlabel("$t$ [s]") 
    ax1.set_ylabel(r"$\phi_1$ [rad]") 
 
    ax2.plot(t_list, phi2, color=c2)
    ax2.set_xlabel("$t$ [s]") 
    ax2.set_ylabel(r"$\phi_2$ [rad]") 

    ax3.plot(t_list, phi3, color=c3) 
    ax3.set_xlabel("$t$ [s]") 
    ax3.set_ylabel(r"$\phi_3$ [rad]") 

    # ω1, ω2, ω3の時系列変化
    ax4.plot(t_list, omega1, label=r'$\omega_1$', color=c1)
    ax4.plot(t_list, omega2, label=r'$\omega_2$', color=c2)
    ax4.plot(t_list, omega3, label=r'$\omega_3$', color=c3)
    ax4.set_xlabel("$t$ [s]")
    ax4.set_ylabel(r"$\omega$ [rad]")
    ax4.legend()
    # after plotting your 3 lines on ax (or ax5)
    ax4.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.)
    # make room on the right (pick a value that fits your legend)
    ax4.figure.subplots_adjust(right=0.78)

    # 角速度ベクトルの成分表示
    ax6.plot(omega1, omega2, omega3)
    ax6.set_xlabel(r"$\omega_1$")
    ax6.set_ylabel(r"$\omega_2$")
    ax6.set_zlabel(r"$\omega_3$")  
    
    
    fig.tight_layout()  
    plt.show()
    #plt.savefig("spin_attitude.png")
    
    
    # 結果の入った構造体の保存（別のファイルで使用する時用)
    # filename = 'ode_result.pkl'
    # # Save the `OdeResult` object to a file using pickle
    # with open(filename, 'wb') as file:
    #     pickle.dump(x_result, file)