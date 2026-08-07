"""
内部にばね・ダンパを持つスピン衛星の運動を計算する。
ばね・ダンパ系は1つ。位置rに質点mが取り付けられ、動作方向はnに拘束される。
衛星本体の剛体の運動方程式と、ばね・ダンパの運動方程式を連立して解くことで、
減衰のあるスピン衛星の運動の時刻歴を計算する。
計算結果として機体座標系でみた角速度、質点の位置、角運動量のグラフを表示する。
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# -----------------------------------------------------------------------------
# 衛星と内部質点の定数
# -----------------------------------------------------------------------------
J1 = 2.0  # x 軸まわりの慣性モーメント [kg m^2]
J2 = 2.0  # y 軸まわりの慣性モーメント [kg m^2]
J3 = 3.0  # z 軸まわりの慣性モーメント [kg m^2]

m = 0.1    # 内部質点の質量 [kg]
d = 0.02   # ダンパの減衰係数 [N s/m]
k = 0.1    # ばね定数 [N/m]

# 内部質点の基準位置 r と、質点が移動する方向 n
r = np.array([0.0, 1.0, 0.0])
n = np.array([0.0, 0.0, 1.0])

# r の外積行列 r~
# r~ @ n はベクトルの外積 r × n と等しい。
r1 = r[0]
r2 = r[1]
r3 = r[2]
r_tilde = np.array(
    [
        [0.0, -r3, r2],
        [r3, 0.0, -r1],
        [-r2, r1, 0.0],
    ]
)

r_cross_n = r_tilde @ n

# 変位 z = 0 のときの慣性テンソル J*
j_star = np.diag([J1, J2, J3])


def calculate_inertia(z):
    """r, n と内部質点の変位 z から慣性テンソル J(z) を計算する。"""
    r_n_term = np.outer(r, n) + np.outer(n, r)
    n_plane_term = np.dot(n, n) * np.eye(3) - np.outer(n, n)

    return j_star - m * z * r_n_term + m * z**2 * n_plane_term


def calculate_inertia_rate(z, dz):
    """慣性テンソルの時間微分 dJ/dt を r, n から計算する。"""
    r_n_term = np.outer(r, n) + np.outer(n, r)
    n_plane_term = np.dot(n, n) * np.eye(3) - np.outer(n, n)

    return -m * dz * r_n_term + 2.0 * m * z * dz * n_plane_term


# -----------------------------------------------------------------------------
# 運動方程式 (状態方程式として記述)
# state = [omega1, omega2, omega3, z, dz]
# -----------------------------------------------------------------------------
def equations_of_motion(time, state):
    """角速度、内部質点の変位、変位速度の時間微分を返す。"""
    # solve_ivp が求める形式に合わせて time を受け取るが、
    # この運動方程式は時刻を直接には使わない。
    del time

    omega1 = state[0]
    omega2 = state[1]
    omega3 = state[2]
    z = state[3]
    dz = state[4]

    omega = np.array([omega1, omega2, omega3])

    # omega~ : omega と任意のベクトルの外積を作る行列
    omega_tilde = np.array(
        [
            [0.0, -omega3, omega2],
            [omega3, 0.0, -omega1],
            [-omega2, omega1, 0.0],
        ]
    )

    # 内部質点の変位 z を含む慣性テンソル J(z)
    inertia = calculate_inertia(z)

    # 慣性テンソルの時間微分 dJ/dt
    inertia_rate = calculate_inertia_rate(z, dz)

    # 運動方程式を次の 4 元連立一次方程式にする。
    #
    # coefficient @ [domega1, domega2, domega3, ddz] = right_hand_side
    #
    # r~ n から角運動と並進運動の結合項を作る。
    coefficient = np.array(
        [
            [inertia[0, 0], inertia[0, 1], inertia[0, 2], m * r_cross_n[0]],
            [inertia[1, 0], inertia[1, 1], inertia[1, 2], m * r_cross_n[1]],
            [inertia[2, 0], inertia[2, 1], inertia[2, 2], m * r_cross_n[2]],
            [m * r_cross_n[0], m * r_cross_n[1], m * r_cross_n[2], m],
        ]
    )

    # 角運動方程式の右辺
    angular_right_hand_side = (
        -omega_tilde @ inertia @ omega
        -inertia_rate @ omega
        -m * dz * (omega_tilde @ r_cross_n)
    )

    # 内部質点の並進運動方程式の右辺
    position_vector = r + z * n
    linear_right_hand_side = (
        -d * dz
        -k * z
        -m * n @ omega_tilde @ omega_tilde @ position_vector
    )

    right_hand_side = np.array(
        [
            angular_right_hand_side[0],
            angular_right_hand_side[1],
            angular_right_hand_side[2],
            linear_right_hand_side,
        ]
    )

    acceleration = np.linalg.solve(coefficient, right_hand_side)

    domega1 = acceleration[0]
    domega2 = acceleration[1]
    domega3 = acceleration[2]
    ddz = acceleration[3]

    return np.array([domega1, domega2, domega3, dz, ddz])


def calculate_angular_momentum(states):
    """各時刻の機体座標系における角運動量を計算する。"""
    angular_momentum = []

    for state in states:
        omega = state[0:3]
        z = state[3]
        dz = state[4]

        inertia = calculate_inertia(z)

        # H = J omega + m dz (r~ n)
        momentum = inertia @ omega + m * dz * r_cross_n
        angular_momentum.append(momentum)

    return np.array(angular_momentum)


def main():
    # -------------------------------------------------------------------------
    # 計算条件。変更したい場合はここを書き換える。
    # -------------------------------------------------------------------------
    duration = 500.0       # 解析時間 [s]
    time_step = 0.1        # 結果を取得する時間間隔 [s]

    omega1_initial = 0.2   # [rad/s]
    omega2_initial = 0.0   # [rad/s]
    omega3_initial = 1.0   # [rad/s]
    z_initial = 0.0       # [m]
    dz_initial = 0.0       # [m/s]

    initial_state = np.array(
        [
            omega1_initial,
            omega2_initial,
            omega3_initial,
            z_initial,
            dz_initial,
        ]
    )

    # duration を変えても、結果を取得する時間間隔は time_step に保つ。
    times = np.arange(0.0, duration + 0.5 * time_step, time_step)
    times = times[times <= duration]
    solution = solve_ivp(
        equations_of_motion,
        [0.0, duration],
        initial_state,
        t_eval=times,
        method="DOP853",
        rtol=1.0e-9,
        atol=1.0e-11,
    )

    if not solution.success:
        raise RuntimeError("Numerical integration failed: " + solution.message)

    # solve_ivp の出力を「1 行が 1 時刻」の配列に変換する。
    states = solution.y.T

    # -------------------------------------------------------------------------
    # 計算結果のプロット
    # -------------------------------------------------------------------------
    figure = plt.figure(figsize=(13.0, 10.0))
    figure.suptitle("Case 4: Original Damper, r=1.0  (J1=2, J2=2, J3=3)", fontsize=16)

    angular_velocity_plot = figure.add_subplot(2, 2, 1)
    angular_velocity_plot.plot(times, states[:, 0], label=r"$\omega_1$")
    angular_velocity_plot.plot(times, states[:, 1], label=r"$\omega_2$")
    angular_velocity_plot.plot(times, states[:, 2], label=r"$\omega_3$")
    angular_velocity_plot.set_xlabel(r"Time $t$ [s]")
    angular_velocity_plot.set_ylabel("Angular velocity [rad/s]")
    angular_velocity_plot.grid(True, alpha=0.3)
    angular_velocity_plot.legend()

    displacement_plot = figure.add_subplot(2, 2, 2)
    displacement_plot.plot(times, states[:, 3], color="tab:purple")
    displacement_plot.set_xlabel(r"Time $t$ [s]")
    displacement_plot.set_ylabel(r"Displacement $z$ [m]")
    displacement_plot.grid(True, alpha=0.3)

    momentum = calculate_angular_momentum(states)
    momentum_plot = figure.add_subplot(2, 2, 3)
    momentum_plot.plot(times, momentum[:, 0], label=r"$H_1$")
    momentum_plot.plot(times, momentum[:, 1], label=r"$H_2$")
    momentum_plot.plot(times, momentum[:, 2], label=r"$H_3$")
    momentum_plot.plot(
        times,
        np.linalg.norm(momentum, axis=1),
        color="black",
        linestyle="--",
        label=r"$\|H\|$",
    )
    momentum_plot.set_xlabel(r"Time $t$ [s]")
    momentum_plot.set_ylabel("Angular momentum [N m s]")
    momentum_plot.grid(True, alpha=0.3)
    momentum_plot.legend()

    phase_plot = figure.add_subplot(2, 2, 4, projection="3d")
    phase_plot.plot(states[:, 0], states[:, 1], states[:, 2])
    phase_plot.scatter(
        states[0, 0], states[0, 1], states[0, 2],
        color="tab:green", label="start",
    )
    phase_plot.scatter(
        states[-1, 0], states[-1, 1], states[-1, 2],
        color="tab:red", label="end",
    )
    phase_plot.set_xlabel(r"$\omega_1$ [rad/s]")
    phase_plot.set_ylabel(r"$\omega_2$ [rad/s]")
    phase_plot.set_zlabel(r"$\omega_3$ [rad/s]")
    phase_plot.legend()

    figure.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    figure.savefig("case4_damper_original.png", dpi=150) 
    plt.show()


if __name__ == "__main__":
    main()

