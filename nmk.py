import stiffness
import numpy as np
from scipy.interpolate import interp1d
import sys

def run(wavePath: str, wave_dt: float, step: float, pga: float,
        m1: float, m2: float, T1: float, T2: float,
        khh_e: float, khh_p: float, xy: float, C: float, alpha: float,
        TT0: float, vzeta: float):

    gamma = 0.5
    beta = 0.25
    # beta = 1.0 / 6.0

    factor = 0.125
    # mass parameters
    m3 = m2 * factor
    m4 = m2 * (1. - factor)

    w = np.genfromtxt(wavePath, delimiter=',')

    t = np.linspace(0.0,wave_dt*(len(w)-1), len(w))

    amp = np.max(np.abs(w))
    wf = interp1d(t, w, kind='linear')  # 可选择其他插值方法如'cubic'

    # stiffness parameters

    k1h_e = 4 * np.pi **2 / T1 ** 2 * m1
    k2h_e = 4 * np.pi **2 / T1 ** 2 * m2
    k2v_e = 4 * np.pi **2 / T2 ** 2 * m2
    kvv_e = 4 * np.pi **2 / (T2 * TT0) ** 2 * m2
    # damping parameters
    c1h = 2 * 0.03 * np.sqrt(k1h_e * m1)
    c2h = 2 * 0.03 * np.sqrt(k2h_e * m2)
    chh = c2h * (k1h_e / (k2h_e + k1h_e)) ** 2\
        + c1h * (k2h_e / (k2h_e + k1h_e)) ** 2
    cvv = 2 * vzeta * np.sqrt(kvv_e * m2)
    c2v = 2 * 0.03 * np.sqrt(k2v_e * m2)


    # define stiffness functions
    khh = stiffness.DoubleLinear(khh_e, khh_p, xy)
    k12 = stiffness.Linear((k2h_e * k1h_e) / (k1h_e + k2h_e))

    kvv = stiffness.Linear(kvv_e)
    k2v = stiffness.Linear(k2v_e)


    # define mass matrix
    m = np.array([[m1, 0., 0., 0.],\
                [0., m2, 0., 0.],\
                [0., 0., m3, 0.],\
                [0., 0., 0., m4]])

    # define damping matrix
    c = np.array([[chh, - chh,  0., 0.],\
                [- chh, chh,  0., 0.],\
                [0., 0., cvv + c2v, - c2v],\
                [0., 0., - c2v, c2v]])

    # 计算迭代的恢复力和刚度矩阵
    def getKbyDx(x: np.ndarray, r: np.ndarray, dx: np.ndarray) -> tuple:

        x_r = np.array([[x[0,0], x[1,0] - x[0,0], x[2,0], x[3,0] - x[2,0]]]).T
        dx_r = np.array([[dx[0,0], dx[1,0] - dx[0,0], dx[2,0], dx[3,0] - dx[2,0]]]).T

        khh_r, khh_k = khh(x_r[0,0], r[0,0], dx_r[0,0])
        k12_r, k12_k = k12(x_r[1,0], r[1,0], dx_r[1,0])
        kvv_r, kvv_k = kvv(x_r[2,0], r[2,0], dx_r[2,0])
        k2v_r, k2v_k = k2v(x_r[3,0], r[3,0], dx_r[3,0])

        k = np.array([[khh_k + k12_k, - k12_k, 0., 0.],
                    [- k12_k, k12_k, 0., 0.],
                    [0., 0., kvv_k + k2v_k, -k2v_k],\
                    [0., 0., -k2v_k, k2v_k]])
        
        rn = np.array([[khh_r, k12_r, kvv_r, k2v_r]]).T

        fs = np.array([[khh_r - k12_r, k12_r, kvv_r - k2v_r, k2v_r]]).T

        return k, fs, rn


    # 每一步迭代前得到初始刚度矩阵
    def getKbyV(x: np.ndarray, r: np.ndarray, v: np.ndarray) -> tuple:
        k, _, _ = getKbyDx(x, r, v)
        fs = np.array([[r[0,0] - r[1,0], r[1,0], r[2,0] - r[3,0], r[3,0]]]).T
        return k, fs

    # newmark 系数
    def a1(dt: float) -> np.ndarray:
       return 1.0 / beta / dt / dt * m + gamma / beta / dt * c

    # newmark 系数
    def a2(dt: float) -> np.ndarray:
        return 1.0 / beta / dt  * m + (gamma / beta - 1.0) * c

    # newmark 系数
    def a3(dt: float) -> np.ndarray:
        return (0.5 / beta - 1.0) * m + dt * (0.5 * gamma / beta - 1.0) * c

    # 进行一步积分
    def integrateOneStep(x: np.ndarray, r: np.ndarray, v: np.ndarray,\
                        p_: np.ndarray, dt: float) -> tuple:
        # 拿到初始状态
        xj = x.copy()
        rj = r.copy()
        kj, fsj = getKbyV(xj, rj, v)
        j = 0

        # 迭代求解
        while True:
            Rj_ = p_ - fsj - a1(dt) @ xj
            if np.linalg.norm(Rj_) < 1e-4 or j >= 1000:
                break
            kj_ = kj + a1(dt)
            dxj = np.linalg.inv(kj_) @ Rj_
            kj, fsj, rj = getKbyDx(xj, rj, dxj)
            xj += dxj
            j += 1
        
        # 因为迭代次数达限制退出
        if j >= 1000:
            return x, r, False
        # 正常退出
        else:
            return xj, rj, True

    # 计算外力中的惯性项等
    def getPiLeft(x: np.ndarray, v: np.ndarray, a: np.ndarray, dt: float) -> np.ndarray:
        return a1(dt) @ x + a2(dt) @ v + a3(dt) @ a

    # 更新状态
    def update(x1: np.ndarray, x2: np.ndarray, v1: np.ndarray,\
                a1: np.ndarray, r2: np.ndarray, dt: float) -> tuple:
        v2 = gamma / beta / dt * (x2 - x1) + (1.0 - gamma / beta) * v1\
            + dt * (1.0 - 0.5 * gamma / beta) * a1
        a2 = 1.0 / beta / dt / dt * (x2 - x1) - 1.0 / beta / dt * v1\
            - (0.5 / beta - 1.0) * a1
        return x2, v2, a2, r2
    
    # 初始状态
    x = np.array([[0., 0., 0., 0.]]).T
    v = np.array([[0., 0., 0., 0.]]).T
    p = m @ np.array([[1.0, 1.0, -0.65, -0.65]]).T * pga / amp * 9806

    r = np.array([[0., 0., 0., 0.]]).T
    _, fs = getKbyV(x, r, v)
    a = np.linalg.inv(m) @ (p * wf(0.0) - c @ v - fs)

    record = []
    record.append([0.0, x[0, 0], x[2, 0]])

    # 主循环
    for i in range((len(w) - 2) * int(wave_dt/step)):
        # 开始按整个时间步长尝试积分
        begin = i * step
        end = (i + 1) * step
        dt = step
        times = 0
        while True:
            pi_ = p * wf(begin + dt) + getPiLeft(x, v, a, dt)\
                + np.array([[-np.sign(v[0,0]) * C * abs(v[0,0]) ** alpha, 0., 0., 0.]]).T
            xi, ri, converge = integrateOneStep(x, r, v, pi_, dt)
            # 如果迭代结果收敛
            if converge:
                x, v, a, r = update(x, xi, v, a, ri, dt)
                begin += dt
                # 如果总积分时间达到步长，退出循环
                if abs(end - begin) < 1e-8:
                    break
                # 如果没达到步长，继续积分
                else:
                    dt = end - begin
                    times = 0
            else:
                # 积分步长缩减多次仍不收敛，直接终止
                if times > 5:
                    print('convergence failed')
                    sys.exit(1)
                # 步长缩减
                dt /= 2
            times += 1
        record.append([i * step, x[0, 0], x[2, 0]])
    record = np.array(record)
    return record

if __name__ == '__main__':
    wavePath = sys.argv[1]
    wave_dt = float(sys.argv[2])
    step = float(sys.argv[3])
    pga = float(sys.argv[4])
    m1 = float(sys.argv[5])
    m2 = float(sys.argv[6])
    T1 = float(sys.argv[7])
    T2 = float(sys.argv[8])
    khh_e = float(sys.argv[9])
    khh_p = float(sys.argv[10])
    xy = float(sys.argv[11])
    C = float(sys.argv[12])
    alpha = float(sys.argv[13])
    TT0 = float(sys.argv[14])
    vzeta = float(sys.argv[15])
    record = run(wavePath, wave_dt, step, pga,
                 m1, m2, T1, T2,
                 khh_e, khh_p, xy, C, alpha,
                 TT0, vzeta)
    t = record[:, 0]
    np.savetxt('time_temp', t, delimiter=',')
    h = record[:, 1]
    np.savetxt('hori_temp', h, delimiter=',')
    v = record[:, 2]
    np.savetxt('vert_temp', v, delimiter=',')