import numpy as np
from scipy.optimize import minimize
import sys
import json



def optimization(K_t: float, f_t: float,
                 dmin: float, dmax: float,
                 Dmin: float, Dmax: float,
                 nmin: float, nmax: float,
                 fbmin: float, fbmax: float):
    tau = 740
    G = 79000
    def zx(x):
        return x[3] * x[2] + 2 * x[0]

    def ix1(x):
        C = x[1] / x[0]
        parm = (4 * C - 1) / (4 * C - 4) + 0.615 / C
        return tau * np.pi * x[1] ** 2 * x[2] / parm / G / x[0] - f_t

    def ix2(x):
        C = x[1] / x[0]
        parm = (4 * C - 1) / (4 * C - 4) + 0.615 / C
        fn = tau * np.pi * x[1] ** 2 * x[2] / parm / G / x[0]
        return x[3] - fn / x[2] - 1.2 * x[0]

    def ix3(x):
        return x[1] - x[0]


    def ex1(x):
        return G * x[0] ** 4 / 8 / (x[1] ** 3)/ x[2] - K_t

    def fn(x):
        C = x[1] / x[0]
        parm = (4 * C - 1) / (4 * C - 4) + 0.615 / C
        return tau * np.pi * x[1] ** 2 * x[2] / parm / G / x[0]
    # 设置初始点
    x0 = np.array([200,2000,3,1000])

    # 设置求解参数
    bnds = ((dmin, dmax),
            (Dmin, Dmax),
            (nmin, nmax),
            (fbmin, fbmax))

    cons = ({'type': 'ineq', 'fun': ix1},
                {'type': 'ineq', 'fun': ix2},
                {'type': 'ineq', 'fun': ix3},
                {'type': 'eq', 'fun': ex1},
                )

    # 调用minimize函数求解
    sol = minimize(zx, x0, bounds=bnds, constraints=cons, method='SLSQP')
    return sol


if __name__ == '__main__':
    K_t = float(sys.argv[1])
    f_t = float(sys.argv[2])
    dmin = float(sys.argv[3])
    dmax = float(sys.argv[4])
    Dmin = float(sys.argv[5])
    Dmax = float(sys.argv[6])
    nmin = float(sys.argv[7])
    nmax = float(sys.argv[8])
    fbmin = float(sys.argv[9])
    fbmax = float(sys.argv[10])

    sol = optimization(K_t, f_t,
                        dmin, dmax,
                        Dmin, Dmax,
                        nmin, nmax,
                        fbmin, fbmax)
    
    to_save = {}
    to_save['message'] = sol.message
    to_save['success'] = int(sol.success)
    to_save['status'] = sol.status
    to_save['fun'] = sol.fun
    to_save['x'] = sol.x.tolist()
    to_save['nit'] = sol.nit
    to_save['jac'] = sol.jac.tolist()
    to_save['nfev'] = sol.nfev
    to_save['njev'] = sol.njev
    f = open(r"support_design.json", "w")
    f.write(json.dumps(to_save, indent=4))
    f.close()
    
