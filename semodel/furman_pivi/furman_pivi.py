from dataclasses import dataclass, field
import numpy as np
import scipy.special as sp
import math

@dataclass
class SecondaryEmissionParams:
    # 真電子（backscattered, rediffused それぞれで扱うパラメータ）および
    # 真二次電子（true secondary）に関するパラメータをまとめたクラス
    # --- 真電子 (backscattered) 用パラメータ ---
    E_hat_e: float                  # 基準エネルギー
    W: float                        # 幅（エネルギースペクトルのスケールパラメータ）
    p: float                        # 指数パラメータ
    P_infty_e: float                # δₑ(∞) = P_{1,e}(infty)
    P_hat_e: float                  # δₑのピーク値
    e1: float                       # 角度依存性の係数1
    e2: float                       # 角度依存性の指数
    
    # --- 反拡散 (rediffused) 用パラメータ ---
    P_infty_r: float                # δᵣ(∞) = P_{1,r}(infty)
    Er: float                       # スケールパラメータ（エネルギー）
    r: float                        # 指数パラメータ
    r1: float                       # 角度依存性の係数1
    r2: float                       # 角度依存性の指数
    
    # --- 真二次電子 (true secondary) 用パラメータ ---
    delta_hat_ts: float             # 真二次電子のピーク値
    s: float                        # エネルギースペクトルの補正パラメータ
    E_hat_ts: float                 # 真二次電子の基準エネルギー
    t1: float                       # 角度依存性（第1項）係数
    t2: float                       # 角度依存性（第1項）指数
    t3: float                       # 角度依存性（第2項）係数
    t4: float                       # 角度依存性（第2項）指数
    
    # --- エネルギースペクトルの補正用パラメータ ---
    sigma_e: float                  # 真電子エネルギースペクトルの分散（調整される）
    q: float                        # 反拡散におけるエネルギー分布の指数
    pn: list = field(default_factory=list)         # 真二次電子の形状パラメータリスト（各nに対して）
    epsilon_n: list = field(default_factory=list)    # 真二次電子のエネルギースケールパラメータリスト
        

# -------------------------------
# 各プロセスの計算関数（パラメータを明示的に渡す）
# -------------------------------

def delta_e(params: SecondaryEmissionParams, E0, theta0):
    """
    真電子（backscattered）に対するSEY
    
    E_hat_e,
    W,
    p,
    P_infty_e,
    P_hat_e,
    e1, e2
    """
    # 指数関数的なエネルギー依存性
    t = -((np.abs(E0 - params.E_hat_e) / params.W) ** params.p) / params.p
    delta_e_0 = params.P_infty_e + (params.P_hat_e - params.P_infty_e) * np.exp(t)
    
    # 角度依存性の補正
    angular_factor = 1 + params.e1 * (1 - np.cos(theta0) ** params.e2)
    
    return delta_e_0 * angular_factor

def delta_r(params: SecondaryEmissionParams, E0, theta0):
    """
    反拡散 (rediffused) に対するSEY
    
    Er,
    r,
    P_infty_r,
    r1, r2
    """
    t = -((E0 / params.Er) ** params.r)
    delta_r_0 = params.P_infty_r * (1.0 - np.exp(t))
    
    angular_factor = 1 + params.r1 * (1 - np.cos(theta0) ** params.r2)
    
    return delta_r_0 * angular_factor

def delta_ts(params: SecondaryEmissionParams, E0, theta0):
    """
    真二次電子 (true secondary) に対するSEY
    
    t1, t2, t3, t4,
    E_hat_ts,
    s,
    
    """
    # 角度依存性（2番目の補正項）
    angular_factor34 = 1 + params.t3 * (1 - np.cos(theta0) ** params.t4)
    E_hat = params.E_hat_ts * angular_factor34

    def D(x):
        return params.s * x / (params.s - 1 + x**params.s)
    
    delta_ts_0 = params.delta_hat_ts * D(E0 / E_hat)
    
    # 角度依存性（1番目の補正項）
    angular_factor12 = 1 + params.t1 * (1 - np.cos(theta0) ** params.t2)
    
    return delta_ts_0 * angular_factor12

# -------------------------------
# エネルギースペクトルの関数（f_{1,e}, f_{1,r}, f_{1,ts}）
# -------------------------------

def energy_e(params: SecondaryEmissionParams, E, E0, theta0):
    """
    真電子 (backscattered) に対するエネルギースペクトル f_{1,e}
    
    E: 出力エネルギー
    E0: 入射エネルギー
    theta0: 入射角
    
    sigma_e
    
    delta_e
    """
    # 定義域：Eが負またはE0より大きい場合は0を返す
    if E < 0 or E > E0:
        return 0.0

    # sigma_e の補正（元コード中の定数 aa, bb, cc, dd を用いた補正）
    aa, bb, cc, dd = 1.88, 2.5, 1e-2, 1.5e2
    sigma_modified = (params.sigma_e - aa) + bb * (1 + np.tanh(cc * (E0 - dd)))
    
    delta0 = delta_e(params, E0, theta0)
    
    a = 2 * np.exp(-((E - E0) ** 2) / (2 * sigma_modified ** 2))
    b = np.sqrt(2 * np.pi) * sigma_modified * sp.erf(E0 / (np.sqrt(2) * sigma_modified))
    
    return delta0 * a / b

def energy_r(params: SecondaryEmissionParams, E, E0, theta0):
    """
    反拡散 (rediffused) に対するエネルギースペクトル f_{1,r}
    
    q
    
    delta_r
    """
    if E < 0 or E > E0:
        return 0.0
    
    delta0 = delta_r(params, E0, theta0)
    
    a = (params.q + 1) * E ** params.q
    b = E0 ** (params.q + 1)
    
    return delta0 * a / b

def energy_ts(params: SecondaryEmissionParams, E, E0, theta0):
    """
    真二次電子 (true secondary) に対するエネルギースペクトル f_{ts}
    
    複数の寄与（n=1,...,M; M=10）を足し合わせたものを返す。
    
    epsilon_n
    pn
    
    delta_ts
    """
    if E < 0:
        return 0.0
    if E > E0:
        return 0.0

    M = 10
    f_total = 0.0
    delta0 = delta_ts(params, E0, theta0)
    p_val = delta0 / M
    
    # n=1～Mについて各項を足し合わせる（n=0は除外）
    for n in range(1, M + 1):
        Pn = sp.binom(M, n) * p_val ** n * (1 - p_val) ** (M - n)
        # pn と epsilon_n はリストの各要素
        a = Pn * (E / params.epsilon_n[n - 1]) ** (params.pn[n - 1] - 1) * np.exp(-E / params.epsilon_n[n - 1])
        b = params.epsilon_n[n - 1] * sp.gamma(params.pn[n - 1]) * sp.gammainc(n * params.pn[n - 1], E0 / params.epsilon_n[n - 1])
        c = sp.gammainc((n - 1) * params.pn[n - 1], (E0 - E) / params.epsilon_n[n - 1])
        
        fn = a / b * c
        f_total += fn * n
    
    return f_total

# -------------------------------
# サンプル（モンテカルロ法によるエネルギーサンプル生成）の関数（フィッティングには直接不要な場合もある）
# -------------------------------

def sample_e(params: SecondaryEmissionParams, E0):
    """
    真電子に対応するエネルギーサンプル生成（モンテカルロ用）
    """
    aa, bb, cc, dd = 1.88, 2.5, 1e-2, 1.5e2
    sigma_modified = (params.sigma_e - aa) + bb * (1 + np.tanh(cc * (E0 - dd)))
    u = np.random.rand()
    return E0 + np.sqrt(2) * sigma_modified * sp.erfinv((u - 1) * sp.erf(E0 / (np.sqrt(2) * sigma_modified)))

def sample_r(params: SecondaryEmissionParams, E0):
    """
    反拡散に対応するエネルギーサンプル生成
    """
    u = np.random.rand()
    return u ** (1 / (params.q + 1)) * E0

def sample_ts(params: SecondaryEmissionParams, E0, n):
    """
    真二次電子の n 番目の寄与に対するエネルギーサンプル生成
    """
    pn_val = params.pn[n - 1]
    epsilon = params.epsilon_n[n - 1]
    cdf = sp.gammainc(pn_val, E0 / epsilon)
    u = np.random.rand()
    x = u * cdf
    if x < 1e-12:
        x = 0.0
    return epsilon * sp.gammaincinv(pn_val, x)
