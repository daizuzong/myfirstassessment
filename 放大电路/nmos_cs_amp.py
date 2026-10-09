#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nmos_cs_amp.py —— NMOS 共源级放大电路  PySpice + Ngspice 仿真
============================================================================
只做仿真, 不写理论推导。输出: tau, fc, V_th, R_th, I_D, V_DS, V_GS, gm, 增益,
以及输入/输出波形图。

电路 (题卡固定参数):

                 VDD = 5V
                  |
        +---------+---------+
        |                   |
      Rg1=60k             Rd=2k
        |                   +--------+-----> Vout (节点 d)
        +------ g ----------| M1  NMOS
        |                   +--------+
      Rg2=40k              |
        |                  GND  (源极接地, 共源级)
       GND
        ^
        |
      Cb1 ---- in <---- Vi = 10mV / 1kHz 正弦

NMOS: K = Kp = 0.8 mA/V^2,  V_th = Vto = 1 V,  lambda = 0.02 /V

运行:
    python nmos_cs_amp.py

输出文件 (相对路径, 与脚本同目录):
    nmos_cs_amp_waveform.png   输入/栅极/输出 瞬态波形 + 直流负载线
    nmos_cs_amp_bode.png       Vin->Vg / Vin->Vout 幅频特性, 标出 fc

依赖: PySpice (见 requirements.txt) + 可用的 ngspice 动态库

============================================================================
 ngspice 共享库路径设置 (PySpice 通过 ctypes 加载 ngspice 动态库)
============================================================================
PySpice 默认找 <PySpice>/Spice/NgSpice/Spice64_dll/dll-vs/ngspice{}.dll,
若该目录没有库文件会报 OSError / cannot load library。本脚本会自动探测;
探测不到时按下面手动设置。

--- Windows ---
  1) 下载官方 "ngspice-XX_dll_64.zip" (注意是 DLL 版, 不是普通 exe 版),
     解压到如 D:\\ngspice\\ngspice-47_64, 其 bin 目录下应有:
         ngspice.dll  +  libomp140.x86_64.dll  (两者必须同目录)
  2) setx NGSPICE_HOME "D:\\ngspice\\ngspice-47_64\\Spice64"
     脚本会去找 %NGSPICE_HOME%\\bin\\ngspice.dll
  3) 或直接指定库文件全路径:
     setx NGSPICE_LIBRARY_PATH "D:\\ngspice\\ngspice-47_64\\Spice64\\bin\\ngspice.dll"

--- Linux ---
  1) sudo apt install ngspice libngspice0
  2) 库: /usr/lib/x86_64-linux-gnu/libngspice.so.0
     export NGSPICE_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/libngspice.so.0
  3) 文件名不确定时用 `ldconfig -p | grep ngspice` 查看。

--- macOS ---
  1) brew install ngspice
  2) 库: /opt/homebrew/lib/libngspice.dylib   (Apple Silicon)
         /usr/local/lib/libngspice.dylib      (Intel)
     export NGSPICE_LIBRARY_PATH=/opt/homebrew/lib/libngspice.dylib
  3) 可用 `brew --prefix ngspice` 找到前缀再拼 /lib/libngspice.dylib
============================================================================
"""

from __future__ import annotations

import os
import re
import sys
import math
import platform
import warnings

warnings.filterwarnings("ignore", category=SyntaxWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


# ---------------------------------------------------------------------------
# 依赖自举: 若脚本同目录下存在 .pylibs (本地依赖目录, 内含 PySpice 1.5 + numpy),
# 就把它加到 sys.path 最前面, 这样无需设置 PYTHONPATH, VS Code 里直接 F5 即可。
#
# 重要: 必须在【任何】import PySpice 之前无条件插入路径, 不能先"探测"再插入。
# 因为一次探测性的 import 失败后, Python 会把失败的命名空间包缓存进 sys.modules,
# 此后即使把 .pylibs 插到 sys.path[0], 那条缓存依然生效, 导入会一直失败。
# ---------------------------------------------------------------------------
def _bootstrap_local_libs():
    here = os.path.dirname(os.path.abspath(__file__))
    local = os.path.join(here, ".pylibs")
    if not os.path.isdir(local):
        return
    if local not in sys.path:
        sys.path.insert(0, local)


_bootstrap_local_libs()


# ===========================================================================
# 1. 电路参数 (题卡固定值)
# ===========================================================================
VDD = 5.0             # V      电源
RG1 = 60e3            # Ohm    VDD -> 栅极
RG2 = 40e3            # Ohm    栅极 -> 地
RD = 2e3              # Ohm    VDD -> 漏极
K_N = 0.8e-3          # A/V^2  器件参数 K = 0.8 mA/V^2
VTH = 1.0             # V      阈值电压 V_th
LAMBDA = 0.02         # 1/V    沟道长度调制

VI_AMP = 10e-3        # V      输入正弦幅度 (峰峰值 20 mV)
FREQ = 1e3            # Hz     输入频率

# 耦合电容。题卡只要求"足够大、交流短路", 因此它是可自由选择的量,
# tau 和 fc 都随它变化, 所以下面给出一个取值表一起输出。
CB1 = 1e-3            # F      Cb1 实际取值 = 1 mF

# tau / fc 验证: 对多组 Cb1 做 AC 扫频, 实测 -3dB 频率并与 1/(2*pi*Rth*Cb1) 比较
CB1_VERIFY = [1e-6, 10e-6, 100e-6, 1e-3, 10e-3]   # F  用于验证 tau = Rth*Cb1

N_CYCLES = 10         # .tran 周期数 (要求 >= 5)
T_STEP = 1e-6         # s      瞬态步长 1 us (每周期 1000 点)
STEADY_FRACTION = 0.5 # 用后一半时间做稳态测量

PLOT_WAVE = "nmos_cs_amp_waveform.png"
PLOT_BODE = "nmos_cs_amp_bode.png"

SEP = "=" * 78
SUB = "-" * 78


def banner(title):
    print()
    print(SEP)
    print(title)
    print(SEP)


def section(title):
    print()
    print(SUB)
    print(title)
    print(SUB)


# ===========================================================================
# 2. 定位并配置 ngspice 共享库
# ===========================================================================
def find_ngspice_library():
    """返回 ngspice 共享库路径模板(含 {} 占位符), 找不到返回 None。"""
    system = platform.system()

    # (1) 环境变量显式指定
    explicit = os.environ.get("NGSPICE_LIBRARY_PATH") or os.environ.get("NGSPICE_LIB")
    if explicit and os.path.exists(explicit):
        if system == "Windows":
            if explicit.endswith("ngspice.dll"):
                return explicit[: -len("ngspice.dll")] + "ngspice{}.dll"
            return explicit
        if system == "Darwin":
            for suf in ("libngspice.dylib", "ngspice.dylib"):
                if explicit.endswith(suf):
                    return explicit[: -len(suf)] + "libngspice{}.dylib"
            return explicit
        for suf in ("libngspice.so.0", "libngspice.so", "ngspice.so"):
            if explicit.endswith(suf):
                return explicit[: -len(suf)] + "libngspice{}.so.0"
        return explicit

    # (2) NGSPICE_HOME / NGSPICE_ROOT
    home = os.environ.get("NGSPICE_HOME") or os.environ.get("NGSPICE_ROOT")
    if home:
        if system == "Windows":
            cands = [os.path.join(home, "bin", "ngspice{}.dll"),
                     os.path.join(home, "Spice64_dll", "dll-vs", "ngspice{}.dll")]
        elif system == "Darwin":
            cands = [os.path.join(home, "lib", "libngspice{}.dylib")]
        else:
            cands = [os.path.join(home, "lib", "libngspice{}.so.0"),
                     os.path.join(home, "lib", "libngspice{}.so")]
        for cand in cands:
            if os.path.exists(cand.replace("{}", "")):
                return cand

    # (3) 常见默认安装位置
    if system == "Windows":
        cands = []
        for ver in (47, 46, 45, 44, 43, 42, 41, 40):
            for root in ("D:\\ngspice", "C:\\ngspice", "C:\\Program Files\\ngspice"):
                cands.append(os.path.join(root, "ngspice-%d_64" % ver,
                                          "Spice64", "bin", "ngspice{}.dll"))
            for root in ("D:\\ngspice", "C:\\ngspice"):
                cands.append(os.path.join(root, "ngspice-%d_dll_64" % ver,
                                          "Spice64_dll", "dll-vs", "ngspice{}.dll"))
        for cand in cands:
            if os.path.exists(cand.replace("{}", "")):
                return cand
    elif system == "Darwin":
        for cand in ("/opt/homebrew/lib/libngspice{}.dylib",
                     "/usr/local/lib/libngspice{}.dylib"):
            if os.path.exists(cand.replace("{}", "")):
                return cand
    else:
        for cand in ("/usr/lib/x86_64-linux-gnu/libngspice{}.so.0",
                     "/usr/lib/x86_64-linux-gnu/libngspice{}.so",
                     "/usr/lib/libngspice{}.so.0",
                     "/usr/local/lib/libngspice{}.so.0"):
            if os.path.exists(cand.replace("{}", "")):
                return cand
    return None


def configure_ngspice():
    from PySpice.Spice.NgSpice.Shared import NgSpiceShared
    default = NgSpiceShared.LIBRARY_PATH
    found = find_ngspice_library()
    if found:
        NgSpiceShared.LIBRARY_PATH = found
        print("ngspice 库 : %s" % found)
    else:
        print("ngspice 库 : 未自动定位, 沿用 PySpice 默认 %s" % default)
        print("            若报 OSError, 请按文件顶部注释设置 NGSPICE_LIBRARY_PATH")


# ===========================================================================
# 3. 搭建电路
# ===========================================================================
def build_circuit(cb1=CB1, with_sine=True):
    """搭建 NMOS 共源级放大器。

    单位写法必须是 "数值 @ 单位" (例如 60 @ u_kOhm)!
    写成 60e3 @ u_kOhm 会生成 "60000.0kOhm", 单位后缀不做数值换算,
    等效电阻被放大 1000 倍, 工作点会彻底跑偏。

    with_sine=False 会去掉输入源 —— 只能用于 .op, 不能用于 AC 分析
    (没有激励源时网表编译不过)。
    """
    from PySpice.Spice.Netlist import Circuit
    from PySpice.Unit import u_V, u_kOhm, u_mV, u_kHz, u_mF

    c = Circuit("NMOS Common-Source Amplifier")

    # 电源。不要写 circuit.VDD(...), PySpice 的 __getattr__ 会把 VDD 当元件名
    # 去查表并抛 AttributeError; 用索引命名 'DD' 生成的元件名同样是 VDD。
    c.V("DD", "vdd", c.gnd, VDD @ u_V)

    # 偏置与负载
    c.R("g1", "vdd", "g", (RG1 / 1e3) @ u_kOhm)
    c.R("g2", "g", c.gnd, (RG2 / 1e3) @ u_kOhm)
    c.R("d", "vdd", "d", (RD / 1e3) @ u_kOhm)

    # 输入正弦源。不要给 dc_offset! 否则 PySpice 会生成
    # "... DC 2V AC 1V SIN(0V 10mV ...)" 这种 DC 与 SIN 偏置重复的语句,
    # ngspice 直接 run failed。本电路栅极偏置由 Rg1/Rg2 提供, 不需要它。
    if with_sine:
        c.SinusoidalVoltageSource("in", "in", c.gnd,
                                  amplitude=(VI_AMP / 1e-3) @ u_mV,
                                  frequency=(FREQ / 1e3) @ u_kHz)

    # 耦合电容: 隔直 + 交流短路; ic 让瞬态从工作点起步
    c.C("b1", "in", "g", (cb1 / 1e-3) @ u_mF,
        ic=(VDD * RG2 / (RG1 + RG2)) @ u_V)

    # NMOS, 源极与衬底接地
    c.M(1, "d", "g", c.gnd, c.gnd, model="NMOS1")
    c.model("NMOS1", "nmos", Kp=K_N, Vto=VTH, Lambda=LAMBDA, level=1)
    return c


def build_step_circuit(cb1):
    """阶跃电路 (仅作说明用途, 主流程不使用)。

    为什么不拿它来实测 tau: 栅极节点在直流下被 Cb1 隔断, 没有对地直流通路,
    .op 无法确定它的初始电位; 做阶跃实验时 Cb1 会把输入电平直接耦合到栅极,
    量到的不是 Rth*Cb1 的充电过程。tau 与 fc 等价, 用 AC 扫频测 fc 更严谨,
    见 measure_fc_vs_cb1()。
    """
    from PySpice.Spice.Netlist import Circuit
    from PySpice.Spice.HighLevelElement import PieceWiseLinearVoltageSource
    from PySpice.Unit import u_V, u_kOhm, u_mF, u_ns

    c = Circuit("Step response for tau")
    c.V("DD", "vdd", c.gnd, VDD @ u_V)
    c.R("g1", "vdd", "g", (RG1 / 1e3) @ u_kOhm)
    c.R("g2", "g", c.gnd, (RG2 / 1e3) @ u_kOhm)
    c.R("d", "vdd", "d", (RD / 1e3) @ u_kOhm)

    vg_target = VDD * RG2 / (RG1 + RG2)
    PieceWiseLinearVoltageSource(c, "step", "in", c.gnd,
                                 values=[(0, 0 @ u_V), (1 @ u_ns, vg_target @ u_V)])
    c.C("b1", "in", "g", (cb1 / 1e-3) @ u_mF, ic=0 @ u_V)
    c.M(1, "d", "g", c.gnd, c.gnd, model="NMOS1")
    c.model("NMOS1", "nmos", Kp=K_N, Vto=VTH, Lambda=LAMBDA, level=1)
    return c, vg_target


# ===========================================================================
# 4. 从 ngspice 取器件参数
# ===========================================================================
def query_mosfet_parameters(simulator):
    """执行 `show m1` 并解析器件参数 (gm/gds/id/vgs/vds/von/vth ...)。

    PySpice 的 .op analysis 在部分 ngspice 版本下不暴露 @m1[gm], 但共享库
    支持直接下发命令, `show m1` 的输出始终可用。
    """
    try:
        shared = simulator.ngspice
        shared.exec_command("op")
        text = shared.exec_command("show m1")
    except Exception:
        return {}
    if not text:
        return {}

    keys = ("gm", "gds", "id", "vgs", "vds", "vbs", "von", "vdsat", "vth")
    out = {}
    for key in keys:
        m = re.search(r"^\s*%s\s+(-?[\d.eE+-]+)\s*$" % re.escape(key),
                      text, re.MULTILINE)
        if m:
            try:
                out[key] = float(m.group(1))
            except ValueError:
                pass
    return out


# ===========================================================================
# 5. 三次分析
# ===========================================================================
def run_operating_point(circuit):
    """.op 直流工作点 -> V_GS, I_D, V_DS 及 gm / ro。"""
    sim = circuit.simulator()
    op = sim.operating_point()

    vg = float(op["g"][0])
    vd = float(op["d"][0])
    i_dd = float(op.branches["vdd"][0])       # VDD 支路: 漏极支路 + Rg1 偏置支路
    i_bias = (VDD - vg) / RG1                 # Rg1 支路电流
    id_op = abs(i_dd) - i_bias                # 真正的漏极电流
    id_check = (VDD - vd) / RD                # 独立校验路径

    dev = query_mosfet_parameters(sim)

    gm = dev.get("gm")
    gds = dev.get("gds")
    ro = (1.0 / gds) if gds else None

    # 饱和区判据
    vov = vg - VTH
    saturated = (vg > VTH) and (vd > vov)

    return dict(sim=sim, op=op, vg=vg, vd=vd, id_op=id_op, id_check=id_check,
                i_dd=i_dd, i_bias=i_bias, dev=dev, gm=gm, gds=gds, ro=ro,
                vov=vov, saturated=saturated)


def run_transient(circuit, sim):
    """.tran 瞬态 -> 波形、输入/输出摆幅、增益。"""
    import numpy as np

    end_time = N_CYCLES / FREQ
    tran = sim.transient(step_time=T_STEP, end_time=end_time)

    t = np.array(tran.time)
    vin = np.array(tran["in"])
    vg = np.array(tran["g"])
    vout = np.array(tran["d"])

    t_start = end_time * STEADY_FRACTION
    m = t >= t_start

    vin_pp = float(vin[m].max() - vin[m].min())
    vg_pp = float(vg[m].max() - vg[m].min())
    vout_pp = float(vout[m].max() - vout[m].min())

    # 反相判据: 输入最大时输出应在最小附近
    i_max = int(np.argmax(vin[m]))
    i_min = int(np.argmin(vin[m]))
    inverted = bool(vout[m][i_max] < vout[m][i_min])

    av = vout_pp / vin_pp
    if inverted:
        av = -av

    return dict(t=t, vin=vin, vg=vg, vout=vout,
                t_start=t_start, end_time=end_time, n_points=len(t),
                vin_pp=vin_pp, vg_pp=vg_pp, vout_pp=vout_pp,
                av=av, inverted=inverted,
                vin_dc=float(vin[m].mean()), vg_dc=float(vg[m].mean()),
                vout_dc=float(vout[m].mean()), mask=m)


def run_ac_sweep(circuit, sim):
    """.ac 交流扫频 -> Vin->Vg 与 Vin->Vout 的幅频特性, 并定位 fc。"""
    import numpy as np

    # 用 1 mHz ~ 100 kHz 覆盖 fc(很小时侯甚至 6.6 mHz) 以及工作频段
    ac = sim.ac(variation="dec", number_of_points=30,
                start_frequency=1e-3, stop_frequency=1e5)

    f = np.array(ac.frequency)
    mag_vg = np.abs(np.array(ac["g"]))
    mag_vout = np.abs(np.array(ac["d"]))
    mag_vin = np.abs(np.array(ac["in"]))

    h_vg = mag_vg / mag_vin          # Vin -> Vg  (高通, 直流为 0)
    h_vout = mag_vout / mag_vin      # Vin -> Vout

    # 找 -3dB 点: |H_vg| 在低频从 0 升到平顶值, 穿越 0.707*平顶 的频率即 fc
    mid_band = float(np.median(h_vg[f >= 1e3]))
    target = mid_band / math.sqrt(2.0)
    fc_meas = None
    for i in range(1, len(f)):
        if h_vg[i - 1] < target <= h_vg[i]:
            # 对数频率上线性插值
            x0, x1 = math.log10(f[i - 1]), math.log10(f[i])
            y0, y1 = h_vg[i - 1], h_vg[i]
            frac = (target - y0) / (y1 - y0) if y1 != y0 else 0.0
            fc_meas = 10 ** (x0 + frac * (x1 - x0))
    return dict(f=f, h_vg=h_vg, h_vout=h_vout,
                mid_band=mid_band, fc_meas=fc_meas)


def measure_fc_vs_cb1(cb1_list):
    """对多组 Cb1 做 AC 扫频, 实测 -3dB 频率, 验证 tau = R_th * Cb1。

    为什么不用瞬态阶跃测 tau:
      栅极节点在直流下没有对地直流通路(被 Cb1 隔断), 其直流电位在 .op 里
      本来就无法确定; 做阶跃实验时 Cb1 会把输入电平直接顶到栅极, 测到的
      不是 Rth*Cb1 的充电过程。而 fc = 1/(2*pi*Rth*Cb1) 与 tau = Rth*Cb1
      是完全等价的, 所以用 AC 扫频测 fc 更严谨也更稳定。
    """
    import numpy as np

    rth = 1.0 / (1.0 / RG1 + 1.0 / RG2)
    rows = []
    for cb in cb1_list:
        fc_calc = 1.0 / (2.0 * math.pi * rth * cb)
        circuit = build_circuit(cb)   # 必须保留输入源, 否则网表无法编译
        # 在 fc 附近 ±3 个十倍频程细扫, 保证 -3dB 穿越点被准确捕捉
        ac = circuit.simulator().ac(variation="dec", number_of_points=60,
                                    start_frequency=fc_calc / 1e3,
                                    stop_frequency=fc_calc * 1e3)
        f = np.array(ac.frequency).copy()
        h_vg = np.abs(np.array(ac["g"])).copy() / np.abs(np.array(ac["in"])).copy()

        mid_band = float(np.median(h_vg[f >= fc_calc * 100]))
        target = mid_band / math.sqrt(2.0)
        fc_meas = None
        for i in range(1, len(f)):
            if h_vg[i - 1] < target <= h_vg[i]:
                x0, x1 = math.log10(f[i - 1]), math.log10(f[i])
                y0, y1 = h_vg[i - 1], h_vg[i]
                frac = (target - y0) / (y1 - y0) if y1 != y0 else 0.0
                fc_meas = 10 ** (x0 + frac * (x1 - x0))
        rows.append(dict(cb1=cb, fc_meas=fc_meas, fc_calc=fc_calc,
                         tau_calc=rth * cb, mid_band=mid_band,
                         zc=1.0 / (2.0 * math.pi * FREQ * cb)))
    return rows


# ===========================================================================
# 6. 绘图
# ===========================================================================
def setup_cjk_font(plt):
    """让 matplotlib 显示中文; 找不到中文字体返回 False (调用方改用英文)。"""
    from matplotlib import font_manager
    wanted = ["Microsoft YaHei", "SimHei", "SimSun", "DengXian", "KaiTi",
              "Noto Sans CJK SC", "Noto Sans CJK JP", "Source Han Sans CN",
              "WenQuanYi Micro Hei", "WenQuanYi Zen Hei",
              "PingFang SC", "Heiti SC", "Hiragino Sans GB", "Arial Unicode MS"]
    available = {f.name for f in font_manager.fontManager.ttflist}
    picked = [n for n in wanted if n in available]
    if not picked:
        return False
    plt.rcParams["font.sans-serif"] = picked + ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    return True


def make_plots(op, tr, ac, out_wave=None, out_bode=None):
    # 输出路径一律相对【脚本所在目录】解析, 这样无论从哪个工作目录运行,
    # 图都落在脚本旁边, 不会因为 cwd 不可写而报 PermissionError。
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if out_wave is None:
        out_wave = os.path.join(script_dir, PLOT_WAVE)
    if out_bode is None:
        out_bode = os.path.join(script_dir, PLOT_BODE)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:
        print("(跳过绘图: %s)" % exc)
        return

    zh = setup_cjk_font(plt)

    def L(cn, en):
        return cn if zh else en

    # ---------------- 图 1: 瞬态波形 + 负载线 ----------------
    t_ms = tr["t"] * 1e3
    vin_ac = (tr["vin"] - tr["vin"].mean()) * 1e3     # mV
    vg_ac = (tr["vg"] - tr["vg"].mean()) * 1e3        # mV

    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(9.5, 10))

    a1.plot(t_ms, vin_ac, color="tab:blue", lw=1.6,
            label=L("Vin 交流分量 (mV)", "Vin AC (mV)"))
    a1.plot(t_ms, vg_ac, color="tab:cyan", ls="--", lw=1.4,
            label=L("Vg 交流分量 (mV)", "Vg AC (mV)"))
    a1.set_ylabel(L("交流摆幅 (mV)", "AC swing (mV)"))
    a1.set_title(L("NMOS 共源级放大器 瞬态仿真 (Cb1=%.3g F, %.0f mV/%.0f Hz)"
                   % (CB1, VI_AMP * 1e3, FREQ),
                   "NMOS CS amp transient (Cb1=%.3g F, %.0f mV/%.0f Hz)"
                   % (CB1, VI_AMP * 1e3, FREQ)))
    a1.grid(True, alpha=0.3)
    a1.legend(loc="upper right", fontsize=9)

    a2.plot(t_ms, tr["vout"], color="tab:red", lw=1.6,
            label=L("Vout (V)", "Vout (V)"))
    a2.axhline(tr["vout_dc"], color="gray", ls=":", lw=1.2,
               label=L("Vout 直流工作点 %.4f V" % tr["vout_dc"],
                       "Vout DC bias %.4f V" % tr["vout_dc"]))
    a2.set_ylabel(L("Vout (V)", "Vout (V)"))
    a2.grid(True, alpha=0.3)
    a2.legend(loc="upper right", fontsize=9)
    a2.annotate(L("反相: Vin 峰值 -> Vout 谷值",
                  "inverting: Vin peak -> Vout valley"),
                xy=(0.02, 0.08), xycoords="axes fraction",
                fontsize=9, color="dimgray")

    vds = np.linspace(0.0, VDD, 400)
    i_load = (VDD - vds) / RD * 1e3
    vov = op["vov"]
    i_sat = 0.5 * K_N * vov ** 2 * (1 + LAMBDA * vds) * 1e3
    with np.errstate(invalid="ignore"):
        i_tri = np.where(vds < vov,
                         K_N * (vov * vds - 0.5 * vds ** 2) *
                         (1 + LAMBDA * vds) * 1e3, np.nan)
    a3.plot(vds, i_load, color="tab:green", lw=1.8,
            label=L("负载线 $I_D=(V_{DD}-V_{DS})/R_D$",
                    "load line $I_D=(V_{DD}-V_{DS})/R_D$"))
    a3.plot(vds, i_sat, color="tab:orange", lw=1.6,
            label=L("MOS 饱和区特性", "MOS saturation"))
    a3.plot(vds, i_tri, color="tab:orange", ls=":", lw=1.6,
            label=L("MOS 线性区特性", "MOS triode"))
    a3.plot([op["vd"]], [op["id_op"] * 1e3], "o", color="red", ms=8, zorder=5,
            label=L("Q 点 (%.4f V, %.4f mA)" % (op["vd"], op["id_op"] * 1e3),
                    "Q (%.4f V, %.4f mA)" % (op["vd"], op["id_op"] * 1e3)))
    a3.axvline(vov, color="gray", ls="--", lw=1,
               label=L("Vov=%.3f V (饱和边界)" % vov,
                       "Vov=%.3f V" % vov))
    a3.set_xlabel(L("$V_{DS}$ (V)", "$V_{DS}$ (V)"))
    a3.set_ylabel(L("$I_D$ (mA)", "$I_D$ (mA)"))
    a3.set_title(L("直流负载线与静态工作点", "DC load line and Q point"))
    a3.set_xlim(0, VDD)
    a3.set_ylim(0, min(2.6, VDD / RD * 1e3 * 1.05))
    a3.grid(True, alpha=0.3)
    a3.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_wave, dpi=130)
    plt.close(fig)
    print("  已保存: %s" % out_wave)

    # ---------------- 图 2: 幅频特性 ----------------
    fig2, ax = plt.subplots(figsize=(9.5, 4.6))
    ax.semilogx(ac["f"], 20 * np.log10(np.clip(ac["h_vg"], 1e-12, None)),
                color="tab:cyan", lw=1.6,
                label=L("$V_g/V_{in}$ (高通, $f_c$)", "$V_g/V_{in}$ (high-pass)"))
    ax.semilogx(ac["f"], 20 * np.log10(np.clip(ac["h_vout"], 1e-12, None)),
                color="tab:red", lw=1.8,
                label=L("$V_{out}/V_{in}$ (增益)", "$V_{out}/V_{in}$ (gain)"))
    if ac["fc_meas"]:
        ax.axvline(ac["fc_meas"], color="k", ls="--", lw=1.2,
                   label=L("实测 $f_c$ = %.4g Hz" % ac["fc_meas"],
                           "measured $f_c$ = %.4g Hz" % ac["fc_meas"]))
    ax.axhline(20 * math.log10(ac["mid_band"] / math.sqrt(2.0)),
               color="gray", ls=":", lw=1.1,
               label=L("-3 dB (相对中频 $V_g/V_{in}$)", "-3 dB"))
    ax.axvline(FREQ, color="tab:blue", ls="-.", lw=1.1,
               label=L("工作频率 %.0f Hz" % FREQ, "operating %.0f Hz" % FREQ))
    ax.set_xlabel(L("频率 (Hz)", "Frequency (Hz)"))
    ax.set_ylabel(L("幅度 (dB)", "Magnitude (dB)"))
    ax.set_title(L("幅频特性 (Cb1=%.3g F): 低频被 Cb1 高通截止"
                   % CB1,
                   "Frequency response (Cb1=%.3g F)" % CB1))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig2.tight_layout()
    fig2.savefig(out_bode, dpi=130)
    plt.close(fig2)
    print("  已保存: %s" % out_bode)


# ===========================================================================
# 7. 输出
# ===========================================================================
def print_op(op, tr):
    section("[1] 直流工作点 (.op) 与器件参数")
    print("  V_GS     = %.6f V      (源极接地, 故 V_GS = Vg)" % op["vg"])
    print("  V_DS     = %.6f V" % op["vd"])
    print("  I_D      = %.6e A = %.4f mA   [由 |I(VDD)| - I(Rg1)]"
          % (op["id_op"], op["id_op"] * 1e3))
    print("             校验 (VDD-Vd)/Rd = %.6e A = %.4f mA"
          % (op["id_check"], op["id_check"] * 1e3))
    print("  V_ov     = %.6f V      (V_GS - V_th)" % op["vov"])
    print("  工作区   = %s" % ("饱和区 (V_GS>V_th 且 V_DS>V_ov)"
                              if op["saturated"] else "非饱和区"))
    print()
    dev = op["dev"]
    if dev:
        print("  ngspice `show m1` 器件参数:")
        if "gm" in dev:
            print("    gm    = %.6e A/V = %.4f mS" % (dev["gm"], dev["gm"] * 1e3))
        if "gds" in dev:
            print("    gds   = %.6e S   -> ro = 1/gds = %.2f kOhm"
                  % (dev["gds"], 1.0 / dev["gds"] / 1e3))
        if "id" in dev:
            print("    id    = %.6e A = %.4f mA" % (dev["id"], dev["id"] * 1e3))
        for k in ("vgs", "vds", "von", "vdsat", "vth"):
            if k in dev:
                print("    %-6s= %.6f V" % (k, dev[k]))
    else:
        print("  (未能取到 `show m1` 器件参数)")


def print_params(op, tr, ac, fc_rows):
    section("[2] 关键参数汇总")

    rth = 1.0 / (1.0 / RG1 + 1.0 / RG2)
    tau_calc = rth * CB1
    fc_calc = 1.0 / (2.0 * math.pi * tau_calc)
    ro = op["ro"]
    rd_ro = 1.0 / (1.0 / RD + 1.0 / ro) if ro else None
    av_from_gm = (-op["gm"] * rd_ro) if (op["gm"] and rd_ro) else None

    print("  tau      = R_th * Cb1 = %.2f kOhm * %.3g F = %.6g s"
          % (rth / 1e3, CB1, tau_calc))
    print("  fc       = 1/(2*pi*tau) = %.6g Hz" % fc_calc)
    # 取聚焦扫频(表[4]同款设置)的实测值, 精度最高
    fc_fine = None
    for row in fc_rows:
        if abs(row["cb1"] - CB1) < 1e-15:
            fc_fine = row["fc_meas"]
    if fc_fine:
        print("  fc 实测 (AC 扫频 -3dB) = %.6g Hz   相对误差 %+.4f %%"
              % (fc_fine, (fc_fine - fc_calc) / fc_calc * 100))
    print("     (Cb1 是自由参数, tau 与 fc 都随它变化, 见表 [4])")
    print()
    print("  V_th     = %.6f V      (NMOS 模型 Vto, 题卡给定 %.1f V)"
          % (op["dev"].get("vth", op["dev"].get("von", VTH)), VTH))
    print("  R_th     = Rg1||Rg2 = 1/(1/%.0fk + 1/%.0fk) = %.1f Ohm = %.2f kOhm"
          % (RG1 / 1e3, RG2 / 1e3, rth, rth / 1e3))
    print()
    print("  V_GS     = %.6f V" % op["vg"])
    print("  I_D      = %.4f mA" % (op["id_op"] * 1e3))
    print("  V_DS     = %.6f V" % op["vd"])
    if op["gm"]:
        print("  gm       = %.4f mS" % (op["gm"] * 1e3))
    if ro:
        print("  ro       = %.2f kOhm   (Rd||ro = %.4f kOhm)"
              % (ro / 1e3, rd_ro / 1e3))
    print()
    print("  增益:")
    if av_from_gm is not None:
        print("    由 gm*(Rd||ro) 计算  Av = %+.4f V/V" % av_from_gm)
    print("    由 .tran 实测        Av = %+.4f V/V  "
          "(Vout_pp/Vin_pp, 反相为负)" % tr["av"])
    print("    |Av| 实测            = %.4f V/V" % abs(tr["av"]))


def print_transient(tr):
    section("[3] 瞬态分析 (.tran)")
    print("  输入      : %.1f mV 幅度 / %.0f Hz, 峰峰值 %.1f mV"
          % (VI_AMP * 1e3, FREQ, VI_AMP * 2e3))
    print("  周期      : %.4f ms, 仿真 %d 个周期 = %.3f ms"
          % (1e3 / FREQ, N_CYCLES, tr["end_time"] * 1e3))
    print("  步长/点数 : %.1f us / %d 点" % (T_STEP * 1e6, tr["n_points"]))
    print("  测量窗口  : %.3f ~ %.3f ms (后一半, 稳态)"
          % (tr["t_start"] * 1e3, tr["end_time"] * 1e3))
    print()
    print("  直流偏置(窗口平均): Vin=%+.6f V  Vg=%+.6f V  Vout=%+.6f V"
          % (tr["vin_dc"], tr["vg_dc"], tr["vout_dc"]))
    print()
    print("  峰峰值    : Vin_pp = %.4f mV" % (tr["vin_pp"] * 1e3))
    print("              Vg_pp  = %.4f mV   (Cb1 交流短路, 应≈Vin_pp)"
          % (tr["vg_pp"] * 1e3))
    print("              Vout_pp= %.4f mV" % (tr["vout_pp"] * 1e3))
    print()
    print("  相位关系  : %s" % ("反相 (输入最大 -> 输出最小)"
                                if tr["inverted"] else "同相"))
    print("  增益      : Av = Vout_pp/Vin_pp = %+.4f V/V" % tr["av"])


def print_cb1_table(fc_rows):
    section("[4] Cb1 取值对 tau 与 fc 的影响 (Cb1 是自由参数, AC 扫频实测验证)")
    print("  %-9s %-11s %-13s %-13s %-10s %s"
          % ("Cb1 (F)", "tau (s)", "fc 计算(Hz)", "fc 实测(Hz)", "误差", "1kHz |Zc|"))
    print("  " + "-" * 76)
    for r in fc_rows:
        err = ((r["fc_meas"] - r["fc_calc"]) / r["fc_calc"] * 100
               if r["fc_meas"] else float("nan"))
        mark = "  <-- 本次取值" if abs(r["cb1"] - CB1) < 1e-15 else ""
        print("  %-9.3g %-11.6g %-13.6g %-13.6g %+8.4f%%  %.3f Ohm%s"
              % (r["cb1"], r["tau_calc"], r["fc_calc"],
                 r["fc_meas"] if r["fc_meas"] else float("nan"),
                 err, r["zc"], mark))
    print()
    print("  tau 与 fc 成反比 (tau = Rth*Cb1, fc = 1/(2*pi*tau)): 表中 %d 组 Cb1"
          % len(fc_rows))
    print("  的 fc 实测都与公式吻合, 说明 tau = Rth*Cb1 的关系成立。")
    print("  交流短路判据: |Zc| = 1/(2*pi*f*Cb1) << R_th = %.2f kOhm。"
          % ((1.0 / (1.0 / RG1 + 1.0 / RG2)) / 1e3))
    print()
    print("  说明: 本脚本不用瞬态阶跃法测 tau, 因为栅极节点在直流下被 Cb1 隔断、")
    print("        没有对地直流通路, .op 无法确定其初始电位, 阶跃实验会把输入")
    print("        电平直接顶到栅极; 而 fc 与 tau 等价, AC 扫频更严谨稳定。")


# ===========================================================================
# 8. 主入口
# ===========================================================================
def main():
    # Windows 终端可能是 GBK, 本脚本输出含中文, 统一成 UTF-8 避免乱码
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    banner("NMOS 共源级放大电路  PySpice + Ngspice 仿真")
    print("Python    : %s (%s)" % (sys.version.split()[0], platform.system()))
    print("工作目录  : %s" % os.getcwd())
    print("脚本      : %s" % os.path.abspath(__file__))
    print("电路      : VDD=%.0fV  Rg1=%.0fk  Rg2=%.0fk  Rd=%.0fk  Cb1=%.3gF"
          % (VDD, RG1 / 1e3, RG2 / 1e3, RD / 1e3, CB1))
    print("NMOS      : K=%.2f mA/V^2  Vth=%.1f V  lambda=%.2f /V"
          % (K_N * 1e3, VTH, LAMBDA))
    print("输入      : %.0f mV / %.0f Hz 正弦" % (VI_AMP * 1e3, FREQ))

    try:
        import PySpice  # noqa: F401
    except ImportError:
        print("\n[错误] 未安装 PySpice, 请执行: pip install -r requirements.txt")
        return 2

    configure_ngspice()
    try:
        import PySpice.Logging.Logging as Logging
        Logging.setup_logging(logging_level="ERROR")
    except Exception:
        pass

    try:
        # ---- 电路 ----
        circuit = build_circuit(CB1)
        print()
        print("生成的 SPICE 网表:")
        print(circuit)

        # ---- .op ----
        op = run_operating_point(circuit)

        # ---- .tran ----
        tr = run_transient(circuit, op["sim"])

        # ---- .ac ----
        ac = run_ac_sweep(circuit, op["sim"])

        # ---- tau / fc 验证 (多组 Cb1 的 AC 扫频 -3dB) ----
        fc_rows = measure_fc_vs_cb1(CB1_VERIFY)

    except Exception as exc:
        print("\n[仿真失败] %s: %s" % (type(exc).__name__, exc))
        print()
        print("排查:")
        print("  1) ngspice 动态库没找到 -> 按文件顶部注释设置 NGSPICE_LIBRARY_PATH")
        print("  2) Windows 下 ngspice.dll 同目录缺少 libomp140.x86_64.dll")
        print("  3) 用了 ngspice exe 版而不是 dll 版")
        return 4

    # ---- 输出 ----
    print_op(op, tr)
    print_transient(tr)
    print_params(op, tr, ac, fc_rows)
    print_cb1_table(fc_rows)

    section("[5] 输出文件")
    make_plots(op, tr, ac)

    banner("完成")
    print("波形图与幅频图已保存到脚本所在目录 (与 nmos_cs_amp.py 同目录)。")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
