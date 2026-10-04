import os
os.environ.pop('NGSPICE_LIBRARY_PATH', None)

from PySpice.Spice.NgSpice.Shared import NgSpiceShared
NgSpiceShared.NGSPICE_PATH = r'C:\Users\daizuzong\miniconda3\envs\pyspice_env'
NgSpiceShared.LIBRARY_PATH = r'C:\Users\daizuzong\miniconda3\envs\pyspice_env\Library\bin\ngspice.dll'
"""
RC 低通滤波器 —— PySpice + Ngspice 交流小信号(AC)分析仿真
========================================================

电路结构:
    Vin —— R(1 kΩ) —— Vout —— C(1 μF) —— GND

分析内容:
    1) 交流小信号 AC 扫频: 10 Hz ~ 100 kHz (每十倍频程 20 个点, 对数扫频)
    2) 绘制幅频响应曲线 (增益 dB - 频率, 对数坐标)
    3) 计算理论截止频率  fc = 1 / (2 * pi * R * C) ≈ 159.2 Hz
       (该频率处增益下降 3 dB, 即 -3.01 dB, 对应 |H| = 1/sqrt(2) ≈ 0.707)

运行前提:
    pip install PySpice numpy matplotlib
    并安装 Ngspice 仿真器 (https://ngspice.sourceforge.io/),
    确保其 DLL / 可执行文件能被 PySpice 找到(通常加入 PATH 即可)。
"""

import math

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 0) matplotlib 中文字体设置(否则坐标轴中文会显示为方块)
#    按本机可用字体从上到下依次回退, Windows 通常装有微软雅黑/黑体
# ---------------------------------------------------------------------------
plt.rcParams["font.sans-serif"] = [
    "Microsoft YaHei",      # 微软雅黑
    "SimHei",               # 黑体
    "Noto Sans CJK SC",     # Linux 常用
    "WenQuanYi Zen Hei",    # Linux 备用
    "Arial Unicode MS",     # macOS 备用
]
plt.rcParams["axes.unicode_minus"] = False   # 正常显示负号

# ---------------------------------------------------------------------------
# 1) PySpice 相关导入
#    u_V / u_kΩ / u_uF / u_Hz / u_kHz 等是带量纲的物理量, 便于自动换算单位
# ---------------------------------------------------------------------------
from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import *


# ---------------------------------------------------------------------------
# 2) 电路参数
# ---------------------------------------------------------------------------
R = float(1 @ u_kΩ)   # 电阻 1 kΩ = 1000 Ω
C = float(1 @ u_uF)   # 电容 1 μF = 1e-6 F


# ---------------------------------------------------------------------------
# 3) 理论截止频率: fc = 1 / (2 * pi * R * C)
#    tau = R*C 称为时间常数; 输出滞后于输入 45°、增益降到 -3.01 dB 的频率即 fc
# ---------------------------------------------------------------------------
fc = 1.0 / (2.0 * math.pi * R * C)
print("理论截止频率: fc = 1/(2*pi*R*C) = {:.2f} Hz".format(fc))


# ---------------------------------------------------------------------------
# 4) 搭建电路
#    节点命名: 'vin'(输入) -> R -> 'out'(输出) -> C -> 地(gnd)
#   注意: 不要用 'in' 作为节点名 —— 'in' 是 Python 关键字, PySpice 会打警告,
#        而且无法用 analysis.in 这种属性写法读取, 只能写成 analysis['in'];
#        命名为 'vin' 后两种写法都可以用。
#    SinusoidalVoltageSource 交流源幅值取 1 V, 方便按 |Vout| 直接读增益
# ---------------------------------------------------------------------------
circuit = Circuit("RC 低通滤波器 (1 kΩ / 1 μF)")

circuit.SinusoidalVoltageSource(
    "input",            # 电源名称(SPICE 中写作 Vin)
    "vin",              # 正节点
    circuit.gnd,        # 负节点(接地)
    amplitude=1 @ u_V,  # 交流幅值 1 V
)
circuit.R(1, "vin", "out", R)               # R1: vin --(1 kΩ)-- out, 数值按欧姆计
circuit.C(1, "out", circuit.gnd, C)          # C1: out --(1 μF)-- 地, 数值按法拉计

# ---------------------------------------------------------------------------
# 5) 交流小信号 (AC) 扫频分析
#    从 10 Hz 到 100 kHz, 对数坐标每十倍频程取 20 个采样点
#    结果中的电压为复数(同时含幅度和相位信息)
# ---------------------------------------------------------------------------
simulator = circuit.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.ac(
    start_frequency=10 @ u_Hz,      # 起始频率
    stop_frequency=100 @ u_kHz,     # 截止(终止)频率
    number_of_points=20,            # 每十倍频程点数
    variation="dec",                # 对数(decade)扫频
)

# 提取数据并转为普通 numpy 数组(去掉量纲)
frequency = np.asarray(analysis.frequency)          # 频率, 单位 Hz
v_in = np.asarray(analysis["vin"])              # 输入节点复电压
v_out = np.asarray(analysis.out)                    # 输出节点复电压

# 幅频响应: 增益(线性) = |Vout| / |Vin|; 增益(dB) = 20*lg(增益线性)
gain_linear = np.abs(v_out) / np.abs(v_in)
gain_db = 20.0 * np.log10(gain_linear)

# ---------------------------------------------------------------------------
# 6) 在仿真结果里找到最接近理论截止频率 fc 的点, 校验其增益是否为 -3.01 dB
# ---------------------------------------------------------------------------
idx_fc = int(np.argmin(np.abs(frequency - fc)))
print("仿真中离 fc 最近的频率点: {:.2f} Hz -> 增益 {:.2f} dB (理论 -3.01 dB)".format(
    frequency[idx_fc], gain_db[idx_fc]))
print("低频段(10 Hz)增益: {:.2f} dB | 高频段(100 kHz)增益: {:.2f} dB".format(
    gain_db[0], gain_db[-1]))

# ---------------------------------------------------------------------------
# 7) 绘制幅频响应曲线
#    理论公式对照: |H(jf)| = 1 / sqrt(1 + (f/fc)^2)
# ---------------------------------------------------------------------------
# 理论曲线用的密集频率轴
f_theory = np.logspace(np.log10(10.0), np.log10(100e3), 1000)
h_theory = 1.0 / np.sqrt(1.0 + (f_theory / fc) ** 2)      # 理论增益(线性)
gain_theory_db = 20.0 * np.log10(h_theory)

fig, ax = plt.subplots(figsize=(10, 6))

# 仿真数据(圆点连线)
ax.semilogx(frequency, gain_db, "o-", markersize=4,
            linewidth=1.5, color="blue", label="PySpice 仿真结果")
# 理论曲线(灰色虚线, 用于对照)
ax.semilogx(f_theory, gain_theory_db, "--", color="gray",
            linewidth=1.5, label="理论曲线 20lg|H(jf)|")

# 标出截止频率 fc 与 -3.01 dB 参考线
ax.axvline(fc, color="red", linestyle=":", linewidth=1.5,
           label="fc = {:.1f} Hz".format(fc))
ax.axhline(-3.0103, color="green", linestyle=":", linewidth=1.0)
ax.annotate("-3.01 dB", xy=(1.05 * fc, -3.0103),
            xytext=(2.5 * fc, -3.5), color="green",
            arrowprops=dict(arrowstyle="->", color="green"))

# 坐标轴与标题
ax.set_xlabel("频率 f (Hz)")
ax.set_ylabel("增益 (dB)")
ax.set_title("RC 低通滤波器幅频响应 (R = 1 kΩ, C = 1 μF)")
ax.set_xlim(10, 100e3)
ax.set_ylim(-50, 5)
ax.grid(True, which="both", linestyle="--", alpha=0.5)
ax.legend(loc="upper right")

plt.tight_layout()

# 保存图片: 即使在没有图形界面的环境下也能看到结果
output_png = "rc_lowpass_magnitude.png"
plt.savefig(output_png, dpi=150, bbox_inches="tight")
print("幅频响应曲线已保存为: {}".format(output_png))

# 仅在有图形界面的交互后端下才弹出窗口; Agg 等无界面后端会自动跳过
if "agg" not in matplotlib.get_backend().lower():
    plt.show()
