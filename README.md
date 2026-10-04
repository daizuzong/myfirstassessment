# myfirstassessment
网址：https://daizuzong.github.io/myfirstassessment/
任务要求2：分支：原稿之外的起草 合并: 将写好内容整合到原稿
- `cd`：进入文件夹
-  `ls`：列出所有文件和子文件夹
-  `mkdir`：新建一个文件夹
-  `git init`：在当前文件夹初始化 Git 仓库
-  `git add .`：把当前所有改动添加到“暂存区” 
- AI代码部分:第一次提示词:"请帮我写一份基于html的贪吃蛇游戏"结果:为我生成了一份基于html的代码，并且附带暂停和计分功能
-  第二次:"请为我注释每一步代码的意义"
-  第三次:"我需要添加一个ai自己玩的模式，而且至少要能连吃15个食物"结果:ai可以独立吃到100分以上，甚至能通关
-  第四次:"设计一个计分板，可以记录最近10局的分数，并且可以自动追加新一局的功能，同时对于历史分数可以通过点击主动查看完整列表"
-  第五次:"请给贪吃蛇设立一个新的浅色页面，可以一键切换，而且刷新后可以保存配置"
-  以上过程提示词皆无调整，AI模型理解力在初阶游戏设计上非常强
-  对于第一阶段的实现采用了VScode和deepseek harness作为辅助工具，去编译贪吃蛇核心逻辑,如if (nh.x < 0 || nh.y < 0 || nh.x >= GRID || nh.y >= GRID) {
        return gameOver("撞到墙了! 😵");   
      }当在坐标系上到达零点（墙上）位置则结束游束
-  对于第二阶段的AI智能化，将蛇头与食物坐标输入，并将身体设为障碍，使ai避障逻辑躲开身体输出路线去寻找食物同时加分
-  对于验证“我通过以下方式验证阶段二功能：

连续运行AI模式3次，每次AI均成功吃完15个食物，并弹出‘AI挑战成功！’弹窗，全程无人工干预。

在AI运行过程中，观察蛇的移动轨迹，确认其能有效避开墙壁和自身，未出现卡死或撞墙现象。

切换手动模式测试，确认手动控制不受AI模式影响，两者互不干扰。
测试结果表明AI策略稳定可靠，达到了‘连续吃15个食物不死’的硬指标要求。”
-  对于戴维南定律的验证，给ai提示词:“请帮我写个基于PySpice的二端网络，其中两个电源均与一个R=0.2Ω的电子串联，其中一个电源为7V，另一个为6.2V。同时他们并联与一个负载电阻R3=3.2Ω。并给出R3支路的电流”
验证思路:通过戴维南定律计算等效电源和电阻后通过负载的电流大小，再用pyspice模拟原电路的负载电流。最后发现结果一致
- 以下是模拟结果与手算:<img width="2009" height="738" alt="two_source_network" src="https://github.com/user-attachments/assets/26d8ccb8-ebd4-41b3-bcbd-a61831e5ac42" />
- 图二:<img width="3072" height="4096" alt="手算" src="https://github.com/user-attachments/assets/9611cb34-7b5e-45d3-a64f-3448e3344d39" />

| 参数 | 理论值 | 仿真值 | 误差 (%) | 备注 |
|------|--------|--------|----------|------|
| U | 6.6 | 6.6 | 0.00 | 无 |
| I | 2.0 | 2.0 | 0.00 | 无 |

- NMOS:<img width="1235" height="1300" alt="nmos_cs_amp_waveform" src="https://github.com/user-attachments/assets/0e60c691-10d4-4a8d-8d8d-3068f6b1e18d" />
- 图二:<img width="1235" height="598" alt="nmos_cs_amp_bode" src="https://github.com/user-attachments/assets/15b1bc9d-c82f-428d-a2fc-b5929a128abe" />
-手算:
- 以下是放大电路的对比表：

| 参数 | 理论值 | 仿真值 | 误差 (%) | 备注 |
|------|--------|--------|----------|------|
| V_GS (V) | 2.00 | 2.00 | 0.00 | 分压固定 |
| I_D (mA) | 4.3 | 4.331 | 0.01 | 饱和区(V_GS>V_th 且 V_DS>V_ov) |
| V_DS (V) | 4.1 | 4.133858 | 0.01 | > V_GS-Vth |
| gm (mS) | 8.66 | 8.66 | 0.00 | sqrt(2*K*I_D) |
| Av (V/V) | -1.7 | -1.7050 | 0.00 | 反相 |

- 低通滤波器:<img width="6144" height="8192" alt="手算 png(1)" src="https://github.com/user-attachments/assets/b3c79523-9f18-436f-ab30-390ab12fd506" />
- 图像:<img width="1482" height="880" alt="rc_lowpass_magnitude" src="https://github.com/user-attachments/assets/d7770a7f-d8c1-48ae-ab08-9d0179498de3" />

- 以下是低通滤波器对比表:

| 参数 | 理论值 | 仿真值 | 误差 (%) | 备注 |
|------|--------|--------|----------|------|
| f | 159 Hz | 159 Hz | 0.00 | 无 |
| τ | 1 ms | 1 ms | 0.01 | 无 |
