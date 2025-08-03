# AstroLibrary - 天體力學 Python 實作

[![CI](https://github.com/ewdlop/AstroLibrary/actions/workflows/ci.yml/badge.svg)](https://github.com/ewdlop/AstroLibrary/actions/workflows/ci.yml)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## 專案簡介

本專案示範如何使用 Python 及開源套件（Astropy、SciPy 等）實作天體力學（Celestial Mechanics），涵蓋：
- **現代經典力學**：行星運動、開普勒問題  
- **廣義相對論修正**：近星行星進動、時鐘延遲  

範例程式碼基於 [Astropy](https://github.com/astropy/astropy) 與 SciPy 撰寫，並搭配 NumPy、Matplotlib 做數值計算與可視化。

*Not scientifically tested - for educational and demonstration purposes.*

---

## 環境安裝

```bash
# 建議使用 virtualenv / conda 建立隔離環境
pip install -r requirements.txt
```

### 相依套件

- `astropy>=5.0` - 天文計算核心
- `scipy>=1.9.0` - 科學計算與數值方法
- `numpy>=1.21.0` - 數值陣列計算
- `matplotlib>=3.5.0` - 資料視覺化
- `jupyter>=1.0.0` - 互動式筆記本
- `pytest>=7.0.0` - 單元測試框架

---

## 理論背景

### 1. 經典力學（Newtonian Mechanics）

- 牛頓萬有引力定律：

  $$F = G \frac{m_1 m_2}{r^2}$$

- 行星的二體問題、開普勒三定律
- 開普勒方程式數值求解

### 2. 廣義相對論（General Relativity）修正

- 施瓦茲席爾德 (Schwarzschild) 度規下的行星軌道進動

  $$\Delta\omega \approx \frac{6\pi GM}{c^2 a(1 - e^2)}$$

- 時間延遲（Shapiro Delay）

---

## 目錄結構

```
AstroLibrary/
│
├── src/                    # 核心原始碼
│   ├── __init__.py
│   ├── classical.py       # 經典力學相關函式
│   ├── relativity.py      # 相對論修正函式
│   └── utils.py           # 公用工具函式
│
├── notebooks/             # Jupyter 筆記本示例
│   ├── kepler.ipynb       # 開普勒問題演示
│   └── precession.ipynb   # 軌道進動分析
│
├── tests/                 # 單元測試
│   ├── test_classical.py
│   ├── test_relativity.py
│   └── test_utils.py
│
├── data/                  # 觀測或模擬資料
├── requirements.txt       # Python 相依套件
├── .github/workflows/     # CI/CD 設定
└── README.md
```

---

## 快速開始

### 基本使用範例

```python
import numpy as np
from src.classical import solve_kepler, orbital_elements_to_state, orbital_period
from src.relativity import pericenter_precession, mercury_precession_test
from src.utils import deg_to_rad, AU, solar_mass, G

# 地球軌道參數
a_earth = 1.0 * AU  # 半長軸
e_earth = 0.0167    # 離心率
mu_sun = G * solar_mass

# 計算軌道週期
T = orbital_period(a_earth, mu_sun)
print(f"地球軌道週期: {T / (365.25 * 24 * 3600):.3f} 年")

# 水星近日點進動測試
mercury_test = mercury_precession_test()
print(f"水星GR進動: {mercury_test['calculated_gr_precession_per_year']*100:.1f} 角秒/世紀")
```

### 範例一：經典二體問題與開普勒方程式

```python
# src/classical.py 主要功能演示
import numpy as np
from src.classical import solve_kepler, orbital_elements_to_state

# 求解開普勒方程式
M = np.pi / 3  # 平均近點角 (60度)
e = 0.5        # 離心率
E = solve_kepler(M, e)
print(f"偏近點角 E = {E:.4f} rad")

# 軌道要素轉狀態向量
a = 1.0 * AU
mu = G * solar_mass
pos, vel = orbital_elements_to_state(a, e, M, mu)
print(f"位置: {pos/AU} AU")
print(f"速度: {vel/1000} km/s")
```

### 範例二：廣義相對論進動修正

```python
# src/relativity.py 主要功能演示
from src.relativity import pericenter_precession, shapiro_delay

# 水星軌道參數
a_mercury = 5.791e10  # m
e_mercury = 0.2056
M_sun = solar_mass

# 計算近日點進動
precession = pericenter_precession(a_mercury, e_mercury, M_sun)
print(f"每軌道進動角: {np.degrees(precession)*3600:.6f} 角秒")

# Shapiro 時延示例
r1, r2 = 1.0*AU, 0.7*AU  # 地球-金星距離
r_closest = 7e8  # 最接近太陽距離
delay = shapiro_delay(r1, r2, r_closest, M_sun)
print(f"Shapiro 延遲: {delay*1e6:.1f} 微秒")
```

---

## Jupyter 筆記本示例

### 1. kepler.ipynb - 開普勒問題演示

- 載入 `src/classical.py`
- 設定地球–太陽系參數
- 演示開普勒方程式求解
- 畫出橢圓軌道位置
- 太陽系行星軌道比較
- 霍曼轉移軌道計算

### 2. precession.ipynb - 進動效應分析

- 載入 `src/relativity.py`
- 以水星軌道為例，計算近日點進動
- 與觀測值比較（~43角秒/世紀）
- 可視化進動角隨時間之變化
- 光線偏折與 Shapiro 延遲計算

---

## 測試與品質保證

### 執行測試

```bash
# 執行所有測試
pytest tests/ -v

# 執行特定測試模組
pytest tests/test_classical.py -v

# 產生覆蓋率報告
pytest --cov=src --cov-report=html
```

### 測試範圍

- **經典力學測試** (`test_classical.py`)
  - 開普勒方程式求解精度
  - 軌道要素與狀態向量轉換
  - 能量與角動量守恆
  - 開普勒第三定律驗證

- **相對論測試** (`test_relativity.py`)
  - 水星進動計算準確性
  - 光線偏折公式驗證
  - Schwarzschild 半徑計算
  - 物理一致性檢查

- **工具函式測試** (`test_utils.py`)
  - 單位轉換精度
  - 物理常數正確性
  - 數學工具函式
  - 天文計算驗證

---

## 可視化範例

### 軌道可視化

```python
import matplotlib.pyplot as plt
from src.classical import orbital_elements_to_state

# 繪製不同離心率的軌道
eccentricities = [0.0, 0.2, 0.5, 0.8]
plt.figure(figsize=(10, 8))

for e in eccentricities:
    M_vals = np.linspace(0, 2*np.pi, 200)
    x_coords, y_coords = [], []
    
    for M in M_vals:
        pos, _ = orbital_elements_to_state(AU, e, M, G*solar_mass)
        x_coords.append(pos[0]/AU)
        y_coords.append(pos[1]/AU)
    
    plt.plot(x_coords, y_coords, label=f'e = {e}', linewidth=2)

plt.scatter([0], [0], color='orange', s=200, label='太陽')
plt.axis('equal')
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.title('不同離心率的橢圓軌道')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### 進動效應演示

```python
from src.relativity import pericenter_precession_per_year

# 比較不同行星的GR進動
planets = {
    '水星': {'a': 0.387*AU, 'e': 0.2056},
    '金星': {'a': 0.723*AU, 'e': 0.0067},
    '地球': {'a': 1.000*AU, 'e': 0.0167},
    '火星': {'a': 1.524*AU, 'e': 0.0934},
}

for name, data in planets.items():
    T = orbital_period(data['a'], G*solar_mass)
    precession = pericenter_precession_per_year(data['a'], data['e'], solar_mass, T)
    print(f"{name}: {precession*100:.3f} 角秒/世紀")
```

---

## 持續整合與部署

GitHub Actions 自動化流程：

```yaml
# .github/workflows/ci.yml
name: CI
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ['3.8', '3.9', '3.10', '3.11', '3.12']
    steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v4
      with: python-version: ${{ matrix.python-version }}
    - run: pip install -r requirements.txt
    - run: pytest --maxfail=1 --disable-warnings -q --cov=src
```

---

## 未來拓展

- **多體問題**：N-體模擬、Barnes–Hut 或螺旋法加速
- **數值積分器**：Runge–Kutta、Symplectic integrator 比較
- **更深入相對論**：利用 `EinsteinPy` 或自行實作 GR 數值积分
- **資料驅動**：整合實際天文觀測資料，進行參數擬合
- **互動式介面**：Web 應用程式或 GUI 工具

---

## 參考資料

- **Astropy 官方文件**：[https://docs.astropy.org/](https://docs.astropy.org)
- **SciPy Optimize 模組**：[https://docs.scipy.org/doc/scipy/reference/optimize.html](https://docs.scipy.org/doc/scipy/reference/optimize.html)
- **經典力學與廣義相對論教材**：
  - [Physics Books Collection](https://github.com/manjunath5496/Physics-Books)
  - [Physics Books Repository](https://github.com/rgjha/PhysicsBooks)

---

## 授權與貢獻

本專案採用開源授權，歡迎提交 issue 和 pull request 改進程式碼。

**注意**：本程式碼僅供教育與演示用途，未經過嚴格的科學驗證，請勿用於實際的太空任務計算。
