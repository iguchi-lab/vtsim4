# 温熱計算プログラムvtsim3の使い方

---
## 1.はじめに

  vtsimは、建築環境工学・建築設備を学ぶ大学生、研究者向けに開発された温熱計算プログラムで下記の特徴がある。

*   ブラウザ上（Google Clab）で動作する。環境設定が不要。
*   自由度の高いプログラム言語pythonで記述。
*   速度が求められる計算部分はc++で記述。
*   ユーザーは、プログラム言語pythonに加え、pandas、numpyなどの知識が必要。
*   熱・換気回路網による計算をベースとし、節点（ノード）と回路網（ネットワーク）、各種条件の設定で動作。

---
## 2.インストール方法

  Google Colabで下記のコードを実行して、インストールする。

```
!pip install git+https://github.com/iguchi-lab/vtsim3
from vtsim import vtsim as vt
```
  あわせて下記も実行しておくことを便利である。

  ※matplotlibは、グラフをの描画に必要。japanize-matplotlibで、日本語にも対応できる。

  ※numpyとpandasは、行列などを扱う際に必要。

```
!pip install japanize-matplotlib
import matplotlib.pyplot as plt
import japanize_matplotlib

import numpy as np
import pandas as pd
```

---
## 3.vtsim関連モジュールの構成

  vtsimは、下記のモジュールで構成される。

- vtsim

  vtsim本体の計算機能
- archenvlib

  各種物理定数の他、湿り空気、太陽位置、直散分離、快適性などを計算
- fan_sts

  送風ファンのPQ特性
- cof_ground

  地盤の応答係数

---
## 4.vtsimでできること

ノードとネットワークの構成により、下記のような様々な計算が可能。

- 換気回路網計算による、多数室、部位間の空気の移動
- 熱回路網計算による、多数室、部位間の温度や熱流
- 空気と共に移動する物質の室内濃度

![ノードとネットワークの設定例](sample01.png)
