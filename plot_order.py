import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import json
import numpy as np

with open("results_iTEBD_chi70.json", 'r', encoding='utf-8') as file:
    # JSONを辞書としてロード
    data = json.load(file)

# 各行の最後の値を取り出す
last_values = []
for key, value in data.items():
    # スペース区切りで分割し、最後の値を取得
    last_value = float(value.split()[-1])
    last_values.append(last_value)

# 結果を表示
# for value in last_values:
#     print(value)
print(last_value)
# print(file)
magnetic_field = []
for i in range(51):
    magnetic_field.append(i/10)


plt.plot(magnetic_field,last_values, marker="o", linestyle="dashed")
plt.xlabel("B")
plt.ylabel("|S^z|^2")
plt.legend()
plt.savefig("S^z_itebd_chi70.png")
plt.show()
