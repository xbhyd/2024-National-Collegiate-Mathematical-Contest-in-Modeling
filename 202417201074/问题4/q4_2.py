# 运行此代码前，在X.xlsx和Y.xlsx中第一行前插一行，用(-100,100)，步长0.1填充
import pandas as pd
import sympy as sm


# 计算两点距离
def length(x1, y1, x2, y2):
    return sm.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


file_path_1 = "X.xlsx"
file_path_2 = "Y.xlsx"
file_path = "result4.xlsx"
df_1 = pd.read_excel(file_path_1)
df_2 = pd.read_excel(file_path_2)
df_3 = pd.read_excel(file_path, index_col=0, sheet_name="位置")
df_4 = pd.read_excel(file_path, index_col=0, sheet_name="速度")

print(df_1)
print(df_2)
print(df_3)
print(df_4)

for t in range(-100, 101):
    for i in range(0, 224):
        if i == 0:
            df_3.loc['龙头x (m)', str(t) + ' s'] = df_1.iat[i, 10 * t + 1000]
            df_3.loc['龙头y (m)', str(t) + ' s'] = df_2.iat[i, 10 * t + 1000]
            df_4.loc['龙头 (m/s)', str(t) + ' s'] = 1.000000
        elif i == 222:
            df_3.loc['龙尾x (m)', str(t) + ' s'] = df_1.iat[i, 10 * t + 1000]
            df_3.loc['龙尾y (m)', str(t) + ' s'] = df_2.iat[i, 10 * t + 1000]
            if t == -100:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 10 * t + 1001]
                y_b = df_2.iat[i, 10 * t + 1001]
                D1 = length(x_a, y_a, x_b, y_b)
                v = D1 / 0.1
                df_4.loc['龙尾  (m/s)', str(t) + ' s'] = v
            elif t == 100:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 10 * t + 999]
                y_b = df_2.iat[i, 10 * t + 999]
                D1 = length(x_a, y_a, x_b, y_b)
                v = D1 / 0.1
                df_4.loc['龙尾  (m/s)', str(t) + ' s'] = v
            else:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 10 * t + 999]
                y_b = df_2.iat[i, 10 * t + 999]
                x_c = df_1.iat[i, 10 * t + 1001]
                y_c = df_2.iat[i, 10 * t + 1001]
                D1 = length(x_a, y_a, x_b, y_b)
                D2 = length(x_a, y_a, x_c, y_c)
                v = (D2 + D1) / 0.2  # 某时刻后龙尾的速度
                df_4.loc['龙尾  (m/s)', str(t) + ' s'] = v
        elif i == 223:
            df_3.loc['龙尾（后）x (m)', str(t) + ' s'] = df_1.iat[i, 10 * t + 1000]
            df_3.loc['龙尾（后）y (m)', str(t) + ' s'] = df_2.iat[i, 10 * t + 1000]
            if t == -100:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 1]
                y_b = df_2.iat[i, 1]
                D1 = length(x_a, y_a, x_b, y_b)
                v = D1 / 0.1
                df_4.loc['龙尾（后） (m/s)', str(t) + ' s'] = v
            elif t == 100:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 10 * t + 999]
                y_b = df_2.iat[i, 10 * t + 999]
                D1 = length(x_a, y_a, x_b, y_b)
                v = D1 / 0.1
                df_4.loc['龙尾（后） (m/s)', str(t) + ' s'] = v
            else:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 10 * t + 999]
                y_b = df_2.iat[i, 10 * t + 999]
                x_c = df_1.iat[i, 10 * t + 1001]
                y_c = df_2.iat[i, 10 * t + 1001]
                D1 = length(x_a, y_a, x_b, y_b)
                D2 = length(x_a, y_a, x_c, y_c)
                v = (D2 + D1) / 0.2
                df_4.loc['龙尾（后） (m/s)', str(t) + ' s'] = v
        else:
            df_3.loc['第' + str(i) + '节龙身x (m)', str(t) + ' s'] = df_1.iat[i, 10 * t + 1000]
            df_3.loc['第' + str(i) + '节龙身y (m)', str(t) + ' s'] = df_2.iat[i, 10 * t + 1000]
            if t == -100:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 1]
                y_b = df_2.iat[i, 1]
                D1 = length(x_a, y_a, x_b, y_b)
                v = D1 / 0.1
                df_4.loc['第' + str(i) + '节龙身  (m/s)', str(t) + ' s'] = v
            elif t == 100:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 10 * t + 999]
                y_b = df_2.iat[i, 10 * t + 999]
                D1 = length(x_a, y_a, x_b, y_b)
                v = D1 / 0.1
                df_4.loc['第' + str(i) + '节龙身  (m/s)', str(t) + ' s'] = v
            else:
                x_a = df_1.iat[i, 10 * t + 1000]
                y_a = df_2.iat[i, 10 * t + 1000]
                x_b = df_1.iat[i, 10 * t + 999]
                y_b = df_2.iat[i, 10 * t + 999]
                x_c = df_1.iat[i, 10 * t + 1001]
                y_c = df_2.iat[i, 10 * t + 1001]
                D1 = length(x_a, y_a, x_b, y_b)
                D2 = length(x_a, y_a, x_c, y_c)
                v = (D2 + D1) / 0.2
                df_4.loc['第' + str(i) + '节龙身  (m/s)', str(t) + ' s'] = v

# 写入文件
with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
    df_3.to_excel(writer, sheet_name="位置", float_format="%.6f")
    df_4.to_excel(writer, sheet_name="速度", float_format="%.6f")
