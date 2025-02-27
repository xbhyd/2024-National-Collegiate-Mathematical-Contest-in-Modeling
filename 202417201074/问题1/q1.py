import pandas as pd
import sympy as sm
import math


def length(x1, y1, x2, y2):
    return sm.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def f(theta, X, Y, D):
    return sm.sqrt((a[0] * theta * sm.cos(theta) - X) ** 2 + (a[0] * theta * sm.sin(theta) - Y) ** 2) - D


# 定义符号
x = sm.symbols("x")  # 方位角
b = sm.symbols("b")  # 积分下限，第t秒的方位角
a = sm.symbols("a")  # 螺线参数

# 求解螺线参数 a
a = sm.solve(2 * sm.pi * a - 0.55, a)

# 指定 Excel 文件的路径
file_path = 'result1.xlsx'

# 使用 pandas 读取 Excel 文件
df_1 = pd.read_excel(file_path, index_col=0, sheet_name="位置")
df_2 = pd.read_excel(file_path, index_col=0, sheet_name="速度")
df_3 = df_1.copy()  # 存放第-0.0001秒的运动坐标
df_4 = df_1.copy()  # 存放第0.0001秒的运动坐标

# 先计算0-300s的龙头位置
# 循环遍历 t 值，并计算对应的积分下限 b
for t in range(0, 301):
    print('第' + str(t) + 's计算')
    # 进行积分，传递正确的积分上下限
    fx = sm.integrate(a[0] * sm.sqrt(x ** 2 + 1), (x, b, sm.pi * 32)) - t
    # 解方程，求解 b
    res = sm.nsolve(fx, b, 16)
    # 求解 r x y
    r = res * a[0]
    x_t = r * sm.cos(res)
    y_t = r * sm.sin(res)
    df_1.loc['龙头x (m)', str(t) + ' s'] = x_t
    df_1.loc['龙头y (m)', str(t) + ' s'] = y_t
    flag = 0  # 判断是不是龙头的标志
    theta0 = res
    while r < a[0] * 32 * sm.pi:
        # 定义距离差
        if flag == 0:
            delta_d = 2.86
            theta1 = theta0  # 初始角度，可以根据需要调整
            flag += 1
            x0 = x_t
            y0 = y_t
        else:
            delta_d = 1.65
            theta1 = theta0
            flag += 1
        delta = 0.5
        # 牛顿法求解
        while abs(f(theta1, x0, y0, delta_d)) > 1e-8:
            if f(theta1, x0, y0, delta_d) > 0:
                theta1 -= delta
                delta /= 2
            else:
                theta1 += delta
        # 计算相邻点的极坐标
        r1 = a[0] * theta1
        x1 = r1 * sm.cos(theta1)
        y1 = r1 * sm.sin(theta1)
        if r1 <= a[0] * 32 * sm.pi:
            if flag <= 221:
                df_1.loc['第' + str(flag) + '节龙身x (m)', str(t) + ' s'] = x1
                df_1.loc['第' + str(flag) + '节龙身y (m)', str(t) + ' s'] = y1
            elif flag == 222:
                df_1.loc['龙尾x (m)', str(t) + ' s'] = x1
                df_1.loc['龙尾y (m)', str(t) + ' s'] = y1
            elif flag == 223:
                df_1.loc['龙尾（后）x (m)', str(t) + ' s'] = x1
                df_1.loc['龙尾（后）y (m)', str(t) + ' s'] = y1
            else:
                break
        x0 = r1 * sm.cos(theta1)
        y0 = r1 * sm.sin(theta1)
        r = r1
        theta0 = theta1

df_2.loc[['龙头 (m/s)'], :] = 1.000000

# 下面开始计算-0.0001s的（点是否在螺线内不管，在前面df_1中已经控制）
# 这里只是为了计算-0.0001s时各点位置便于速度计算，最终某时刻某点的速度是否存在以df_1中对应时间的该点是否有坐标来判断
for t in range(0, 301):
    print('第' + str(t) + 's计算')
    # 进行积分，传递正确的积分上下限
    fx = sm.integrate(a[0] * sm.sqrt(x ** 2 + 1), (x, b, sm.pi * 32)) - t + 0.0001  # 0.0001s为设定的间隔时间
    # 解方程，求解 b
    res = sm.nsolve(fx, b, 16)
    # 求解 r x y
    r = res * a[0]
    x_t = r * sm.cos(res)
    y_t = r * sm.sin(res)
    df_3.loc['龙头x (m)', str(t) + ' s'] = x_t
    df_3.loc['龙头y (m)', str(t) + ' s'] = y_t
    flag = 0  # 判断是不是龙头的标志
    theta0 = res
    while r < a[0] * 32 * sm.pi + 0.55:  # 这里0.55为手动添加的微小扰动，防止某点在正秒时在螺线范围内，但0.0001s前在螺线外
        # 定义距离差
        if flag == 0:
            delta_d = 2.86
            theta1 = theta0  # 初始角度，可以根据需要调整
            flag += 1
            x0 = x_t
            y0 = y_t
        else:
            delta_d = 1.65
            theta1 = theta0
            flag += 1
        delta = 0.5
        # 牛顿法求解
        while abs(f(theta1, x0, y0, delta_d)) > 1e-8:
            if f(theta1, x0, y0, delta_d) > 0:
                theta1 -= delta
                delta /= 2
            else:
                theta1 += delta
        # 计算相邻点的极坐标
        r1 = a[0] * theta1
        x1 = r1 * sm.cos(theta1)
        y1 = r1 * sm.sin(theta1)
        if r1 <= a[0] * 32 * sm.pi + 0.55:  # 这里的微小扰动原理同上
            if flag <= 221:
                df_3.loc['第' + str(flag) + '节龙身x (m)', str(t) + ' s'] = x1
                df_3.loc['第' + str(flag) + '节龙身y (m)', str(t) + ' s'] = y1
            elif flag == 222:
                df_3.loc['龙尾x (m)', str(t) + ' s'] = x1
                df_3.loc['龙尾y (m)', str(t) + ' s'] = y1
            elif flag == 223:
                df_3.loc['龙尾（后）x (m)', str(t) + ' s'] = x1
                df_3.loc['龙尾（后）y (m)', str(t) + ' s'] = y1
            else:
                break
        x0 = r1 * sm.cos(theta1)
        y0 = r1 * sm.sin(theta1)
        r = r1
        theta0 = theta1

# 下面开始计算0.0001s的
for t in range(0, 301):
    print('第' + str(t) + 's计算')
    # 进行积分，传递正确的积分上下限
    fx = sm.integrate(a[0] * sm.sqrt(x ** 2 + 1), (x, b, sm.pi * 32)) - t - 0.0001  # 0.0001s为设定的间隔时间
    # 解方程，求解 b
    res = sm.nsolve(fx, b, 16)
    # 求解 r x y
    r = res * a[0]
    x_t = r * sm.cos(res)
    y_t = r * sm.sin(res)
    df_4.loc['龙头x (m)', str(t) + ' s'] = x_t
    df_4.loc['龙头y (m)', str(t) + ' s'] = y_t
    flag = 0  # 判断是不是龙头的标志
    theta0 = res
    while r < a[0] * 32 * sm.pi + 0.55:  # 这里0.55为手动添加的微小扰动，防止某点在正秒时在螺线范围内，但0.0001s前在螺线外
        # 定义距离差
        if flag == 0:
            delta_d = 2.86
            theta1 = theta0  # 初始角度，可以根据需要调整
            flag += 1
            x0 = x_t
            y0 = y_t
        else:
            delta_d = 1.65
            theta1 = theta0
            flag += 1
        delta = 0.5
        # 牛顿法求解
        while abs(f(theta1, x0, y0, delta_d)) > 1e-8:
            if f(theta1, x0, y0, delta_d) > 0:
                theta1 -= delta
                delta /= 2
            else:
                theta1 += delta
        # 计算相邻点的极坐标
        r1 = a[0] * theta1
        x1 = r1 * sm.cos(theta1)
        y1 = r1 * sm.sin(theta1)
        if r1 <= a[0] * 32 * sm.pi + 0.55:  # 这里的微小扰动原理同上
            if flag <= 221:
                df_4.loc['第' + str(flag) + '节龙身x (m)', str(t) + ' s'] = x1
                df_4.loc['第' + str(flag) + '节龙身y (m)', str(t) + ' s'] = y1
            elif flag == 222:
                df_4.loc['龙尾x (m)', str(t) + ' s'] = x1
                df_4.loc['龙尾y (m)', str(t) + ' s'] = y1
            elif flag == 223:
                df_4.loc['龙尾（后）x (m)', str(t) + ' s'] = x1
                df_4.loc['龙尾（后）y (m)', str(t) + ' s'] = y1
            else:
                break
        x0 = r1 * sm.cos(theta1)
        y0 = r1 * sm.sin(theta1)
        r = r1
        theta0 = theta1

# 开始用df_1、df_3、df_4对应点的坐标距离差分计算第1节龙身到221节龙身速度
for i in range(1, 222):
    for t in range(0, 301):
        if math.isnan(df_1.at['第' + str(i) + '节龙身x (m)', str(t) + ' s']) is not True:
            x_a = df_1.at['第' + str(i) + '节龙身x (m)', str(t) + ' s']  # 某时刻某点的x坐标
            y_a = df_1.at['第' + str(i) + '节龙身y (m)', str(t) + ' s']  # 某时刻某点的y坐标
            x_b = df_3.at['第' + str(i) + '节龙身x (m)', str(t) + ' s']  # 某时刻某点前0.0001s的x坐标
            y_b = df_3.at['第' + str(i) + '节龙身y (m)', str(t) + ' s']  # 某时刻某点前0.0001s的y坐标
            x_c = df_4.at['第' + str(i) + '节龙身x (m)', str(t) + ' s']  # 某时刻某点后0.0001s的x坐标
            y_c = df_4.at['第' + str(i) + '节龙身y (m)', str(t) + ' s']  # 某时刻某点后0.0001s的y坐标
            D1 = length(x_a, y_a, x_b, y_b)
            D2 = length(x_a, y_a, x_c, y_c)
            v = (D2 + D1) / 0.0002  # 某时刻某点的速度
            df_2.loc['第' + str(i) + '节龙身  (m/s)', str(t) + ' s'] = v

# 计算两个龙尾速度
for t in range(0, 301):
    if math.isnan(df_1.at['龙尾x (m)', str(t) + ' s']) is not True:
        x_a = df_1.at['龙尾x (m)', str(t) + ' s']  # 某时刻前龙尾的x坐标
        y_a = df_1.at['龙尾y (m)', str(t) + ' s']  # 某时刻前龙尾的y坐标
        x_b = df_3.at['龙尾x (m)', str(t) + ' s']  # 某时刻前龙尾前0.0001s的x坐标
        y_b = df_3.at['龙尾y (m)', str(t) + ' s']  # 某时刻前龙尾前0.0001s的y坐标
        x_c = df_4.at['龙尾x (m)', str(t) + ' s']  # 某时刻前龙尾后0.0001s的x坐标
        y_c = df_4.at['龙尾y (m)', str(t) + ' s']  # 某时刻前龙尾后0.0001s的y坐标
        D1 = length(x_a, y_a, x_b, y_b)
        D2 = length(x_a, y_a, x_c, y_c)
        v = (D2 + D1) / 0.0002  # 某时刻前龙尾的速度
        df_2.loc['龙尾  (m/s)', str(t) + ' s'] = v
for t in range(0, 301):
    if math.isnan(df_1.at['龙尾（后）x (m)', str(t) + ' s']) is not True:
        x_a = df_1.at['龙尾（后）x (m)', str(t) + ' s']  # 某时刻后龙尾的x坐标
        y_a = df_1.at['龙尾（后）y (m)', str(t) + ' s']  # 某时刻后龙尾的y坐标
        x_b = df_3.at['龙尾（后）x (m)', str(t) + ' s']  # 某时刻后龙尾前0.0001s的x坐标
        y_b = df_3.at['龙尾（后）y (m)', str(t) + ' s']  # 某时刻后龙尾前0.0001s的y坐标
        x_c = df_4.at['龙尾（后）x (m)', str(t) + ' s']  # 某时刻后龙尾后0.0001s的x坐标
        y_c = df_4.at['龙尾（后）y (m)', str(t) + ' s']  # 某时刻后龙尾后0.0001s的y坐标
        D1 = length(x_a, y_a, x_b, y_b)
        D2 = length(x_a, y_a, x_c, y_c)
        v = (D2 + D1) / 0.0002  # 某时刻后龙尾的速度
        df_2.loc['龙尾（后） (m/s)', str(t) + ' s'] = v

# 写入文件
with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
    df_1.to_excel(writer, sheet_name="位置", float_format="%.6f")
    df_2.to_excel(writer, sheet_name="速度", float_format="%.6f")
