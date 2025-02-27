import pandas as pd
import sympy as sm
import math


def length(x1, y1, x2, y2):
    return sm.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def f(theta, X, Y, D):
    return sm.sqrt((a[0] * theta * sm.cos(theta) - X) ** 2 + (a[0] * theta * sm.sin(theta) - Y) ** 2) - D


file_path = 'result2.xlsx'
time = 412.47383765876293  # q2(1).py运行出的盘入终止时间

# 定义符号
x = sm.symbols("x")  # 方位角
b = sm.symbols("b")  # 积分下限，第t秒的方位角
a = sm.symbols("a")  # 螺线参数

# 求解螺线参数 a
a = sm.solve(2 * sm.pi * a - 0.55, a)

df = pd.read_excel(file_path, index_col=0, sheet_name='Sheet1')
df_2 = df.copy()  # 存放第-0.0001秒的运动坐标
df_3 = df.copy()  # 存放第0.0001秒的运动坐标

print('第' + str(time) + 's计算')
# 进行积分，传递正确的积分上下限
fx = sm.integrate(a[0] * sm.sqrt(x ** 2 + 1), (x, b, sm.pi * 32)) - time
# 解方程，求解 b
res = sm.nsolve(fx, b, 16)
# 求解 r x y
r = res * a[0]
print('半径：', r)
x_t = r * sm.cos(res)
y_t = r * sm.sin(res)
df.loc['龙头', '横坐标x (m)'] = x_t
df.loc['龙头', '纵坐标y (m)'] = y_t
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
        print(flag)
        if flag <= 221:
            df.loc['第' + str(flag) + '节龙身', '横坐标x (m)'] = x1
            df.loc['第' + str(flag) + '节龙身', '纵坐标y (m)'] = y1
        elif flag == 222:
            df.loc['龙尾', '横坐标x (m)'] = x1
            df.loc['龙尾', '纵坐标y (m)'] = y1
        elif flag == 223:
            df.loc['龙尾（后）', '横坐标x (m)'] = x1
            df.loc['龙尾（后）', '纵坐标y (m)'] = y1
        else:
            break
    x0 = r1 * sm.cos(theta1)
    y0 = r1 * sm.sin(theta1)
    r = r1
    theta0 = theta1

# 显示数据框内容
print(df)

# 下面开始计算-0.0001s的（点是否在螺线内不管，在前面df中已经控制）
# 这里只是为了计算-0.0001s时各点位置便于速度计算，最终某时刻某点的速度是否存在以df中对应时间的该点是否有坐标来判断
print('第' + str(time) + 's计算')
# 进行积分，传递正确的积分上下限
fx = sm.integrate(a[0] * sm.sqrt(x ** 2 + 1), (x, b, sm.pi * 32)) - time + 0.0001  # 0.0001s为设定的间隔时间
# 解方程，求解 b
res = sm.nsolve(fx, b, 16)
# 求解 r x y
r = res * a[0]
x_t = r * sm.cos(res)
y_t = r * sm.sin(res)
df_2.loc['龙头', '横坐标x (m)'] = x_t
df_2.loc['龙头', '纵坐标y (m)'] = y_t
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
        print(flag)
        if flag <= 221:
            df_2.loc['第' + str(flag) + '节龙身', '横坐标x (m)'] = x1
            df_2.loc['第' + str(flag) + '节龙身', '纵坐标y (m)'] = y1
        elif flag == 222:
            df_2.loc['龙尾', '横坐标x (m)'] = x1
            df_2.loc['龙尾', '纵坐标y (m)'] = y1
        elif flag == 223:
            df_2.loc['龙尾（后）', '横坐标x (m)'] = x1
            df_2.loc['龙尾（后）', '纵坐标y (m)'] = y1
        else:
            break
    x0 = r1 * sm.cos(theta1)
    y0 = r1 * sm.sin(theta1)
    r = r1
    theta0 = theta1

# 下面开始计算0.0001s的（点是否在螺线内不管，在前面df中已经控制）
# 这里只是为了计算0.0001s时各点位置便于速度计算，最终某时刻某点的速度是否存在以df中对应时间的该点是否有坐标来判断
print('第' + str(time) + 's计算')
# 进行积分，传递正确的积分上下限
fx = sm.integrate(a[0] * sm.sqrt(x ** 2 + 1), (x, b, sm.pi * 32)) - time - 0.0001  # 0.0001s为设定的间隔时间
# 解方程，求解 b
res = sm.nsolve(fx, b, 16)
# 求解 r x y
r = res * a[0]
x_t = r * sm.cos(res)
y_t = r * sm.sin(res)
df_3.loc['龙头', '横坐标x (m)'] = x_t
df_3.loc['龙头', '纵坐标y (m)'] = y_t
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
        print(flag)
        if flag <= 221:
            df_3.loc['第' + str(flag) + '节龙身', '横坐标x (m)'] = x1
            df_3.loc['第' + str(flag) + '节龙身', '纵坐标y (m)'] = y1
        elif flag == 222:
            df_3.loc['龙尾', '横坐标x (m)'] = x1
            df_3.loc['龙尾', '纵坐标y (m)'] = y1
        elif flag == 223:
            df_3.loc['龙尾（后）', '横坐标x (m)'] = x1
            df_3.loc['龙尾（后）', '纵坐标y (m)'] = y1
        else:
            break
    x0 = r1 * sm.cos(theta1)
    y0 = r1 * sm.sin(theta1)
    r = r1
    theta0 = theta1

df.loc['龙头', '速度 (m/s)'] = 1.000000
# 开始用df、df_2、df_3对应点的坐标距离差分计算第1节龙身到221节龙身速度
for i in range(1, 222):
    if math.isnan(df.at['第' + str(i) + '节龙身', '横坐标x (m)']) is not True:
        x_a = df.at['第' + str(i) + '节龙身', '横坐标x (m)']  # 某点的x坐标
        y_a = df.at['第' + str(i) + '节龙身', '纵坐标y (m)']  # 某点的y坐标
        x_b = df_2.at['第' + str(i) + '节龙身', '横坐标x (m)']  # 某点前0.0001s的x坐标
        y_b = df_2.at['第' + str(i) + '节龙身', '纵坐标y (m)']  # 某点前0.0001s的y坐标
        x_c = df_3.at['第' + str(i) + '节龙身', '横坐标x (m)']  # 某点后0.0001s的x坐标
        y_c = df_3.at['第' + str(i) + '节龙身', '纵坐标y (m)']  # 某点后0.0001s的y坐标
        D1 = length(x_a, y_a, x_b, y_b)
        D2 = length(x_a, y_a, x_c, y_c)
        v = (D2 + D1) / 0.0002  # 某点的速度
        df.loc['第' + str(i) + '节龙身', '速度 (m/s)'] = v

# 计算两个龙尾速度
# 前龙尾
if math.isnan(df.at['龙尾', '横坐标x (m)']) is not True:
    x_a = df.at['龙尾', '横坐标x (m)']  # 前龙尾的x坐标
    y_a = df.at['龙尾', '纵坐标y (m)']  # 前龙尾的y坐标
    x_b = df_2.at['龙尾', '横坐标x (m)']  # 前龙尾前0.0001s的x坐标
    y_b = df_2.at['龙尾', '纵坐标y (m)']  # 前龙尾前0.0001s的y坐标
    x_c = df_3.at['龙尾', '横坐标x (m)']  # 前龙尾后0.0001s的x坐标
    y_c = df_3.at['龙尾', '纵坐标y (m)']  # 前龙尾后0.0001s的y坐标
    D1 = length(x_a, y_a, x_b, y_b)
    D2 = length(x_a, y_a, x_c, y_c)
    v = (D2 + D1) / 0.0002  # 某时刻前龙尾的速度
    df.loc['龙尾', '速度 (m/s)'] = v

# 后龙尾
if math.isnan(df.at['龙尾（后）', '横坐标x (m)']) is not True:
    x_a = df.at['龙尾（后）', '横坐标x (m)']  # 后龙尾的x坐标
    y_a = df.at['龙尾（后）', '纵坐标y (m)']  # 后龙尾的y坐标
    x_b = df_2.at['龙尾（后）', '横坐标x (m)']  # 后龙尾前0.0001s的x坐标
    y_b = df_2.at['龙尾（后）', '纵坐标y (m)']  # 后龙尾前0.0001s的y坐标
    x_c = df_3.at['龙尾（后）', '横坐标x (m)']  # 后龙尾后0.0001s的x坐标
    y_c = df_3.at['龙尾（后）', '纵坐标y (m)']  # 后龙尾后0.0001s的y坐标
    D1 = length(x_a, y_a, x_b, y_b)
    D2 = length(x_a, y_a, x_c, y_c)
    v = (D2 + D1) / 0.0002  # 某时刻后龙尾的速度
    df.loc['龙尾（后）', '速度 (m/s)'] = v

# 保存数据
df.to_excel(file_path, sheet_name="Sheet1")
