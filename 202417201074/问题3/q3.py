import sympy as sm


def sum_up(l):
    add_all = 0
    for i in l:
        add_all += i
    if add_all > 0:
        return False
    else:
        return True


# 根据板凳两个把手获取该板凳的四个顶点以便于后续相交分析
def get_point_of_rectangle(x0, y0, x1, y1, delta):
    p1 = (x0 + 0.15 * (y1 - y0) / delta + 0.275 * (x0 - x1) / delta,
          y0 + 0.15 * (x0 - x1) / delta + 0.275 * (y0 - y1) / delta)
    p2 = (x1 + 0.15 * (y1 - y0) / delta - 0.275 * (x0 - x1) / delta,
          y1 + 0.15 * (x0 - x1) / delta - 0.275 * (y0 - y1) / delta)
    p3 = (x1 - 0.15 * (y1 - y0) / delta - 0.275 * (x0 - x1) / delta,
          y1 - 0.15 * (x0 - x1) / delta - 0.275 * (y0 - y1) / delta)
    p4 = (x0 - 0.15 * (y1 - y0) / delta + 0.275 * (x0 - x1) / delta,
          y0 - 0.15 * (x0 - x1) / delta + 0.275 * (y0 - y1) / delta)
    return [p1, p2, p3, p4]


def f(theta, X, Y, D):
    return sm.sqrt((a * theta * sm.cos(theta) - X) ** 2 + (a * theta * sm.sin(theta) - Y) ** 2) - D


# 判断点是否在四边形内（射线相交法）
def point_in_irregular_quadrilateral(point, vertices):
    x, y = point
    num_vertices = len(vertices)
    num_crossings = 0
    for i in range(num_vertices):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % num_vertices]
        if (y1 > y) != (y2 > y) and x < ((x2 - x1) * (y - y1) / (y2 - y1) + x1):
            num_crossings += 1
    return num_crossings % 2 == 1


def calculate_intersection_true_or_false(a):
    time = 300  # 初始化时间
    delta_time = 50  # 二分法时间变化
    x = sm.symbols("x")  # 方位角
    b = sm.symbols("b")  # 积分下限，第t秒的方位角
    while abs(delta_time >= 1e-8):
        results = []
        # 进行积分，传递正确的积分上下限
        fx = sm.integrate(a * sm.sqrt(x ** 2 + 1), (x, b, sm.pi * 32)) - time
        # 解方程，求解 b
        theta0 = sm.nsolve(fx, b, 16)
        # 求解 r x y
        r = theta0 * a
        x_t1 = r * sm.cos(theta0)
        y_t1 = r * sm.sin(theta0)
        flag = 0
        while flag <= 100:
            theta1 = theta0  # 初始角度，可以根据需要调整
            # 定义距离差
            if flag == 0:
                delta_d = 2.86
                x0 = x_t1
                y0 = y_t1
            else:
                delta_d = 1.65
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
            r1 = a * theta1
            x1 = r1 * sm.cos(theta1)
            y1 = r1 * sm.sin(theta1)
            if flag == 1:
                x_t2 = x1
                y_t2 = y1
                [p1, p2, p3, p4] = get_point_of_rectangle(x_t1, y_t1, x_t2, y_t2, delta_d)
            else:
                [pt1, pt2, pt3, pt4] = get_point_of_rectangle(x0, y0, x1, y1, delta_d)
                result1 = point_in_irregular_quadrilateral(p1, [pt1, pt2, pt3, pt4])
                result2 = point_in_irregular_quadrilateral(p2, [pt1, pt2, pt3, pt4])
                result3 = point_in_irregular_quadrilateral(p3, [pt1, pt2, pt3, pt4])
                result4 = point_in_irregular_quadrilateral(p4, [pt1, pt2, pt3, pt4])
                if result1 is not True and result2 is not True and result3 is not True and result4 is not True:
                    results.append(0)
                else:
                    if flag == 2:
                        results.append(0)  # 龙头和第一节龙身一定相交，为方便记作0
                    else:
                        results.append(1)  # 相交为1
            x0 = r1 * sm.cos(theta1)
            y0 = r1 * sm.sin(theta1)
            theta0 = theta1
        print(delta_time)
        if sum_up(results):
            time += delta_time
        else:
            time -= delta_time
            delta_time /= 2

    print("time:", time)
    # 定义符号
    xp = sm.symbols("x")  # 方位角
    bp = sm.symbols("b")  # 积分下限，第t秒的方位角
    print('第' + str(time) + 's计算')
    # 进行积分，传递正确的积分上下限
    fxp = sm.integrate(a * sm.sqrt(xp ** 2 + 1), (xp, bp, sm.pi * 32)) - time
    # 解方程，求解 b
    resp = sm.nsolve(fxp, bp, 16)
    # 求解 r x y
    rp = resp * a
    print("极轴:", rp)
    return rp


p = 0.55
delta_p = 0.05
# 定义符号螺线参数
a = sm.symbols("a")
# 求解螺线参数 a
a = sm.solve(2 * sm.pi * a - p, a)[0]
distance = calculate_intersection_true_or_false(a)
while abs(distance - 4.5) > 1e-8:
    print("p:", p)
    if distance <= 4.5:
        p -= delta_p
        # 定义符号螺线参数
        a = sm.symbols("a")
        # 求解螺线参数 a
        a = sm.solve(2 * sm.pi * a - p, a)[0]
        distance = calculate_intersection_true_or_false(a)
    else:
        p += delta_p
        delta_p /= 2
        # 定义符号螺线参数
        a = sm.symbols("a")
        # 求解螺线参数 a
        a = sm.solve(2 * sm.pi * a - p, a)[0]
        distance = calculate_intersection_true_or_false(a)

print(p)
