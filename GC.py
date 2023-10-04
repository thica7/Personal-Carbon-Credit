#from scipy.optimize import linprog
import numpy as np
from sympy import Symbol


def find_all_paths(graph, start, end, visited=[], path=[], distance=0):
    visited.append(start)
    path.append(start)
    
    if start == end:
        #print("Path:", path, "Distance:", distance)
        distanceset.append(distance)
        pathset.append(path)
    
    else:
        for i in range(len(graph[start])):
            if graph[start][i] != 0 and i not in visited:
                new_distance = distance + graph[start][i]
                find_all_paths(graph, i, end, visited[:], path[:], new_distance)
                
    #visited.remove(start)
    #path.remove(start)
    return pathset, distanceset


def main_path(start, end):
    # adjacency matrix representing the weighted graph
    adjacency_matrix = [
            [  0,   3,   0,   0,   0,   0,   0,   0,   0,   0],
            [  3,   0,   5,   0,   0,   0,   0,   0,   0,   0],
            [  0,   5,   0,   4,   0,   6,   0,   0,   0,   0],
            [  0,   0,   5,   0,   8,   0,   0,   0,   0,   0],
            [  0,   0,   0,   8,   0,   8,   0,   0,   7,   0],
            [  0,   0,   6,   0,   8,   0,   2,   0,   0,   0],
            [  0,   0,   0,   0,   0,   2,   0,   9,   0,   0],
            [  0,   0,   0,   0,   0,   0,   9,   0,   0,   0],
            [  0,   0,   0,   0,   7,   0,   0,   0,   0,  10],
            [  0,   0,   0,   0,   0,   0,   0,   0,  10,   0]]

    start_point = start  # Index of the start point
    end_point = end    # Index of the end point

    print("All possible paths and their corresponding distances:")
    return find_all_paths(adjacency_matrix, start_point, end_point)


#假设k为交通方式：0为私家车，1为出租车，2为公交车，3为步行

def beta(n):
    beta = [[0.031, 0.056, -3.34, -0.197],
            [0.0066, 0.113, -0.031, -0.144], 
            [0.011, 0.057, 2.521, 1.02]]
    income = [20000, 40000, 60000, 60000, 80000, 300000]
    
    if income[n] < 45999:
        user_beta = beta[0]
    elif income[n] > 68999:
        user_beta = beta[2]
    else:
        user_beta = beta[1]
    return user_beta

    

def CC(price_c, pathset, distanceset):
    CC_set = []
    User_type = ["Car", "Taxi", "Bus", "Walk"]

    #Emission CO2
    per_emission = [4.0, 2, 0.1, 0]

    for k in range(len(per_emission)):  #交通类型
        for i in range(len(pathset)):   #路径
            for n in range(6):          #用户1-6
                emission = per_emission[k] * distanceset[i]
                CC_mnpt = emission * price_c
                beta_n = beta(n)
                CC_set.append((beta_n[1]*CC_mnpt, 0, User_type[k], n, pathset[i]))
                #print("Carbon Cost for", User_type[k], "User", n, "path", pathset[i], "Distance:", distanceset[i], "is", CC_mnpt)
    
    return CC_set


#假定k为交通方式：0为私家车，1为出租车，2为公交车，3为步行

def TC(pathset, distanceset, v):
    TC_set = []
    for k in [TCs, TCc, TCg, TCb]: 
        TC_set += k(pathset, distanceset, v)
                         # TC_set = TC_set[:][:2] + (k[1],) + TC_set[:][2:]
    return TC_set    

def TCs(pathset, distanceset, v):
    TCs_set = []

    F_i = 8 #假设单位油价
    parking_fee = 2 #假设单位停车费用
    d_n = 15 #假设单位油量行驶里程，km/L
    v_f = 60 #假设畅通的速度，km/h

    for i in range(len(pathset)):
        for n in range(6):
            TCs = parking_fee + F_i * distanceset[i]/d_n + VOT(n) * (1/v + t_f(v)-t_r(v) + sigma()+(1/v-1/v_f))
            beta_n = beta(n)
            TCs_set.append((beta_n[0]*TCs, 1, "Car", n, pathset[i]))
            #print("Time Cost for Car User", n, "path", pathset[i], "is", TCs)

    return TCs_set


def TCc(pathset, distanceset, v):
    TCc_set = []

    F_i = 2.4 #假设单位公里油价
    F_0 =10 #假设起步价10元2.5km
    v_f = 60 #假设畅通的速度，km/h

    for i in range(len(pathset)):
        for n in range(6):
            TCc = F_0 + F_i * max(distanceset[i]-2.5, 0) + float(VOT(n) * (1/v + t_f(v)-t_r(v) + sigma()+(1/v-1/v_f)))
            beta_n = beta(n)
            TCc_set.append((beta_n[0]*TCc, 1, "Taxi", n, pathset[i]))
            #print("Time Cost for Taxi User", n, "path", pathset[i], "is", TCc)

    return TCc_set


def TCg(pathset, distanceset, v):
    TCg_set = []

    F_0 = 2 #假设票价2元
    v_f = 75 #假设畅通的速度，km/h

    for i in range(len(pathset)):
        for n in range(6):
            TCg = F_0 + VOT(n) * (1/v + t_f(v)-t_r(v) + sigma()+(1/v-1/v_f))
            beta_n = beta(n)
            TCg_set.append((beta_n[0]*TCg, 1, "Bus", n, pathset[i]))
            #print("Time Cost for Bus User", n, "path", pathset[i], "is", TCg)

    return TCg_set


def TCb(pathset, distanceset, v):
    TCb_set = []

    v_f = 60 #假设畅通的速度，km/h

    for i in range(len(pathset)):
        for n in range(6):
            TCb = VOT(n) * (1/v + t_f(v)-t_r(v) + sigma()+(1/v-1/v_f))
            beta_n = beta(n)
            TCb_set.append((beta_n[0]*TCb, 1, "Walk", n, pathset[i]))
            #print("Time Cost for Walk User", n, "path", pathset[i], "is", TCb)

    return TCb_set


def t_f(v):
    return 1 / v


def t_r(v):
    return 1 / v


def VOT(n):
    I = [20000, 40000, 60000, 60000, 80000, 300000]
    dy = 250 / 12
    hd = 9.16
    return I[n]/(dy * hd)


def sigma():
    return 1


def GC(start, end, price_c, v):
    #GC_i = betaT_m * TC_mnpt + betaC_m * CC_mnpt + betaR_m * RB_i + betaH_m * H_i

    pathset, distanceset = main_path(start, end)

    cc = CC(price_c, pathset, distanceset)
    tc = TC(pathset, distanceset, v)
    

    GC = {
        "TC": tc,
        "CC": cc,
        "GC": [1*t1[0]+0.6*t2[0] for t1, t2 in zip(tc, cc)]
    }

    #GC = CC(price_c, pathset, distanceset) + TCs(pathset, distanceset, v)
    #print(GC["TC"])
    # for i in range(len(GC["TC"])):
    #     print( "\n", GC["GC"][i])

    return GC


def departure_rate(GC, demand_tot, pathset):
    r_count = len(pathset)
    r = [[[1 for n in range(4)] for path in range(r_count)]for m in range(6)]
    #r = [30 for i in range(r_count * 6 * 4)]


def GC_min(start, end, v, price):
    price_c = price
    GC_idx = GC(start, end, price_c, v)
    GC_mini = []
    method_list = []
    for m in range(6):
        min = 10000
        i_CC = 0
        for i in range(8):
            if min >= GC_idx["GC"][6 * i + m]:
                min = GC_idx["GC"][6 * i + m]
                i_CC = i
        GC_mini.append((GC_idx["GC"][6 * i_CC + m],
                        GC_idx["CC"][6 * i_CC + m][2],
                        GC_idx["CC"][6 * i_CC + m][0],
                        m))
    print(GC_mini)
    return GC_mini

def price(start, end, v, amount, price_c):

    user_credit = GC_min(start, end, v, price_c)
    credit_used = 30 * sum(x[2] for x in user_credit)
    diff = amount - credit_used
    return diff

def bisection(start, end, v, amount, tol=0.0000001, max_iter=1000):
    a = -100
    b = 100

    if price(start, end, v, amount, a) * price(start, end, v, amount, b) > 0:
        print(price(start, end, v, amount, a), price(start, end, v, amount, b))
        raise ValueError("f(a) and f(b) must have different signs")
    
    iter_count = 0
    while abs(b - a) / 2.0 > tol and iter_count < max_iter:
        c = (a+b) / 2.0
        if price(start, end, v, amount, c) == 0:
            return c
        if price(start, end, v, amount, c) * price(start, end, v, amount, a) < 0:
            b = c
        else:
            a = c
        iter_count += 1
        print(iter_count, b, a)
    print(price(start, end, v, amount, (a + b) / 2.0))
    return (a + b) / 2.0

if __name__ == "__main__":
    pathset = []
    distanceset = []
    Credit_sum = 10000
    #demand_tot = 4000
    price_c = bisection(0, 7, 35, Credit_sum)
    print(price_c)

    #GC = GC(1, 7, price_c, 35)
    #for i in range(len(GC["TC"])):
        #print(GC["TC"][i], "\n", GC["CC"][i])