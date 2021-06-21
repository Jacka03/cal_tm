import matplotlib.pyplot as plt
import numpy as np
import math


def choose(tem_list, cou=1):
    """
    根据tm标准差降序排序，返回前count个标准差小的分割方法
    :param tem_list:
    :param cou:
    :return:
    """
    tem_list = np.array(tem_list)
    tem_list = tem_list[np.lexsort(tem_list.T)]
    # count = 10  # 在二维数组中获取前count个tm标准差最小的
    # TODO 获取前cou个、还是前cou种、还是：
    # print(tem_list[:cou])
    # tem_list = np.delete(tem_list, -1, axis=1)  # 删除最后一列
    # print(tem_list[:cou])
    if len(tem_list) > cou:
        return tem_list[:cou]

    return tem_list


def show_w(x, y, head):
    plt.scatter(x, y)
    tem_str = "min:{:3f},max:{:3f},max-min={:3f},std={:4f}".format(min(y), max(y), max(y) - min(y), np.std(y))
    plt.legend([tem_str])
    plt.title(head)
    # plt.ylim(80, 90)
    plt.show()


def get_gene(path):
    data = ''
    with open(path, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip()
            data += line
    if data.islower():
        data = data.upper()
    return data


def cal_tm(temp_gene):
    """
    计算一小段基因(temp_gene)的tm
    :param temp_gene:
    :return: 这段基因的tm
    """
    AATT = ATTA = TAAT = CAGT = GTCA = CTGA = GACT = CGGC = GCCG = GGCC = 0
    for i in range(len(temp_gene) - 1):
        if (temp_gene[i:i + 2] == 'AA') | (temp_gene[i:i + 2] == 'TT'):
            AATT += 1
        elif temp_gene[i:i + 2] == 'AT':
            ATTA += 1
        elif temp_gene[i:i + 2] == 'TA':
            TAAT += 1
        elif (temp_gene[i:i + 2] == 'CA') | (temp_gene[i:i + 2] == 'TG'):
            CAGT += 1
        elif (temp_gene[i:i + 2] == 'GT') | (temp_gene[i:i + 2] == 'AC'):
            GTCA += 1
        elif (temp_gene[i:i + 2] == 'CT') | (temp_gene[i:i + 2] == 'AG'):
            CTGA += 1
        elif (temp_gene[i:i + 2] == 'GA') | (temp_gene[i:i + 2] == 'TC'):
            GACT += 1
        elif temp_gene[i:i + 2] == 'CG':
            CGGC += 1
        elif temp_gene[i:i + 2] == 'GC':
            GCCG += 1
        elif (temp_gene[i:i + 2] == 'GG') | (temp_gene[i:i + 2] == 'CC'):
            GGCC += 1

    H = AATT * (-7.9) + ATTA * (-7.2) + TAAT * (-7.2) + CAGT * (-8.5) + GTCA * (-8.4) + CTGA * (-7.8) + GACT * (
        -8.2) + CGGC * (-10.6) + GCCG * (-9.8) + GGCC * (-8) + 0.1 + 2.3
    S = AATT * (-22.2) + ATTA * (-20.4) + TAAT * (-21.3) + CAGT * (-22.7) + GTCA * (-22.4) + CTGA * (-21) + GACT * (
        -22.2) + CGGC * (-27.2) + GCCG * (-24.4) + GGCC * (-19.9) - 2.8 + 4.1 - 1.4
    # TODO 钠离子浓度需要重新设置
    return (H * 1000) / (S + 1.987 * math.log10((10 ** -4) / 4)) - 273.15 + 16.6 * math.log10(1.)

