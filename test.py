import math
import numpy as np
import pandas as pd
from main import show_tm, show_w, cal_all_tm, find_top_k


def get_gene(path):
    data = ''
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
    return (H * 1000) / (S + 1.987 * math.log10((10 ** -4) / 4)) - 273.15 + 16.6 * math.log10(1.2)


def cal_first_tm():
    """
    计算第一段与第二段的tm值并且计算其标准差
    :return:
    """
    tem_res = []
    for i in range(max_len - min_len):
        mid_cut = min_len + i
        fir_gene = gene[:mid_cut]  # 第一段基因
        fir_tm = cal_tm(fir_gene)
        for j in range(max_len - min_len):
            end_cut = mid_cut + min_len + j
            sec_gene = gene[mid_cut:end_cut]  # 第二段基因
            sec_tm = cal_tm(sec_gene)
            tem_res.append([mid_cut, fir_tm, end_cut, sec_tm, np.std([fir_tm, sec_tm])])
    return tem_res


def choose(tem_list, cou):
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
    tem_list = np.delete(tem_list, -1, axis=1)  # 删除最后一列
    # print(tem_list[:cou])
    if len(tem_list) > cou:
        return tem_list[:cou]

    return tem_list


def cal_next_tm(ans, answer):
    tem_res = []
    for i in range(len(ans)):  # 遍历上一轮选择到的最优的
        sta = int(ans[i, -2])  # 这段gene开始
        for j in range(max_len - min_len):  #
            end = sta + min_len + j  # 这段gene结束
            if end >= len(gene) - 1:
                end = len(gene) - 1
            tem_tm = cal_tm(gene[sta: end])  # 计算这段gene的tm
            bef_tm = ans[i, 1::2]  # 取出前面所有tm
            bef_tm = np.append(bef_tm, tem_tm)  # 将这段gene的tm添加到之前中
            tm_std = np.std(bef_tm)  # 计算标准差
            bef_arr = ans[i, :]  # 获取数组，转化为列表
            bef_arr = bef_arr.tolist()
            tem_gene_tm = [end, tem_tm, tm_std]
            tem_list = bef_arr + tem_gene_tm

            if sta + min_len > len(gene) - 1:
                answer1.append(tem_list)  # TODO最后一段是独立好还是分开好
                break
            elif end == len(gene) - 1:
                answer.append(tem_list)
                break
            else:
                tem_res.append(tem_list)
    return tem_res


def optimize(gene, index_t, tm_t, lasts_t):
    avg = np.mean(tm_t)
    tem = abs(tm_t - avg)
    max_index = np.argmax(tem)  # 最大值的下标

    if index_t[0] != 0:
        index_t = np.insert(index_t, 0, [0])  # 在数组开头添加0，方便计算
    else:
        print("func optimize insert error")

    left = index_t[max_index]  # 最大值左部    切割位点
    right = index_t[max_index + 1]  # 最大值有部
    # print("left:{0}, right:{1}, lasts:{2}".format(left, right, lasts_t))

    if [left, right] == lasts_t:
        second = find_top_k(tm_t, 3)  # 选择第二大
        tm_t = list(tm_t)
        max_index = tm_t.index(second)
        left = index_t[max_index]  # 最大值左部    切割位点
        right = index_t[max_index + 1]  # 最大值有部
        # print(left, right, second, max_index)

    fre = 5  # 偏移的最大的个数
    tem_list = []
    for i in range(fre):
        # 左半部，向右
        tem_index = index_t.copy()
        tem_index[max_index] = left + i
        tem_std, _ = cal_all_tm(tem_index)
        tem_list.append(([max_index, left + i, tem_std]))

        if left - i >= 0:
            # 左半部向左
            tem_index = index_t.copy()
            tem_index[max_index] = left - i
            tem_std, _ = cal_all_tm(tem_index)
            tem_list.append([max_index, left - i, tem_std])

        if right + i >= len(gene) - 1:
            tem_index = index_t.copy()
            tem_index[max_index + i] = right + i
            tem_std, _ = cal_all_tm(tem_index)
            tem_list.append([max_index + 1, right - i, tem_std])

        tem_index = index_t.copy()
        tem_index[max_index + 1] = right - i
        tem_std, _ = cal_all_tm(tem_index)
        tem_list.append([max_index + 1, right - i, tem_std])

    tem_list = np.array(tem_list)
    tem_list = list(tem_list[np.lexsort(tem_list.T)])

    index_t[int(tem_list[0][0])] = tem_list[0][1]
    _, tm_t = cal_all_tm(index_t)

    index_t = list(index_t)
    # print(index_t)
    if index_t[0] == 0:
        del (index_t[0])

    return index_t, tm_t, [index_t[max_index - 1], index_t[max_index]]


def over_lap(test_gene, index_t, tm_t):
    print(np.std(tm_t))
    if index_t[0] > 10:
        index_t = np.insert(index_t, 0, [0])  # 在数组开头添加0，方便计算
    gene_list = []
    for i in range(len(tm_t)):  # 将gene截取出来存放在一个list中
        gene_list.append(
            [test_gene[int(index_t[i]):int(index_t[i + 1])], index_t[i], index_t[i], index_t[i + 1], index_t[i + 1],
             tm_t[i]])

    over_size = 8
    temp_avg_tm = [0, 1]
    while temp_avg_tm[-1] != temp_avg_tm[-2]:
        # 不用对标准差影响最大的，只是找最大的
        # avg = np.mean(tm_t)
        # tem = abs(tm_t - avg)
        # max_index = np.argmax(tem)  # 最大值的下标
        max_index = int(np.argmax(tm_t))  # 需要截取的基因的下标
        print(max_index, tm_t[max_index])
        tem_gene = gene_list[max_index][0]  # 需要更改的基因片段
        if (gene_list[max_index][1] != gene_list[max_index][2]) | (
                gene_list[max_index][3] != gene_list[max_index][4]):  # 说明这片段被处理了
            v = [x[-1] for x in gene_list]
            v = np.array(v)
            ind = np.array(index_t)
            ind = ind[1:]
            show_w(ind, v, "end1")

            min_ind = np.argmin(v)
            print(min_ind, ind[min_ind], v[min_ind])
            print(np.std(v))
            print(v)
            return 0
        start = gene_list[max_index][2]
        end = gene_list[max_index][3]
        new_tm_list = tm_t.copy()  # 复制所有的片段的tm方便计算标准差
        new_cut = []
        for i in range(over_size):
            for j in range(over_size):
                gene_len = len(tem_gene)
                new_tm = cal_tm(tem_gene[i:gene_len - j])  # 新片段的tm
                # 顺便计算一下标准差
                new_tm_list[max_index] = new_tm
                tm_std = np.std(new_tm_list)
                # print(new_tm, tm_std)
                tt = [start + i, end - j, new_tm, tm_std]
                new_cut.append(tt)
        tem_list = np.array(new_cut)
        tem_list = tem_list[np.lexsort(tem_list.T)]

        gene_list[max_index][2] = int(tem_list[0, 0])
        gene_list[max_index][3] = int(tem_list[0, 1])
        gene_list[max_index][5] = tem_list[0, 2]
        tm_t[max_index] = tem_list[0, 2]
        temp_avg_tm.append(tm_std)  # 用于结束判断
    print(gene_list)
    # 将大的往下降


def cal_first_tm2(tm_mean):
    res = []
    for i in range(max_len - min_len):
        fir_gene = gene[0:min_len + i]
        fir_tm = cal_tm(fir_gene)
        std = np.abs(tm_mean - fir_tm)
        res.append([min_len + i, fir_tm, std])
    return res


if __name__ == '__main__':
    gene = get_gene('test_gene.txt')
    min_len, max_len = 20, 30
    res = cal_first_tm()
    count = 20  # 每一代取标准差最小的前count个

    answer = []
    # 尝试控制answer不为0
    while len(answer) == 0:
        res = choose(res, count)
        res = cal_next_tm(res, answer)

    answer
    ##

    tm = np.array(answer[0][1:-1:2])
    # tm_std = np.mean(tm)  # 选择tm平均值作为贪心算法初始化值
    tm_std = np.min(tm)  # 选择tm最小值作为贪心算法初始化值

    print(tm_std)

    ans1 = cal_first_tm2(tm_std)
    answer = []
    answer1 = []
    # 尝试控制answer不为0
    count = 5  # 每一代取标准差最小的前count个
    while len(ans1) > 0:
        ans1 = choose(ans1, count)
        ans1 = cal_next_tm(gene, ans1, answer, answer1)
    answer
    tm = np.array(answer[0][1:-1:2])
    tm_std = np.mean(tm)
    print(tm_std)

    dataframe = pd.DataFrame(answer)
    dataframe.to_csv('answer.csv')

    dataframe1 = pd.DataFrame(answer1)
    dataframe1.to_csv('answer1.csv')
    # show_tm(answer)
    """从answer中选取一个较好的作为下面迭代的开始"""

    index = np.array(answer[0][:-1:2])
    tm = np.array(answer[0][1:-1:2])
    show_w(index, tm, "init")

    lasts = [0, 0]
    for i in range(10):
        index, tm, lasts = optimize(gene, index, tm, lasts)
        # show_w(index, tm, "{0}".format(i))

    print(len(index), len(tm))
    if index[0] < min_len:
        index = index[1:]
        tm = tm[1:]

    over_lap(gene, index, tm)
