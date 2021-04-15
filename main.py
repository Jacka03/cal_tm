#!/usr/bin/env python
# -*- coding: UTF-8 -*-


'''
@Project ：PycharmProjects
@File    ：test3.py
@IDE     ：PyCharm
@Author  ：jacka
@Date    ：2020/12/22 23:20
'''

##
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import random


def quick_sort(data_list):
    length = len(data_list)
    quick_sort_c(data_list, 0, length - 1)


def quick_sort_c(data_list, begin, end):
    """
    可以递归的函数调用
    """
    if begin >= end:
        return
    else:
        # 获取分区数据partition_data最后的下标
        index = partition(data_list, begin, end)
        # print(data_list)
        quick_sort_c(data_list, begin, index - 1)
        quick_sort_c(data_list, index + 1, end)


def partition(data_list, begin, end):
    # 选择最后一个元素作为分区键
    partition_key = data_list[end]

    # index为分区键的最终位置
    # 比partition_key小的放左边，比partition_key 大的放右边
    index = begin
    for i in range(begin, end):
        if data_list[i] < partition_key:
            data_list[i], data_list[index] = data_list[index], data_list[i]
            index += 1

    data_list[index], data_list[end] = data_list[end], data_list[index]
    return index


def find_top_k(data_list, K):
    length = len(data_list)
    begin = 0
    end = length - 1
    index = partition(data_list, begin, end)  # 这里的partition函数就是上面快排用到的函数
    while index != length - K:
        if index > length - K:
            end = index - 1
            index = partition(data_list, begin, index - 1)
        else:
            begin = index + 1
            index = partition(data_list, index + 1, end)
    return data_list[index]


def calTM(fasta):
    AATT = ATTA = TAAT = CAGT = GTCA = CTGA = GACT = CGGC = GCCG = GGCC = 0

    for i in range(len(fasta) - 1):
        if (fasta[i:i + 2] == 'AA') | (fasta[i:i + 2] == 'TT'):
            AATT += 1
        elif fasta[i:i + 2] == 'AT':
            ATTA += 1
        elif fasta[i:i + 2] == 'TA':
            TAAT += 1
        elif (fasta[i:i + 2] == 'CA') | (fasta[i:i + 2] == 'TG'):
            CAGT += 1
        elif (fasta[i:i + 2] == 'GT') | (fasta[i:i + 2] == 'AC'):
            GTCA += 1
        elif (fasta[i:i + 2] == 'CT') | (fasta[i:i + 2] == 'AG'):
            CTGA += 1
        elif (fasta[i:i + 2] == 'GA') | (fasta[i:i + 2] == 'TC'):
            GACT += 1
        elif fasta[i:i + 2] == 'CG':
            CGGC += 1
        elif fasta[i:i + 2] == 'GC':
            GCCG += 1
        elif (fasta[i:i + 2] == 'GG') | (fasta[i:i + 2] == 'CC'):
            GGCC += 1

    H = AATT * (-7.9) + ATTA * (-7.2) + TAAT * (-7.2) + CAGT * (-8.5) + GTCA * (-8.4) + CTGA * (-7.8) + GACT * (
        -8.2) + CGGC * (-10.6) + GCCG * (-9.8) + GGCC * (-8) + 0.1 + 2.3
    S = AATT * (-22.2) + ATTA * (-20.4) + TAAT * (-21.3) + CAGT * (-22.7) + GTCA * (-22.4) + CTGA * (-21) + GACT * (
        -22.2) + CGGC * (-27.2) + GCCG * (-24.4) + GGCC * (-19.9) - 2.8 + 4.1 - 1.4

    return (H * 1000) / (S + 1.987 * math.log10((10 ** -4) / 4)) - 273.15 + 16.6 * math.log10(1.2)


gene = "taagcacctgtaggatcgtacaggtttacgcaagaaaatggtttgttatagtcgaataacaccgtgcgtgttgactatttt" \
       "acctctggcggtgatatactagagaaagaggagaaatactagatgaccatgattacgccaagcgcgcaattaaccctcact" \
       "aaagggaacaaaagctggagctccaccgcggtggcggcagcactagagctagtggatcccccgggctgtagaaattcgata" \
       "tcaagcttatcgataccgtcgacctcgagggggggcccggtacccaattcgccctatagtgagtcgtattacgcgcgctca" \
       "ctggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttc" \
       "gccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaataataacgctgatagt" \
       "gctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctg" \
       "ttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata"


gene = gene.upper()

gene = "CGTTTTAAAGGGCCCGCGCGTTGCCGCCCCCTCGGCCCGCCATGCTGCTATCCGTGCCGCTGCTGCTCGGCCTCCTCGGCCTGGCCGTCGCCGAGCCTGCCGTCTACTTCAAGGAGCAGTTTCTGGACGGAGACGGGTGGACTTCCCGCTGGATCGAATCCAAACACAAGTCAGATTTTGGCAAATTCGTTCTCAGTTCCGGCAAGTTCTACGGTGACGAGGAGAAAGATAAAGGTTTGCAGACAAGCCAGGATGCACGCTTTTATGCTCTGTCGGCCAGTTTCGAGCCTTTCAGCAACAAAGGCCAGACGCTGGTGGTGCAGTTCACGGTGAAACATGAGCAGAACATCGACTGTGGGGGCGGCTATGTGAAGCTGTTTCCTAATAGTTTGGACCAGACAGACATGCACGGAGACTCAGAATACAACATCATGTTTGGTCCCGACATCTGTGGCCCTGGCACCAAGAAGGTTCATGTCATCTTCAACTACAAGGGCAAGAACGTGCTGATCAACAAGGACATCCGTTGCAAGGATGATGAGTTTACACACCTGTACACACTGATTGTGCGGCCAGACAACACCTATGAGGTGAAGATTGACAACAGCCAGGTGGAGTCCGGCTCCTTGGAAGACGATTGGGACTTCCTGCCACC"
gene = "GGCAGATGCGATCCAGCGGCTCTGGGGGCGGCAGCGGTGGTAGCAGCTGGTACCTCCCGCCGCCTCTGTTCGGAGGGTCGCGGGGCACCGAGGTGCTTTCCGGCCGCCCTCTGGTCGGCCACCCAAAGCCGCGGGCGCTGATGATGGGTGAGGAGGGGGCGGCAAGATTTCGGGCGCCCCTGCCCTGAACGCCCTCAGCTGCTGCCGCCGGGGCCGCTCCAGTGCCTGCGAACTCTGAGGAGCCGAGGCGCCGGTGAGAGCAAGGACGCTGCAAACTTGCGCAGCGCGGGGGCTGGGATTCACGCCCAGAAGTTCAGCAGGCAGACAGTCCGAAGCCTTCCCGCAGCGGAGAGATAGCTTGAGGGTGCGCAAGACGGCAGCCTCCGCCCTCGGTTCCCGCCCAGACCGGGCAGAAGAGCTTGGAGGAGCCAAAAGGAACGCAAAAGGCGGCCAGGACAGCGTGCAGCAGCTGGGAGCCGCCGTTCTCAGCCTTAAAAGTT"

min_len = 15  # 每段序列最小长度
max_len = 45  # 每段序列最大长度


def calSumTM(test_gene):
    res = []
    for i in range(max_len - min_len):
        mid = min_len + i
        tem1 = test_gene[:mid]
        sum1 = calTM(tem1)
        for j in range(max_len - min_len):
            end = mid + min_len + j
            tem2 = test_gene[mid: end]
            sum2 = calTM(tem2)
            res.append([mid, sum1, end, sum2, np.std([sum1, sum2])])
    return res


##
def choose(res):
    res = np.array(res)
    res = res[np.lexsort(res.T)]
    # print(res)
    dif = max_len - min_len  # 每段基因最大长度与最小长度的差值为选择的个数
    dif = 10  # TODO 先简单的选择dif个
    """
    x_len = len(res[0]) - 1
    if res[8, x_len] == res[9, x_len]:  # TODO 确定前面四个
        cou = 0
        for i in range(1, len(res)):
            if cou == 4:
                return res[:i, :]
            if res[i - 1, x_len] != res[i, x_len]:
                if i > 10:
                    return res[:i, :]
                cou = cou + 1
    """
    if len(res) >= dif:
        return res[:dif, :]
    else:
        return res


##
def cal(test_gene, ans, answer):
    """
    根据前面的切割位点以及tm找出下一个最合适的切割位点
    :param test_gene: 测试基因
    :param ans: 前面分割好得到的切割位点以及每段的tm
    :param answer: 存放由贪心算法得到的结果
    :return:
    """
    res = []
    for i in range(len(ans)):  # 遍历上一轮选择到的最优的
        sta = int(ans[i, -2])  # 这段gene开始

        if sta < len(test_gene):
            for j in range(max_len - min_len):  #
                end = sta + min_len + j  # 这段gene结束

                if end >= len(test_gene):
                    end = len(test_gene)
                tem_gene = test_gene[sta: end]  # 获取这段gene
                # print(tem_gene)
                tem_tm = calTM(tem_gene)  # 计算这段gene的tm

                bef_tm = ans[i, 1::2]  # 取出前面所有tm
                bef_tm = np.append(bef_tm, tem_tm)  # 将这段gene的tm添加到之前中

                tm_std = np.std(bef_tm)  # 计算标准差
                bef_arr = ans[i, :]  # 获取数组，转化为列表
                bef_arr = bef_arr.tolist()
                tem_gene = [end, tem_tm, tm_std]
                tem_list = bef_arr + tem_gene
                if end >= len(test_gene):
                    answer.append(tem_list)
                    # print("v", sta)
                    break
                else:
                    res.append(tem_list)
    return res


def show_tm(answer):
    """
    输出种情况的最大值和最小值
    :param answer:
    :return:
    """
    for i in range(len(answer)):
        temp = np.array(answer[i])
        x = temp[:-1:2]
        y = temp[1:-1:2]
        plt.scatter(x, y)
        str = "min:{:4f},max:{:4f},max-min={:4f}".format(min(temp[1:-1:2]), max(temp[1:-1:2]),
                                                         max(temp[1:-1:2]) - min(temp[1:-1:2]))
        plt.legend([str])
        str = "min-len:{},max-len:{}".format(min_len, max_len)
        plt.title(str)
        plt.show()


times = 8


def cal_all_tm(arr):
    """
    求这种切割位点的tm的标准差
    :param arr: 切割位点
    :return: np.std(tm_list)标准差， tm_list：每段的tm组成的list
    """
    tm_list = []
    arr = arr.astype(int)
    for i in range(1, len(arr)):
        # print(arr[i - 1], arr[i])
        # print(gene[arr[i-1]: arr[i]])
        tm = calTM(gene[arr[i - 1]: arr[i]])
        tm_list.append(tm)

    return np.std(tm_list), tm_list


def show_w(x, y, head):
    # print(len(x), len(y))
    # tem = [for i in range(len(x))]
    # plt.title("{0}, {1}".format(max(x), min()))
    plt.scatter(x, y)
    tem_str = "min:{:3f},max:{:3f},max-min={:3f},std={:4f}".format(min(y), max(y), max(y) - min(y), np.std(y))
    plt.legend([tem_str])
    plt.title(head)
    plt.show()


"""
使用误差最大的数会使得陷入局部最优解
"""


def incr(index, tm, lasts):
    """
    通过左右移动切割位点使得tm标准差下降
    :param index: 切割位点
    :param tm: 每段的Tm
    :return:
    """

    mean = np.mean(tm)  # 均值
    # print(max(tm))
    # print(np.argmax(tm))
    t_tm = abs(tm - mean)
    max_index = np.argmax(t_tm)  # 最大值的下标
    # print("max_index", max_index)

    if index[0] != 0:
        index = np.insert(index, 0, [0])  # 在数组开头添加0，方便计算
    else:
        print("todo:", index)

    left = index[max_index]  # 最大值左部    切割位点
    right = index[max_index + 1]  # 最大值有部
    # print(left, right)
    # print([left, right], lasts[-2])

    # 当出现在最大值之间振动的时候，选择第二大
    if [left, right] == lasts[-2]:
        # 找出对标准差影响前n大的tm
        tem = find_top_k(t_tm, 2)  # 选择第二大
        t_tm = list(t_tm)
        max_index = t_tm.index(tem)

        left = index[max_index]  # 最大值左部    切割位点
        right = index[max_index + 1]  # 最大值有部
        print("更改", max_index, left, right)

    # 对于左部
    tem_list = []
    for i in range(1, times):
        # 左半部
        # 向右
        tem_index = index.copy()
        tem_index[max_index] = left + i
        tem_std, _ = cal_all_tm(tem_index)
        tem_list.append([max_index, left + i, tem_std])

        if left - i >= 0:
            # 向坐
            tem_index = index.copy()
            tem_index[max_index] = left - i
            tem_std, _ = cal_all_tm(tem_index)
            tem_list.append([max_index, left - i, tem_std])

        # 右半部
        # 向右
        if right + i >= len(gene):
            tem_index = index.copy()

            tem_index[max_index + 1] = right + i
            tem_std, _ = cal_all_tm(tem_index)
            tem_list.append([max_index + 1, right + i, tem_std])

        tem_index = index.copy()
        tem_index[max_index + 1] = right - i
        tem_std, _ = cal_all_tm(tem_index)
        tem_list.append([max_index + 1, right - i, tem_std])

    tem_list = np.array(tem_list)
    res = list(tem_list[np.lexsort(tem_list.T)])

    #
    index[int(res[0][0])] = res[0][1]
    # print(index)
    # index[max_index] = res[0, 0]

    _, tm = cal_all_tm(index)
    index = list(index)
    if index[0] == 0:
        del (index[0])

    print(index)
    print(tm)
    return index, tm, [left, right]


def main():
    # print(len(gene))
    res = calSumTM(gene)

    answer = []
    while len(res) > 0:
        res = choose(res)
        res = np.delete(res, -1, axis=1)  # 删除最后一列
        res = cal(gene, res, answer)
    # show_tm(answer)
    tem = pd.DataFrame(answer)
    tem.to_csv('answer.csv')

    """最后得到的res的长度有可能不一样，如何处理？？？"""
    # TODE 是否应给删除最后一列,维度不一样能不能
    index = np.array(answer[0][:-1:2])
    tm = np.array(answer[0][1:-1:2])

    show_w(index, tm, "init")
    lasts = [[123, 456], [23, 57], [246, 456]]
    for i in range(100):
        index1, tm1, last = incr(index, tm, lasts)
        lasts.append(last)
        index = index1
        tm = tm1
        stri = "{}".format(i)
        # if i % 10 == 0 :
        show_w(index, tm, stri)


if __name__ == '__main__':
    main()

##
