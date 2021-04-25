import math
import numpy as np
import matplotlib.pyplot as plt


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


def show_w(x, y, head):
    # print(len(x), len(y))
    # tem = [for i in range(len(x))]
    # plt.title("{0}, {1}".format(max(x), min()))
    plt.scatter(x, y)
    tem_str = "min:{:3f},max:{:3f},max-min={:3f},std={:4f}".format(min(y), max(y), max(y) - min(y), np.std(y))
    plt.legend([tem_str])
    plt.title(head)
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
        tm = cal_tm(gene[arr[i - 1]: arr[i]])
        tm_list.append(tm)

    return np.std(tm_list), tm_list


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


def cal_next_tm():
    """
    计算第三段到最后的切割位点
    :return:
    """
    tem_res = []
    for i in range(len(res)):  # 遍历上一轮选择到的最优的
        fir_cut = int(res[i, -2])  # 这段gene开始

        for j in range(max_len - min_len):  #
            sec_cut = fir_cut + min_len + j  # 这段gene结束
            if sec_cut > len(gene) - 1:
                sec_cut = len(gene) - 1
            tem_tm = cal_tm(gene[fir_cut: sec_cut])  # 计算这段gene的tm
            bef_tm = res[i, 1::2]  # 取出前面所有tm
            bef_tm = np.append(bef_tm, tem_tm)  # 将这段gene的tm添加到之前中
            tm_std = np.std(bef_tm)  # 计算标准差
            bef_arr = res[i, :]  # 获取数组，转化为列表
            bef_arr = bef_arr.tolist()
            tem_gene_tm = [sec_cut, tem_tm, tm_std]
            tem_list = bef_arr + tem_gene_tm

            if fir_cut + min_len > len(gene) - 1:
                answer1.append(tem_list)  # TODO最后一段是独立好还是分开好
                break
            elif sec_cut == len(gene) - 1:
                fir_ans.append(tem_list)
                break
            else:
                tem_res.append(tem_list)
    return tem_res


def cal_first_tm2(tm_mean):
    tem_res = []
    for i in range(max_len - min_len):
        fir_gene = gene[0:min_len + i]
        fir_tm = cal_tm(fir_gene)
        tem_res.append([min_len + i, fir_tm, np.std([tm_mean, fir_tm])])
    return tem_res


def optimize_f1(index_list, tm_list):  # 全局迭代，从左到右
    # print(index_list, tm_list)
    print(np.std(tm_list))
    """
    全局迭代，
    :param index_list:上次得到最好的剪切位置
    :param tm_list:每个剪切片段的tm
    :return:
    """
    flag = [-1, -2, -3]  # 提前停止遍历的终止判断
    tem_max_len = max_len + 5
    tem_min_len = min_len - 5
    bias_len = 5
    while flag[-1] != flag[-3]:
        for i in range(len(tm_list)):  # 对于整条
            left = index_list[i]
            right = index_list[i + 1]
            # 遍历
            tem_result = []
            for j in range(1 - bias_len, bias_len):  # 对于每个片段,
                tem_left = left + j  # 左边新的切割位点
                if tem_left < 0:
                    # 当第一个位点小于0时，舍弃该情况
                    continue
                # 左边移动后前一个片段的长度判断：
                if i >= 1 and (
                        tem_left - index_list[i - 1] < tem_min_len or tem_left - index_list[i - 1] > tem_max_len):
                    # 当左边位点移动导致左边片段过长或者过短时，舍弃该情况
                    continue

                for k in range(1 - bias_len, bias_len):
                    tem_right = right + k  # 右边新的切割位点
                    if tem_right >= len(gene):
                        # 右边切割位点超越右边界
                        continue
                    # 右节点移动后下一个片段长度判断
                    if i + 1 < len(tm_list) and (
                            index_list[i + 2] - tem_right > tem_max_len or index_list[i + 2] - tem_right < tem_min_len):
                        # 当右边位点移动导致右边片段过长或者过短时，舍弃该情况
                        continue
                    # 当前片段长度是否正确
                    if tem_min_len > (tem_right - tem_left) or (tem_right - tem_left) > tem_max_len:
                        # 当前片段过长或者过短时
                        continue

                    tem_index_list = index_list.copy()
                    tem_index_list[i] = tem_left
                    tem_index_list[i + 1] = tem_right
                    tem_std, erere = cal_all_tm(tem_index_list)
                    tem_result.append([tem_left, tem_right, tem_std])
                    # print(erere)
            tem_result = choose(tem_result, cou=1)
            index_list[i] = tem_result[0, 0]
            index_list[i + 1] = tem_result[0, 1]
            best_std, tm_list = cal_all_tm(index_list)  # 本次迭代得到最好的std和tm_list
            print(best_std)

        flag.append(best_std)
        show_w(index_list[1:], tm_list, "d")
        # print(best_std)
    return index_list, tm_list

"""
def over_lap(test_gene, index_t, tm_t):
    # print(np.std(tm_t))
    if index_t[0] > 10:
        index_t = np.insert(index_t, 0, [0])  # 在数组开头添加0，方便计算
    gene_list = []
    for i in range(len(tm_t)):  # 将gene截取出来存放在一个二维list中
        gene_list.append(
            [test_gene[int(index_t[i]):int(index_t[i + 1])], index_t[i], index_t[i], index_t[i + 1], index_t[i + 1],
             tm_t[i]])

    over_size = 5
    temp_avg_tm = [0, 1]
    while temp_avg_tm[-1] != temp_avg_tm[-2]:  # 迭代，防止只是剪切那一段
        # 不用对标准差影响最大的，只是找最大的
        # avg = np.mean(tm_t)
        # tem = abs(tm_t - avg)
        # max_index = np.argmax(tem)  # 最大值的下标
        max_index = int(np.argmax(tm_t))  # 需要截取的基因的下标
        # print(max_index, tm_t[max_index])
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
        tm_std = 0
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


def over_lap1(tem_index, tem_tm):
    gene_list = []
    for i in range(len(tem_tm)):  # 将gene截取出来存放在一个二维list中
        gene_list.append(
            [gene[int(tem_index[i]):int(tem_index[i + 1])], tem_index[i], tem_index[i], tem_index[i + 1],
             tem_index[i + 1],
             tem_tm[i]])

    over_size = 5
    temp_avg_tm = [0, 1]
    while temp_avg_tm[-1] != temp_avg_tm[-2]:  # 迭代，防止只是剪切那一段
        # 不用对标准差影响最大的，只是找最大的
        # avg = np.mean(tm_t)
        # tem = abs(tm_t - avg)
        # max_index = np.argmax(tem)  # 最大值的下标
        max_index = int(np.argmax(tem_tm))  # 需要截取的基因的下标
        # print(max_index, tm_t[max_index])
        tem_gene = gene_list[max_index][0]  # 需要更改的基因片段
        if (gene_list[max_index][1] != gene_list[max_index][2]) | (
                gene_list[max_index][3] != gene_list[max_index][4]):  # 说明这片段被处理了
            v = [x[-1] for x in gene_list]
            v = np.array(v)
            ind = np.array(tem_index)
            ind = ind[1:]
            show_w(ind, v, "end1")

            min_ind = np.argmin(v)
            print(min_ind, ind[min_ind], v[min_ind])
            print(np.std(v))
            print(v)
            return 0
        start = gene_list[max_index][2]
        end = gene_list[max_index][3]
        new_tm_list = tem_tm.copy()  # 复制所有的片段的tm方便计算标准差
        new_cut = []
        tm_std = 0
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
        tem_tm[max_index] = tem_list[0, 2]
        temp_avg_tm.append(tm_std)  # 用于结束判断
    print(gene_list)
    # 将大的往下降


def over_lap2(index_list, tm_list):
    gene_list = []
    for i in range(len(tm_list)):  # 将gene截取出来存放在一个二维list中
        # 【gene片段，原来第一个切割位点，修改后第一个切割位点，原来第二个切割位点，修改后第二个切割位点，修改后片段tm】
        gene_list.append(
            [gene[int(index_list[i]):int(index_list[i + 1])], index_list[i], index_list[i],
             index_list[i + 1], index_list[i + 1],
             tm_list[i]])

    over_size = 5
    temp_avg_tm = [0, 1]

    tem_max_op = 10  # 相邻两个片段间隔最大值

    tem_max_len = max_len  # 切割后每个片段长度范围
    tem_min_len = min_len - 10
    x = 0
    while temp_avg_tm[-1] != temp_avg_tm[-2]:  # 迭代，终止条件

        for i in range(len(gene_list)):
            new_tm_list = tm_list.copy()
            tem_result = []
            for j in range(over_size):
                for k in range(over_size):
                    # TODO 判断越界
                    # 取片段， 计算片段长度。 对于该片段，获取本次切割的左右位点
                    tem_gene = gene_list[i][0]  # TODO 片段处理出错
                    gene_len = len(tem_gene) - 1
                    left = gene_list[i][2] - gene_list[i][1] + j
                    right = gene_len - (gene_list[i][4] - gene_list[i][3]) - k
                    # 越界判断

                    # print(left, right)
                    if right - left < tem_min_len:
                        # print(left, right, right - left)
                        continue
                    if i > 0 and gene_list[i][2] + j - gene_list[i - 1][3] > tem_max_op:
                        continue
                    if i + 1 < len(gene_list) and gene_list[i + 1][2] - (gene_list[i][3] - k) > tem_max_op:
                        continue
                    sss = tem_gene[int(left): int(right)]
                    tem_tm = cal_tm(sss)
                    new_tm_list[i] = tem_tm
                    tm_std = np.std(new_tm_list)
                    tem_result.append([gene_list[i][2] + j, gene_list[i][3] - k, tem_tm, tm_std])
            # if len(tem_result) == 0:
            #     break
            tem_list = choose(tem_result, 1)

            gene_list[i][2] = int(tem_list[0, 0])
            gene_list[i][3] = int(tem_list[0, 1])
            gene_list[i][5] = tem_list[0, 2]
            tm_list[i] = tem_list[0, 2]
            temp_avg_tm.append(tem_list[0, 3])
            # index_list[i] = int(tem_list[0, 0])
            # index_list[i+1] = int(tem_list[0, 1])
        x = x + 1
        show_w(index_list[1:], tm_list, x)

    a = np.argsort(tm_list)
    print(a)

    for i in range(len(a)):
        if a[i] < 1 or a[i] + 2 > len(a):  # 先不管第一段和最后一段
            continue
        test_tm_list = tm_list.copy()
        test_result = []
        for j in range(int(gene_list[a[i]][2] - gene_list[a[i] - 1][3])):
            print(gene_list[a[i] + 1][2], gene_list[a[i]][3], gene_list[a[i] + 1][2] - gene_list[a[i]][3])
            for k in range(int(gene_list[a[i] + 1][2] - gene_list[a[i]][3])):
                print(gene_list[a[i]][2], gene_list[a[i]][2] - j)
                print(gene_list[a[i]][3], gene_list[a[i]][3] + k)

                test_gene = gene[gene_list[a[i]][2] - j: gene_list[a[i]][3] + k]
                test_tm = cal_tm(test_gene)
                print(test_tm, gene_list[a[i]][5])
                test_tm_list[a[i]] = test_tm
                test_std = np.std(test_tm_list)
                test_result.append([gene_list[a[i]][2] - j, gene_list[a[i]][3] + k, test_tm, test_std])

        test_result = choose(test_result, 1)
        gene_list[a[i]][2] = test_result[0, 0]
        gene_list[a[i]][3] = test_result[0, 1]
        gene_list[a[i]][5] = test_result[0, 2]
        tm_list[a[i]] = test_result[0, 2]
        show_w(index_list[1:], tm_list, "test")

    for i in range(len(gene_list)):
        print("原来+{0}，更改{1}".format(gene_list[i][4] - gene_list[i][1], gene_list[i][3] - gene_list[i][2]))

"""
def over_lap3(index_list, tm_list):
    index_list = index_list.astype(int)
    gene_list = []
    for i in range(len(tm_list)):  # 将gene截取出来存放在一个二维list中
        # 【原来第一个切割位点，修改后第一个切割位点，原来第二个切割位点，修改后第二个切割位点，修改后片段tm】
        gene_list.append([index_list[i], index_list[i], index_list[i + 1], index_list[i + 1], tm_list[i]])

    over_size = 6
    temp_avg_tm = [0, 1]  # flag 结束迭代的标志

    tem_max_op = 10  # 相邻两个片段间隔最大值

    tem_min_len = min_len - 10  # 切割后每个片段最小值
    x = 0
    while temp_avg_tm[-1] != temp_avg_tm[-2]:  # 迭代，终止条件

        for i in range(len(gene_list)):
            new_tm_list = tm_list.copy()
            tem_result = []
            for j in range(over_size):
                for k in range(over_size):
                    left = gene_list[i][1] + j
                    right = gene_list[i][2] - k

                    if right - left < tem_min_len:  # 长度小于限定值
                        continue
                    if (i == 0 and gene_list[i][1] - 0 > tem_max_op) or (
                            i > 0 and left - gene_list[i - 1][2] > tem_max_op):  # 这段与前一段的距离
                        continue
                    if (i + 1 == len(gene_list) and len(gene) - gene_list[i][2] - 1 > tem_max_op) or (
                            i + 1 < len(gene_list) and gene_list[i + 1][1] - right > tem_max_op):
                        continue

                    tem_tm = cal_tm(gene[left: right])
                    new_tm_list[i] = tem_tm
                    tm_std = np.std(new_tm_list)
                    tem_result.append([left, right, tem_tm, tm_std])
            if len(tem_result) == 0:
                continue
            tem_result = choose(tem_result, 1)
            gene_list[i][1] = int(tem_result[0, 0])
            gene_list[i][2] = int(tem_result[0, 1])
            tm_list[i] = tem_result[0, 2]
            temp_avg_tm.append(tem_result[0, 3])
        x = x + 1
        show_w(index_list[1:], tm_list, x)

    a = np.argsort(tm_list)
    print(a)

    # 经过剪切后，在迭代一次，进行扩展
    for i in range(len(a)):
        if a[i] < 1 or a[i] + 2 > len(a):  # 先不管第一段和最后一段
            continue
        test_tm_list = tm_list.copy()
        test_result = []
        for j in range(int(gene_list[a[i]][1] - gene_list[a[i] - 1][2])):
            print(gene_list[a[i] + 1][1], gene_list[a[i]][2], gene_list[a[i] + 1][1] - gene_list[a[i]][2])
            for k in range(int(gene_list[a[i] + 1][1] - gene_list[a[i]][2])):

                print(gene_list[a[i]][1], gene_list[a[i]][1] - j)
                print(gene_list[a[i]][2], gene_list[a[i]][2] + k)

                test_gene = gene[gene_list[a[i]][1] - j: gene_list[a[i]][2] + k]
                test_tm = cal_tm(test_gene)
                print(test_tm, gene_list[a[i]][4])
                test_tm_list[a[i]] = test_tm
                test_std = np.std(test_tm_list)
                test_result.append([gene_list[a[i]][1] - j, gene_list[a[i]][2] + k, test_tm, test_std])
        if len(test_result) == 0:
            continue
        # 下表从0开始
        test_result = choose(test_result, 1)
        gene_list[a[i]][1] = test_result[0, 0]
        gene_list[a[i]][2] = test_result[0, 1]
        gene_list[a[i]][4] = test_result[0, 2]
        tm_list[a[i]] = test_result[0, 2]
        show_w(index_list[1:], tm_list, "test")

    for i in range(len(gene_list)):
        print("原来+{0}，更改{1}".format(gene_list[i][3] - gene_list[i][0], gene_list[i][2] - gene_list[i][1]))


if __name__ == '__main__':
    gene = get_gene('test_gene/test_gene.txt')
    gene = "CGTTTTAAAGGGCCCGCGCGTTGCCGCCCCCTCGGCCCGCCATGCTGCTATCCGTGCCGCTGCTGCTCGGCCTCCTCGGCCTGGCCGTCGCCGAGCCTGCCGTCTACTTCAAGGAGCAGTTTCTGGACGGAGACGGGTGGACTTCCCGCTGGATCGAATCCAAACACAAGTCAGATTTTGGCAAATTCGTTCTCAGTTCCGGCAAGTTCTACGGTGACGAGGAGAAAGATAAAGGTTTGCAGACAAGCCAGGATGCACGCTTTTATGCTCTGTCGGCCAGTTTCGAGCCTTTCAGCAACAAAGGCCAGACGCTGGTGGTGCAGTTCACGGTGAAACATGAGCAGAACATCGACTGTGGGGGCGGCTATGTGAAGCTGTTTCCTAATAGTTTGGACCAGACAGACATGCACGGAGACTCAGAATACAACATCATGTTTGGTCCCGACATCTGTGGCCCTGGCACCAAGAAGGTTCATGTCATCTTCAACTACAAGGGCAAGAACGTGCTGATCAACAAGGACATCCGTTGCAAGGATGATGAGTTTACACACCTGTACACACTGATTGTGCGGCCAGACAACACCTATGAGGTGAAGATTGACAACAGCCAGGTGGAGTCCGGCTCCTTGGAAGACGATTGGGACTTCCTGCCACC"
    gene = "GGCAGATGCGATCCAGCGGCTCTGGGGGCGGCAGCGGTGGTAGCAGCTGGTACCTCCCGCCGCCTCTGTTCGGAGGGTCGCGGGGCACCGAGGTGCTTTCCGGCCGCCCTCTGGTCGGCCACCCAAAGCCGCGGGCGCTGATGATGGGTGAGGAGGGGGCGGCAAGATTTCGGGCGCCCCTGCCCTGAACGCCCTCAGCTGCTGCCGCCGGGGCCGCTCCAGTGCCTGCGAACTCTGAGGAGCCGAGGCGCCGGTGAGAGCAAGGACGCTGCAAACTTGCGCAGCGCGGGGGCTGGGATTCACGCCCAGAAGTTCAGCAGGCAGACAGTCCGAAGCCTTCCCGCAGCGGAGAGATAGCTTGAGGGTGCGCAAGACGGCAGCCTCCGCCCTCGGTTCCCGCCCAGACCGGGCAGAAGAGCTTGGAGGAGCCAAAAGGAACGCAAAAGGCGGCCAGGACAGCGTGCAGCAGCTGGGAGCCGCCGTTCTCAGCCTTAAAAGTT"
    gene = "ATGAGATTTAGTTCAACGGATATGCAATACCAAAAGATGCTATTTGCTGCTATTCTATTTATTTGTGCATTAAGTTCGAAGAAGATCTCAATCTATAATGAAGAAATGATAGTAGCTGGTTGTTTTATAGGCTTTCTCATATTCAGTCGGAAGAGTTTAGGTAAGACTTTCCAAGCCACTCTCGACGGGAGAATCGAGTCTATTCAGGAAGAATCGCAGCAATTCTCCAATCCTAACGAAGTCCTTCCTCCGGAATCCAATGAACAACAACGATTACTTAGGATCAGCTTGCAAATTTGCGGCACCGTAGTAGAATCATTACCAACGGCACGCTGTGCGCCTAAGTGCGAAAAGACAGTGCAAGCTTTGTTATGCCGAAACCTAAATGTTAAGTCAGAAACACTTCTAAATGCCACTTCTTCCCGTCGCATCCGTCTTCAGGCCGATATAGTCACAGGGTTTAACTTTGGGGTGAGTGAAAGTGGGTGTACGTTGAAAACTTCTATCGTAGAACTAATTCGAGAGGGCTTGGTAGTCTTAAAAATAGCCTAA"
    print(len(gene))
    min_len, max_len = 15, 35
    res = cal_first_tm()
    count = 20  # 每一代取标准差最小的前count个

    fir_ans = []  # 初步切割结果
    answer1 = []
    # 尝试控制answer不为0
    while len(fir_ans) == 0:
        res = choose(res, count)
        res = np.delete(res, -1, axis=1)  # 删除最后一列
        res = cal_next_tm()
    fir_ans = choose(fir_ans)
    if len(answer1) > 0 and len(answer1[0]) == len(answer1[-1]):
        answer1 = choose(answer1)
        if fir_ans[0, -1] > answer1[0, -1]:
            fir_ans = answer1
    index = np.array(fir_ans[0][:-1:2])
    tm = np.array(fir_ans[0][1::2])
    show_w(index, tm, "init")

    tm_init = np.mean(tm)  # 选择tm最小值作为贪心算法初始化值
    res = cal_first_tm2(tm_init)
    fir_ans = []
    answer1 = []
    count = 5  # 每一代取标准差最小的前count个
    while len(fir_ans) == 0:
        res = choose(res, count)
        res = np.delete(res, -1, axis=1)  # 删除最后一列
        res = cal_next_tm()

    fir_ans = choose(fir_ans)
    if len(answer1) > 0 and len(answer1[0]) == len(answer1[-1]):
        answer1 = choose(answer1)
        if fir_ans[0, -1] > answer1[0, -1]:
            fir_ans = answer1
    index = np.array(fir_ans[0][:-1:2])
    tm = np.array(fir_ans[0][1::2])
    show_w(index, tm, "init1")

    # print(index, tm)
    print(np.std(tm))
    """从answer中选取一个较好的作为下面迭代的开始"""
    # 对整体遍历
    index = np.insert(index, 0, [0])
    index, tm = optimize_f1(index, tm)
    # print()
    # for i in range(len(index) - 1):
    #     print(index[i + 1] - index[i], end=' ')
    # print()
    #print(index)
    over_lap3(index, tm)
