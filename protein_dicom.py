import os
import numpy as np
import pandas as pd
import time
import cv2
from itertools import chain

# 获取所有图片
def process_file_name(file_dir):
    L = []  # 路径
    jpg_files = []  # 文件名称
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            # 判断是否是jpg文件
            if os.path.splitext(file)[1] == '.jpg':
                L.append(os.path.join(root, file))
                jpg_files.append(file)
    return L, jpg_files


# 自动分割
def grouping(filename, img_name, protein, protein_type):
    """
    :param filename: 文件路径
    :param protein: 蛋白质浓度列表
    :param my_range: TE序列
    :return: 候选框位置、TE序列、蛋白质浓度、0.2-0.8平均灰度值
    """
    row = len(protein)

    # cv2.imdecode 和 cv2.imencode避免中文路径的干扰
    # 以灰度图的形式读
    im_gray = cv2.imdecode(np.fromfile(filename, dtype=np.uint8), cv2.IMREAD_GRAYSCALE)

    # 高斯滤波并找到边框
    result = blur_OTSU(im_gray)
    # 所有的候选框
    selected_rec = Contours(result)
    TE = img_name.split('.')[0]
    selected_rec.insert(4, 'TE', TE[2:])
    selected_rec.insert(4, 'protein', 0)
    selected_rec.insert(0, 'row', 'row1')

    selected_rec = selected_rec.sort_values(by='y', axis=0, ascending=True)

    # 总共需要检测的个数
    total_row = selected_rec.shape[0]
    # 一层中需要检测的个数
    group_row = int(total_row / row)

    # 标记框在第几层
    for i in range(row):
        selected_rec.iloc[i * group_row:i * group_row + group_row + 1, 0] = i + 1
        # img1所需要的参数设置
        # if i == 0:
        #     selected_rec.iloc[i * group_row:i * group_row + group_row + 1, 0] = i + 1
        # elif i == 1:
        #     selected_rec.iloc[1 + i * group_row:i * group_row + group_row + 2, 0] = i + 1
        # else:
        #     selected_rec.iloc[2 + i * group_row:2 + i * group_row + group_row, 0] = i + 1

    selected_avg = insert_parameter(im_gray, selected_rec, row, protein, protein_type)
    #
    return selected_avg


# 插入标识变量
def insert_parameter(im_gray, selected_rec, row, protein, protein_type):
    # 注意：reset_index不改变自身
    selected_rec = selected_rec.reset_index(drop=True)
    # 将层数设为索引
    selected_rec.set_index('row', inplace=True)

    # 获得对应灰度值DataFrame
    selected_avg = avg_gray(im_gray, selected_rec)

    # 填入protein值
    selected_rec_protein = pd.DataFrame(columns=selected_avg.columns)

    # print(selected_avg)
    for i in range(1, row + 1):
        # 按照横坐标进行排序
        temp_rows = selected_avg.loc[str(i)]
        temp_rows = temp_rows.sort_values(by='x', axis=0, ascending=True)
        temp_rows['protein'] = protein[i - 1]
        selected_rec_protein = selected_rec_protein.append(temp_rows)

    selected_rec_protein.insert(5, 'protein_type', list(chain.from_iterable(protein_type)))
    # print(selected_rec_protein)

    return selected_rec_protein


# 滤波
def blur_OTSU(gray_img):
    blur = cv2.GaussianBlur(gray_img, (5, 5), 0)
    # blur = cv2.bilateralFilter(gray_img, 9, 75, 75)
    ret1, result1 = cv2.threshold(blur, 0, 255, cv2.THRESH_TOZERO + cv2.THRESH_OTSU)
    # cv2.imshow("result1", result1)
    # cv2.waitKey()
    # print(result1)

    return result1


# 找边框
def Contours(blur_result):
    contours, hierarchy = cv2.findContours(blur_result, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

    # 画出所有的待检测图片
    # contours_result = cv2.drawContours(blur_result, contours, -1, (255, 255, 0), 1)
    # cv2.imshow('contours', contours_result)
    # cv2.waitKey()

    selected_rect = pd.DataFrame(columns=('x', 'y', 'w', 'h'))

    # 画出候选框
    for i in range(len(contours)):
        # 不旋转矩形候选框
        # x,y为矩形左上角坐标；w,h为矩形长宽
        x, y, w, h = cv2.boundingRect(contours[i])
        selected_rect.loc[i] = {'x': x, 'y': y, 'w': w, 'h': h}
        contour_result = cv2.rectangle(blur_result, (x, y), (x + w, y + h), (255, 0, 0), 1)

        # 旋转矩形候选框
        # rect = cv2.minAreaRect(contours[i])
        # box = cv2.boxPoints(rect)
        # box = np.int0(box)
        # contour_result = cv2.drawContours(blur_result, [box], 0, (255, 0, 0), 1)
    # print(selected_rect.sort_values(by=['y', 'x'], axis=0, ascending=True))
    # cv2.imshow("contour_result", contour_result)
    # cv2.waitKey()
    return selected_rect


# 求平均灰度
def avg_gray(im_gray, selected_rec_TE):
    gray_list = pd.DataFrame(im_gray)
    avg_value = np.zeros(1)

    value_25 = np.zeros(1)
    value_50 = np.zeros(1)
    value_75 = np.zeros(1)

    # 遍历每一个待选框，找出对应的灰度值信息
    for rec in selected_rec_TE.values:
        little_rec = gray_list.iloc[rec[1]:rec[1] + rec[3], rec[0]:rec[0] + rec[2]]  # [y:y+h,x:x+w]
        little_rec_list = blur_OTSU(little_rec.values)

        # 筛选所有非零值
        mask_0 = little_rec_list != 0
        # list_0 = little_rec_list[mask_0]
        # 对原图进行筛选，而不是对高斯滤波之后的图进行筛选
        list_0 = little_rec.values[mask_0]

        # 获取 25% 50% 75%
        gray_25 = int(np.percentile(list_0, 25))
        gray_50 = int(np.percentile(list_0, 50))
        gray_75 = int(np.percentile(list_0, 75))

        # 取值在20%-80%范围内的所有数
        mask_p20 = list_0 > np.percentile(list_0, 20)
        list_p20 = list_0[mask_p20]

        mask_p80 = list_p20 < np.percentile(list_p20, 80)
        list_p80 = list_p20[mask_p80]

        # 计算平均灰度值
        avg_gray = list_p80.sum() / list_p80.size

        # 记录所有的灰度值
        avg_value = np.append(avg_value, avg_gray)
        value_25 = np.append(value_25, gray_25)
        value_50 = np.append(value_50, gray_50)
        value_75 = np.append(value_75, gray_75)

    avg_value = np.delete(avg_value, 0)
    value_25 = np.delete(value_25, 0)
    value_50 = np.delete(value_50, 0)
    value_75 = np.delete(value_75, 0)

    # 将灰度值添加到信息DataFrame中
    selected_rec_TE['avg_gray'] = avg_value
    selected_rec_TE['25_gray'] = value_25
    selected_rec_TE['50_gray'] = value_50
    selected_rec_TE['75_gray'] = value_75

    return selected_rec_TE


# 对图片进行处理
def img_processing(
        filename='D:/Users/zcz/PycharmProjects/cv_processing/protein_process/img/img2_51_75/img2_51_75_new_name'):
    # 获取文件路径及文件名称
    imgs_dir, img_name = process_file_name(filename)

    img_lower_root = []

    root_list = filename.split('/')
    root = '/'.join(root_list[:-1])
    print(root)

    # 新建输出文件夹根目录 按照output+当前时间的方法命名，避免重复覆盖的问题
    final_root_dir = root + '/output_' + time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()) + '/'
    for dirs in imgs_dir:
        # 去除根目录名称
        dirs = dirs.split('\\')[1:]
        img_lower_root.append('/'.join(dirs))
        # 去除文件名称
        dirs = dirs[:-1]

        # 新建对应输出路径
        final_dirs = final_root_dir + './'.join(dirs)
        # 判断是否已存在相关路径
        isExists = os.path.exists(final_dirs)
        # 若不存在就新建，否则跳过
        if not isExists:
            pass
            os.makedirs(final_root_dir + './'.join(dirs))

    all_img = pd.DataFrame()

    # 蛋白质类型声明
    A = "clIscA1/clCry4"
    B = "clIscA1/clCry4(LB +FAC)"
    D = "zero"
    C = "clIscA1"
    E = "clIscA1(LB +FAC)"
    # 蛋白质浓度
    protein_img1 = [[0.09, 0.17, 0.15, 0.17, 0.15, 0],
                    [0.07, 0.15, 0.13, 0.15, 0.13, 0],
                    [0.08, 0.14, 0.13, 0.14, 0.13],
                    [0.07, 0.13, 0.12, 0.13, 0.12],
                    [0.07, 0.12, 0.11, 0.12, 0.11],
                    [0.06, 0.11, 0.1, 0.11, 0.1],
                    [0.07, 0.1, 0.1, 0.1, 0.1],
                    [0.06, 0.09, 0.09, 0.09, 0.09]]

    protein_img1_change = [[0.0901, 0.1701, 0.1501, 0.1702, 0.1502, 0.0001],
                           [0.0701, 0.1501, 0.1301, 0.1502, 0.1302, 0.0002],
                           [0.0801, 0.1401, 0.1303, 0.1402, 0.1304],
                           [0.0702, 0.1301, 0.1201, 0.1302, 0.1202],
                           [0.0701, 0.1201, 0.1101, 0.1202, 0.1102],
                           [0.0601, 0.1101, 0.1001, 0.1102, 0.1002],
                           [0.0702, 0.1001, 0.1003, 0.1002, 0.1004],
                           [0.0602, 0.0901, 0.0902, 0.0902, 0.0903]]

    protein_img2 = [[1.74, 0.53, 1.24, 3.48, 0.75, 1.1],
                    [1.45, 7, 3.05, 3.29, 0.9, 2.1],
                    [0.62, 3.14, 2.68, 0.2, 0.71, 4],
                    [0.57, 1.78, 1.74, 4.46, 0.66, 3.48],
                    [0.44, 2.87, 0.89, 3.25, 0.4, 2],
                    [0.33, 1.62, 7.18, 5.43, 0.2, 1.43],
                    [0.0001, 3.55, 8.27, 4.46, 0.26, 1.13],
                    [0.0002, 1.92, 1.12, 3.19, 1.74, 0.9]]

    protein_img3 = [[0.96, 0.23, 0.7, 2.35, 0.70],
                    [0.58, 5.27, 0.12, 1.56, 0.68],
                    [0.39, 0.18, 9.07, 1.9, 2.51],
                    [0.32, 4.11, 0.42, 1.06, 4.28],
                    [0.19, 2.1, 0.35, 11.25, 2.36],
                    [0.0001, 1.9, 0.4, 11.3, 9.1],
                    [0.0002, 1.5, 5.95, 10.96, 10.51]]

    # 蛋白质类型
    protein_type_img1 = [[A, B, A, B, A, D],
                         [A, B, A, B, A, D],
                         [A, B, A, B, A],
                         [A, B, A, B, A],
                         [B, B, A, B, A],
                         [B, B, A, B, A],
                         [B, B, A, B, A],
                         [B, B, A, B, A]]

    protein_type_img2 = [[C, C, C, C, C, C],
                         [C, C, C, C, C, C],
                         [C, C, C, C, C, C],
                         [C, C, C, C, C, C],
                         [C, C, C, C, C, C],
                         [C, C, C, C, C, C],
                         [C, C, C, C, C, C],
                         [C, C, C, C, C, C]]

    protein_type_img3 = [[E, E, E, E, E],
                         [E, E, E, E, E],
                         [E, E, E, E, E],
                         [E, E, E, E, E],
                         [E, E, E, E, E],
                         [E, E, E, E, E],
                         [E, E, E, E, E]]

    for img, output_name, lower_root in zip(imgs_dir, img_name, img_lower_root):
        # print(img)

        # cv2.imdecode 和 cv2.imencode避免中文路径的干扰
        # 以灰度图的形式读

        temp_img = grouping(img, output_name, protein_img2, protein_type_img2)
        all_img = all_img.append(temp_img)

    print(all_img)
    all_img.to_csv('img/img2.csv')
    pass


if __name__ == '__main__':
    img_processing()
    pass
