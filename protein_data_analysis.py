import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import os


def gray_bar(raw_img, protein, TE_range, y_label):
    def polygon_under_graph(xlist, ylist, y_start):
        """
        Construct the vertex list which defines the polygon filling the space under
        the (xlist, ylist) line graph.  Assumes the xs are in ascending order.
        """
        # print(xlist)
        # print(ylist)
        return [(xlist[0], y_start), *zip(xlist, ylist), (xlist[-1], y_start)]

    cnames = [
        '#FF8C00', '#F0F8FF',
        '#FAEBD7',
        '#00FFFF',
        '#BDB76B',
        '#8B008B',
        '#556B2F',
        '#9932CC',
        '#7FFFD4',
        '#F0FFFF',
        '#F5F5DC',
        '#FFE4C4',
        '#000000',
        '#FFEBCD',
        '#0000FF',
        '#8A2BE2',
        '#A52A2A',
        '#DEB887',
        '#5F9EA0',
        '#7FFF00',
        '#D2691E',
        '#FF7F50',
        '#6495ED',
        '#FFF8DC',
        '#DC143C',
        '#00FFFF',
        '#00008B',
        '#008B8B',
        '#B8860B',
        '#A9A9A9',
        '#006400',

        '#8B0000',
    ]
    # 画数据结果图
    ax = plt.axes(projection='3d')
    min_gray = np.min(raw_img['avg_gray'])
    max_gray = np.max(raw_img['avg_gray'])

    min_TE = np.min(TE_range)
    max_TE = np.max(TE_range)

    min_protein = np.min(protein)
    max_protein = np.max(protein)
    # 柱状图相关参数
    # ax.set_zlim(min_gray - 5, max_gray + 5)
    # ax.set_xlabel('TE')
    # ax.set_ylabel('protein')
    # ax.set_zlabel('avg_gray')
    # bar_width = bar_length = 0.2

    verts = []
    for i in range(len(protein)):
        raw_img_temp = raw_img.iloc[i * 25:(i + 1) * 25]
        raw_img_temp = raw_img_temp.sort_values(by='TE', axis=0, ascending=True)
        verts.append(
            polygon_under_graph(raw_img_temp['TE'].values, 255 - raw_img_temp['avg_gray'].values, 30))

    poly = PolyCollection(verts, facecolors=cnames[:len(protein)], alpha=.7)
    ax.add_collection3d(poly, zs=protein, zdir='y')

    ax.set_xlabel('TE (ms)')

    # 'Protein:clIscA1(mg/mL)'
    # 'Protein:clIscA1(mg/mL)\nLB+FAC'
    # 'Protein:clIscA1/clCry4(mg/mL)'
    # 'Protein:clIscA1/clCry4(mg/mL)\nLB+FAC'

    # E.coli-clIscA1 (LB+FAC) Cell Density(OD)
    # E.coli-clIscA1/clCry4 (LB+FAC) Cell Density(OD)
    ax.set_ylabel(y_label)
    ax.set_zlabel('Mean gray value')
    ax.set_xlim(min_TE - 10, max_TE + 10)
    ax.set_ylim(min_protein, max_protein)
    # img_group1 30,90
    # img_group2 70,250
    # img_group3 55,105
    ax.set_zlim(30, 95)

    plt.show()
    pass


# 对整体文件进行分类
def classif(protein_data):
    all_typeA = protein_data[protein_data['protein_type'] == A]
    all_typeB = protein_data[protein_data['protein_type'] == B]
    all_typeC = protein_data[protein_data['protein_type'] == C]
    # all_typeA.to_csv('img/img1_typeA.csv')
    all_typeA = all_typeA.sort_values(by='protein', axis=0, ascending=True)
    all_typeB = all_typeB.sort_values(by='protein', axis=0, ascending=True)
    all_typeC = all_typeC.sort_values(by='y', axis=0, ascending=True)
    all_typeA.to_csv('img/img1_typeA.csv')
    all_typeB.to_csv('img/img1_typeB.csv')
    all_typeC.to_csv('img/img1_typeC.csv')


def output_gray(raw_img, protein, TE_range):
    three_gray = raw_img.set_index('protein')
    three_gray = three_gray.set_index('TE', append=True)
    three_gray = three_gray.drop(['x', 'y', 'w', 'h', 'avg_gray', 'protein_type'], axis=1)
    three_gray['25_gray'] = 255 - three_gray['25_gray']
    three_gray['50_gray'] = 255 - three_gray['50_gray']
    three_gray['75_gray'] = 255 - three_gray['75_gray']
    print(three_gray)
    three_gray.to_csv('img/img1/img1_typeA/img1_typeA_three_gray.csv')

    avg_gray_list = pd.DataFrame(index=TE_range)
    for i in range(len(protein)):
        temp = raw_img.iloc[i * 25:(i + 1) * 25]
        temp = temp.sort_values(by='TE', axis=0, ascending=True)
        temp_gray = temp.iloc[:, 9]
        avg_gray_list[str(protein[i])] = np.ones_like(temp_gray) * 255 - temp_gray.values

    avg_gray_list.to_csv('img/img1/img1_typeA/img1_typeA_avg_gray.csv')


if __name__ == '__main__':
    # 蛋白质类型声明
    A = "clIscA1/clCry4"
    B = "clIscA1/clCry4(LB +FAC)"
    C = "zero"
    # TE范围
    TE_range = range(9, 254, 9)

    protein_data = pd.read_csv('img/img1/img1.csv')
    typeA_data = pd.read_csv('img/img1/img1_typeA/img1_typeA.csv')
    typeB_data = pd.read_csv('img/img1/img1_typeB/img1_typeB.csv')
    typeC_data = pd.read_csv('img/img1/img1_typeC/img1_typeC.csv')

    selected_data = typeA_data
    protein_density = selected_data['protein'].unique()
    TE_range = selected_data['TE'].unique()
    TE_range.sort()
    y_label = selected_data['protein_type'][0]
    print(selected_data)

    gray_bar(selected_data, protein_density, TE_range, y_label)
    output_gray(selected_data, protein_density, TE_range)
    pass
