import numpy as np
import pandas as pd
from numpy.linalg import matrix_rank
from scipy.stats import f
import matplotlib.pyplot as plt
from tabulate import tabulate

def Raw_Gene_Data(Path):
    data = pd.read_csv(Path, sep=" ").values
    g_data = []
    g_symbol = []
    for i in data:
        y = []
        g = i[0].split('\t')
        g_symbol.append(g[49])
        g.remove(g[0])
        for i in range(len(g)):
            if (i < 48):
                y.append(float(g[i]))
        g_data.append(y)
    return g_data,g_symbol


def Cancer_Gene_Data(Path_1,Path_2,Path_3,Path_4):
    list_1=[]
    list_2=[]
    list_3=[]
    list_4=[]
    
    f = open(Path_1, "r")
    k=0
    for line in f:
        if(k>=2):
            list_1.append(line.strip())
        k=k+1
    f.close()
    f = open(Path_4, "r")
    k = 0
    for line in f:
        if (k >= 2):
            list_4.append(line.strip())
        k = k + 1
    f.close()
    f = open(Path_2, "r")
    k = 0
    for line in f:
        if (k >= 2):
            list_2.append(line.strip())
        k = k + 1
    f.close()
    f = open(Path_3, "r")
    k = 0
    for line in f:
        if (k >= 2):
            list_3.append(line.strip())
        k = k + 1
    f.close()
    return list_1,list_2,list_3,list_4


def Model_Matrix():
    M = []
    M_hat = []
    for i in range(4):
        for j in range(12):
            if i == 0:
                M.append([0, 0, 1, 0])
                M_hat.append([0, 1, 1, 0])
            elif i == 1:
                M.append([1, 0, 0, 0])
                M_hat.append([1, 0, 1, 0])
            elif i == 2:
                M.append([0, 0, 0, 1])
                M_hat.append([0, 1, 0, 1])
            elif i == 3:
                M.append([0, 1, 0, 0])
                M_hat.append([1, 0, 0, 1])
    return M,M_hat

def Calculate_P_value(g_data,M,M_hat):
    P_val = []
    for i in g_data:
        temp1 = np.matmul(np.matmul(M,np.linalg.pinv(np.matmul(np.transpose(M),M))),np.transpose(M))
        temp2 = np.matmul(np.matmul(M_hat, np.linalg.pinv(np.matmul(np.transpose(M_hat), M_hat))), np.transpose(M_hat))
        stat=(np.matmul(np.matmul(np.transpose(i),np.subtract(temp1,temp2)),i))/(np.matmul(np.matmul(np.transpose(i),np.subtract(np.identity(48),temp1)),i))
        dfn=48-matrix_rank(M)
        dfd=matrix_rank(M)-matrix_rank(M_hat)
        fstat=stat*(dfn)/(dfd)
        pvalue=f.cdf(fstat,dfd,dfn)
        P_val.append(1-pvalue)
    return P_val

def Short_Gene_Symbols(P_val,g_symbol):
    g_Short_List = []
    for i in range(len(P_val)):
        g = []
        if (P_val[i] <= 0.05):
            g.append(g_symbol[i])
            g.append(P_val[i])
            g_Short_List.append(g)
    return g_Short_List


def Intersection_List(list_1, list_2, list_3, list_4,g_Short_List):
    list_1_intersec = []
    list_2_intersec = []
    list_3_intersec = []
    list_4_intersec = []
    for gene in g_Short_List:
        if gene[0] in list_1:
            list_1_intersec.append(gene)
        if gene[0] in list_2:
            list_2_intersec.append(gene)
        if gene[0] in list_3:
            list_3_intersec.append(gene)
        if gene[0] in list_4:
            list_4_intersec.append(gene)
    return list_1_intersec,list_2_intersec,list_3_intersec,list_4_intersec

def Gene_Avg_Value(List):
    g_Effect = []
    i = 0
    for gene in List:
        g_id = gene[0]
        for j in range(len(g_symbol)):
            if(g_symbol[j]==g_id):
                g_Avg = [g_id]
                gdata = g_data[j]
                g_Avg.append(np.mean(gdata[0:11]))
                g_Avg.append(np.mean(gdata[12:23]))
                g_Avg.append(np.mean(gdata[24:35]))
                g_Avg.append(np.mean(gdata[36:47]))
                g_Effect.append(g_Avg)
        i = i + 1
    return g_Effect

def Gene_Effect(List):
    W_SM_vs_NSM_up = []
    W_SM_vs_NSM_down = []
    M_SM_vs_NSM_up = []
    M_SM_vs_NSM_down = []
    for i in range(len(List)):
        gene=List[i]
        g_id=gene[0]
        if(gene[1]>gene[2]):
            M_SM_vs_NSM_down.append(g_id)
        else:
            M_SM_vs_NSM_up.append(g_id)
        if (gene[3] > gene[4]):
            W_SM_vs_NSM_down.append(g_id)
        else:
            W_SM_vs_NSM_up.append(g_id)
    return W_SM_vs_NSM_up,W_SM_vs_NSM_down,M_SM_vs_NSM_up,M_SM_vs_NSM_down

def Plot(P_val):
    plt.hist(P_val,bins=25)
    plt.xlabel("P_values")
    plt.ylabel("Frequency")
    plt.savefig('hist.png')
    plt.show()


if __name__ == "__main__":
    Path = "../data/Raw Data_GeneSpring.txt"
    Path_1 = "../data/XenobioticMetabolism1.txt"
    Path_2 = "../data/FreeRadicalResponse.txt"
    Path_3 = "../data/DNARepair1.txt"
    Path_4 = "../data/NKCellCytotoxicity.txt"

    g_data, g_symbol=Raw_Gene_Data(Path)
    M, M_hat=Model_Matrix()
    P_val = Calculate_P_value(g_data,M,M_hat)
    ## Question 1 :Generated P values  using 2-way ANOVA framework:

    print("\n")
    
    ## Question 2:Draw the Histogram of p-values:

    Plot(P_val)

    ## Question 3:  n0 is taken same as n because of the histogram observed:

    ## Question 4: it is not possible to sortlist based on FDR as n0=n:

    ## Question 5: shortlist rows:

    g_Short_List =Short_Gene_Symbols(P_val,g_symbol)
    

    ## Question 6: Intersect with the gene lists:   
    
    list_1, list_2, list_3, list_4=Cancer_Gene_Data(Path_1,Path_2,Path_3,Path_4)
    # list_1, list_2, list_3, list_4=Cancer_Gene_Data()
     
    list_1_intersec,list_2_intersec,list_3_intersec,list_4_intersec=Intersection_List(list_1, list_2, list_3, list_4,g_Short_List)

    print("Intersection with Xenobiotic metabolism Genes :")
    print(list_1_intersec)
    print("\n")
    
    print("Intersection with Free Radical Response Genes :")
    print(list_2_intersec)
    print("\n")
    
    print("Intersection with DNA Repair Genes :")
    print(list_3_intersec)
    print("\n")
  
    print("Intersection with Natural Killer Cell Cytotoxicity Genes :")
    print(list_4_intersec)
    print("\n")

    ## Question 7: increase/decrease Genes in Women Smokers vs non-Smokers and Men Smokers vs non-Smokers

    Effect_list_1 = Gene_Avg_Value(list_1_intersec)
    Effect_list_3 = Gene_Avg_Value(list_3_intersec)
    Effect_list_4 = Gene_Avg_Value(list_4_intersec)

    print("********************For Xenobiotic metabolism Genes:**************************")
    print(tabulate(Effect_list_1,headers=['Geneid', 'M_Non_Smoker', 'M_Smoker', 'W_Non_Smoker', 'W_Smoker']))

    W_SM_vs_NSM_up, W_SM_vs_NSM_down, M_SM_vs_NSM_up, M_SM_vs_NSM_down = Gene_Effect(Effect_list_1)
    print("Up Genes Smokers vs non-Smokers Women:", list(dict.fromkeys(W_SM_vs_NSM_up)))
    print("Down Genes Smokers vs non-Smokers Women:", list(dict.fromkeys(W_SM_vs_NSM_down)))
    print("Up Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_SM_vs_NSM_up)))
    print("Down Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_SM_vs_NSM_down)))

    print("\n")
    W_Up_List=W_SM_vs_NSM_up
    W_Down_List=W_SM_vs_NSM_down
    M_Up_List =M_SM_vs_NSM_up
    M_Down_List =M_SM_vs_NSM_down

    print("*********************For DNA Repair Genes:*********************")
    print(tabulate(Effect_list_3,headers=['Geneid', 'M_Non_Smoker', 'M_Smoker', 'W_Non_Smoker', 'W_Smoker']))
    W_SM_vs_NSM_up, W_SM_vs_NSM_down, M_SM_vs_NSM_up, M_SM_vs_NSM_down = Gene_Effect(Effect_list_3)
    print("Up Genes Smokers vs non-Smokers Women:",list(dict.fromkeys(W_SM_vs_NSM_up)))
    print("Down Genes Smokers vs non-Smokers Women:",list(dict.fromkeys(W_SM_vs_NSM_down)))
    print("Up Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_SM_vs_NSM_up)))
    print("Down Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_SM_vs_NSM_down)))


    print("\n")
    W_Up_List.extend(W_SM_vs_NSM_up)
    W_Down_List.extend(W_SM_vs_NSM_down)
    M_Up_List.extend(M_SM_vs_NSM_up)
    M_Down_List.extend(M_SM_vs_NSM_down)

    print("**********************For Natural Killer Cell Cytotoxicity Genes:*************************")
    print(tabulate(Effect_list_4,headers=['Geneid', 'M_Non_Smoker', 'M_Smoker', 'W_Non_Smoker', 'W_Smoker']))
    W_SM_vs_NSM_up, W_SM_vs_NSM_down, M_SM_vs_NSM_up, M_SM_vs_NSM_down = Gene_Effect(Effect_list_4)
    print("Up Genes Smokers vs non-Smokers Women:", list(dict.fromkeys(W_SM_vs_NSM_up)))
    print("Down Genes Smokers vs non-Smokers Women:", list(dict.fromkeys(W_SM_vs_NSM_down)))
    print("Up Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_SM_vs_NSM_up)))
    print("Down Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_SM_vs_NSM_down)))

    print("\n")
    W_Up_List.extend(W_SM_vs_NSM_up)
    W_Down_List.extend(W_SM_vs_NSM_down)
    M_Up_List.extend(M_SM_vs_NSM_up)
    M_Down_List.extend(M_SM_vs_NSM_down)
    
    print("\n")
    print("***********************For OverAll types of Genes ********************")
    print("Up Genes Smokers vs non-Smokers Women:", list(dict.fromkeys(W_Up_List)))
    print("Down Genes Smokers vs non-Smokers Women:", list(dict.fromkeys(W_Down_List)))
    print("Up Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_Up_List)))
    print("Down Genes Smokers vs non-Smokers  Men :", list(dict.fromkeys(M_Down_List)))