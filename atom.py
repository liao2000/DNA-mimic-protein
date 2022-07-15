import math
import matplotlib.pyplot as plt
from io import TextIOWrapper
import numpy as np


class Atom:
    """
    ATOM        記錄類型
    serial      原子序號
    name        原子名稱
    resName     殘基名稱
    chainID     鏈標記號
    resSeq      residue sequence number 殘基序列號
    x           直角 x 坐標 (埃)
    y           直角 y 坐標 (埃)
    z           直角 z 坐標 (埃)
    occupancy   占有率
    tempFactor  溫度因子
    segID       segment identifier
    """

    def __init__(self, ATOM, serial, name, resName, chainID, resSeq, x, y, z, occupancy, tempFactor, segID):
        self.ATOM = ATOM
        self.serial = int(serial)
        self.name = name
        self.resName = resName
        self.chainID = chainID
        self.resSeq = int(resSeq)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.occupancy = float(occupancy)
        self.tempFactor = float(tempFactor)
        self.segID = segID

    # 將每個座標減掉自己的平均達到總和為零的效果
    def centroid(list: list):
        sumX, sumY, sumZ = 0, 0, 0
        for e in list:
            sumX += e.x
            sumY += e.y
            sumZ += e.z
        avgX, avgY, avgZ = sumX / len(list), sumY / len(list), sumZ / len(list)
        for e in list:
            e.x -= avgX
            e.y -= avgY
            e.z -= avgZ

    # 將原子輸出成 PDB 格式
    # PDB格式參考
    # https://jerkwin.github.io/2015/06/05/PDB%E6%96%87%E4%BB%B6%E6%A0%BC%E5%BC%8F%E8%AF%B4%E6%98%8E/
    def output(self, index=None):
        if index is None:
            index = self.serial

        res = "%4s  %5d  %-4s%-3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %s" % (
            self.ATOM,
            index,
            self.name,
            self.resName,
            self.chainID,
            self.resSeq,
            self.x,
            self.y,
            self.z,
            self.occupancy,
            self.tempFactor,
            self.segID
        )
        return res

    def outputln(self, index=None):
        return self.output(index)+'\n'

    # 對x軸旋轉 (參數為 radius)
    def rotateX(self, rad: float):
        # 1   0     0
        # 0  cos  -sin
        # 0  sin   cos
        y = self.y*math.cos(rad) - self.z*math.sin(rad)
        self.z = self.y*math.sin(rad) + self.z*math.cos(rad)
        self.y = y

    # 對y軸旋轉 (參數為 radius)
    def rotateY(self, rad: float):
        #  cos  0  sin
        #   0   1   0
        # -sin  0  cos
        x = self.x*math.cos(rad)+self.z*math.sin(rad)
        self.z = - self.x*math.sin(rad)+self.z*math.cos(rad)
        self.x = x

    # 對z軸旋轉 (參數為 radius)
    def rotateZ(self, rad: float):
        #  cos  -sin   0
        #  sin   cos   0
        #   0     0    1
        x = self.x*math.cos(rad)-self.y*math.sin(rad)
        self.y = self.x*math.sin(rad)+self.y*math.cos(rad)
        self.x = x

    def copy(self):
        return Atom(self.ATOM, self.serial, self.name, self.resName, self.chainID, self.resSeq, self.x, self.y, self.z, self.occupancy, self.tempFactor, self.segID)


class Protein:
    def __init__(self, list):
        self.list = list

    # 從 PDB 檔讀入檔案
    def fromFile(f: TextIOWrapper):
        lines = f.readlines()
        results = []
        for line in lines:
            x = line.split()
            if len(x) > 0 and x[0] == 'END':
                continue
            if len(x) > 0 and x[0] != 'ATOM':
                continue
            if len(x) > 0 and x[0] == 'TER':
                continue

            if len(x) == 12:
                atom = Atom(x[0], x[1], x[2], x[3], x[4],
                            x[5], x[6], x[7], x[8], x[9], x[10], x[11])
                results.append(atom)
        return Protein(results)

    def centroid(self):
        sumX, sumY, sumZ = 0, 0, 0
        for e in self.list:
            sumX += e.x
            sumY += e.y
            sumZ += e.z
        avgX, avgY, avgZ = sumX / \
            len(self.list), sumY / len(self.list), sumZ / len(self.list)
        for e in self.list:
            e.x -= avgX
            e.y -= avgY
            e.z -= avgZ

    # 對x軸旋轉 (參數為 radius)
    def rotateX(self, rad: float):
        for e in self.list:
            e.rotateX(rad)

    # 對y軸旋轉 (參數為 radius)
    def rotateY(self, rad: float):
        for e in self.list:
            e.rotateY(rad)

    # 對z軸旋轉 (參數為 radius)
    def rotateZ(self, rad: float):
        for e in self.list:
            e.rotateZ(rad)

    def output(self):
        res = ""
        for e in self.list:
            res += e.outputln()
        return res

    # 對xy平面做投影
    def projectToXYScatter(self, fileName: str = "title"):
        x = []
        y = []
        for e in self.list:
            x.append(e.x)
            y.append(e.y)

        plt.figure(figsize=(7, 5))          # 顯示圖框架大小
        plt.style.use("ggplot")             # 使用ggplot主題樣式
        plt.xlabel("x")                     # 設定x座標標題及粗體
        plt.ylabel("y")                     # 設定y座標標題及粗體
        plt.title(fileName,
                  fontsize=15,
                  fontweight="bold")        # 設定標題、字大小及粗體
        plt.scatter(x,                      # x軸資料
                    y,                      # y軸資料
                    c="r",                  # 點顏色
                    s=10,                   # 點大小
                    marker="D")             # 點樣式
        plt.savefig("%s.jpg" % (fileName,))  # 儲存圖檔
        plt.close()

    # 保留特定原子
    def filter(self, target: list = ['ASP', 'GLU']):
        newList = []
        for e in self.list:
            if (e.resName in target) and (e.name in ['OD1', 'OD2', 'OE1', 'OD2']):
                newList.append(e)
        self.list = newList

    # 保留 x,y,z 座標，生成 numpy 格式
    def toNumpy(self):
        ret = np.zeros((len(self.list), 3), dtype=float)
        for i, e in enumerate(self.list):
            ret[i][0] = e.x
            ret[i][1] = e.y
            ret[i][2] = e.z
        return ret

    # 透過 numpy 調整座標值
    def fromNumpy(self, matrix):
        i = 0
        for rows in matrix:
            self.list[i].x = rows[0]
            self.list[i].y = rows[1]
            self.list[i].z = rows[2]
            i += 1

    def copy(self):
        newList = [protein.copy() for protein in self.list]
        return Protein(newList)


class RMSD:
    """
    modified from: https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
    """
    def rmsd(V, W):
        diff = np.array(V) - np.array(W)
        N = len(V)
        return np.sqrt((diff*diff).sum()/N)

    def kabsch_rmsd(P, Q):
        # without translate(centroid), please centroid in advance.
        return RMSD.rmsd(RMSD.kabsch_rotate(P, Q), Q)

    def kabsch_rotate(P, Q):
        U = RMSD.kabsch(P, Q)
        return np.dot(P, U)

    def kabsch(P, Q):
        """
        Using the Kabsch algorithm with two sets of paired point P and Q, centered
        around the centroid. Each vector set is represented as an NxD
        matrix, where D is the the dimension of the space.
        The algorithm works in three steps:
        - a centroid translation of P and Q (assumed done before this function
        call)
        - the computation of a covariance matrix C
        - computation of the optimal rotation matrix U
        For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm
        Parameters
        ----------
        P : array
            (N,D) matrix, where N is points and D is dimension.
        Q : array
            (N,D) matrix, where N is points and D is dimension.
        Returns
        -------
        U : matrix
            Rotation matrix (D,D)
        """
        # Computation of the covariance matrix
        C = np.dot(np.transpose(P), Q)

        # Computation of the optimal rotation matrix
        # This can be done using singular value decomposition (SVD)
        # Getting the sign of the det(V)*(W) to decide
        # whether we need to correct our rotation matrix to ensure a
        # right-handed coordinate system.
        # And finally calculating the optimal rotation matrix U
        # see http://en.wikipedia.org/wiki/Kabsch_algorithm
        V, S, W = np.linalg.svd(C)
        d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

        if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # Create Rotation matrix U
        U = np.dot(V, W)

        return U
