# -*- coding: utf-8 -*-

"""
「新規性」と書いてある部分がこの研究の新規性である。
つまりこの部分を消せば既存の研究になる。
"""
import numpy as np
import pandas as pd
import re
import sys
import gc

import preprocessor as pp
import algorism as alg
import parsing

"""Class Node"""
class Node:
  def __init__(self, no, name, child, child_num):
    # nodesリストのインデックスがorderになっているのでorderは属性に入れない
    self.no = no # position of the text glycan data
    self.name = name # 結合情報を入れる

    self.elder = None # left (immediately elder sibling)
    self.elder_num = None

    self.younger = None # right (immediately younger sibling)
    self.younger_num = None

    self.parent = None
    self.parent_num = None
    
    self.child = [] # down
    self.child.append(child)

    self.child_num = [] # to from OTMM structure
    self.child_num.append(child_num)

  def add_order(self, order): # 必要かわからないので放置
    self.order = order

  def add_elder(self, elder, elder_num):
    self.elder = elder
    self.elder_num = elder_num

  def add_younger(self, younger, younger_num):
    self.younger = younger
    self.younger_num = younger_num

  def add_parent(self, parent, parent_num):
    self.parent = parent
    self.parent_num = parent_num

  def add_child(self, child, child_num):
    self.child.append(child)
    self.child_num.append(child_num)

"""前処理を行う関数"""
def make_OTMM(row):
  glycan = row["IUPAC Condensed"]
  glycan = pp.get_structure(glycan)
  glycan = pp.create_siblings(glycan)
  glycan = pp.create_parent(glycan)
  glycan = pp.set_instance(glycan)

  row["IUPAC Condensed"] = glycan

  return row

"""データの表示"""
def print_OTMM(glycan):
  print("This is a glycan for OTMM\n")
  for i in range(len(glycan)):
    print("ID:", glycan[i].no)
    print("name:", glycan[i].name)

    if glycan[i].parent == None:
      print("parent:", glycan[i].parent)
    else:
      print("parent:", glycan[i].parent.name, glycan[i].parent.no)

    if glycan[i].child == None:
      print("child:", glycan[i].child)
    else:
      print("child:", glycan[i].child.name, glycan[i].child.no)

    if glycan[i].elder == None:
      print("elder sibling:", glycan[i].elder)
    else:
      print("elder sibling:", glycan[i].elder.name, glycan[i].elder.no)
    
    if glycan[i].younger == None:
      print("younger sibling:", glycan[i].younger)
    else:
      print("younger sibling:", glycan[i].younger.name, glycan[i].younger.no)
    print()

def main(argv):
  """データの前処理"""
  # データをdfに読み込み
  # 必ずGoogleColabにglycan_data.csvをアップロードすること
  df = pd.read_csv("glycan_data.csv") # 細田先生から頂いたデータ

  # データが無い行を消去
  df = df.dropna(subset=['IUPAC Condensed'])

  drop_index = [] #消去する糖鎖のindexを格納する

  # ?が入っている糖鎖は構造が不明な糖鎖
  for row in df.itertuples(): # itertuples()で行ごとにタプルで取り出す
    if "?" in row[2]:
      drop_index.append(row[0])

  # ?が含まれる糖鎖構造を消去
  df = df.drop(drop_index)

  # indexを振りなおす
  df = df.reset_index(drop=True)

  # タンパク質を消去
  def delete_PROT(row):
    if re.fullmatch(r'\(a[0-9]-|\(b[0-9]-$', row["IUPAC Condensed"][-4:]): # [-4:]は文字列の後ろから4文字まで
      row["IUPAC Condensed"] = row["IUPAC Condensed"][:-4] # [:-4]は後ろから4文字を省いたテキスト

    return row

  #1行ごとにdelete_PROT()を呼び出す
  df = df.apply(delete_PROT, axis=1) # axis=1とすることで1行ずつの処理になる


  #1行ごとにseparate_structure()を呼び出す
  df = df.apply(pp.separate_structure, axis=1) # axis=1とすることで1行ずつの処理になる

  """
  データ量の制限
  とりあえず正しく実行できれば良いので、データを101個に減らす。
  """
  df = df.head(9) # データを101個にしてみる
  # df = df

  # プログラム的にIUPAC Condensedが2番目に来ること！（row[2]と記述している部分があるから）
  df = df.drop(['Motifs', 'Subsumption Level', 'ChEBI', 'Monoisotopic Mass', 'Species'], axis=1)

  df = df.apply(make_OTMM, axis=1)

  def find_mark(row):
    i = 0
    glycan = row["IUPAC Condensed"]
    for node in glycan:
      # if node.child == None and node.elder == None: # 厳密にはyoungerがいるのでOTMMでは末端ではない
      #   node.leaf_order = i
      #   i += 1
      if node.child == None and node.younger == None: # これが厳密な末端の葉
        node.leaf_order = i
        i += 1
    row["IUPAC Condensed"] = glycan
    return row

  df = df.apply(find_mark, axis=1)

  """アルゴリズム"""
  label_set = set()

  for row in df.itertuples():
    glycan = row[2]
    for node in glycan:
      label_set.add(node.name)

  # pandas1.5以上はlistに
  label_set = list(label_set)

  print("Number of labels:", len(label_set))
  # print(label_set)

  n = 5 # とりあえず5個にした（状態の数については先生と要相談）
  state_set = set()
  for i in range(n):
    state_set.add(i)

  # pandas1.5以上はlistに
  state_set = list(state_set)

  print("Number of states:", len(state_set))
  # print(state_set)

  """パラメータのinitialize"""

  """初期状態確率分布 π"""
  np.random.seed(seed=0)
  pi = np.random.rand(len(state_set))
  for i in range(len(pi)):
    pi = pi/pi.sum()
  """対数変換"""
  pi = np.log(pi)

  """状態遷移確率 A
  parent - node(me) ・・・ A_a
  """
  np.random.seed(seed=1)
  a_a = np.random.rand(len(state_set), len(state_set))
  a_a = a_a/a_a.sum(axis=1).reshape(-1, 1)
  """対数変換"""
  a_a = np.log(a_a)

  """elder - node(me) ・・・A_b"""
  np.random.seed(seed=2)
  a_b = np.random.rand(len(state_set), len(state_set))
  a_b = a_b/a_b.sum(axis=1).reshape(-1, 1)
  """対数変換"""
  a_b = np.log(a_b)

  """ラベル出力確率分布 B"""
  np.random.seed(seed=3)
  matrix_B = np.random.rand(len(state_set), len(label_set))
  matrix_B = matrix_B/matrix_B.sum(axis=1).reshape(-1, 1)
  """対数変換"""
  matrix_B = np.log(matrix_B)
  b = pd.DataFrame(matrix_B, index=state_set, columns=label_set)

  """初期の尤度"""
  likelihood = alg.calc_likelihood(df, pi, a_a, a_b, b, state_set)
  L = sum(likelihood) # π（全てを掛け合わせる）なのでsumでよい

  """しきい値"""
  epsilon = 100 # とりあえずこのくらい

  gc.collect() # メモリの開放（メモリリーク対策）

  new_pi, new_a_a, new_a_b, new_b = alg.EM(df, pi, a_a, a_b, b, state_set, label_set, L, epsilon)
  print("pi", pi)
  print("pi updated", new_pi)

  """
  Parsing
  viterbiアルゴリズムで経路探索
  """
  print("\nParsing")
  glycan = df["IUPAC Condensed"][0] # glycan[0]
  parsing.parse_glycan(glycan, state_set, new_pi, new_a_a, new_a_b, new_b)

  return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))