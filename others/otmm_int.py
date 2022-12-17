# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import re
import sys
import gc
import time
import os
import csv

import preprocessor as pp
# import algorism as alg
# import algorism_decimal as alg
# import algorism_round as alg
import algorism_int as alg
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
  num_data = input('Number of data. Input the number or "max":')
  epsilon = input("Epsilon:")
  n = input("Number of states:")

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
  とりあえず正しく実行できれば良いので、データを201個に減らす。
  """
  if num_data == "max":
    df = df # データ全部
  else:
    df = df.head(int(num_data)) # データを201個にしてみる
  print("Number of data:", len(df))

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

  n = int(n) # とりあえず5個にした（状態の数については先生と要相談）
  state_set = set()
  for i in range(n):
    state_set.add(i)

  # pandas1.5以上はlistに
  state_set = list(state_set)

  print("Number of states:", len(state_set))
  # print(state_set)

  """パラメータのinitialize"""
  digit = 5 # 5桁の整数で計算

  """初期状態確率分布 π"""
  np.random.seed(seed=0)
  pi = np.random.rand(len(state_set))
  pi = pi/pi.sum()
  """対数変換"""
  pi = np.log(pi)
  pi = pi * (10**digit)
  pi = pi.astype('int64')

  """状態遷移確率 A
  parent - node(me) ・・・ A_a
  """
  np.random.seed(seed=1)
  a_a = np.random.rand(len(state_set), len(state_set))
  a_a = a_a/a_a.sum(axis=1).reshape(-1, 1)
  """対数変換"""
  a_a = np.log(a_a)
  a_a = a_a * (10**digit)
  a_a = a_a.astype('int64')

  """elder - node(me) ・・・A_b"""
  np.random.seed(seed=2)
  a_b = np.random.rand(len(state_set), len(state_set))
  a_b = a_b/a_b.sum(axis=1).reshape(-1, 1)
  """対数変換"""
  a_b = np.log(a_b)
  a_b = a_b * (10**digit)
  a_b = a_b.astype('int64')

  """ラベル出力確率分布 B"""
  np.random.seed(seed=3)
  matrix_B = np.random.rand(len(state_set), len(label_set))
  matrix_B = matrix_B/matrix_B.sum(axis=1).reshape(-1, 1)
  """対数変換"""
  matrix_B = np.log(matrix_B)
  b = pd.DataFrame(matrix_B, index=state_set, columns=label_set)
  b = b * (10**digit)
  b = b.astype('int64')

  """初期の尤度"""
  print("\nStart time.perf_counter()")
  start = time.perf_counter() # time.time()よりtime.perf_counter()の方が精度が高い
  likelihood = alg.calc_likelihood(df, pi, a_a, a_b, b, state_set)
  L = sum(likelihood) # π（全てを掛け合わせる）なのでsumでよい
  L = L / (10**digit)

  """しきい値"""
  epsilon = int(epsilon) # とりあえずこのくらい

  gc.collect() # メモリの開放（メモリリーク対策）

  print("\nLearning")
  new_pi, new_a_a, new_a_b, new_b, L_all = alg.EM(df, pi, a_a, a_b, b, state_set, label_set, L, epsilon)
  end = time.perf_counter()

  new_pi = new_pi / (10**digit)
  new_a_a = new_a_a / (10**digit)
  new_a_b = new_a_b / (10**digit)
  new_b = new_b / (10**digit)
  print("\nEnd time.perf_counter()\n")
  print("pi updated\n", new_pi)
  print("α updated\n", new_a_a)
  print("β updated\n", new_a_b)
  print("b updated\n", new_b)

  print("\nThe processing time in learning:", end-start)

  # 結果を格納するディレクトリを作成
  dir_path =  'result'+'_'+'novelty'+'_'+str(len(df))+'_'+str(epsilon)+'_'+str(len(state_set))+'_'+str(len(label_set))
  os.makedirs(dir_path, exist_ok=True)
  os.chdir(dir_path)  # 相対パス

  # output (新規性)
  time_name = 'time'+'_'+str(len(df))+'_'+str(epsilon)+'_'+str(len(state_set))+'_'+str(len(label_set))+'.txt'
  with open(time_name, 'w') as f:
    f.write("Time in learning: "+str(end-start)+ " seconds")

  pi_name = 'pi'+'_'+str(len(df))+'_'+str(epsilon)+'_'+str(len(state_set))+'.csv'
  np.savetxt(pi_name, new_pi, delimiter=',')
  a_a_name = 'a_a'+'_'+str(len(df))+'_'+str(epsilon)+'_'+str(len(state_set))+'_'+str(len(state_set))+'.csv'
  np.savetxt(a_a_name, new_a_a, delimiter=',')
  a_b_name = 'a_b'+'_'+str(len(df))+'_'+str(epsilon)+'_'+str(len(state_set))+'_'+str(len(state_set))+'.csv'
  np.savetxt(a_b_name, new_a_b, delimiter=',')
  b_name = 'b'+'_'+str(len(df))+'_'+str(epsilon)+'_'+str(len(state_set))+'_'+str(len(label_set))+'.csv'
  new_b.to_csv(b_name)

  # Likelihoodを出力
  likelihood_name = 'likelihood'+'_'+str(len(df))+'_'+str(epsilon)+'_'+str(len(state_set))+'_'+str(len(label_set))+'.csv'
  with open(likelihood_name, 'wt', encoding='utf-8') as f:
    # ライター（書き込み者）を作成
    writer = csv.writer(f)
    # ライターでデータ（リスト）をファイルに出力
    writer.writerow(L_all)

  os.chdir("..")
  """
  Parsing
  viterbiアルゴリズムで経路探索
  """
  print("\nParsing")
  # 解析する糖鎖を1つ選ぶ
  glycan = df["IUPAC Condensed"][0] # glycan[0]
  # ParsinにDecimal型が扱えない処理があるのでfloatに変換
  new_pi = new_pi.astype(np.float64)
  new_a_a = new_a_a.astype(np.float64)
  new_a_b = new_a_b.astype(np.float64)
  new_b = new_b.astype('float64')
  parsing.parse_glycan(glycan, state_set, new_pi, new_a_a, new_a_b, new_b)

  return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))