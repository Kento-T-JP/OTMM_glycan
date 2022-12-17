"""OTMMの前処理"""

"""
(新規性の無い)既存の研究に対応したプログラムである
"""

import re
import copy

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

# タンパク質を消去
def delete_PROT(row):
  if re.fullmatch(r'\(a[0-9]-|\(b[0-9]-$', row["IUPAC Condensed"][-4:]): # [-4:]は文字列の後ろから4文字まで
    row["IUPAC Condensed"] = row["IUPAC Condensed"][:-4] # [:-4]は後ろから4文字を省いたテキスト
  return row

# 糖鎖のテキストデータを行ごとに分ける
def separate_structure(row):
  text = row["IUPAC Condensed"]

  stack = [] # stackを利用
  result = []
  first_flag = False # 最初の単糖を処理したかどうか
  for i, word in enumerate(text):
    if i == 0: # 最初の単糖
      stack.append(0)
      first_flag = True

    if word == '(':
      if first_flag == True: #最初の単糖
        result.append(text[0:i])
        first_flag = False
        stack.append(i)
      else:
        result.append(text[stack.pop():i]) # 単糖
        stack.append(i) # 括弧

    elif word == ')':
      result.append(text[stack.pop():i + 1]) # 括弧
      if text[i+1] != '[':
        stack.append(i + 1) # 単糖

    elif word == '[':
      result.append(text[i]) # 角括弧
      stack.append(i + 1) # 単糖

    elif word == ']':
      result.append(text[i]) # 角括弧
      if text[i+1] != '[':
        stack.append(i + 1) # 単糖
    
    # 途中で括弧などが終わっていてもスタックをすべて使いきれるようにする
    if i == len(text)-1: # "文字数-1"が最後の文字のインデックス
      result.append(text[stack.pop():i + 1])

  row["IUPAC Condensed"] = result
  return row

# 親子関係を抽出（この時点では1人の親が複数の子供を持つ）
def get_structure(result):
  nodes = [] # add nodes of glycan (tree)
  parent_index = set()

  i = len(result)-1
  while i >= 0:
    # ルートのみ例外処理（親がいないので結合情報を入れない）
    if i == len(result)-1:
      if re.fullmatch(r'\([a-z][0-9]-[0-9]\)', result[i-1]) or re.fullmatch(r'\([a-z][0-9]-[0-9]\/[0-9]\)', result[i-1]) or re.fullmatch(r'\([a-z][0-9]-[a-zA-Z]-[0-9]\)', result[i-1]): # 結合様式→子供が1人だけの時
        # print("OK00")
        # print(i)
        nodes.append(Node(i, result[i], result[i-2], i-2)) # nodes[0]==root
        parent_index.add(i)
        # print("root:", result[i-2], i-2, "←", result[i], i, "☆bond:", result[i-1])

      elif result[i-1] == "]": # 分岐の時
        # print("OK01")
        # print(i)
        nodes.append(Node(i, result[i], result[i-3], i-3)) # nodes[0]==root
        parent_index.add(i)
        # print("root:", result[i-3], i-3, "←", result[i], i, "☆bond:", result[i-2])

    # 通常の処理
    elif re.fullmatch(r'[0-9a-zA-Z]+', result[i]) or re.fullmatch(r'[a-zA-z]-[0-9a-zA-Z]+', result[i]) or re.fullmatch(r'[0-9a-zA-Z]+[0-9]\/[0-9][a-zA-Z]', result[i]): # 単糖の場合
      # 1:1の親子関係を認識
      if re.fullmatch(r'\([a-z][0-9]-[0-9]\)', result[i-1]) or re.fullmatch(r'\([a-z][0-9]-[0-9]\/[0-9]\)', result[i-1]) or re.fullmatch(r'\([a-z][0-9]-[a-zA-Z]-[0-9]\)', result[i-1]): # 結合様式
        # print("OK1")
        # print(i)
        nodes.append(Node(i, result[i], result[i-2], i-2)) # 既存の研究
        parent_index.add(i)
        # print(result[i-2], i-2, "←", result[i], i, "☆bond:", result[i-1])
      # 分岐の始まりの場合
      elif result[i-1] == "]":
        parent = result[i] # 分岐元の親を保存（最初の分岐の親）。必ず"["の直後に対応する"]"が来るので1つの変数に適時入れるだけで問題ない
        i_parent = i
        # print(i)
        i, parent_index = branch(result, i-1, nodes, parent_index, parent, i_parent) # iは"["に対応する"]の1つ前のインデックス
        # print("after", i)
        # print("＜branch関数終わり＞")

      elif i == 0: #末尾（葉）の場合
        # print(i)
        nodes.append(Node(i, result[i], None, None)) # 既存の研究
        # 親はいない

    # 分岐の始まりの場合
    elif result[i] == "]":
      if result[i+1] == "[":
        i, parent_index = branch(result, i, nodes, parent_index, parent, i_parent)
      else:
        parent = result[i+1]
        i_parent = i+1
        i, parent_index = branch(result, i, nodes, parent_index, parent, i_parent)
    
    # 分岐の終わりの場合
    elif result[i] == "[":
      # 分岐が連続する場合
      if result[i-1] == "]":
        # 最初の分岐の親がこの分岐の親
        # print("＜新たなbranch＞")
        # print("OK2")
        for node in nodes:
          if node.no == i_parent:
            node.add_child(result[i-3], i-3)
        # print(result[i-3], i-3, "←", parent, i_parent, "☆bond:", result[i-2]) # その他1:1の親子関係については通常の処理で対応可能

      # その階層の分岐が終わる場合（最後の分岐）
      elif re.fullmatch(r'\([a-z][0-9]-[0-9]\)', result[i-1]) or re.fullmatch(r'\([a-z][0-9]-[0-9]\/[0-9]\)', result[i-1]) or re.fullmatch(r'\([a-z][0-9]-[a-zA-Z]-[0-9]\)', result[i-1]): # 結合様式
        # print("OK3")
        # print("＜新たなbranch＞")
        for node in nodes:
          if node.no == i_parent:
            node.add_child(result[i-2], i-2)
            # print(node.no, node.child)
        # print(result[i-2], i-2, "←", parent,i_parent,  "☆bond:", result[i-1])
        # print("＜branch終わり＞\n")
    i -= 1

  # print("\nnumber of nodes:",len(nodes), "\n")
  # for node in nodes:
    # print(node.no, node.name, node.child)
  # print()
  # print(parent_index)

  return nodes

# get_structureで使う関数
def branch(result, j, nodes, parent_index, parent, j_parent):
  # print("\n＜branch始まり＞")
  start = j
  # print(j)
  # print(j_parent)
  while result[j] != "[":
    # 分岐中の親子関係を記述
    if j == start: # ]のとき
      if j+1 in parent_index: # 一つ右のノードが既出の親の場合
        # print(j+1)
        if result[j+1] != result[-1]: # ルートノードではないとき（ルートノードは既に子どもを登録済み）
          for node in nodes:
            if node.no == j+1: # j-2の親が見つかったら
              node.add_child(result[j-2], j-2)
      elif result[j+1] == "[": # 1つ右が"["の場合
          # print(j)
          for node in nodes:
            if node.no == j_parent:
              node.add_child(result[j-2], j-2)
          # print(result[j-2], j-2, "←", parent, j_parent, "☆bond:", result[j-1]) # その他1:1の親子関係については通常の処理で対応可能
      else: # 初めて親が出た時
        # print(j+1)
        nodes.append(Node(j+1, result[j+1], result[j-2], j-2)) # 既存の研究
        parent_index.add(j+1)
        # print(result[j-2], j-2, "←", result[j+1], j+1, "☆bond:", result[j-1])

    # 分岐中の1:1の親子関係を認識
    elif re.fullmatch(r'[0-9a-zA-Z]+', result[j]) or re.fullmatch(r'[a-zA-z]-[0-9a-zA-Z]+', result[j]) or re.fullmatch(r'[0-9a-zA-Z]+[0-9]\/[0-9][a-zA-Z]', result[j]): # 単糖
      if re.fullmatch(r'\([a-z][0-9]-[0-9]\)', result[j-1]) or re.fullmatch(r'\([a-z][0-9]-[0-9]\/[0-9]\)', result[j-1]) or re.fullmatch(r'\([a-z][0-9]-[a-zA-Z]-[0-9]\)', result[j-1]): # 結合様式
        # print(j)
        nodes.append(Node(j, result[j], result[j-2], j-2)) # 既存の研究
        parent_index.add(j)
        # print(result[j-2], j-2, "←", result[j], j, "☆bond:", result[j-1])

    # 分岐の中の分岐に対応
    elif result[j] == "]":
      if result[j+1] == "[": # 3つ以上の分岐のときの2つ目以上の分岐
        parent = copy.deepcopy(parent) # 今の親を維持
        j_parent = copy.deepcopy(j_parent) # # 今の親を維持
      else: # 1つ目の分岐のとき
        parent = result[j+1] # 1つ右が親（参照渡しを利用してglobalに影響を与える）
        j_parent = j+1 # 1つ右が親（参照渡しを利用してglobalに影響を与える）
      # print(j)
      # print("\n＜新たなbranch＞")
      j, parent_index = branch(result, j, nodes, parent_index, parent, j_parent) # 再帰
      # print(j)
      j -= 1 # この処理で再帰後のresult[j]は[になり、次のif文に引っかからなくなる

      # 一つ左が結合様式の時（j_parent_branchの子どもの時）
      if re.fullmatch(r'\([a-z][0-9]-[0-9]\)', result[j-1]) or re.fullmatch(r'\([a-z][0-9]-[0-9]\/[0-9]\)', result[j-1]) or re.fullmatch(r'\([a-z][0-9]-[a-zA-Z]-[0-9]\)', result[j-1]): # 結合様式
        if re.fullmatch(r'[0-9a-zA-Z]+', result[j-2]) or re.fullmatch(r'[a-zA-z]-[0-9a-zA-Z]+', result[j-2]) or re.fullmatch(r'[0-9a-zA-Z]+[0-9]\/[0-9][a-zA-Z]', result[j-2]): # 単糖
          # print(j_parent_branch)
          for node in nodes:
            if node.no == j_parent:
              node.add_child(result[j-2], j-2)

      # 1つ左が"]"のとき（分岐が3つ以上になっているとき or "]["←このようになっているとき）
      if result[j-1] == "]":
        if re.fullmatch(r'\([a-z][0-9]-[0-9]\)', result[j-2]) or re.fullmatch(r'\([a-z][0-9]-[0-9]\/[0-9]\)', result[j-2]) or re.fullmatch(r'\([a-z][0-9]-[a-zA-Z]-[0-9]\)', result[j-2]): # 結合様式
          if re.fullmatch(r'[0-9a-zA-Z]+', result[j-3]) or re.fullmatch(r'[a-zA-z]-[0-9a-zA-Z]+', result[j-3]) or re.fullmatch(r'[0-9a-zA-Z]+[0-9]\/[0-9][a-zA-Z]', result[j-3]): # 単糖
            for node in nodes:
              if node.no == j_parent:
                node.add_child(result[j-3], j-3)
      # print("＜branch関数終わり＞")

    j -= 1 # どんどん前に行く

    # 葉にも対応
    # print(j+1)
    if result[j] == "[":
      # print(j+1)
      nodes.append(Node(j+1, result[j+1], None, None)) # 既存の研究

      if result[j-1] == "]": # 1つ上の階層の親に行かないようにする
        pass

  return j+1, parent_index #"["の1つ前

"""OTMM独自の処理"""
# OTMM用に兄弟関係を抽出（1人の親が1人の子供を持つようにする）
def create_siblings(nodes):
  for node in nodes:
    if len(node.child) >= 2: # 子供が2人以上いるとき
      # print(node.no, node.child)
      parent = node
      for i in range(0,len(parent.child)):
        # print(i)
        if i == 0: # 先頭の子供がi=0（すなわち、eldest（長男）のとき）
          child = parent.child[i] # 子供を入れる
          child_num = parent.child_num[i]
          # print(child_num)
          sibling = parent.child[i+1]
          sibling_num = parent.child_num[i+1]
          # print(sibling_num)
          for node in nodes:
            if node.no == child_num:
              node.add_younger(sibling, sibling_num) # 若い
              node.add_parent(parent.name, parent.no) # 親の処理
            elif node.no == sibling_num:
              node.add_elder(child, child_num) # 年上

        elif i == len(parent.child)-1: # youngestのとき
          youngest = parent.child[i]
          youngest_num = parent.child_num[i]
          #  print(youngest_num)
          elder = parent.child[i-1]
          elder_num = parent.child_num[i-1]
          #  print(elder_num)
          for node in nodes:
            if node.no == youngest_num:
               node.add_elder(elder, elder_num) # 年上のみ

        else:
          elder = parent.child[i]
          elder_num = parent.child_num[i]
          #  print(elder_num)
          younger = parent.child[i+1]
          younger_num = parent.child_num[i+1]
          #  print(younger_num)
          for node in nodes:
            if node.no == elder_num:
               node.add_younger(younger, younger_num) # 若い
            elif node.no == younger_num:
               node.add_elder(elder, elder_num) # 年上
      
      #更新する
      del parent.child[1:len(parent.child)] # 最も近い子供以外を消去
      parent.child = parent.child[0]
      del parent.child_num[1:len(parent.child_num)] # 最も近い子供以外を消去
      parent.child_num = parent.child_num[0]

    elif len(node.child) == 1: # 子供が1人の時
      node.child = node.child[0]
      node.child_num = node.child_num[0]

    elif node.child[0] == None:
      node.child = None
      node.child_num = None

  return nodes

# 親を登録（上記の兄弟の処理が終わっている前提）
def create_parent(nodes):
  for node in nodes:
    if node.child != None:
      parent = node
      child = node.child
      child_num = node.child_num
      for node in nodes:
        if node.no == child_num:
          node.add_parent(parent.name, parent.no)

  return nodes

"""単糖のインスタンスをNodeの属性（parent、child、elder、younger）に代入"""

# 単糖の名前ではなくNodeインスタンスを入れる
def set_instance(nodes):
  for node in nodes:
    # parentをNodeインスタンスに
    if node.parent != None:
      for instance in nodes:
        if instance.no == node.parent_num:
          node.parent = instance # 親のインスタンスを入れる
    # childをNodeインスタンスに
    # if node.child[0] != None: # childのみリスト
    if node.child != None: # childのみリスト
      for instance in nodes:
        # if instance.no == node.child_num[0]:
        if instance.no == node.child_num:
          # node.child[0] = instance # 子のインスタンスを入れる
          node.child = instance
    # youngerをNodeインスタンスに
    if node.younger != None:
      for instance in nodes:
        if instance.no == node.younger_num:
          node.younger = instance # youngerのインスタンスを入れる
    # elderをNodeインスタンスに
    if node.elder != None:
      for instance in nodes:
        if instance.no == node.elder_num:
          node.elder = instance # elderのインスタンスを入れる

  return nodes