"""
主にlikelihoodとlearningで使用するアルゴリズム
新規性のある研究も既存の研究もこの部分は同じ
"""

import math
import pandas as pd
import copy # Pythonでは関数においてミュータブルな変数は参照渡しのため
# import gc
from decimal import Decimal

"""
log(x+y)の近似
sum of logs
"""
# math.exp(0) = 1なので条件分岐する必要ない
# log(1+e^(y-x))をすることで、e^(y-x)が非常に小さい場合log(1+e^(y-x))≒0を出力できる→ほぼxが返される
# とにかく、一番大きいはずの0が出力されない
def smoothmax(x, y):
  try:
    return x + Decimal(str(math.log(Decimal(str(1)) + Decimal(str(math.exp(y-x)))))) # exp(±710)くらいでエラーが発生
  except: 
    # logsumexpと同じ処理
    if y-x > 0: # e^(y-x)が非常に大きくなるので1+ e^(y-x)≒e^(y-x)とする
      return copy.deepcopy(y)
    elif y-x < 0:
      return copy.deepcopy(x)

# math.expの指数の範囲を狭めたが特に変化なし
# def smoothmax(x, y):
#   if abs(y-x) <= 20:
#     return x + math.log(1 + math.exp(y-x))
#   else: 
#     # logsumexpと同じ処理
#     if y-x > 0: # e^(y-x)が非常に大きくなるので1+ e^(y-x)≒e^(y-x)とする
#       return y
#     elif y-x < 0:
#       return x

"""Calculate upward and backward"""
def calc_up(q, p, back, a_a, b, state_set):
  #calculate up
  if p.child == None:
    return b.at[q, p.name]
  else:
    total = Decimal(str(0))
    for m in state_set:
      if back.at[m, p.child.no] == 100: # あり得ない値（100）を発見
        print("100(unexpected value) found at up") # 計算していないbackを使っている状態
      # j is the p's eldest children (thus, j == p.child)
      if total == 0 and m == state_set[0]: # 最初の値はそのまま代入
        total += a_a[q][m] + back.at[m, p.child.no] # 対数なので足し算は掛け算
      else:
        total = copy.deepcopy(smoothmax(total, a_a[q][m] + back.at[m, p.child.no])) #totalは和なのでsmoothmaxにtotalを代入すればよい
    return b.at[q, p.name] + total

def calc_back(m, j, up, back, a_b, state_set):
  #calculate back
  # print(back)
  if j.younger == None: # youngest child and root
    if up.at[m, j.no] == 100:
      print("100(unexpected value) found at back") # 計算していないupを使っている状態
    return up.at[m, j.no]
  else:
    total = Decimal(str(0))
    for l in state_set:
      if total == 0 and l == state_set[0]: # 最初の値はそのまま代入
        total += a_b[m][l] + back.at[l, j.younger.no] # 対数なので足し算は掛け算
      else:
        total = copy.deepcopy(smoothmax(total, a_b[m][l] + back.at[l, j.younger.no])) #totalは和なのでsmoothmaxにtotalを代入すればよい
    return up.at[m, j.no] + total

def calc_likelihood(df, pi, a_a, a_b, b, state_set):
  i = 0
  likelihood = []
  for row in df.itertuples():
    # count how many glycans did we use
    if i == 0:
      print("Number of glycans in likelihood", i)
    elif (i+1)%100 == 0:
      print("Number of glycans in likelihood", i+1) # 0からインデックスが始まるため

    glycan = row[2]

    # 識別子で区別する
    sugar_id = []
    for sugar in glycan:
      sugar_id.append(sugar.no)
    sugar_id = sorted(sugar_id)

    # make upward prob
    up = [[Decimal(str(100))] * len(glycan) for i in [Decimal(str(100))] * len(state_set)]
    up = pd.DataFrame(up, index=state_set, columns=sugar_id)
    
    # make backward prob
    back = [[Decimal(str(100))] * len(glycan) for i in [Decimal(str(100))] * len(state_set)]
    back = pd.DataFrame(back, index=state_set, columns=sugar_id)

    for node in reversed(glycan):
      for state in state_set:
        # print(node.no)
        up.at[state, node.no] = copy.deepcopy(calc_up(state, node, back, a_a, b, state_set))
        # print(type(up.at[state, node.no]))
        back.at[state, node.no] = copy.deepcopy(calc_back(state, node, up, back, a_b, state_set))
        # print(type(back.at[state, node.no]))
    # print(up, "\n")
    
    like = Decimal(str(0))
    for l in state_set:
      if like == 0 and l == state_set[0]:
        like += pi[l] + up.at[l, sugar_id[-1]] # 最もインデックスが大きい数（[-1]）は必ずルート
      else:
        like = copy.deepcopy(smoothmax(like, pi[l] + up.at[l, sugar_id[-1]]))
    likelihood.append(like)
    i += 1

    # メモリリーク（対策）
    # del up, back
    # gc.collect()
  
  return likelihood

"""
Learning(EMアルゴリズム)
"""

"""Forward probability"""
def calc_forward(l, j, down, forward, up, a_a, a_b, b, pi, state_set):
  #calculate forward
  if j.parent == None and j.elder == None and j.younger == None: # when j == root
    # if not re.fullmatch(r'[0-9a-zA-Z]+', j.name): # ルートではないときTrue
    #   print(j.name, j.no)
    # return pi[l]
    return Decimal(str(0)) # 論文で言及されていない部分（とりあえず確率は1でいいらしい）, log(1) = 0
  elif j.elder == None: # eldest children
    total = Decimal(str(0))
    for q in state_set:
      if down.at[q, j.parent.no] == 100:
        print("100(unexpected value) found at forward!")
      if total == 0 and q == state_set[0]:
        total += a_a[q][l] + down.at[q, j.parent.no] + b.at[q, j.parent.name]
      else:
        total = copy.deepcopy(smoothmax(total, a_a[q][l] + down.at[q, j.parent.no] + b.at[q, j.parent.name]))
    return total
  else:
    total = Decimal(str(0))
    for m in state_set:
      if forward.at[m, j.elder.no] == 100:
        print("100(unexpected value) found at forward!!")
      if total == 0 and m == state_set[0]:
        total += a_b[m][l] + forward.at[m, j.elder.no] + up.at[m, j.elder.no] # 対数なので足し算は掛け算
      else:
        total = copy.deepcopy(smoothmax(total, a_b[m][l] + forward.at[m, j.elder.no] + up.at[m, j.elder.no]))
    return total

"""Downward probability"""
def calc_down(l, j, forward, back, pi, a_b, state_set):
  #calculate downward
  if j.parent == None and j.elder == None and j.younger == None: # when j == root
    return pi[l]
  elif j.younger == None:
    if forward.at[l, j.no] == 100:
      print("100(unexpected value) found at down!")
    return forward.at[l, j.no]
  else:
    total = Decimal(str(0))
    for m in state_set:
      if forward.at[l, j.no] == 100:
        print("100(unexpected value) found at down!!")
      if total == 0 and m == state_set[0]:
        total += a_b[l][m] + back.at[m, j.younger.no] # 対数なので足し算は掛け算
      else:
        total = copy.deepcopy(smoothmax(total, a_b[l][m] + back.at[m, j.younger.no]))
    return forward.at[l, j.no] + total

"""Learning（EMアルゴリズム）"""
def EM(df, pi_original, a_a_original, a_b_original, b_original, state_set, label_set, L, epsilon):
  t = 0  
  L_all = [] # 全体の期待値
  print("likelihood before learning", L, "\n")
  L_all.append(L)

  # 参照渡しのため値を別のアドレス（変数）にコピー
  pi = copy.deepcopy(pi_original)
  a_a = copy.deepcopy(a_a_original)
  a_b = copy.deepcopy(a_b_original)
  b = copy.deepcopy(b_original)

  while True:
    print("t(epoc)=", t)
    # L_T = [] # record all likelihoods of T_u

    # 5 of pseudo code
    # initialize mu of T_u
    mu_aa_t = [[Decimal(str(0))] * len(state_set) for i in [Decimal(str(0))] * len(state_set)]

    mu_ab_t = [[Decimal(str(0))] * len(state_set) for i in [Decimal(str(0))] * len(state_set)]

    B_t = [[Decimal(str(0))] * len(label_set) for i in [Decimal(str(0))] * len(state_set)]
    mu_b_t = pd.DataFrame(B_t, index=state_set, columns=label_set)
        
    mu_pi_t = [Decimal(str(0))] * len(state_set)

    # 6 of pseudo code
    count = 0
    for row in df.itertuples():
      if count == 0:
        print("Number of glycans in learning", count)
      elif (count+1)%100 == 0:
        print("Number of glycans in learning", count+1) # 0からインデックスが始まるため

      # 7 of pseudo code
      glycan = row[2]

      # initialize mu of the glycan
      mu_aa_u = [[Decimal(str(0))] * len(state_set) for i in [Decimal(str(0))] * len(state_set)]

      mu_ab_u = [[Decimal(str(0))] * len(state_set) for i in [Decimal(str(0))] * len(state_set)]

      B_u = [[Decimal(str(0))] * len(label_set) for i in [Decimal(str(0))] * len(state_set)]
      mu_b_u = pd.DataFrame(B_u, index=state_set, columns=label_set)
        
      mu_pi_u = [Decimal(str(0))] * len(state_set)

      # 識別子で区別する
      sugar_id = []
      for sugar in glycan:
        sugar_id.append(sugar.no)
      sugar_id = sorted(sugar_id)
      # make upward prob
      up = [[Decimal(str(100))] * len(glycan) for i in [Decimal(str(100))] * len(state_set)]
      up = pd.DataFrame(up, index=state_set, columns=sugar_id)
      
      # make backward prob
      back = [[Decimal(str(100))] * len(glycan) for i in [Decimal(str(100))] * len(state_set)]
      back = pd.DataFrame(back, index=state_set, columns=sugar_id)

      # make forward prob
      forward = [[Decimal(str(100))] * len(glycan) for i in [Decimal(str(100))] * len(state_set)]
      forward = pd.DataFrame(forward, index=state_set, columns=sugar_id)

      # make downward prob
      down = [[Decimal(str(100))] * len(glycan) for i in [Decimal(str(100))] * len(state_set)]
      down = pd.DataFrame(down, index=state_set, columns=sugar_id)

      # U, Bを計算する(7 to 9 of pseudo code)
      for node in reversed(glycan):
        for state in state_set:
          # print(node.no)
          up.at[state, node.no] = copy.deepcopy(calc_up(state, node, back, a_a, b, state_set))
          back.at[state, node.no] = copy.deepcopy(calc_back(state, node, up, back, a_b, state_set))

      # D, Fを計算する(10 to 12 of pseudo code)
      for node in glycan:
        for state in state_set:
          # print(node.no)
          forward.at[state, node.no] = copy.deepcopy(calc_forward(state, node, down, forward, up, a_a, a_b, b, pi, state_set))
          down.at[state, node.no] = copy.deepcopy(calc_down(state, node, forward, back, pi, a_b, state_set))

      # calculate L(T_u), 尤度
      L_u = Decimal(str(0)) # 入力データ1つ1つの尤度を保存
      # L_u2 = 0
      i = glycan[0].no # optional
      # print(i)
      for l in state_set:
        if L_u == 0 and l == state_set[0]:
          L_u += up.at[l, i] + down.at[l, i]
        else:
          L_u = copy.deepcopy(smoothmax(L_u, up.at[l, i] + down.at[l, i]))
        # if L_u2 == 0:
        #   L_u2 = pi_em[l] + up_em.at[l, sugar_id[-1]]
        # else:
        #   L_u2 = smoothmax(L_u2, pi_em[l] + up_em.at[l, sugar_id[-1]])
      # print(L_u, L_u2) # L_u = L_u2だった
      # L_T.append(L_u)

      # 期待値の計算（13 of pseudo code）
      # calculate mu_aa_u[q][l]
      for q in state_set:
        for l in state_set:
          exist_child = False # 該当するchildがいるかどうか
          total = Decimal(str(0))
          for p in glycan:
            if p.child != None:
              exist_child = True
              # p.child.no means j in the thesis of OTMM
              if total == 0:
                total += down.at[q, p.no] + b.at[q, p.name] + a_a[q][l] + back.at[l, p.child.no]
              else:
                total = copy.deepcopy(smoothmax(total, down.at[q, p.no] + b.at[q, p.name] + a_a[q][l] + back.at[l, p.child.no]))
          if exist_child == True or total != 0: # 本来total×L_uなのでtotal=0の時は0（対数変換によって掛け算が足し算になっている）
            mu_aa_u[q][l] = total - L_u # 対数の引き算は割り算
          else:
            mu_aa_u[q][l] = Decimal(str(0)) # 0にしてパラメータの更新に影響が出ないようにする
          # print(total, L_u)
      
      # calculate mu_ab_u[q][l]
      for q in state_set:
        for l in state_set:
          exist_younger = False # younger siblingがいるかどうか
          total = Decimal(str(0))
          for j in glycan:
            if j.younger != None: # X(j) != empty set
              exist_younger = True
              # p.child.no means j in the thesis of OTMM
              if total == 0:
                total += forward.at[q, j.no] + a_b[q][l] + back.at[l, j.younger.no] + up.at[q, j.no]
              else:
                total = copy.deepcopy(smoothmax(total, forward.at[q, j.no] + a_b[q][l] + back.at[l, j.younger.no] + up.at[q, j.no]))
          if exist_younger == True or total != 0:
            mu_ab_u[q][l] = total - L_u # 対数の引き算は割り算
          else:
            mu_ab_u[q][l] = Decimal(str(0))

      # calculate mu_b_u[m][o_h]
      for m in state_set:
        for o_h in label_set:
          exist_oh = False # whether the sugar(o_h) exists or not
          total = Decimal(str(0))
          for i in glycan:
            if i.name == o_h:
              exist_oh = True
              if total == 0:
                total += down.at[m, i.no] + up.at[m, i.no]
              else:
                total = copy.deepcopy(smoothmax(total, down.at[m, i.no] + up.at[m, i.no]))
          if exist_oh == True or total != 0:
            mu_b_u.at[m, o_h] = total - L_u # 対数の引き算は割り算
          else:
            mu_b_u.at[m, o_h] = Decimal(str(0)) # 入力データには単糖o_hが存在しないので期待値は0

      # calculate  mu_pi_u[m]
      for m in state_set:
        mu_pi_u[m] = (pi[m] + up.at[m, sugar_id[-1]]) - L_u # 最も大きい識別子がルートノード

      # for each param, do mu_t = mu_t + mu_u (14 of pseudo code)
      # mu_aa_t = mu_aa_t + mu_aa_u
      for q in state_set:
        for l in state_set:
          if mu_aa_u[q][l] == 0:
            mu_aa_t[q][l] += Decimal(str(0))
          else:
            if mu_aa_t[q][l] == 0:
              mu_aa_t[q][l] += copy.deepcopy(mu_aa_u[q][l])
            else:
              mu_aa_t[q][l] = copy.deepcopy(smoothmax(mu_aa_t[q][l], mu_aa_u[q][l]))

      # mu_ab_t = mu_ab_t + mu_ab_u
      for q in state_set:
        for l in state_set:
          if mu_ab_u[q][l] == 0:
            mu_ab_t[q][l] += Decimal(str(0))
          else:
            if  mu_ab_t[q][l] == 0:
              mu_ab_t[q][l] += copy.deepcopy(mu_ab_u[q][l])
            else:
              mu_ab_t[q][l] = copy.deepcopy(smoothmax(mu_ab_t[q][l], mu_ab_u[q][l]))

      # mu_b_t = mu_b_t + mu_b_u
      for m in state_set:
        for o_h in label_set:
          if mu_b_u.at[m, o_h] == 0:
            mu_b_t.at[m, o_h] += Decimal(str(0))
          else:
            if mu_b_t.at[m, o_h] == 0:
              mu_b_t.at[m, o_h] += copy.deepcopy(mu_b_u.at[m, o_h])
            else:
              mu_b_t.at[m, o_h] = copy.deepcopy(smoothmax(mu_b_t.at[m, o_h],  mu_b_u.at[m, o_h]))

      # mu_pi_t = mu_pi_t + mu_pi_u
      for m in state_set:
        if mu_pi_u[m] == 0:
          mu_pi_t[m] += Decimal(str(0))
        else:
          if mu_pi_t[m] == 0:
            mu_pi_t[m] += copy.deepcopy(mu_pi_u[m])
          else:
            mu_pi_t[m] = copy.deepcopy(smoothmax(mu_pi_t[m], mu_pi_u[m]))
      
      # メモリの開放（メモリリーク対策）
      # del up, back, forward, down
      # del mu_pi_u, mu_aa_u, mu_ab_u, mu_b_u
      # gc.collect()
      count += 1

    # update a_a
    # print(mu_aa_t)
    for q in state_set:
      denominator = Decimal(str(0))
      for l_dash in state_set:
        if denominator == 0 and l_dash == state_set[0]:
          denominator += copy.deepcopy(mu_aa_t[q][l_dash])
        else:
          denominator = copy.deepcopy(smoothmax(denominator, mu_aa_t[q][l_dash]))
      for l in state_set:
        numerator = copy.deepcopy(mu_aa_t[q][l])
        a_a[q][l] = numerator - denominator # 対数の引き算は割り算

    # update a_b
    # print(mu_ab_t)
    for q in state_set:
      denominator = Decimal(str(0))
      for l_dash in state_set:
        if denominator == 0 and l_dash == state_set[0]:
          denominator += copy.deepcopy(mu_ab_t[q][l_dash])
        else:
          denominator = copy.deepcopy(smoothmax(denominator, mu_ab_t[q][l_dash]))
      for l in state_set:
        numerator = copy.deepcopy(mu_ab_t[q][l])
        a_b[q][l] = numerator - denominator

    # update b
    for m in state_set:
      denominator = Decimal(str(0))
      for o_i in label_set:
        if denominator == 0 and o_i == label_set[0]:
          denominator += copy.deepcopy(mu_b_t.at[m, o_i])
        else:
          denominator = copy.deepcopy(smoothmax(denominator, mu_b_t.at[m, o_i]))
      for o_h in label_set:
        numerator = copy.deepcopy(mu_b_t.at[m, o_h])
        # print(numerator, denominator)
        b.at[m, o_h] = numerator - denominator

    # update pi
    for m in state_set:
      numerator = copy.deepcopy(mu_pi_t[m])
      denominator = Decimal(str(0))
      for k in state_set:
        if denominator == 0 and k == state_set[0]:
          denominator += copy.deepcopy(mu_pi_t[k])
        else:
          denominator = copy.deepcopy(smoothmax(denominator, mu_pi_t[k]))
      pi[m] = numerator - denominator

    t += 1
    # calculate L_t(T) using the param updated (17 of the pseudo code)
    likelihoods = calc_likelihood(df, pi, a_a, a_b, b, state_set)
    # print(a_a_em[0])
    # print(a_b_em[0])
    # print(pi_em)
    # print(b_em)
    # print(likelihoods)
    L_all.append(sum(likelihoods)) # 対数（likelihood）の足し算は掛け算
    
    # L_all.append(sum(L_T)) # 現在のパラメータを使って計算をしたいので
    print("Likelihood", L_all[t])
    print("Difference of the likelihoods", L_all[t] - L_all[t-1], "\n")
    if abs(L_all[t] - L_all[t-1]) < epsilon: # or t == 5: # 止まらなかったら困るのでとりあえずt==5で止める
      return pi, a_a, a_b, b, L_all