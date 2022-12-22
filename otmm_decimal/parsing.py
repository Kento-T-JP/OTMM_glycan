"""
Parsing (Retrieve what was learned by retrieving the most likely state paths)
Viterbiアルゴリズムで経路探索
"""

import numpy as np
import pandas as pd

# 識別子で区別する
def parse_glycan(glycan, state_set, pi, a_a, a_b, b):
    sugar_id = []
    for sugar in glycan:
        sugar_id.append(sugar.no)
        sugar_id = sorted(sugar_id)

    """Make data variable for phi"""
    phi_up = np.zeros((len(state_set), len(glycan)))
    phi_up = pd.DataFrame(phi_up, index=state_set, columns=sugar_id)

    phi_back = np.zeros((len(state_set), len(glycan)))
    phi_back = pd.DataFrame(phi_back, index=state_set, columns=sugar_id)

    """Calculate phi"""
    # bottom-up and right-to-left dynamic programming procedure
    for node in reversed(glycan):
        for state in state_set:
            # print(node.no)
            phi_up.at[state, node.no] = calc_phi_up(state, node, phi_back, a_a, b, state_set)
            phi_back.at[state, node.no] = calc_phi_back(state, node, phi_up, phi_back, a_b, state_set)

    p = most_likely_prob(pi, phi_up, state_set, sugar_id)
    print("The probabilitythat all labels are outputted along the most likely state transition:", p)
    z = most_likely_state(glycan, state_set, sugar_id, pi, a_a, a_b, b, phi_up, phi_back)
    for node in glycan:
        print("Node:", node.name, node.no,"Most likely state:", z[node.no])
    print()

"""
Calculate psi
top-down and left-to-right manner
"""
def calc_phi_up(q, p, phi_back, a_a, b, state_set):
  #calculate phi up
  if p.child == None:
    return b.at[q, p.name]
  else:
    prob = []
    for m in state_set:
      if phi_back.at[m, p.child.no] == 0:
        print("Unexpected value at calc phi up") # 計算していないphi_backを使っている状態
      # Σの部分をmaxに
      # j is the p's eldest children (thus, j == p.child)
      prob.append(a_a[q][m] + phi_back.at[m, p.child.no]) # 対数なので足し算は掛け算
    return b.at[q, p.name] + max(prob)

def calc_phi_back(m, j, phi_up, phi_back, a_b, state_set):
  #calculate phi back
  if j.younger == None: # youngest
    if phi_up.at[m, j.no] == 0:
      print("Unexpected value at calc phi back") # 計算していないphi_upを使っている状態
    return phi_up.at[m, j.no]
  else:
    prob = []
    for l in state_set:
      # Σの部分をmaxに
      prob.append(a_b[m][l] + phi_back.at[l, j.younger.no]) # 対数なので足し算は掛け算
    return phi_up.at[m, j.no] + max(prob)

# returns the most likely state for the given node
def psi_up(q, p, phi_back, a_a, state_set):
  if p.child == None: # ψ_upにおいてpは必ず子供を持つのでこの条件は該当しないはず
    print("ERROR! node p shoud have a child!")
  else:
    prob = np.empty(0)
    for m in state_set: # m proceeds from 0(state_set[0]) to state_set[-1]
      if phi_back.at[m, p.child.no] == 0:
        print("Unexpected value at psi up") # 計算していないbackを使っている状態
      # Σの部分をargmaxに
      # j is the p's eldest children (thus, j == p.child)
      # print(q, m, p.child.no)
      prob = np.insert(prob, m, a_a[q][m] + phi_back.at[m, p.child.no]) # 対数なので足し算は掛け算
    return np.argmax(prob) # m represents a state and return it

# returns the most likely state for the given node
def psi_back(m, j, phi_back, a_b, state_set):
  if j.younger == None: # ψ_backにおいてjは必ずyoungerを持つのでこの条件は該当しないはず
    print("ERROR! j shoud have a younger sibling!")
  else:
    prob = np.empty(0)
    for l in state_set: # l proceeds from 0(state_set[0]) to state_set[-1]
      # Σの部分をargmaxに
      prob = np.insert(prob, l, a_b[m][l] + phi_back.at[l, j.younger.no]) # 対数なので足し算は掛け算
    return np.argmax(prob) # l represents a state and return it

def most_likely_prob(pi, phi_up, state_set, sugar_id):
  likely = []
  for l in state_set:
    likely.append(pi[l] + phi_up.at[l, sugar_id[-1]])

  most_likely = max(likely)
  return most_likely

"""Calculate the most likely state for the given node"""
def most_likely_state(glycan, state_set, sugar_id, pi, a_a, a_b, b, phi_up, phi_back):
  z = dict()
  # top-down and left-to-right manner
  for node in glycan:
    # print(node.no)
    if node.no == sugar_id[-1]: # root
      prob = np.empty(0)
      for l in state_set:
        prob = np.insert(prob, l, pi[l] + phi_up.at[l, node.no]) # node.no==sugar_id[-1]
      z[node.no] = np.argmax(prob) # return the most likely state of the root node
    elif node.elder == None and node.parent != None:
      z[node.no] = psi_up(z[node.parent.no], node.parent, phi_back, a_a, state_set)
    else:
      z[node.no] = psi_back(z[node.elder.no], node.elder, phi_back, a_b, state_set)
  return z