# Glycan
This is a repository regarding OTMM(Ordered Tree Markov Model) in glycoinfomatics.

All of the programs are written in Python.

Used for the senior thesis.

# 実行結果の説明
各csvファイルとtxtファイルが実行結果になります。

全て対数変換された値で格納されています。

## 名前
*  piは*Initial state probability* π

*  a_aは*State transition probability(parent→child)* α

*  a_bは*State transition probability(elder sibling→younger sibling)* β

*  bは*Label output probability* b

*  timeは学習に費やした時間

*  likelihoodは入力データ全体に対するモデルの尤度

## ファイル名
*  pi_学習した糖鎖データ数_しきい値(ε)_状態数.csv

*  a_a_学習した糖鎖データ数_しきい値(ε)_状態数(parent)_状態数(child).csv

*  a_b_学習した糖鎖データ数_しきい値(ε)_状態数(elder)_状態数(younger).csv

*  b_学習した糖鎖データ数_しきい値(ε)_状態数_ラベル（単糖）の数.csv

*  time_学習した糖鎖データ数_しきい値(ε)_状態数_ラベルの数.txt

*  likelihood_学習した糖鎖データ数_しきい値(ε)_状態数_ラベルの数.csv

## 内容の説明
piは、列番号が状態の名前です。

a_aとa_bは、行番号と列番号が状態の名前です。行番号から列番号へ遷移します（例：0行1列なら状態0から状態1への遷移）。

bはインデックスに状態名が、カラムに単糖名が入っています。
