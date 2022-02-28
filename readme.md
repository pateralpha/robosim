# File List

- simulation.cpp: シミュレーションのメインプログラム
- machine.h: マシン関連のクラスとか変数が押し込められたファイル
- ip.h, ip.cpp: 補間とかの関数がつまっているかもしれない
- la.h: Linear Algebraという名のベクトルと行列の簡易ライブラリ

<br>

- view_cir.plt: 直線・円弧補間での描画
- view_lag.plt: ラグランジュ補間での描画
- その他 *.plt: なんかグラフ出る

<br>

# compile and run

g++でコンパイルする場合は，cppファイルが2つなので，

```bash
g++ simulation.cpp ip.cpp
```
とかするとa.outが生成されると思います．\
環境によっては-std=c++17か-std=c++14みたいなオプションがいるかもしれない．

コンパイルがでけたら
```bash
./a.out
```
で実行してあげるとa.datとc.datが出てくる予定です．

この状態でgnuplot開いて"view_cir.plt"とか"view_lag.plt"とかをロードしてもらうと描画される...といいなと思います．

<br>

因みにテストした環境はWSL2のdeiban 11.2だと思います（確認方法に自信がない）． \
あとgnuplotは5.2以降じゃないと動かないです．

<br>

Linux全然わからないので依存関係がよくわかってません．許してください．
