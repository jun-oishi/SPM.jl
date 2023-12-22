# SPM

## できること

- HDRファイルの読み込み/書き出し
- 表面形状のプロット
- 探針形状の再構成

## 環境構築

WinsdowsならWSLを使うのがおすすめ

### juliaのインストール

https://julialang.org/downloads/ を参照

### juliaのパッケージのインストール

このディレクトリでjuliaのREPLを起動し、以下のコマンドを実行する。
`Project.toml`に記述された依存パッケージがインストールされる。
この中に含まれる`MDToolbox`がMatsunagaによる探針再構成のプログラム。

```bash
julia> ]                    # `]`キーを押せばpkgモードに入る
(@v1.9) pkg> activate .     # このディレクトリで有効な環境を作成
(SPM) pkg> instantiate      # Project.tomlに記述された依存パッケージをインストール
```

完了すると`Manifest.toml`が作成され実際にインストールされたパッケージの情報が詳細に記述される。
この環境のためにインストールされたパッケージを使用するにはソースコード中で
```julia
import Pkg          # Pkgモジュールをインポート
Pkg.activate(".")   # このディレクトリの環境を有効化
import Plots        # Plotsパッケージをインポート
```
などとする。REPLではパッケージモードに入って操作してもよい。

#### 各パッケージの説明

- `MDToolbox` : Matsunagaによる探針再構成のプログラム
  - `cuDNN` : GPUを使うためのパッケージ(MDToolboxの依存パッケージ)
- `Plots`, `PyPlot` : プロット用
- `IJulia` : jupyter notebookでJuliaを使うためのパッケージ
- `Statistics` : 統計量の計算
- `Flux` : 機械学習のためのパッケージ(differentiable BTRに必要)


### jupyter notebookのインストールと設定

pythonをインストールして

```bash
pip install jupyter
```

とかすればいいです( https://jupyter.org/install も参照)。Anacondaを使っている場合は少し違うらしいです。

jupyterが使えるようになったら自動的にPythonに加えてJuliaのカーネルも使えるようになっているはず(juliaの`IJulia`パッケージが橋渡しになるらしい)。

## 使い方

`test.ipynb`も参照のこと。

### モジュールの構成

`src/SPM.jl`内で`SPM`モジュールを定義しているので、これを`include`すれば`SPM.xxx`でアクセスできる。
他のファイルはサブモジュールを定義していて`src/SPM.jl`はそれらを`include`しているので`SPM.SUBMODULE.関数`で個々のファイルで定義した関数などにアクセスできる。ただし`SPMCore`で定義した関数、構造体は`SPM`モジュールの中で`using`としているので`SPM.関数`でアクセスできる。

#### SPMCore (spm_core.jl)

基本的な機能(表面形状の表現、表面形状に対する基本的な操作)を定義している。

#### HDR (hdr.jl)

HDRファイルの読み込み/書き出しのための関数などを定義している。

#### SPMPlots (plots.jl)

表面形状の可視化などを行う関数などを定義している。

#### BTR (btr.jl)

探針形状の再構成を行う関数などを定義している。

## 不具合など

### 再構成が遅い

探針サイズが大きくなると飛躍的に計算量が増えてなかなか終わらなくなる。
10x10くらいで動作確認をして計算時間を見ながら大きくしていくのが良いと思います。

[マルチスレッド化](https://docs.julialang.org/en/v1/manual/multi-threading/)や[GPUの使用](https://fluxml.ai/Flux.jl/stable/gpu/)(differentiableの場合)で高速化できるかもしれないので、興味があれば調べてみてください。

### 再構成した探針の保存

再構成の条件を`remark`などに入れるようにしているはずなのに、なぜか保存されません。