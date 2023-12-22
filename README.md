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

```julia
julia> ]                    # `]`キーを押せばpkgモードに入る
(@v1.9) pkg> activate .     # このディレクトリで有効な環境を作成
(SPM) pkg> instantiate    # Project.tomlに記述された依存パッケージをインストール
```

### jupyter notebookのインストールと設定

pythonをインストールして

```bash
pip install jupyter
```

とかすればいいです( https://jupyter.org/install も参照)。Anacondaを使っている場合は少し違うらしいですね。

jupyterが使えるようになったら自動的にPythonに加えてJuliaのカーネルも使えるようになっているはず(juliaの`IJulia`パッケージが橋渡しになるらしい)。

## 使い方

`test.ipynb`を参照のこと。

## 不具合など

### 再構成が遅い

探針サイズが大きくなると飛躍的に計算量が増えてなかなか終わらなくなる。
10x10くらいで動作確認をして計算時間を見ながら大きくしていくのが良いと思います。

[マルチスレッド化](https://docs.julialang.org/en/v1/manual/multi-threading/)や[GPUの使用](https://fluxml.ai/Flux.jl/stable/gpu/)(differentiableの場合)で高速化できるかもしれないので、興味があれば調べてみてください。

### 再構成した探針の保存

再構成の条件を`remark`などに入れるようにしているはずなのに、なぜか保存されません。