# 物理ベースアニメーション

![ウブンツ](https://github.com/poly-bear/PBA/actions/workflows/ubuntu.yml/badge.svg)
![マック](https://github.com/poly-bear/PBA/actions/workflows/mac.yml/badge.svg)

https://www.youtube.com/channel/UCsTJawTc_A88lmAUU2DJ1XQ

## ディレクトリ構成

* 3rdparty:  
    サードパーティライブラリ
* opengl:  
    レンダラ with OpenGL
* project2d:  
    物理シミュレーション等の実装(二次元)
* project3d:  
    物理シミュレーション等の実装(三次元)
* projectcuda:  
    物理シミュレーション等の実装 using CUDA
* resource:  
    メッシュ等のリソース
* utils:  
    数学関数等のユーティリティ
* utilscuda:  
    数学関数等のユーティリティ using CUDA
* test:  
    テスト用

## ビルド・実行に必要なもの

- GLFW
- OpenMP
- Eigen (Global Solverを扱う場合)

例:project3d/Demoを実行する
```
$ cd project3d/Demo
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./Program
```

projectcuda以下

- CUDA ツールキット(11.4)

## 操作方法

右マウスボタンでドラッグ:カメラの位置を変更  
スクロール:拡大・縮小

## リポジトリに含まれるサードパーティ

### ライブラリ

* glad
* imgui

### 形状データ

[README.md](/resource/README.md)

を確認されたい。
