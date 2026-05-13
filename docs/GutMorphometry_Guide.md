# GutMorphometry — インストール・ユーザーガイド

マウス腸管（スイスロール標本）H&E切片の形態計測ツール。
Fiji マクロ + Python スクリプトで構成され、絨毛高・絨毛幅・陰窩深・筋層厚を半自動計測する。

---

## 目次

1. [必要環境](#1-必要環境)
2. [インストール手順](#2-インストール手順)
3. [初回動作確認](#3-初回動作確認)
4. [計測ワークフロー](#4-計測ワークフロー)
5. [絨毛の計測（villus）](#5-絨毛の計測-villus)
6. [陰窩の計測（crypt_depth）](#6-陰窩の計測-crypt_depth)
7. [筋層計測の詳細（muscularis）](#7-筋層計測の詳細-muscularis)
8. [キャリブレーション GUI の使い方](#8-キャリブレーション-gui-の使い方)
9. [パラメータ一覧](#9-パラメータ一覧)
10. [アルゴリズム解説](#10-アルゴリズム解説)
11. [出力ファイルの説明](#11-出力ファイルの説明)
12. [トラブルシューティング](#12-トラブルシューティング)

---

## 1. 必要環境

| ソフトウェア | バージョン | 備考 |
|---|---|---|
| Windows | 10 / 11 | macOS 非対応 |
| Fiji (ImageJ2) | 最新版推奨 | https://fiji.sc |
| Anaconda / Miniconda | 任意の最新版 | Python 3.8 以上が |
| Python パッケージ | 下記参照 | base 環境か専用環境にインストール |

**必要な Python パッケージ（標準的な Anaconda base 環境に含まれる）：**
```
numpy
scipy
scikit-image
matplotlib
```

---

## 2. インストール手順

### Step 1 — Fiji のインストール

1. https://fiji.sc から Windows 版をダウンロード・展開する
2. `Fiji.app` フォルダを任意の場所（例：`C:\Users\ユーザー名\Fiji.app`）に置く
3. `Fiji.app\ImageJ-win64.exe` を実行して起動確認する

### Step 2 — Python パッケージの確認・インストール

Anaconda Prompt を開き、必要なパッケージが入っているか確認する：

```bash
python -c "import numpy, scipy, skimage, matplotlib; print('OK')"
```

`OK` と表示されれば問題なし。エラーが出た場合はインストールする：

```bash
pip install numpy scipy scikit-image matplotlib
```

### Step 3 — ツールファイルのコピー

以下の 4 ファイルを対象 PC にコピーする：

| ファイル | コピー先 |
|---|---|
| `Gut_Morphometry.ijm` | `[Fiji.app]\plugins\` |
| `musc_thickness.py` | `[Fiji.app]\scripts\Plugins\AutoRun\` |
| `tune_muscularis.py` | `[Fiji.app]\scripts\Plugins\AutoRun\` |
| `villus_width.py` | `[Fiji.app]\scripts\Plugins\AutoRun\` |

> **AutoRun フォルダについて**  
> `[Fiji.app]\scripts\Plugins\AutoRun\` が存在しない場合は作成する。  
> Fiji は起動時にこのフォルダ内の `.py` ファイルを自動実行しようとするが、  
> 各 `.py` ファイルには Jython ガードが組み込まれており、  
> Fiji 起動時に実行されても何もしない（エラーにならない）。

### Step 4 — Python コマンドの設定

`Gut_Morphometry.ijm` の先頭近くにある以下の変数を、対象 PC の Python パスに合わせて編集する：

```javascript
var g_python_cmd = "python";
// Anaconda のフルパスを指定する場合（推奨）:
// var g_python_cmd = "C:\\Users\\ユーザー名\\anaconda3\\python.exe";
// 仮想環境を使う場合:
// var g_python_cmd = "C:\\Users\\ユーザー名\\anaconda3\\envs\\gut\\python.exe";
```

`python` のまま動く環境もあるが、Fiji が `PATH` を正しく引き継がない場合があるため、  
フルパスの指定が最も確実。

**フルパスの調べ方（Anaconda Prompt で実行）：**
```bash
where python
```

### Step 5 — 設定ファイルフォルダの確認（自動生成）

初回計測後、Fiji のインストールフォルダ直下に以下のフォルダが自動作成される：

```
[Fiji.app]\GutMorphometry_config\
    musc_params.json   ← 筋層キャリブレーション設定（参照色・e_factor）
    pixel_scale.txt    ← スケール保存値（µm/px）
```

例：`C:\Users\ユーザー名\Fiji.app\GutMorphometry_config\`

手動で作成する必要はない。このフォルダは同じ PC 上の全セッションで共有される。

> **なぜ AutoRun フォルダに置かないのか**  
> Fiji は起動時に `scripts\Plugins\AutoRun\` 内の全ファイルをスクリプトとして実行しようとする。  
> `.txt` ファイルもその対象になるため、`pixel_scale.txt` を AutoRun に置くと  
> 起動のたびに Log ウィンドウに数値が表示されてしまう。  
> `GutMorphometry_config\` はその問題を避けるための専用フォルダ。

---

## 3. 初回動作確認

1. Fiji を起動する
2. `Plugins > Gut Morphometry [G]` を選択（または `G` キー）
3. H&E 画像を開く（または開いている状態で起動）
4. セッションセットアップダイアログが表示されれば起動成功

---

## 4. 計測ワークフロー

### 画像ファイル名の規則

マクロはファイル名からメタデータを自動解析する：

```
[動物ID]_[部位コード]_[切片番号].tif

例: 622_U_1.tif
    └─ 動物622、上部小腸（U）、切片1
```

**部位コード：**

| コード | 部位 |
|---|---|
| `U` | 上部小腸（Upper small intestine） |
| `M` | 中部小腸（Mid small intestine） |
| `L` | 下部小腸（Lower small intestine） |
| `C` | 大腸（Colon） |

部位コードの判定は**前方一致**で行われる。`M1`、`U2` のようにサフィックス付きのコードも正しく認識される。  
例：`281_M1_3.tif` → 部位コード `M1` → 中部小腸として扱われる（`villus` 計測が有効）。  
ファイル名が規則に合わない場合は、起動時にダイアログで手動入力できる。

### セッションセットアップ

マクロ起動時に以下を設定する：

| 項目 | 説明 |
|---|---|
| スケールキャリブレーション | スケールバー（500 µm）に合わせて線を引き、実寸を設定する。スケールが既に設定済みの場合はスキップ可 |
| ループ位置 | スイスロールの巻き位置（outer / middle / inner）。セッション中に変更可能 |

### メインメニュー

計測ループのたびにダイアログが表示される：

```
┌──────────────────────────────────────────────────────────┐
│  Measurement type:                     [villus        ▼] │
│  Loop position:                        [outer         ▼] │
│  Measurements per session (villus / crypt): [10]         │
│  ☐ Pick background colour (villus width segmentation)   │
│  Notes (optional): [_________________________________]   │
│  ☑ Calibrate muscularis threshold　　　　　　　　　　　　　│
│                                          [OK]  [Cancel]  │
└──────────────────────────────────────────────────────────┘
```

| 項目 | 説明 |
|---|---|
| Measurement type | 計測タイプを選ぶ（下表参照） |
| Loop position | `outer` / `middle` / `inner` — 同一画像内でいつでも切り替え可能 |
| Measurements per session | 1 セッションで連続計測する本数（デフォルト 10）。`villus` と `crypt_depth` に適用される |
| Pick background colour | ON にするとバックグラウンドピッカーが起動する（後述）。デフォルト OFF |
| Notes | 任意のメモ。CSV の Notes カラムに追記される |
| Calibrate muscularis threshold | 筋層計測時のみ有効。ON のとき調整 GUI を起動する |

**計測タイプ一覧（小腸）：**

| タイプ | ツール | 内容 |
|---|---|---|
| `villus` | ポリライン → ポリゴン | 絨毛高（中心線）＋絨毛幅（境界線）を一括計測 |
| `crypt_depth` | 直線 | 陰窩の開口から底部まで |
| `muscularis` | ポリゴン（自動） | 筋層を囲んで Python で自動厚み計算 |
| `--- Save and Exit ---` | — | オーバーレイ画像を保存してセッション終了 |

大腸（C）では `villus` は表示されない。

---

## 5. 絨毛の計測（villus）

`villus` を選ぶと**連続計測モード**に入る。1 回の操作で複数の絨毛を続けて計測できる。  
各絨毛について **villus_height**（高さ）・**villus_width**（幅）・**villus_area**（面積）が CSV に記録され、  
共通の番号（VH-N / VW-N）でペアを識別できる。

計測本数とバックグラウンドカラーピッカーはメインメニューで設定する（セクション 4 参照）。  
指定本数に達したらループが自動終了し、メインメニューに戻る。

### Step 1 — 中心線（高さ）

1. **セグメント化ポリライン**ツールが自動選択される
2. 絨毛の**基部から先端**に向かって中心線をトレースする
3. `T` キーで ROI マネージャーに登録する
4. OK をクリックする → Step 2 へ進む

- ダイアログタイトルに進捗が表示される（例：`Villus-3  Step 1 — Draw centerline  (3/10)`）
- **目標本数に達したら自動終了**（Step 1 ダイアログは表示されない）
- **早期終了**: Step 1 で何も描かずに OK をクリックする → その時点でループ終了

### Step 2 — アウトライン（幅）

1. **ポリゴン**ツールが自動選択される
2. 絨毛の外形に沿ってポリゴンを描く（T キー不要、OK で確定）
3. Python が幅プロファイルを計算し、QC 画像を表示する
4. OK をクリックすると `villus_width` が記録され、**自動的に次の絨毛の Step 1 へ**移る

> 幅をスキップしたい場合は、Step 2 でポリゴンを描かずに OK をクリックする。  
> `villus_height` は既に記録済みなので失われない。次の絨毛の Step 1 に移る。

### バックグラウンドカラーピッカー

セッション開始時に「Pick background colour」にチェックを入れると起動する。

1. 画像内のバックグラウンド領域（組織のない白い部分）をクリックする
2. OK をクリックする

クリックした点の RGB 値が `musc_params.json` の `bg_rgb` フィールドに保存される。  
保存後、`villus_width.py` が自動的にその色を使って絨毛ポリゴン内の組織と背景を分離するため、  
幅計算の精度が向上する。

> 既に `bg_rgb` が設定済みの場合はスキップ（チェック OFF）してよい。  
> 別の染色バッチや動物に変わったときに再設定することを推奨。

### QC 画像の見方

```
villus QC overlay
─────────────────────────────
  シアン   : ポリゴン輪郭（ROI 境界）
  黄色     : 中心線（スムージング後）
  白破線   : 20% / 80% 位置（mean_mid の計算範囲の境界）
  緑       : 25% / 50% / 75% 位置の幅コード
  赤       : 無効点（マスク外にはみ出たレイ）
```

白い破線 2 本が QC 画像に表示される。これらは中心線の **20% および 80% の弧長位置** を法線方向に横断する線で、  
`mean_mid`（代表幅）の算出範囲（先端・基部 20% を除いた中央域）を示す。

### 出力値の説明

| 値 | 意味 |
|---|---|
| `villus_height` | 中心線の弧長（µm）。Notes に `label=VH-N` が付く |
| `villus_width`（= `mean_mid`） | 中心線の 20〜80% 区間における幅の平均（µm）。先端・基部のノイズを除外した代表値 |
| `villus_area` | ポリゴン頂点の靴紐（Shoelace）公式による面積（µm²）。Notes に `pair=VH-N` が付く |
| `w25` / `w50` / `w75` | 中心線の 25% / 50% / 75% の位置での幅（パーセンタイルではなく位置） |
| `wmax` | プロファイル全体の最大幅 |
| `n` | 有効な幅サンプル点の数 |
| `pair=VH-N` | 対応する `villus_height` エントリの番号 |

> **注意**: `w25 < w50 > w75` となる場合は中腹が最も広い（先細り・根元細り型）。  
> 解剖学的に正常な絨毛では typical なパターン。

---

## 6. 陰窩の計測（crypt_depth）

`crypt_depth` を選ぶと**連続計測モード**に入る。計測本数はメインメニューの「Measurements per session」で設定する（デフォルト 10）。

1. **直線**ツールが自動選択される（折れ線ではなく直線）
2. 陰窩の開口から底部まで直線を引く
3. `T` キーで ROI マネージャーに登録する
4. OK をクリックすると CSV に記録され、次の陰窩のダイアログに自動的に移る
5. **目標本数に達したら自動終了**、またはダイアログで何も描かずに OK をクリックすることで早期終了できる

ダイアログのタイトルに進捗が表示される（例：`CD-3  (3/10)`）。

> **T キー押し忘れ救済機能**：T を押さずに OK をクリックした場合、描線が検出されると  
> 「Forgot to press T?」ダイアログが表示される。その場で T を押して登録できる。  
> それでも登録しない場合はキャンセルでループを終了できる。

---

## 7. 筋層計測の詳細（muscularis）

筋層（muscularis propria）の平均厚みを、中心軸 EDT アルゴリズムで自動計測する。

### ワークフロー

```
Step 1: （キャリブレーション有効の場合）tune_muscularis GUI で参照色を設定
   ↓
Step 2: 500 µm グリッドを表示（位置感覚の参考用）
   ↓
Step 3: ポリゴン ROI を描画 → Python で厚み計測
   ↓
  （キャリブレーション無効のとき）4パネル検証画像を表示
   ↓
Step 4: 結果を CSV に追記 + QC 画像を保存
```

### ポリゴン ROI の描画

- シアンの格子線（500 µm 間隔）がガイドとして表示される
- ポリゴンツールで筋層を囲む
  - クリックで頂点を配置
  - ダブルクリックまたは右クリックでポリゴンを閉じる
  - 筋層の上下に少し余白を含めて囲む

---

## 8. キャリブレーション GUI の使い方

`tune_muscularis.py` が起動し、3 パネルの GUI が開く。

```
┌──────────────────┬──────────────────┬──────────────────┐
│   Panel 1        │   Panel 2        │   Panel 3        │
│  オリジナル画像  │  検出マスク      │  E-H 特徴マップ  │
│  （ここでクリック）│  （緑オーバーレイ）│  + 境界線（黄）  │
└──────────────────┴──────────────────┴──────────────────┘
```

### 推奨手順

**1. 参照色の設定（2〜3クリック）**

| ボタン | 操作 |
|---|---|
| `1. Pick Pink (muscle)` | ボタン押下 → Panel 1 のピンク（筋層）部分をクリック |
| `2. Pick Purple (crypt)` | ボタン押下 → Panel 1 の紫（陰窩）部分をクリック |
| `3. Pick Background` | ボタン押下 → Panel 1 の白い背景部分をクリック（任意） |

> **ヒント**：ピンクの筋層帯が小さい場合は、matplotlib ツールバーの虫眼鏡ボタンで  
> Panel 1 をズームしてから参照色を選ぶと精度が上がる。

**2. ファクタースライダーで微調整**

スライダーを動かすと Panel 2 のマスクがリアルタイムで更新される。

- スライダー右（値大）→ より厳密（ピンクに近い画素のみ選択）
- スライダー左（値小）→ より包括的（薄いピンクも含める）

**3. 保存**

`Save & Close` ボタンで設定を保存して閉じる。  
設定は `[Fiji.app]\GutMorphometry_config\musc_params.json` に保存され、  
次回の計測から自動的に読み込まれる。

**4. リセット**

別の動物や染色条件が変わったときは `Clear refs` で参照色をリセットする。

---

## 9. パラメータ一覧

### musc_params.json の内容

```json
{
  "e_factor": 1.0,
  "pink_rgb": [210, 160, 175],
  "purple_rgb": [140, 100, 160],
  "bg_rgb": [250, 248, 248]
}
```

| パラメータ | 型 | 意味 |
|---|---|---|
| `e_factor` | float | スライダー値。デフォルト 1.0（中立）。正方向で厳密化、負方向で緩和 |
| `pink_rgb` | [R, G, B] | ユーザーが選んだ筋層（ピンク）の参照画素の RGB 値 |
| `purple_rgb` | [R, G, B] | ユーザーが選んだ陰窩（紫）の参照画素の RGB 値 |
| `bg_rgb` | [R, G, B] | ユーザーが選んだ背景の参照画素の RGB 値。筋層計測と絨毛幅計測の両方で使用される |

> `bg_rgb` は筋層マスクの精度向上に使われるだけでなく、`villus_width.py` でも  
> 絨毛ポリゴン内部の組織領域（スライドガラスとの境界）を分離するために使われる。

### e_factor の効果

**参照色あり（OD 最近傍分類モード）：**
```
margin = (e_factor - 1.0) × OD_scale

e_factor = 1.0  → margin = 0（純粋な最近傍分類）
e_factor > 1.0  → ピンクが紫・背景よりも「margin」以上近くないと筋層と判定しない（厳格）
e_factor < 1.0  → margin がマイナスになり、紫・背景寄りでも筋層と判定（緩和）
```

**参照色なし（Otsu フォールバックモード）：**
```
threshold = Otsu 閾値 + (e_factor - 1.0) × IQR_scale

e_factor > 1.0  → 閾値が上がり、よりピンクに偏った画素のみ検出
e_factor < 1.0  → 閾値が下がり、より多くの画素を検出
```

### マスク処理パラメータ（固定値）

| パラメータ | 値 | 意味 |
|---|---|---|
| `r_bridge` | `8 µm / pixel_width_um`（最小3、最大20 px） | 筋層内の核を橋渡しするクロージング半径 |
| 最小面積 | `400 µm²` | これより小さい領域を除去するサイズフィルタ |
| 平滑化 sigma | 2.0 px | 境界スムージングの Gaussian 半径（固定） |

---

## 10. アルゴリズム解説

### H&E 色分離（Ruifrok & Johnston 2001）

H&E 染色の RGB 画像を光学密度（OD）空間に変換し、Haematoxylin（H）と Eosin（E）の濃度マップに分解する。

```
OD = -log10(RGB / 255)
[H, E, Residual] = OD × (染色ベクトル行列)^-1
```

### 組織マスク

OD の合計が 0.15 未満の画素を背景として除外する（ポリゴン外の白も除外される）。

### 筋層マスクの生成

**参照色あり（推奨）— OD 最近傍分類：**

```
d_pink   = ||OD_pixel - OD_pink||
d_purple = ||OD_pixel - OD_purple||
d_bg     = ||OD_pixel - OD_bg||   （bg_rgb が設定されている場合）

筋層 = 組織 AND (d_purple - d_pink > margin)
                AND (d_bg - d_pink > margin)  （bg_rgb がある場合）
```

**参照色なし — Otsu フォールバック：**

```
feature = E - H
threshold = Otsu(feature) + (e_factor - 1.0) × scale
筋層 = 組織 AND (feature > threshold)
```

### 絨毛ポリゴン内の組織分離（villus_width.py）

ポリゴン輪郭をそのまま境界として使うと、スライドガラス（無染色部分）も含まれてしまう。  
`bg_rgb` が設定されている場合、以下の方法で組織領域を抽出する：

```
d_bg = ||OD_pixel - OD_background||   （各画素から背景色までの OD 空間距離）
threshold = Otsu(d_bg) for pixels inside polygon
tissue = polygon AND (d_bg > threshold)
```

背景と組織を OD 空間で分離することで、ポリゴン内の染色されていない領域（スライドガラス）を除外できる。

### 形態処理パイプライン（筋層マスク）

```
初期マスク
  → binary_closing(r_bridge)    # 筋細胞核を橋渡し
  → binary_fill_holes           # 内部の穴を埋める
  → binary_closing(disk(3))     # 小さな欠けを修復
  → remove_small_objects        # (20µm)² 未満の断片を除去
  → 最大コンポーネントのみ保持  # 筋層固有層・人工物を排除
  → binary_fill_holes           # 最終コンポーネントの穴を埋める
  → Gaussian平滑化(sigma=2px)   # 境界のギザギザを除去
  → binary_fill_holes           # 平滑化後の穴を埋める
```

### 中心軸 EDT による筋層厚み計測

従来の全画素 EDT 平均ではなく、**中心軸（medial axis）上の画素のみ**でサンプリングする。

```
skel, dist = medial_axis(mask, return_distance=True)
  skel : ボロノイ境界による幾何学的中心軸
  dist : 各中心軸画素から最近傍境界までの距離 r（内接円半径）

# スパー除去（端点を繰り返し削除）
skel = prune_spurs(skel, n_iter=10)

mean_thickness = mean(2 × r) for all axis pixels
effective_length = area / mean_thickness
```

**直感的な理解：**
```
■■■■■■■■■■■■
■ r=1       ■
■   r=3       ■      r = 境界までの最短距離（内接円半径）
■     r=5       ■    local_thickness = 2r
■   r=3       ■      中央の中心軸画素が最も大きな r を持つ
■ r=1       ■
■■■■■■■■■■■■
```

均一な帯状構造では `mean(2r) = W`（帯幅の正確な推定）になる。  
全画素 EDT 平均では `mean(EDT) = W/2` になるため、係数 2 がずれる問題があった。

**スパー除去の意義：**  
マスク端部でのフリンジ（枝分かれ）を除去し、主軸のみを残すことで、端部バイアスを排除する。

### 絨毛幅の計測手順（villus_width.py）

```
1. ポリゴン → ラスタライズ → 粗マスク
2. （bg_rgb がある場合）OD-distance 組織分離 → 精細マスク
3. 中心線をアーク長 1px 間隔で再サンプリング
4. Gaussian 平滑化（sigma = max(2, total_length × 0.03)）
5. 中央差分で接線・法線ベクトルを計算
6. ±法線方向にレイを投射 → マスク境界までの距離を双方向で計算
7. 先端・基部の各 10% を除外 → 有効なプロファイルを収集
8. mean_mid（20〜80%）、w25/50/75、wmax を計算
```

---

## 11. 出力ファイルの説明

### フォルダ構成

計測結果は、オリジナル画像と同じディレクトリに **画像ファイル名と同名のサブフォルダ** が自動作成され、その中にまとめて保存される：

```
[画像フォルダ]/
├── 622_U_1.tif          ← オリジナル画像（変更なし）
└── 622_U_1/             ← 自動作成されるサブフォルダ（画像名と同じ）
    ├── 622_U_1_results.csv              ← 計測データ（全タイプ共通）
    ├── 622_U_1_coords.jsonl             ← ROI 座標データ（JSON Lines）
    ├── 622_U_1_outer_MU1_muscularis.png ← 筋層 QC 画像
    ├── 622_U_1_outer_V1_villus.png      ← 絨毛幅 QC 画像
    └── 622_U_1_overlay.jpg              ← セッション終了時のオーバーレイ
```

複数の画像を同じフォルダに置く場合も、それぞれの画像に対応するサブフォルダが個別に作成されるため、ファイルが混在しない：

```
[画像フォルダ]/
├── 622_U_1.tif
├── 622_U_1/  ← 画像1の全出力
├── 622_U_2.tif
├── 622_U_2/  ← 画像2の全出力
└── ...
```

---

### results/ フォルダ

#### CSV ファイル（例：`622_U_1_results.csv`）

```csv
AnimalID,Region,Section,LoopPosition,MeasurementType,Value_um,Notes
622,U,1,outer,villus_height,312.50,label=VH-1
622,U,1,outer,villus_width,98.30,w25=73.6;w50=100.6;w75=95.7;wmax=107.7;n=99;pair=VH-1
622,U,1,outer,villus_area,28450.20,pair=VH-1
622,U,1,outer,crypt_depth,180.20,
622,U,1,outer,muscularis,85.30,medial_axis_EDT;eff_len=2180um
622,U,1,outer,VH_CD_ratio,1.73,nVH=10;nCD=10;meanVH=312.5um;meanCD=180.2um
```

| カラム | 内容 |
|---|---|
| `AnimalID` | 動物 ID |
| `Region` | 部位コード（U/M/L/C） |
| `Section` | 切片番号 |
| `LoopPosition` | スイスロール位置（outer/middle/inner） |
| `MeasurementType` | 計測タイプ |
| `Value_um` | 計測値（µm） |
| `Notes` | 備考（自動付記あり） |

> データは常に**追記**される。既存行は上書きされない。

**Notes カラムの自動付記：**

| タイプ | 自動付記の内容 |
|---|---|
| `villus_height` | `label=VH-N`（絨毛番号） |
| `villus_width` | `w25=X;w50=X;w75=X;wmax=X;n=X;pair=VH-N` |
| `villus_area` | `pair=VH-N`（対応する絨毛高さの番号） |
| `muscularis` | `medial_axis_EDT;eff_len=Xum`（有効長） |
| `VH_CD_ratio` | `nVH=N;nCD=N;meanVH=Xum;meanCD=Xum`（Save & Exit 時に自動追記） |

> **VH/CD ratio について**: `--- Save and Exit ---` を選択してセッションを終了すると、  
> ループ位置ごとの絨毛高さ平均 ÷ 陰窩深さ平均の比が自動計算され CSV に追記される。  
> 計測数が 0 の場合は追記されない。

#### JSONL ファイル（例：`622_U_1_coords.jsonl`）

ROI の座標情報を 1 行 1 エントリの JSON Lines 形式で記録する。  
後処理（Python / pandas など）での再解析や空間的な検証に使用できる。

```json
{"animal":"622","region":"U","section":"1","loop":"outer","type":"villus_height","label":"VH-1","length_um":312.5,"notes":"label=VH-1","coords_px":[[120.0,450.0],[125.0,380.0]],"coords_um":[[...],[...]]}
```

---

### qc/ フォルダ

#### 筋層 QC 画像（例：`622_U_1_outer_MU1_muscularis.png`）

2 パネル画像。

```
左パネル: 元画像 + シアンの境界線
右パネル: 厚さマップ（hot カラーマップ）
```

**ファイル名の読み方：**
```
622_U_1_outer_MU1_muscularis.png
 │   │  │  │     │
 │   │  │  │     └── MU1 = 同セッション内の筋層計測 1 回目
 │   │  │  └──────── outer = スイスロールの外層
 │   │  └─────────── 1 = 切片番号
 │   └────────────── U = 上部小腸
 └────────────────── 622 = 動物 ID
```

**厚さマップの色の見方（hot カラーマップ）：**

| 色 | 意味 |
|---|---|
| 黒 | 筋層の端（局所的に薄い） |
| 赤〜橙 | 中程度の厚み |
| 黄〜白 | 最も厚い部分（中央軸付近） |

カラーバーの最大値は全筋層画素の 95 パーセンタイル値に設定される。

#### 絨毛 QC 画像（例：`622_U_1_outer_V1_villus.png`）

```
シアン : ポリゴン輪郭（= 組織マスクの境界）
黄色  : 平滑化された中心線
緑    : 25% / 50% / 75% 位置の幅コード
赤    : 無効点（レイがマスク外にはみ出た位置）
タイトル: 組織分離モード（bg_rgb 使用 or ポリゴン境界のみ）
```

#### オーバーレイ JPEG（例：`622_U_1_overlay.jpg`）

セッション終了時（`--- Save and Exit ---` 選択時）に保存される。  
全 ROI オーバーレイが焼き込まれた圧縮画像（JPEG 65% 品質）。  
計測の全体像の確認・記録用。

---

### GutMorphometry_config/（共有設定フォルダ）

```
[Fiji.app]\GutMorphometry_config\
    musc_params.json   ← 筋層キャリブレーション設定
    pixel_scale.txt    ← スケール保存値（µm/px）
```

同じ PC 上の全セッションで共有される。別の PC にツールをコピーする場合は、このフォルダごとコピーすると設定を引き継げる。

- `musc_params.json`：動物や染色ロットが変わったら `Clear refs` でリセットして再キャリブレーションすることを推奨
- `pixel_scale.txt`：スケールバー設定を保存。次回セッションで「Use saved」として再利用される

---

## 12. トラブルシューティング

### "Python script produced no output"

**原因と対処：**

| 原因 | 確認方法 | 対処 |
|---|---|---|
| Python コマンドが見つからない | `g_python_cmd` の値を確認 | `where python` でフルパスを調べて設定 |
| パッケージ未インストール | `python -c "import skimage"` | `pip install scikit-image` |
| スクリプトファイルが見つからない | `py_script` パスを確認 | AutoRun フォルダに `.py` ファイルを配置 |

### キャリブレーション GUI が表示されない

- `tune_muscularis.py` が `musc_thickness.py` と同じフォルダにあるか確認
- matplotlib のバックエンドが GUI を表示できる環境か確認（headless サーバー環境では動作しない）

### 筋層が検出されない / 別の構造物が検出される

1. `Clear refs` で参照色をリセット
2. ツールバーのズームで筋層（ピンク帯）を拡大表示
3. `Pick Pink` → ピンク部分をクリック、`Pick Purple` → 紫部分をクリック
4. 必要なら `Pick Background` → 背景（白）部分をクリック
5. スライダーで微調整 → `Save & Close`

### 絨毛幅の QC 画像で赤点（invalid）が多い

- Step 2 のポリゴンが絨毛の輪郭に密着しすぎている → 少し余白を持たせて描く
- `bg_rgb` が未設定のため組織分離ができず、ポリゴン境界内の空白部分がマスクに含まれている → キャリブレーション GUI で `Pick Background` を実行する

### スケールが正しく設定されない

- `Set scale from 500 µm scale bar?` チェックボックスを ON にしてセッションを開始
- スケールバーの両端に正確に線を引く（拡大して描くと精度が上がる）
- 設定後、`Image > Properties` でスケールを確認できる

### サブフォルダが作成されない

- 画像フォルダへの書き込み権限を確認する（Program Files 内など権限が制限されている場所に画像を置かない）
- `File.makeDirectory` は Fiji macro の組み込み関数。Fiji のバージョンが古い場合は更新を検討する

### バックグラウンドカラーを選んだが幅計算に反映されない

- `musc_params.json` に `bg_rgb` が正しく書き込まれているか確認する（Fiji のステータスバーに `R=... G=... B=... saved` が表示されれば OK）
- `villus_width.py` を AutoRun フォルダに配置しているか確認する

---

*GutMorphometry — 開発者: ksmhp*
