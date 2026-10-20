# 生成出力の読み方

このページは公開コードを実行したときに生成される主な列・尺度を説明します。生成結果一式を配布しているという意味ではありません。公開物には固定入力・説明書・コードが含まれ、実行結果は指定した出力先に作成されます。点推定だけの実行と全TOST実行は区別してください。

## 1. 都道府県コードと支持集合内クラスタID

`common_input_1573.csv` の `prefecture_code` は全国01〜47の都道府県コードです。

一方、`sensitivity/inputs/supports/*/INPUT_METADATA.json` の `source_references[].prefecture_code` は、各支持集合内で1〜Gに再符号化されたクラスタIDです。生成入力 `model_stdin.csv` の `pref`、固定効果列 `pref_feXX`、生成出力 `LEAVE_ONE_PREFECTURE_OUT.csv` の `removed_prefecture_code` も、この支持集合内IDを使います。全国コードとの結合には [集合別対応表](sensitivity/support_cluster_code_crosswalk.csv) を用いてください。例えばP2_k1のクラスタ12は全国コード13（東京都）です。集合をまたいでクラスタ番号だけを結合しないでください。

対応表は既存研究IDによる固定入力との照合から作成した説明資料です。支持集合、クラスタ分け、元metadata、計算コードは変更していません。

## 2. 係数・予測率差・区間

| ファイル／列 | 意味 |
|---|---|
| OLS係数表の `estimate_fraction` | アウトカムが0〜1尺度の回帰係数。 |
| OLS係数表の `estimate_pp` | 上記係数の100倍。連続共変量は保存尺度1単位当たりで、割合共変量の1pp増加当たりではありません。 |
| `FRACTIONAL_COEFFICIENTS.csv` の `estimate` | fractional-logitのlogit尺度係数。率差やppではありません。 |
| `FRACTIONAL_FIT_DIAGNOSTICS.csv` の `delta_pp` | 同じ支持集合でDを全員1／全員0に置いた予測率の差を等重み平均し、100倍した値（D1−D0）。 |
| `THRESHOLDS.csv` の `delta_pp` | 同等性判定の固定幅j/100 pp。推定された率差ではありません。 |
| `SE_fraction`／`SE_pp` | それぞれ比率／pp尺度の標準誤差。 |
| CR2出力の `CI_lower_*`／`CI_upper_*` | 対応する係数・対比の多重性未調整の個別95％参考区間。単位は接尾辞によります。 |

例えばshare変数の0.01増加に対応するOLSの差は、その変数の `estimate_pp` の0.01倍です。対数人数の係数も人数1人増加当たりではありません。全国加法モデルの納期回数は8回、仮算定なしフラグはあり、全国都道府県固定効果は01が基準です。支持集合内モデルの固定効果は各集合のクラスタ1が基準です。

`cr2/TWO_CONTRAST_RESULTS.csv` の `saved_estimate_pp` と `estimate_fraction` は保存精密値による中心、`recomputed_estimate_fraction` は同一入力から再計算した点差です。保存中心との照合を通して区間を接続する仕様は [PROVENANCE_JA.md](PROVENANCE_JA.md) に記載しています。

## 3. 支持集合のD符号と論文での表示方向

| 比較 | 元のD1 | 元のD0 | 論文の表示方向 | 元のD1−D0からの表示変換 |
|---|---|---|---|---|
| P1 | 仮算定なし・6回 | 仮算定あり・6回 | なし6回−あり6回 | そのまま |
| P2 | 仮算定なし・7月開始毎月8回 | 仮算定なし・7月開始毎月9回 | なし9回−なし8回 | 符号反転 |
| P3 | 仮算定あり・6回 | 仮算定あり・12回 | あり12回−あり6回 | 符号反転 |

本文r27の方法Ⅱ5・表1に対応する記載、補足r26のS-M3・表S5注と照合しています。原稿ファイル自体は本リポジトリに含まれません。元入力や係数を反転して再保存するのではなく、表示の際に上表を適用します。P2・P3の下／上側を表示し直す場合は元の上／下側に対応し、tの符号を反転します。p値・MC区間・尾数・診断は対応する側の値を保ち、1−pへ変換しません。

2026年10月20日には、上記の表示方向と側の対応を本文r30の方法Ⅱ5・表1、補足r29のS-M3・表S5注とも照合しました。上のr27／r26は以前の照合時の版番号です。

`FRACTIONAL_DIAGNOSTIC_ARRAYS.csv` の `p_observed` は観測制度下の予測率、`p_D0`・`p_D1` は同じ集合の全員の制度指標だけを0／1に置いた予測率です。実測の所属群別集計ではありません。`eta_*` はlogitの線形予測子です。

## 4. その他の取り違えやすい列

- `ILLUSTRATIVE_MONEY_SCALES.csv` の `point_yen`、`lower_yen`、`upper_yen` は、率差の点推定と参考区間の両端を同じ `fixed_assessment_yen` で換算した円額です。`*_10k_yen` は万円。個別保険者の制度変更後の収入予測や減収上限ではありません。
- `exploration/institution_size` の `ordinary_assessed_sum`、`special_assessed_sum`、`total_assessed_sum` は円、`insured_*` は年度末第1号被保険者数です。`influence_distribution` の `assessed_amount`、`received_amount` も円で、同ファイルの `amount_unit` を併読します。
- `S3_ALL_DELETION_COEFFICIENTS_PP.csv` の `annual_*`、`provisional_none` 列は、その変数の係数（pp）です。入力の0/1フラグではありません。
- 記述分布の `quantile_0.1_percent` は10パーセンタイル値を％尺度で示します。列名内のqは0〜1の確率です。
- `caliper` は背景変数の支持幅（SD尺度）。同等性マージン（pp）とは別です。`neighbor_count` は条件を満たす相手の数で、分析重みやペア割付ではありません。
- `positive_cell`／`negative_cell` は差の引き算で正側／負側に置く制度セルで、効果の正負を意味しません。

## 5. 空欄・NA・0と判定状態

生成出力の空欄・NAには、空セル群で平均が定義されない場合、係数が識別できない場合、診断式を適用できない場合などがあります。`support`、`empty_cell`、`undefined_reason`、`D_estimable` 等を併読し、0に一括置換しないでください。公開入力の月列は全て保存値0/1で、欠測補完を表す0ではありません。

`THRESHOLDS.csv` の `tost_status` は固定規則によるPASS／NONPASS／HOLDです。HOLDは計算診断を優先した状態で、他方の境界の公式p値が通過しない場合もあります。HOLDは同等性成立を示すラベルではありません。また、tost_statusは効果の信頼区間ではありません。詳細は [sensitivity/README.md](sensitivity/README.md) を参照してください。
