# FY2024 介護保険料普通徴収：データと計算コード

2024年度の全国1,573保険者を対象とする固定研究入力、月別納期の閲覧用Excel、変数説明、補足根拠台帳、既存解析の実行用コードです。個人の被保険者データは含みません。保険者名を含む行政集計であり、匿名データではありません。

このリポジトリは著者によるデータ・コードの公開資料です。『社会保障研究』の採択や、同誌の正式なweb appendixとしての認定を意味しません。投稿本文・補足原稿・編集履歴・私的な連絡記録は含めません。

## 最初に見るもの

- data/FY2024_1573保険者_月別納期と研究分類_提出候補.xlsx：保険者ごとの月別納期と研究分類を閲覧するための表
- data/FY2024_月別納期_提出候補.csv：Excelと対応する月別納期。月別収納額ではありません
- data/common_input_1573.csv：数値計算に使う固定15列の入力。十進表記を維持しています
- data/analysis_data_dictionary.csv、data/README_JA.md：変数・単位・分類と出典の説明
- [補足根拠台帳・公開CSV4点](data/evidence/20261020/)：納期確認方法と分類判断の根拠を追跡するための統合更新版
- PROVENANCE_JA.md：再現できる範囲と限界

ExcelでCSVを上書きすると桁・型が変わることがあります。計算には元のCSVを直接使ってください。研究IDは公的保険者番号ではありません。

## 実行

確認環境はLinux（Debian 13）、Python 3.12.14です。Pythonはrequirements.txtの固定版を使います。仮想環境で `python -m pip install -r requirements.txt` を実行してください。Rによるブートストラップまで行う場合は sensitivity/README.md の固定環境も必要です。

Python 3.12.14 は参照・確認環境です。CR2 と固定仕様再構成は Python 3.12 系の別の修正版数でも警告を記録して続行しますが、NumPy 2.3.5・SciPy 1.17.0 の固定版と既存の数値検査は維持します。両モジュールの RUNTIME_ENVIRONMENT.json に実際の版数を記録します。別環境での互換性や末尾桁一致を一律に保証するものではありません。

Python側の全計算：

    python run_all.py --out /任意の新規出力先 --scope deterministic

Rの感度分析を含む全計算：

    python run_all.py --out /任意の新規出力先 --scope all --r-command /固定R実行ファイル

出力先は新しいディレクトリを指定します。既存の結果の読み替え・コピーで成功を表示する処理はありません。Rの全計算には時間がかかります。`deterministic` はPython部分のみで、TOSTを実行したとは扱いません。R実行にはR本体4.6.1、fwildclusterboot 0.14.3、dqrng 0.4.1等の固定版を使います。

## 含む計算

1. Decimal50の全国・制度組合せ別記述集計
2. 全国等重み加法モデル、都道府県クラスタCR2、対比固有Satterthwaite自由度、既存2対比と固定調定額の換算
3. 4重みの点推定・規模分布・既存の影響診断
4. 固定16セルモデルと全国8係数の情報集中・除外診断（記載仕様から再構成したコード）
5. 固定された背景条件による支持集合の選定
6. 固定支持集合の感度分析・全1,600境界のTOST計算

モデル・分類・支持幅・比較方向・seed・反復数を結果に合わせて変更しないでください。個別95％参考区間は探索的・多重性未調整です。関連差と金額換算は制度変更の因果効果や個別自治体の増減収予測ではありません。

## 提出時の扱い

『社会保障研究』の執筆要領3(6)ではweb appendixも査読対象とされています。既公開で分量が多い資料についてはURL記載による代替の定めがありますが、本文での明記・電子媒体での同時提出等の条件があります。正式な掲載先・扱いは編集委員会の判断に従います。
https://www.ipss.go.jp/site-ad/index_japanese/shippitsuyouryou.pdf

公開しただけで投稿条件への適合が認定されるわけではありません。研究データ・独自コードへの新しい包括的な利用許諾は設定していません。第三者コード・依存ライブラリの権利とライセンスは THIRD_PARTY_NOTICES.md を参照してください。

## Column names and output interpretation

See [データの列名と読み方](data/HEADER_GUIDE_JA.md) for the monthly workbook/CSV correspondence, and [生成出力の読み方](OUTPUT_DICTIONARY_JA.md) for coefficient scales, contrast directions and support-local cluster IDs. The workbook also contains a Japanese reading guide. These guides describe generated files; they do not imply that generated result CSVs are included in this repository.

## 2026年10月20日 根拠資料の統合更新

本文 r34・補足資料 r37 に対応する補足根拠台帳の統合更新版v4（Excel、公開CSV4点、説明）を `data/evidence/20261020/` に収録しています。投稿原稿そのものは公開していません。相馬市・小川村の公開条例写しの版識別と引用照合範囲、姫島村を含む電話4件の確認範囲を著者確認に基づき訂正しています。

電話4件（空知中部広域連合・吉野町・新庄村・姫島村）すべてで、納期・年度接続に加えて初回額が仮額かも確認したことを反映しました。初回分析月6月の電話確認は新庄村1件です。旧メール10件はこの電話による性質確認件数に加算しません。個別の質問・回答の逐語記録は未収録で、具体的な回答文言は補作していません。姫島村のFY2024公式原資料そのものの未接続は維持し、電話による性質確認とは区別します。

解析コードと固定分析入力の参照版は [9357842efdb388d267c75f1f35b88fbb1f3fd023](https://github.com/Knowledge373/ltci-fy2024-data-code/tree/9357842efdb388d267c75f1f35b88fbb1f3fd023) です。今回の変更は根拠資料の確認範囲の訂正と説明・ハッシュの更新に限り、同参照版の解析コード、固定入力、研究採用分類、初回分析月、採用分析値は変更していません。補足資料の既存計算コードへの固定リンクは、この解析参照版を示します。根拠追加版を特定する場合は、このリポジトリの当該更新コミットを別に指定してください。

台帳は保存根拠の追跡整理であり、全1,573件の原資料を新しく読み直した記録ではありません。条例・公式公表・電話・メールの方法別件数には重複があり、仮算定性質の直接確認数や誤分類数と同一ではありません。個人の被保険者情報、電話・メールの原文、第三者原本は収録していません。詳細と確認限界は台帳の「読み方」および [台帳の説明](data/evidence/20261020/README_PUBLIC.txt) を参照してください。

変更履歴は [CHANGELOG.md](CHANGELOG.md)、配布ファイルのSHA-256は `MANIFEST.sha256` に記載しています。
