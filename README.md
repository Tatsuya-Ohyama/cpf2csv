# cpf2csv.py

## 概要
ABINIT-MP の .log および .out ファイルから、IFIE の相互作用を出力するプログラム


## 使用方法
```sh
$ cpf2csv.py [-h] -i LOG [-o PREFIX] [-O] [-a] [-t] [-f] [-e] [-s] [-x] [-c] [-d] [-q] [-pc] [--include Frag_No. [Frag_No. ...] | --exclude Frag_No. [Frag_No. ...]]
```

* `-h`, `--help`
	: ヘルプメッセージを表示して終了する。
* `-i LOG`
	: ABINIT-MP の .log および .out ファイル
* `-o PREFIX`
	: 出力ファイルの接頭辞
* `-O`
	: 上書きするプロンプトを表示せずに上書きする (Default: False)。
* `-a, --all`
	: すべての相互作用エネルギーを出力する (`-tfesxcdq` と同じ)。
* `-t, --total`
	: 全エネルギー (kcal/mol) を出力する。
* `-f, --hartree`
	: Hartree-Fock エネルギー (kcal/mol) を出力する。
* `-e, --correlation`
	: 電子相関エネルギー (kcal/mol) を出力する。
* `-s, --electrostatic`
	: 静電エネルギー (ES) (kcal/mol) を出力する。
* `-x, --exchange`
	: 交換反発エネルギー (EX) (kcal/mol) を出力する。
* `-c, --chargetransfer-mix`
	: 電荷移動エネルギー (CT+mix) (kcal/mol) を出力する。
* `-d`, `--dispersion`
	: 分散力エネルギー (DI) (kcal/mol) を出力する。
* `-q`, `--chargetransfer-amount`
	: 電荷移動量 (e; I(row) -> J(col)) を出力する。
* `-pc`, `--partial-charge`
	: 部分電荷を出力する。
* `--include Frag_No. [Frag_No. ...]`
	: 含めるフラグメントを指定する。
* `--exclude Frag_No. [Frag_No. ...]`
	: 含めないフラグメントを指定する。


## 更新履歴
### Ver. 10.3 (2021-06-02)
* PEP8 に合わせてスタイルを変更した。
* 公開した。
