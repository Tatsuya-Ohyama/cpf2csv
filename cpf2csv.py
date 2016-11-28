#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpf2csv.py
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from py_module_basic import basic

# =============== functions =============== #


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "cpf2csv - convert log for ABINIT-MP to CSV", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest = "input", metavar = "LOG", required = True, help = "LOG for ABINIT-MP")
	parser.add_argument("-o", dest = "prefix", help = "prefix for output")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	parser.add_argument("-N", dest = "cancel", action = "store_true", default = False, help = "do NOT cancel the interaction between adjacent fragments")
	args = parser.parse_args()

	basic.check_exist(args.input, 2)

	fragment_atoms = []
	fragment_connects = []
	energy_HFs = []
	energy_MP2s = []
	energy_tots = []

	with open(args.input, "r") as obj_input:
		re_wsp = re.compile(r"[\s\t]+")
		re_fragmentation = re.compile(r"Frag\.   Elec\.   ATOM")
		re_bonded_atom = re.compile(r"Bonded Atom")
		re_IFIE1 = re.compile(r"## ((HF)|(MP2))-IFIE")
		re_IFIE2 = re.compile(r"-{5,}")
		re_empty = re.compile(r"^[\s\t]*\n$")

		flag_read = 0
		for line in obj_input:
			if re_fragmentation.search(line):
				flag_read = 1
			elif re_bonded_atom.search(line):
				flag_read = 2

			elif re_IFIE1.search(line):
				flag_read = 3

				# IFIE 格納用配列を初期化
				for i in range(len(fragment_atoms) + 1):
					if i != 0:
						#  通常の行
						energy_HFs.append([])
						energy_MP2s.append([])
						energy_tots.append([])
					else:
						# 最初の行はラベル
						labels = list(range(1, len(fragment_atoms) + 1))
						labels.insert(0, "")
						energy_HFs.insert(0, labels)
						energy_MP2s.insert(0, labels)
						energy_tots.insert(0, labels)
						continue
					for j in range(len(fragment_atoms) + 1):
						# 通常行の中身
						if i == j:
							energy_HFs[i].append(0.0)
							energy_MP2s[i].append(0.0)
							energy_tots[i].append(0.0)
						elif j != 0:
							energy_HFs[i].append("-")
							energy_MP2s[i].append("-")
							energy_tots[i].append("-")
						else:
							# 最初の列はラベル
							energy_HFs[i].append(i)
							energy_MP2s[i].append(i)
							energy_tots[i].append(i)

			elif flag_read == 3 and re_IFIE2.search(line):
				flag_read = 4

			elif re_empty.search(line) and flag_read != 3:
				flag_read = 0

			elif flag_read == 1:
				# フラグメント構成原子取得
				fragment_num = line[5:13].strip()
				atoms = list(map(lambda x : int(x), re_wsp.split(line[23:].strip())))
				if len(fragment_num) == 0:
					# 継続の場合
					fragment_atoms[len(fragment_atoms) - 1].extend(atoms)

				else:
					# 新たなフラグメントの場合
					fragment_atoms.append(atoms)

			elif flag_read == 2:
				# フラグメント接続情報取得
				datas = list(map(lambda x : int(x), re_wsp.split(line.strip())))
				tmp_connects = []
				for atom in datas[1:3]:
					index = 0
					for info in fragment_atoms:
						index += 1
						if atom in info:
							tmp_connects.append(index)
							break
				fragment_connects.append(tmp_connects)

			elif flag_read == 4:
				# IFIE 取得
				frag_i = int(line[8:13].strip())
				frag_j = int(line[13:21].strip())
				energy_HF = float(line[39:50].strip()) * 627.5095
				energy_MP2 = float(line[50:61].strip()) * 627.5095

				if args.cancel == False:
					if [frag_i, frag_j] in fragment_connects or [frag_j, frag_i] in fragment_connects:
						# フラグメントが共有結合している場合
						energy_HF = 0
						energy_MP2 = 0

				energy_HFs[frag_i][frag_j] = energy_HF
				energy_HFs[frag_j][frag_i] = energy_HF
				energy_MP2s[frag_i][frag_j] = energy_MP2
				energy_MP2s[frag_j][frag_i] = energy_MP2
				energy_tots[frag_i][frag_j] = energy_HF + energy_MP2
				energy_tots[frag_j][frag_i] = energy_HF + energy_MP2

	# 出力ファイル
	prefix = ""
	if args.prefix != None:
		prefix = args.prefix
	else:
		prefix = re.sub(r"\..+$", "", args.input)

	import csv

	output = prefix + "_HF.csv"
	if args.flag_overwrite == False:
		basic.check_overwrite(output)
	with open(output, "w") as obj_output:
		csv_writer = csv.writer(obj_output, lineterminator = "\n")
		csv_writer.writerows(energy_HFs)

	output = prefix + "_MP2.csv"
	if args.flag_overwrite == False:
		basic.check_overwrite(output)
	with open(output, "w") as obj_output:
		csv_writer = csv.writer(obj_output, lineterminator = "\n")
		csv_writer.writerows(energy_MP2s)

	output = prefix + "_tot.csv"
	if args.flag_overwrite == False:
		basic.check_overwrite(output)
	with open(output, "w") as obj_output:
		csv_writer = csv.writer(obj_output, lineterminator = "\n")
		csv_writer.writerows(energy_tots)
