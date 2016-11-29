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
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")
	parser.add_argument("-p", dest = "flag_pieda", action = "store_true", default = False, help = "output with PIEDA results (Default: False)")
	parser.add_argument("-N", dest = "flag_cancel", action = "store_true", default = False, help = "do NOT cancel the interaction between adjacent fragments (default: False)")
	args = parser.parse_args()

	basic.check_exist(args.input, 2)

	fragment_atoms = []
	fragment_connects = []
	energy_HFs = []
	energy_MP2s = []
	energy_tots = []
	energy_ESs = []
	energy_EXs = []
	energy_mixs = []
	energy_DIs = []
	energy_qs = []

	with open(args.input, "r") as obj_input:
		re_wsp = re.compile(r"[\s\t]+")
		re_fragmentation = re.compile(r"Frag\.   Elec\.   ATOM")
		re_bonded_atom = re.compile(r"Bonded Atom")
		re_IFIE = re.compile(r"## ((HF)|(MP2))-IFIE")
		re_separator = re.compile(r"-{5,}")
		re_PIEDA = re.compile(r"## PIEDA")
		re_empty = re.compile(r"^[\s\t]*\n$")

		flag_read = 0
		for line in obj_input:
			if re_fragmentation.search(line):
				flag_read = 1
			elif re_bonded_atom.search(line):
				flag_read = 2

			elif re_IFIE.search(line):
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

			elif flag_read == 3:
				if re_separator.search(line):
					flag_read = 4

			elif flag_read == 5:
				if re_separator.search(line):
					flag_read = 6
					# `PIEDA 格納用配列を初期化
					for i in range(len(fragment_atoms) + 1):
						if i != 0:
							#  通常の行
							energy_ESs.append([])
							energy_EXs.append([])
							energy_mixs.append([])
							energy_DIs.append([])
							energy_qs.append([])
						else:
							# 最初の行はラベル
							labels = list(range(1, len(fragment_atoms) + 1))
							labels.insert(0, "")
							energy_ESs.insert(0, labels)
							energy_EXs.insert(0, labels)
							energy_mixs.insert(0, labels)
							energy_DIs.insert(0, labels)
							energy_qs.insert(0, labels)
							continue
						for j in range(len(fragment_atoms) + 1):
							# 通常行の中身
							if i == j:
								energy_ESs[i].append(0.0)
								energy_EXs[i].append(0.0)
								energy_mixs[i].append(0.0)
								energy_DIs[i].append(0.0)
								energy_qs[i].append(0.0)
							elif j != 0:
								energy_ESs[i].append("-")
								energy_EXs[i].append("-")
								energy_mixs[i].append("-")
								energy_DIs[i].append("-")
								energy_qs[i].append("-")
							else:
								# 最初の列はラベル
								energy_ESs[i].append(i)
								energy_EXs[i].append(i)
								energy_mixs[i].append(i)
								energy_DIs[i].append(i)
								energy_qs[i].append(i)

			elif args.flag_pieda == True and re_PIEDA.search(line):
				flag_read = 5

			elif re_empty.search(line):
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
				frag_j = int(line[13:18].strip())
				energy_HF = float(line[39:50].strip()) * 627.5095
				energy_MP2 = float(line[50:61].strip()) * 627.5095

				if args.flag_cancel == False:
					if [frag_i, frag_j] in fragment_connects or [frag_j, frag_i] in fragment_connects:
						# フラグメントが共有結合している場合
						energy_HF = 0.0
						energy_MP2 = 0.0

				energy_HFs[frag_i][frag_j] = energy_HF
				energy_HFs[frag_j][frag_i] = energy_HF
				energy_MP2s[frag_i][frag_j] = energy_MP2
				energy_MP2s[frag_j][frag_i] = energy_MP2
				energy_tots[frag_i][frag_j] = energy_HF + energy_MP2
				energy_tots[frag_j][frag_i] = energy_HF + energy_MP2

			elif flag_read == 6:
				# PIEDA 取得
				frag_i = int(line[8:13].strip())
				frag_j = int(line[13:18].strip())
				energy_ES = float(line[18:33].strip())
				energy_EX = float(line[33:48].strip())
				energy_mix = float(line[48:63].strip())
				energy_DI = float(line[63:78].strip())
				energy_q = float(line[78:93].strip())

				if args.flag_cancel == False:
					if [frag_i, frag_j] in fragment_connects or [frag_j, frag_i] in fragment_connects:
						# フラグメントが共有結合している場合
						energy_ES = 0.0
						energy_EX = 0.0
						energy_mix = 0.0
						energy_DI = 0.0
						energy_q = 0.0

				energy_ESs[frag_i][frag_j] = energy_ES
				energy_ESs[frag_j][frag_i] = energy_ES
				energy_EXs[frag_i][frag_j] = energy_EX
				energy_EXs[frag_j][frag_i] = energy_EX
				energy_mixs[frag_i][frag_j] = energy_mix
				energy_mixs[frag_j][frag_i] = energy_mix
				energy_DIs[frag_i][frag_j] = energy_DI
				energy_DIs[frag_j][frag_i] = energy_DI
				energy_qs[frag_i][frag_j] = energy_q
				energy_qs[frag_j][frag_i] = energy_q

	# 出力ファイル
	prefix = ""
	if args.prefix != None:
		prefix = args.prefix
	else:
		prefix = re.sub(r"\..{3,4}$", "", args.input)

	import csv

	output = prefix + "_HF.csv"
	if args.flag_overwrite == False:
		basic.check_overwrite(output)
	with open(output, "w") as obj_output:
		csv_writer = csv.writer(obj_output, lineterminator = "\n")
		csv_writer.writerows(energy_HFs)
		sys.stderr.write("%s was created\n" % output)

	output = prefix + "_MP2.csv"
	if args.flag_overwrite == False:
		basic.check_overwrite(output)
	with open(output, "w") as obj_output:
		csv_writer = csv.writer(obj_output, lineterminator = "\n")
		csv_writer.writerows(energy_MP2s)
		sys.stderr.write("%s was created\n" % output)

	output = prefix + "_tot.csv"
	if args.flag_overwrite == False:
		basic.check_overwrite(output)
	with open(output, "w") as obj_output:
		csv_writer = csv.writer(obj_output, lineterminator = "\n")
		csv_writer.writerows(energy_tots)
		sys.stderr.write("%s was created\n" % output)

	if args.flag_pieda == True:
		# PIEDA の結果出力
		output = prefix + "_ES.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_ESs)
			sys.stderr.write("%s was created\n" % output)

		output = prefix + "_EX.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_EXs)
			sys.stderr.write("%s was created\n" % output)

		output = prefix + "_MIX.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_mixs)
			sys.stderr.write("%s was created\n" % output)

		output = prefix + "_DI.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_DIs)
			sys.stderr.write("%s was created\n" % output)

		output = prefix + "_qtrans.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_qs)
			sys.stderr.write("%s was created\n" % output)
