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
	global_option = parser.add_argument_group(title = "global option", description = "")
	global_option.add_argument("-i", dest = "input", metavar = "LOG", required = True, help = "LOG for ABINIT-MP")
	global_option.add_argument("-o", dest = "prefix", help = "prefix for output")
	global_option.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")
	global_option.add_argument("-N", dest = "flag_cancel", action = "store_true", default = False, help = "do NOT cancel the interaction between adjacent fragments (default: False)")

	output_type = parser.add_argument_group(title = "energy type option", description = "energy type for output (default: -t)")
	output_type.add_argument("-a", "--all", dest = "flag_all", action = "store_true", default = False, help = "select all type, same as -tfesxcdq")
	output_type.add_argument("-t", "--total", dest = "flag_total", action = "store_true", default = False, help = "total energy (kcal/mol)")
	output_type.add_argument("-f", "--hartree", dest = "flag_HF", action = "store_true", default = False, help = "Hartree-Fock energy (kcal/mol)")
	output_type.add_argument("-e", "--correlation", dest = "flag_corr", action = "store_true", default = False, help = "electron correlation energy (kcal/mol)")
	output_type.add_argument("-s", "--electrostatic", dest = "flag_ES", action = "store_true", default = False, help = "electrostatic interaction (ES) (kcal/mol)")
	output_type.add_argument("-x", "--exchange", dest = "flag_EX", action = "store_true", default = False, help = "exchange-repulsion energy (EX) (kcal/mol)")
	output_type.add_argument("-c", "--chargetransfer-mix", dest = "flag_CT", action = "store_true", default = False, help = "charge transfer and other interaction energy (CT+mix) (kcal/mol)")
	output_type.add_argument("-d", "--dispersion", dest = "flag_DI", action = "store_true", default = False, help = "dispersion energy (DI) (kcal/mol)")
	output_type.add_argument("-q", "--chargetransfer-amount", dest = "flag_q", action = "store_true", default = False, help = "amount of charge transfer")
	output_type.add_argument("-ac", "--atomic-charge", dest = "flag_ac", action = "store_true", default = False, help = "atomic charge")
	args = parser.parse_args()

	basic.check_exist(args.input, 2)

	flag_total = args.flag_total
	flag_HF = args.flag_HF
	flag_corr = args.flag_corr
	flag_ES = args.flag_ES
	flag_EX = args.flag_EX
	flag_CT = args.flag_CT
	flag_DI = args.flag_DI
	flag_q = args.flag_q

	if args.flag_all:
		flag_total = True
		flag_HF = True
		flag_corr = True
		flag_ES = True
		flag_EX = True
		flag_CT = True
		flag_DI = True
		flag_q = True
	elif flag_total == False:
		if not (flag_HF or flag_corr or flag_ES or flag_EX or flag_CT or flag_DI or flag_q or args.flag_ac):
			# 他のオプションが未指定の場合のみ total オプションを機能させる
			flag_total = True

	flag_pieda = False
	if flag_ES or flag_EX or flag_CT or flag_DI or flag_q:
		flag_pieda = True

	fragment_atoms = []
	fragment_connects = []
	energy_HFs = []
	energy_MP2s = []
	energy_tots = []
	energy_ESs = []
	energy_EXs = []
	energy_CTs = []
	energy_DIs = []
	energy_qs = []
	atomic_charges = [["No.", "Atom", "Atomic pop.", "Net charge"]]

	with open(args.input, "r") as obj_input:
		re_wsp = re.compile(r"[\s\t]+")
		re_fragmentation = re.compile(r"Frag\.   Elec\.   ATOM")
		re_bonded_atom = re.compile(r"Bonded Atom")
		re_IFIE = re.compile(r"## ((HF)|(MP2))-IFIE")
		re_separator = re.compile(r"-{5,}")
		re_PIEDA = re.compile(r"## PIEDA")
		re_charge = re.compile(r"No\. Atom   Atomic pop\.  Net charge")
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
							energy_CTs.append([])
							energy_DIs.append([])
							energy_qs.append([])
						else:
							# 最初の行はラベル
							labels = list(range(1, len(fragment_atoms) + 1))
							labels.insert(0, "")
							energy_ESs.insert(0, labels)
							energy_EXs.insert(0, labels)
							energy_CTs.insert(0, labels)
							energy_DIs.insert(0, labels)
							energy_qs.insert(0, labels)
							continue
						for j in range(len(fragment_atoms) + 1):
							# 通常行の中身
							if i == j:
								energy_ESs[i].append(0.0)
								energy_EXs[i].append(0.0)
								energy_CTs[i].append(0.0)
								energy_DIs[i].append(0.0)
								energy_qs[i].append(0.0)
							elif j != 0:
								energy_ESs[i].append("-")
								energy_EXs[i].append("-")
								energy_CTs[i].append("-")
								energy_DIs[i].append("-")
								energy_qs[i].append("-")
							else:
								# 最初の列はラベル
								energy_ESs[i].append(i)
								energy_EXs[i].append(i)
								energy_CTs[i].append(i)
								energy_DIs[i].append(i)
								energy_qs[i].append(i)

			elif flag_pieda and re_PIEDA.search(line):
				flag_read = 5

			elif re_charge.search(line):
				flag_read = 7

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
				energy_CT = float(line[48:63].strip())
				energy_DI = float(line[63:78].strip())
				energy_q = float(line[78:93].strip())

				if args.flag_cancel == False:
					if [frag_i, frag_j] in fragment_connects or [frag_j, frag_i] in fragment_connects:
						# フラグメントが共有結合している場合
						energy_ES = 0.0
						energy_EX = 0.0
						energy_CT = 0.0
						energy_DI = 0.0
						energy_q = 0.0

				energy_ESs[frag_i][frag_j] = energy_ES
				energy_ESs[frag_j][frag_i] = energy_ES
				energy_EXs[frag_i][frag_j] = energy_EX
				energy_EXs[frag_j][frag_i] = energy_EX
				energy_CTs[frag_i][frag_j] = energy_CT
				energy_CTs[frag_j][frag_i] = energy_CT
				energy_DIs[frag_i][frag_j] = energy_DI
				energy_DIs[frag_j][frag_i] = energy_DI
				energy_qs[frag_i][frag_j] = energy_q
				energy_qs[frag_j][frag_i] = energy_q

			elif flag_read == 7:
				line = line.strip()
				datas = re_wsp.split(line)
				atomic_charges.append(datas)

	# 出力ファイル
	prefix = ""
	if args.prefix != None:
		prefix = args.prefix
	else:
		prefix = re.sub(r"\..{3,4}$", "", args.input)

	import csv

	if flag_total:
		output = prefix + "_tot.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_tots)
			sys.stderr.write("create: %s (total)\n" % output)

	if flag_HF:
		output = prefix + "_HF.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_HFs)
			sys.stderr.write("create: %s (hartree)\n" % output)

	if flag_corr:
		output = prefix + "_MP2.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_MP2s)
			sys.stderr.write("create: %s (correlation)\n" % output)

	if flag_ES:
		output = prefix + "_ES.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_ESs)
			sys.stderr.write("create: %s (electrostatic)\n" % output)

	if flag_EX:
		output = prefix + "_EX.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_EXs)
			sys.stderr.write("create: %s (exchange)\n" % output)

	if flag_CT:
		output = prefix + "_CT.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_CTs)
			sys.stderr.write("create: %s (chargetransfer-mix)\n" % output)

	if flag_DI:
		output = prefix + "_DI.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_DIs)
			sys.stderr.write("create: %s (dispersion)\n" % output)

	if flag_q:
		output = prefix + "_transq.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy_qs)
			sys.stderr.write("create: %s (chargetransfer-amount)\n" % output)

	if args.flag_ac:
		output = prefix + "_ac.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(atomic_charges)
			sys.stderr.write("create: %s (atomic-charge)\n" % output)
