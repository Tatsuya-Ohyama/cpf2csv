#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpf2csv.py
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from py_module_basic import basic
from mods import EnergyData

# =============== functions =============== #
# 二次リストから特定の値を検索し、index を返す
def search_list(query, array, offset = 0):
	index = offset
	for item in array:
		if query in item:
			return index
		index += 1


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
	output_type.add_argument("-q", "--chargetransfer-amount", dest = "flag_q", action = "store_true", default = False, help = "amount of charge transfer (e; I(row) -> J(col))")
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

	# データ読み込み＆解析
	energy = EnergyData.EnergyData(args.input)


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
			csv_writer.writerows(energy.output_energy("Total"))
			sys.stderr.write("create: %s (total)\n" % output)

	if flag_HF:
		output = prefix + "_HF.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_energy("HF"))
			sys.stderr.write("create: %s (hartree)\n" % output)

	if flag_corr:
		output = prefix + "_MP2.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_energy("Total"))
			sys.stderr.write("create: %s (correlation)\n" % output)

	if flag_ES:
		output = prefix + "_ES.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_energy("ES"))
			sys.stderr.write("create: %s (electrostatic)\n" % output)

	if flag_EX:
		output = prefix + "_EX.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_energy("EX"))
			sys.stderr.write("create: %s (exchange)\n" % output)

	if flag_CT:
		output = prefix + "_CT.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_energy("CT"))
			sys.stderr.write("create: %s (chargetransfer-mix)\n" % output)

	if flag_DI:
		output = prefix + "_DI.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_energy("DI"))
			sys.stderr.write("create: %s (dispersion)\n" % output)

	if flag_q:
		output = prefix + "_transq.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_energy("Q"))
			sys.stderr.write("create: %s (chargetransfer-amount)\n" % output)

	if args.flag_ac:
		output = prefix + "_ac.csv"
		if args.flag_overwrite == False:
			basic.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(atomic_charges)
			sys.stderr.write("create: %s (atomic-charge)\n" % output)
