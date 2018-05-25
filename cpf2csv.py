#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpf2csv.py
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import basic_func
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
	output_type.add_argument("-pc", "--partial-charge", dest = "flag_pc", action = "store_true", default = False, help = "partial charge")

	output_range = parser.add_mutually_exclusive_group()
	output_range.add_argument("--include", dest = "include", metavar = "Frag_No.", nargs = "+", help = "")
	output_range.add_argument("--exclude", dest = "exclude", metavar = "Frag_No.", nargs = "+", help = "")

	args = parser.parse_args()

	basic_func.check_exist(args.input, 2)

	output_flag = [
		args.flag_total,
		args.flag_HF,
		args.flag_corr,
		args.flag_ES,
		args.flag_EX,
		args.flag_CT,
		args.flag_DI
	]
	output_suffix = [
		"_Total.csv",
		"_HF.csv",
		"_CR.csv",
		"_ES.csv",
		"_EX.csv",
		"_CT.csv",
		"_DI.csv",
		"_transq.csv",
	]
	output_name = [
		["Total", "Total energy"],
		["HF", "HF energy"],
		["CR", "Correlation energy"],
		["ES", "Electrostatic energy"],
		["EX", "Exchange repulsion energy"],
		["CT", "Charge transfer energy"],
		["DI", "Dispersion force"],
		["Q", "Transfer charge"]
	]

	if args.flag_all:
		output_flag = [True for x in output_flag]
	elif args.flag_total == False:
		if len([True for x in output_flag if x == True]) == 0:
			# 他のオプションが未指定の場合のみ total オプションを機能させる
			output_flag[0] = True


	# データ読み込み＆解析
	energy = EnergyData.EnergyData(args.input)


	# 出力フラグメントの決定
	output_range = []
	if args.include is None and args.exclude is None:
		# 未指定の場合
		output_range = energy.get_label()

	elif args.include is not None:
		# include で指定の場合
		output_range = [int(x) for x in args.include]

	elif args.exclude is not None:
		# exclude で指定の場合
		output_range = list(set(energy.get_label()) - set([int(x) for x in args.exclude]))


	# 出力ファイル
	prefix = ""
	if args.prefix != None:
		prefix = args.prefix
	else:
		prefix = re.sub(r"\..{3,4}$", "", args.input)

	import csv
	for idx, flag in enumerate(output_flag):
		if flag:
			output = prefix + output_suffix[idx]
			if args.flag_overwrite == False:
				basic_func.check_overwrite(output)
			with open(output, "w") as obj_output:
				csv_writer = csv.writer(obj_output, lineterminator = "\n")
				csv_writer.writerows(energy.output_energy(output_name[idx][0], output_range))
				sys.stderr.write("create: {0} ({1})\n".format(output, output_name[idx][1]))

	if args.flag_q:
		output = prefix + "_pc.csv"
		if args.flag_overwrite == False:
			basic_func.check_overwrite(output)
		with open(output, "w") as obj_output:
			csv_writer = csv.writer(obj_output, lineterminator = "\n")
			csv_writer.writerows(energy.output_charge(output_range))
			sys.stderr.write("create: %s (partial charge)\n" % output)
