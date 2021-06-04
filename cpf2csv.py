#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cpf2csv.py
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import os
import csv

from mods.basic_func import *
from mods.FileLogABINITMP import FileLogABINITMP
from mods.FileCpf import FileCpf

from pprint import pprint



# =============== constant =============== #
OUTPUT_SUFFIX = [
	"_Total.csv",
	"_HF.csv",
	"_CR.csv",
	"_ES.csv",
	"_EX.csv",
	"_CT.csv",
	"_DI.csv",
	"_transq.csv",
	"_pc.csv"
]
OUTPUT_NAME = [
	["Total", "Total energy"],
	["HF", "HF energy"],
	["CR", "Correlation energy"],
	["ES", "Electrostatic energy"],
	["EX", "Exchange repulsion energy"],
	["CT", "Charge transfer energy"],
	["DI", "Dispersion force"],
	["Q", "Transfer charge"],
	["PC", "Particle charge"]
]



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="cpf2csv - convert log for ABINIT-MP to CSV", formatter_class=argparse.RawTextHelpFormatter)
	global_option = parser.add_argument_group(title="global option", description="")
	global_option.add_argument("-i", dest="INPUT", metavar="INPUT.(log|out|cpf)", required=True, help=".log, .out or .cpf for ABINIT-MP")
	global_option.add_argument("-o", dest="PREFIX", help="prefix for output")
	global_option.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly (Default: False)")

	output_type = parser.add_argument_group(title="energy type option", description="energy type for output (default: -t)")
	output_type.add_argument("-a", "--all", dest="FLAG_ALL", action="store_true", default=False, help="select all type, same as -tfesxcdq")
	output_type.add_argument("-t", "--total", dest="FLAG_TOTAL", action="store_true", default=False, help="total energy (kcal/mol)")
	output_type.add_argument("-f", "--hartree", dest="FLAG_HF", action="store_true", default=False, help="Hartree-Fock energy (kcal/mol)")
	output_type.add_argument("-e", "--correlation", dest="FLAG_CORR", action="store_true", default=False, help="electron correlation energy (kcal/mol)")
	output_type.add_argument("-s", "--electrostatic", dest="FLAG_ES", action="store_true", default=False, help="electrostatic interaction (ES) (kcal/mol)")
	output_type.add_argument("-x", "--exchange", dest="FLAG_EX", action="store_true", default=False, help="exchange-repulsion energy (EX) (kcal/mol)")
	output_type.add_argument("-c", "--chargetransfer-mix", dest="FLAG_CT", action="store_true", default=False, help="charge transfer and other interaction energy (CT+mix) (kcal/mol)")
	output_type.add_argument("-d", "--dispersion", dest="FLAG_DI", action="store_true", default=False, help="dispersion energy (DI) (kcal/mol)")
	output_type.add_argument("-q", "--chargetransfer-amount", dest="FLAG_Q", action="store_true", default=False, help="amount of charge transfer (e; I(row) -> J(col))")
	output_type.add_argument("-pc", "--partial-charge", dest="FLAG_PC", action="store_true", default=False, help="partial charge")

	output_range = parser.add_mutually_exclusive_group()
	output_range.add_argument("--include", dest="INCLUDE", metavar="Frag_No.", nargs="+", help="")
	output_range.add_argument("--exclude", dest="EXCLUDE", metavar="Frag_No.", nargs="+", help="")

	args = parser.parse_args()

	check_exist(args.INPUT, 2)

	output_flag = [
		args.FLAG_TOTAL,
		args.FLAG_HF,
		args.FLAG_CORR,
		args.FLAG_ES,
		args.FLAG_EX,
		args.FLAG_CT,
		args.FLAG_DI,
		args.FLAG_Q,
		args.FLAG_PC
	]

	if args.FLAG_ALL:
		output_flag = [True for x in output_flag]
	elif args.FLAG_TOTAL == False:
		if len([True for x in output_flag if x == True]) == 0:
			# 他のオプションが未指定の場合のみ total オプションを機能させる
			output_flag[0] = True

	# データ読み込み＆解析
	data_FMO = None
	if os.path.splitext(args.INPUT)[1] == ".cpf":
		data_FMO = FileCpf(args.INPUT)

	else:
		data_FMO = FileLogABINITMP(args.INPUT)

	# 出力フラグメントの決定
	output_range = []
	if args.INCLUDE is None and args.EXCLUDE is None:
		# 未指定の場合
		output_range = data_FMO.get_label()

	elif args.INCLUDE is not None:
		# include で指定の場合
		output_range = [int(x) for x in args.INCLUDE]

	elif args.EXCLUDE is not None:
		# exclude で指定の場合
		output_range = list(set(data_FMO.get_label()) - set([int(x) for x in args.EXCLUDE]))

	# 出力ファイル
	prefix = args.PREFIX
	if prefix is None:
		prefix = os.path.splitext(os.path.basename(args.INPUT))[0]

	for idx, flag in enumerate(output_flag):
		if flag:
			output = prefix + OUTPUT_SUFFIX[idx]
			if args.FLAG_OVERWRITE == False:
				check_overwrite(output)

			if OUTPUT_NAME[idx][0] == "PC":
				with open(output, "w") as obj_output:
					csv_writer = csv.writer(obj_output, lineterminator="\n")
					csv_writer.writerows(data_FMO.output_charge(output_range))
					sys.stderr.write("create: %s (partial charge)\n" % output)
			else:
				with open(output, "w") as obj_output:
					csv_writer = csv.writer(obj_output, lineterminator="\n")
					csv_writer.writerows(data_FMO.output_energy(OUTPUT_NAME[idx][0], output_range))
					sys.stderr.write("create: {0} ({1})\n".format(output, OUTPUT_NAME[idx][1]))
