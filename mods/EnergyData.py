#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EnergyData class
"""

import sys, re
import numpy as np

# =============== const =============== #
au = 627.5095
digit = 4

class EnergyData:
	""" エネルギーデータを扱うクラス """
	def __init__(self, input_file):
		self.__frag_atom = []
		self.__label = []
		self.__energy_HF = None
		self.__energy_CR = None
		self.__energy_ES = None
		self.__energy_EX = None
		self.__energy_CT = None
		self.__energy_DI = None
		self.__energy_Q = None

		self.__charge_atom = []
		self.__charge_frag = []

		self.__load_file(input_file)

	def __load_file(self, input_file):
		""" ファイルを読み込むメソッド """
		re_fragmentation = re.compile(r"Frag\.   Elec\.   ATOM")
		re_IFIE = re.compile(r"## ((HF)|(MP2))-IFIE")
		re_PIEDA = re.compile(r"## PIEDA")
		re_separator = re.compile(r"-{5,}")
		re_charge = re.compile(r"No\. Atom   Atomic pop\.  Net charge")
		distance_idx = []

		flag_read = [0,0]
		with open(input_file, "r") as obj_input:
			for line in obj_input:
				if flag_read[0] == 0:
					# フラグ分類
					if re_fragmentation.search(line):
						flag_read[0] = 1
					elif re_IFIE.search(line):
						flag_read[0] = 2
					elif re_PIEDA.search(line):
						flag_read[0] = 3
					elif re_charge.search(line):
						flag_read[0] = 4

				elif flag_read[0] == 1:
					# フラグメント構成原子の取得
					if len(line.strip()) == 0:
						flag_read = [0,0]
						continue

					label = line[5:13].strip()
					atoms = [int(x) for x in line[23:].strip().split()]
					if label:
						self.__label.append(label)
						self.__frag_atom.append(atoms)
					else:
						self.__frag_atom[-1].extend(atoms)

				elif flag_read[0] == 2:
					# IFIE
					if flag_read[1] == 0:
						# 初期化
						self.__energy_HF = np.zeros((len(self.__frag_atom), len(self.__frag_atom)))
						self.__energy_CR = np.zeros((len(self.__frag_atom), len(self.__frag_atom)))
						self.__charge_frag = [0.0 for i in range(len(self.__frag_atom))]
						flag_read[1] = 1

					elif re_separator.search(line):
						flag_read[1] = 2

					elif flag_read[1] == 2:
						if len(line.strip()) == 0:
							flag_read = [0,0]
							continue

						i = int(line[8:13].strip()) - 1
						j = int(line[13:18].strip()) - 1
						distance = float(line[18:30].strip())

						energies = [float(x.strip()) for x in [line[39:50], line[50:61]]]

						if distance == 0.000000:
							distance_idx.append(line[8:18].strip())
							energies = [0.0 for x in energies]

						self.__energy_HF[i][j] = self.__energy_HF[j][i] = energies[0]
						self.__energy_CR[i][j] = self.__energy_CR[j][i] = energies[1]

				elif flag_read[0] == 3:
					# PIDA
					if flag_read[1] == 0:
						# 初期化
						self.__energy_ES = np.zeros((len(self.__frag_atom), len(self.__frag_atom)))
						self.__energy_EX = np.zeros((len(self.__frag_atom), len(self.__frag_atom)))
						self.__energy_CT = np.zeros((len(self.__frag_atom), len(self.__frag_atom)))
						self.__energy_DI = np.zeros((len(self.__frag_atom), len(self.__frag_atom)))
						self.__energy_Q = np.zeros((len(self.__frag_atom), len(self.__frag_atom)))
						self.__charge_frag = [0.0 for i in range(len(self.__frag_atom))]
						flag_read[1] = 1

					elif re_separator.search(line):
						flag_read[1] = 2

					elif flag_read[1] == 2:
						if len(line.strip()) == 0:
							flag_read = [0,0]
							continue

						i = int(line[8:13].strip()) - 1
						j = int(line[13:18].strip()) - 1

						energies = [float(x.strip()) for x in [line[18:33], line[33:48], line[48:63], line[63:78], line[78:93]]]
						if line[8:18].strip() in distance_idx:
							energies = [0.0 for x in energies]

						self.__energy_ES[i][j] = self.__energy_ES[j][i] = energies[0]
						self.__energy_EX[i][j] = self.__energy_EX[j][i] = energies[1]
						self.__energy_CT[i][j] = self.__energy_CT[j][i] = energies[2]
						self.__energy_DI[i][j] = self.__energy_DI[j][i] = energies[3]
						self.__energy_Q[i][j] = float(line[78:93].strip())
						self.__energy_Q[i][j] = -1 * float(line[78:93].strip())

				elif flag_read[0] == 4:
					# 電荷
					if len(line.strip()) == 0:
						flag_read = [0,0]
						continue

					atom_idx = int(line[:13].strip())
					data_idx = [idx for idx, value in enumerate(self.__frag_atom) if atom_idx in value][0]
					charge = float(line[31:].strip())
					self.__charge_atom.append(charge)
					self.__charge_frag[data_idx] += charge

	def get_label(self, frag_idx = None):
		""" ラベルを返すメソッド """
		if(frag_idx is None):
			return self.__label
		else:
			return self.__label[frag_idx - 1]

	def get_fragment_atom(self, frag_idx = None):
		""" フラグメント構成原子を返すメソッド """
		if(frag_idx is None):
			return self.__frag_atom
		else:
			return self.__frag_atom[frag_idx - 1]

	def get_atom_charge(self, atom_idx = None):
		""" 原子電荷を返すメソッド """
		if atom_idx is None:
			return self.__charge_atom
		else:
			return self.__charge_atom[atom_idx - 1]

	def get_fragment_charge(self, frag_idx = None):
		""" フラグメント電荷を返すメソッド """
		if frag_idx is None:
			return self.__charge_frag
		else:
			return self.__charge_frag[frag_idx - 1]

	def get_energy(self, energy_type = "Total", frag_idx = None):
		""" IFIE エネルギーを返すメソッド """
		energies = None
		if energy_type == "Total":
			energies = (self.__energy_HF + self.__energy_CR) * au
		elif energy_type == "HF":
			energies = self.__energy_HF * au
		elif energy_type == "CR":
			energies = self.__energy_CR * au
		elif energy_type == "ES":
			energies = self.__energy_ES * au
		elif energy_type == "EX":
			energies = self.__energy_EX * au
		elif energy_type == "CT":
			energies = self.__energy_CT * au
		elif energy_type == "DI":
			energies = self.__energy_DI * au
		elif energy_type == "Q":
			energies = self.__energy_Q

		if frag_idx is None:
			return np.round(energies, digit)
		else:
			return np.round(energies[frag_idx[0]][frag_idx[1]], digit)

	def output_energy(self, energy_type = "Total"):
		""" IFIE エネルギーを出力形式で返すメソッド """
		result = [[""] + energy.get_label()]
		row_label = np.matrix(np.array(self.get_label())).T

		if energy_type == "Total":
			result += np.concatenate((row_label, self.get_energy("Total")), axis = 1).tolist()
		if energy_type == "HF":
			result += np.concatenate((row_label, self.get_energy("HF")), axis = 1).tolist()
		elif energy_type == "CR":
			result += np.concatenate((row_label, self.get_energy("CR")), axis = 1).tolist()
		elif energy_type == "ES":
			result += np.concatenate((row_label, self.get_energy("ES")), axis = 1).tolist()
		elif energy_type == "EX":
			result += np.concatenate((row_label, self.get_energy("EX")), axis = 1).tolist()
		elif energy_type == "CT":
			result += np.concatenate((row_label, self.get_energy("CT")), axis = 1).tolist()
		elif energy_type == "DI":
			result += np.concatenate((row_label, self.get_energy("DI")), axis = 1).tolist()
		elif energy_type == "Q":
			result += np.concatenate((row_label, self.get_energy("Q")), axis = 1).tolist()
		return result
