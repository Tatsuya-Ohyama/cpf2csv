#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EnergyData class
"""

import sys, re, copy
import numpy as np

# =============== const =============== #
au = 627.5095
digit = 4
re_atomic_charge = re.compile(r"\d+[\s\t]+\D+(:?[\s\t]+-?\d+\.\d+){2}")


# =============== classes =============== #
class EnergyData:
	""" エネルギーデータを扱うクラス """
	def __init__(self, input_file):
		self._frag_atom = []
		self._label = []
		self._energy_HF = None
		self._energy_CR = None
		self._energy_ES = None
		self._energy_EX = None
		self._energy_CT = None
		self._energy_DI = None
		self._energy_Q = None

		self._charge_atom = []
		self._charge_frag = []

		self._load_file(input_file)

	def _load_file(self, input_file):
		""" ファイルを読み込むメソッド """
		re_fragmentation = re.compile(r"Frag\.   Elec\.   ATOM")
		re_IFIE = re.compile(r"## ((HF)|(MP2))-IFIE")
		re_PIEDA = re.compile(r"## PIEDA")
		re_separator = re.compile(r"-{5,}")
		re_charge = re.compile(r"No\. Atom   Atomic pop\.  Net charge")
		distance_idx = []

		flag_read = [0,0]
		with open(input_file, "r") as obj_input:
			for line_idx, line in enumerate(obj_input):
				if "ERROR" in line:
					sys.stderr.write("ERROR: ERROR in .log at {0}.\n".format(line_idx + 1))
					sys.exit(1)

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
						self._label.append(int(label))
						self._frag_atom.append(atoms)
					else:
						self._frag_atom[-1].extend(atoms)

				elif flag_read[0] == 2:
					# IFIE
					if flag_read[1] == 0:
						# 初期化
						self._energy_HF = np.zeros((len(self._frag_atom), len(self._frag_atom)))
						self._energy_CR = np.zeros((len(self._frag_atom), len(self._frag_atom)))
						self._charge_frag = [0.0 for i in range(len(self._frag_atom))]
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

						self._energy_HF[i][j] = self._energy_HF[j][i] = energies[0]
						self._energy_CR[i][j] = self._energy_CR[j][i] = energies[1]

				elif flag_read[0] == 3:
					# PIDA
					if flag_read[1] == 0:
						# 初期化
						self._energy_ES = np.zeros((len(self._frag_atom), len(self._frag_atom)))
						self._energy_EX = np.zeros((len(self._frag_atom), len(self._frag_atom)))
						self._energy_CT = np.zeros((len(self._frag_atom), len(self._frag_atom)))
						self._energy_DI = np.zeros((len(self._frag_atom), len(self._frag_atom)))
						self._energy_Q = np.zeros((len(self._frag_atom), len(self._frag_atom)))
						self._charge_frag = [0.0 for i in range(len(self._frag_atom))]
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

						self._energy_ES[i][j] = self._energy_ES[j][i] = energies[0]
						self._energy_EX[i][j] = self._energy_EX[j][i] = energies[1]
						self._energy_CT[i][j] = self._energy_CT[j][i] = energies[2]
						self._energy_DI[i][j] = self._energy_DI[j][i] = energies[3]
						self._energy_Q[i][j] = float(line[78:93].strip())
						self._energy_Q[i][j] = -1 * float(line[78:93].strip())

				elif flag_read[0] == 4:
					# 電荷
					if len(line.strip()) == 0:
						flag_read = [0,0]
						continue

					if re_atomic_charge.search(line):
						atom_idx = int(line[:13].strip())
						charge = float(line[31:].strip())
						data_idx = [idx for idx, value in enumerate(self._frag_atom) if atom_idx in value][0]
						self._charge_atom.append([atom_idx, line[14:19].strip(), charge])
						self._charge_frag[data_idx] += charge

	def get_label(self, frag_idx = None):
		""" ラベルを返すメソッド """
		if(frag_idx is None):
			return self._label
		else:
			return self._label[frag_idx - 1]

	def get_fragment_atom(self, frag_idx = None):
		""" フラグメント構成原子を返すメソッド """
		if(frag_idx is None):
			return self._frag_atom
		else:
			return self._frag_atom[frag_idx - 1]

	def get_atom_charge(self, atom_idx = None):
		""" 原子電荷を返すメソッド """
		if atom_idx is None:
			return self._charge_atom
		else:
			return self._charge_atom[atom_idx - 1]

	def get_fragment_charge(self, frag_idx = None):
		""" フラグメント電荷を返すメソッド """
		if frag_idx is None:
			return self._charge_frag
		else:
			return self._charge_frag[frag_idx - 1]

	def get_energy(self, energy_type = "Total", frag_idx = None):
		""" IFIE エネルギーを返すメソッド """
		energies = None
		if energy_type == "Total":
			energies = (self._energy_HF + self._energy_CR) * au
		elif energy_type == "HF":
			energies = self._energy_HF * au
		elif energy_type == "CR":
			energies = self._energy_CR * au
		elif energy_type == "ES":
			energies = self._energy_ES
		elif energy_type == "EX":
			energies = self._energy_EX
		elif energy_type == "CT":
			energies = self._energy_CT
		elif energy_type == "DI":
			energies = self._energy_DI
		elif energy_type == "Q":
			energies = self._energy_Q

		if frag_idx is None:
			return np.round(energies, digit)
		else:
			return np.round(energies[frag_idx[0] - 1][frag_idx[1] - 1], digit)

	def output_energy(self, energy_type = "Total", output_range = None):
		""" IFIE エネルギーを出力形式で返すメソッド """
		result = None

		if energy_type == "Total":
			result = (self.get_energy("Total"))
		elif energy_type == "HF":
			result = self.get_energy("HF")
		elif energy_type == "CR":
			result = self.get_energy("CR")
		elif energy_type == "ES":
			result = self.get_energy("ES")
		elif energy_type == "EX":
			result = self.get_energy("EX")
		elif energy_type == "CT":
			result = self.get_energy("CT")
		elif energy_type == "DI":
			result = self.get_energy("DI")
		elif energy_type == "Q":
			result = self.get_energy("Q")

		label = copy.deepcopy(self.get_label())
		if output_range is not None:
			delete_range = list(set(self.get_label()) - set(output_range))
			for idx in reversed(sorted(delete_range)):
				del(label[idx - 1])
				result = np.delete(result, idx - 1, 0)
				result = np.delete(result, idx - 1, 1)

		result = result.tolist()
		result = [[label[idx]] + value for idx, value in enumerate(result)]
		result = [[""] + label] + result
		return result

	def output_charge(self, output_range = None):
		""" 電荷情報を出力形式で返すメソッド """
		result = []
		result_frag = [["Fragment index", "Fragment charge", ""]]
		result_atom = [["Fragment index", "Atom index", "Atom", "Atomic charge"]]

		cnt_atom = 0
		for frag_idx in range(len(self._frag_atom)):
			if output_range is None or frag_idx + 1 in output_range:
				result_frag.append([frag_idx + 1, self._charge_frag[frag_idx], ""])
				for atom_idx in self._frag_atom[frag_idx]:
					result_atom.append([frag_idx + 1, self._charge_atom[atom_idx - 1][0], self._charge_atom[atom_idx - 1][1], self._charge_atom[atom_idx - 1][2]])

		diff_row = len(result_atom) - len(result_frag)
		if 0 < diff_row:
			result_frag += [["", "", ""] for x in range(diff_row)]
		elif diff_row < 0:
			result_atom += [["", "", "", ""] for x in range(diff_row)]

		return [x + y for x, y in zip(result_frag, result_atom)]
