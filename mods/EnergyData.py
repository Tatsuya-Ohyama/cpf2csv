#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EnergyData class
"""

import sys, re, copy
import numpy as np



# =============== const =============== #
AU = 627.5095
DIGIT = 4
RE_ATOMIC_CHARGE = re.compile(r"\d+[\s\t]+\D+(:?[\s\t]+-?\d+\.\d+){2}")
RE_IFIE = re.compile(r"## ((HF)|(MP2))-IFIE")



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
		"""
		ファイルを読み込むメソッド

		Args:
			input_file (str): ABINIT-MP の .out および .log ファイル

		Returns:
			self
		"""
		distance_idx = []

		flag_read = [0,0]
		with open(input_file, "r") as obj_input:
			for line_idx, line_val in enumerate(obj_input):
				if "ERROR" in line_val:
					sys.stderr.write("ERROR: ERROR in .log at {0}.\n".format(line_idx + 1))
					sys.exit(1)

				if flag_read[0] == 0:
					# フラグ分類
					if "Frag.   Elec.   ATOM" in line_val:
						flag_read[0] = 1
					elif RE_IFIE.search(line_val):
						flag_read[0] = 2
					elif "## PIEDA" in line_val:
						flag_read[0] = 3
					elif "No. Atom   Atomic pop.  Net charge" in line_val:
						flag_read[0] = 4

				elif flag_read[0] == 1:
					# フラグメント構成原子の取得
					if len(line_val.strip()) == 0:
						flag_read = [0,0]
						continue

					label = line_val[5:13].strip()
					atoms = [int(x) for x in line_val[23:].strip().split()]
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

					elif "------" in line_val:
						flag_read[1] = 2

					elif flag_read[1] == 2:
						if len(line_val.strip()) == 0:
							flag_read = [0,0]
							continue

						i = int(line_val[8:13].strip()) - 1
						j = int(line_val[13:18].strip()) - 1
						distance = float(line_val[18:30].strip())

						energies = [float(x.strip()) for x in [line_val[39:50], line_val[50:61]]]

						if distance == 0.000000:
							distance_idx.append(line_val[8:18].strip())
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

					elif "------" in line_val:
						flag_read[1] = 2

					elif flag_read[1] == 2:
						if len(line_val.strip()) == 0:
							flag_read = [0,0]
							continue

						i = int(line_val[8:13].strip()) - 1
						j = int(line_val[13:18].strip()) - 1

						energies = [float(x.strip()) for x in [
							line_val[18:33],
							line_val[33:48],
							line_val[48:63],
							line_val[63:78],
							line_val[78:93]
						]]
						if line_val[8:18].strip() in distance_idx:
							energies = [0.0 for x in energies]

						self._energy_ES[i][j] = self._energy_ES[j][i] = energies[0]
						self._energy_EX[i][j] = self._energy_EX[j][i] = energies[1]
						self._energy_CT[i][j] = self._energy_CT[j][i] = energies[2]
						self._energy_DI[i][j] = self._energy_DI[j][i] = energies[3]
						self._energy_Q[i][j] = float(line_val[78:93].strip())
						self._energy_Q[i][j] = -1 * float(line_val[78:93].strip())

				elif flag_read[0] == 4:
					# 電荷
					if len(line_val.strip()) == 0:
						flag_read = [0,0]
						continue

					if RE_ATOMIC_CHARGE.search(line_val):
						atom_idx = int(line_val[:13].strip())
						charge = float(line_val[31:].strip())
						data_idx = [idx for idx, value in enumerate(self._frag_atom) if atom_idx in value][0]
						self._charge_atom.append([atom_idx, line_val[14:19].strip(), charge])
						self._charge_frag[data_idx] += charge
		return self


	def get_label(self, frag_idx=None):
		"""
		ラベルを返すメソッド

		Args:
			frag_idx (int, optional): フラグメントインデックス (Default: None)

		Returns:
			str: ラベル
		"""
		if(frag_idx is None):
			return self._label
		else:
			return self._label[frag_idx - 1]


	def get_fragment_atom(self, frag_idx=None):
		"""
		フラグメント構成原子を返すメソッド

		Args:
			frag_idx (int, optional): フラグメントインデックス (Default: None)

		Returns:
			list: 構成原子のリスト
		"""
		if(frag_idx is None):
			return self._frag_atom
		else:
			return self._frag_atom[frag_idx - 1]


	def get_atom_charge(self, atom_idx=None):
		"""
		原子電荷を返すメソッド

		Args:
			atom_idx (int, optional): 原子インデックス (Default: None)

		Returns:
			float: 原子電荷
		"""
		if atom_idx is None:
			return self._charge_atom
		else:
			return self._charge_atom[atom_idx - 1]


	def get_fragment_charge(self, frag_idx=None):
		"""
		フラグメント電荷を返すメソッド

		Args:
			frag_idx (int, optional): フラグメントインデックス (Default: None)

		Returns:
			float: フラグメント電荷
		"""
		if frag_idx is None:
			return self._charge_frag
		else:
			return self._charge_frag[frag_idx - 1]


	def get_energy(self, energy_type="Total", frag_idx=None):
		"""
		IFIE エネルギーを返すメソッド

		Args:
			energy_type (str, optional): `Total`, `ES`, `EX`, `CT`, `DI` or `Q` (Default: "Total")
			frag_idx (int, optional): フラグメントインデックス (Default: None)

		Returns:
			list
		"""
		energies = None
		if energy_type == "Total":
			energies = (self._energy_HF + self._energy_CR) * AU
		elif energy_type == "HF":
			energies = self._energy_HF * AU
		elif energy_type == "CR":
			energies = self._energy_CR * AU
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
			return np.round(energies, DIGIT)
		else:
			return np.round(energies[frag_idx[0] - 1][frag_idx[1] - 1], DIGIT)


	def output_energy(self, energy_type="Total", output_range=None):
		"""
		IFIE エネルギーを出力形式で返すメソッド

		Args:
			energy_type (str, optional): `Total`, `ES`, `EX`, `CT`, `DI` or `Q` (Default: "Total")
			output_range (list, optional): 出力するラベルリスト (Default: None)

		Returns:
			list
		"""
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


	def output_charge(self, output_range=None):
		"""
		電荷情報を出力形式で返すメソッド

		Args:
			output_range (list, optional): 出力するラベルリスト (Default: None)

		Returns:
			list
		"""
		result = []
		result_frag = [["Fragment index", "Fragment charge", ""]]
		result_atom = [["Fragment index", "Atom index", "Atom", "Atomic charge"]]

		cnt_atom = 0
		for frag_idx in range(len(self._frag_atom)):
			if output_range is None or frag_idx + 1 in output_range:
				result_frag.append([frag_idx + 1, self._charge_frag[frag_idx], ""])
				for atom_idx in self._frag_atom[frag_idx]:
					result_atom.append([
						frag_idx + 1,
						self._charge_atom[atom_idx - 1][0],
						self._charge_atom[atom_idx - 1][1],
						self._charge_atom[atom_idx - 1][2]
					])

		diff_row = len(result_atom) - len(result_frag)
		if 0 < diff_row:
			result_frag += [["", "", ""] for x in range(diff_row)]
		elif diff_row < 0:
			result_atom += [["", "", "", ""] for x in range(diff_row)]

		return [x + y for x, y in zip(result_frag, result_atom)]
