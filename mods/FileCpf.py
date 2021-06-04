#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CPF File class
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)
import json
import numpy as np
import collections
import copy
from pprint import pprint



# =============== constant =============== #
AU_TO_KCAL = 627.509468804
BOHR_RADIUS = 0.52911772
DIGIT = 4
STRUCTURE_COLUMNS = [
	'Index',
	'Element',
	'AtomName',
	'ResidueName',
	'ResidueNumber',
	'FragmentNumber',
	'CoordinateX',
	'CoordinateY',
	'CoordinateZ',
	'HF_MullikenCharge',
	'MP2_MullikenCharge',
	'HF_NBOCharge',
	'MP2_NBOCharge',
	'HF_ESPCharge',
	'MP2_ESPCharge',
	'ChainID',
	'PDBInsertionCode'
]

CPF_VERSION = {
	'CPF Ver.4.201 (MIZUHO)': 'CPFVersion.Ver4_201_MIZUHO',
	'CPF Ver.7.0 (MIZUHO 4.0) PIEDA': 'CPFVersion.Ver70_MIZUHO40_PIEDA',
	'CPF Open1.0 rev10': 'CPFVersion.Ver_Open1_0_rev10',
}

CPF_FORMAT = {
	"ELECTRON": {
		"length": 5,
		"number": 16
	}
}

IFIE_FORMAT = {
	'CPF Ver.4.2':                    ['Repulsion', 'HF-Electron', 'HF-ES', 'MP2-IFIE', 'SCS-MP2-IFIE', 'MP3-IFIE'    , 'SCS-MP3-IFIE', 'HF-IFIE-BSSE', 'MP2-IFIE-BSSE', 'SCS-IFIE-BSSE', 'MP3-IFIE-BSSE', 'SCS-MP3-IFIE-BSSE'],
	'CPF Ver.4.201':                  ['Repulsion', 'HF-Electron', 'HF-ES', 'MP2-IFIE', 'SCS-MP2-IFIE', 'MP3-IFIE'    , 'SCS-MP3-IFIE', 'HF-IFIE-BSSE', 'MP2-IFIE-BSSE', 'SCS-IFIE-BSSE', 'MP3-IFIE-BSSE', 'SCS-MP3-IFIE-BSSE', 'PIEDA-EX', 'PIEDA-CT', 'PIEDA-dq'],
	'CPF Ver.4.201 (MIZUHO)':         ['Repulsion', 'HF-Electron', 'HF-ES', 'MP2-IFIE', 'SCS-MP2-IFIE', 'MP3-IFIE'    , 'SCS-MP3-IFIE', 'HF-IFIE-BSSE', 'MP2-IFIE-BSSE', 'SCS-IFIE-BSSE', 'MP3-IFIE-BSSE', 'SCS-MP3-IFIE-BSSE', 'PIEDA-EX', 'PIEDA-CT', 'PIEDA-dq'],
	'CPF Ver.7.0 (MIZUHO 4.0) PIEDA': ['Repulsion', 'HF-Electron', 'HF-ES', 'MP2-IFIE', 'SCS-MP2-IFIE', 'MP3-IFIE'    , 'SCS-MP3-IFIE', 'HF-IFIE-BSSE', 'MP2-IFIE-BSSE', 'SCS-IFIE-BSSE', 'MP3-IFIE-BSSE', 'SCS-MP3-IFIE-BSSE', 'Solv-ES', 'Solv-NP', 'PIEDA-EX', 'PIEDA-CT', 'PIEDA-dq'],
	'CPF Ver.7.2':                    ['Repulsion', 'HF-Electron', 'HF-ES', 'MP2-IFIE', 'SCS-MP2-IFIE', 'MP3-IFIE'    , 'SCS-MP3-IFIE', 'HF-IFIE-BSSE', 'MP2-IFIE-BSSE', 'SCS-IFIE-BSSE', 'MP3-IFIE-BSSE', 'SCS-MP3-IFIE-BSSE', 'Solv-ES', 'Solv-NP', 'PIEDA-EX', 'PIEDA-CT', 'PIEDA-dq'],
	'CPF Open1.0 rev10':              ['Repulsion', 'HF-Electron', 'HF-ES', 'MP2-IFIE', 'PR-MP2-IFIE' , 'SCS-MP2-IFIE', 'MP3-IFIE', 'SCS-MP3-IFIE', 'HF-IFIE-BSSE', 'MP2-IFIE-BSSE', 'SCS-MP2-IFIE-BSSE', 'MP3-IFIE-BSSE', 'SCS-MP3-IFIE-BSSE', 'Solv-ES', 'Solv-NP', 'PIEDA-EX', 'PIEDA-CT', 'PIEDA-dq']
}


ENERGY_TYPE = {
	"ES": "HF-ES",
	"EX": "PIEDA-EX",
	"CT": "PIEDA-CT",
	"DI": "MP2-IFIE",
	"Q": "PIEDA-dq"
}



# =============== function =============== #
def parser_split_line_by_length(line, length_value, dtype):
	"""
	function of parser for fixed length

	Args:
		line (str): line
		length_value (int): length
		dtype (str): data type for return values

	Returns:
		list: [val(str), ...]
	"""
	list_length = []
	list_dtype = []
	if type(length_value) == int and type(dtype) == str:
		# 単体指定、繰り返し分割の場合 (同じ長さの固定長の場合)
		if len(line) % length_value == 0:
			list_length = [length_value for _ in range(len(line) // length_value)]
			list_dtype = [dtype for _ in range(len(line) // length_value)]
		else:
			sys.stderr.write("ERROR: the strings in the following lines are not divisible by the specified length.\n")
			sys.stderr.write("       {0}\n".format(line))
			sys.exit(1)

	else:
		# リストで与えられた場合 (様々な長さの固定長の場合)
		if len(length_value) != len(dtype):
			sys.stderr.write("ERROR: length of `length_value` and `` is not matched.\n")
			sys.exit(1)
		list_length = length_value
		list_dtype = dtype

	list_values = []
	pos_start = 0
	pos_end = 0
	for length_elem, dtype_elem in zip(list_length, list_dtype):
		pos_end += length_elem
		value = line[pos_start : pos_end]
		if dtype_elem == "int":
			list_values.append(int(value.strip()))
		elif dtype_elem == "float":
			list_values.append(float(value.strip()))
		else:
			list_values.append(value)
		pos_start = pos_end

	return list_values


def parser_structure(structure_line, version = None):
	"""
	function of parser for structure information in .cpf file

	Args:
		structure_line (str): structure line
		version (str, optional): CPF version (Default: None)

	Returns:
		list
	"""
	if version is None:
		return [
			int(structure_line[0:5].strip()),
			structure_line[6:8],
			structure_line[9:13],
			structure_line[14:17],
			int(structure_line[18:22]),
			int(structure_line[23:27]),
			float(structure_line[28:40].strip()),
			float(structure_line[40:52].strip()),
			float(structure_line[52:64].strip()),
			float(structure_line[64:76].strip()),
			float(structure_line[76:88].strip()),
			float(structure_line[88:100]),
			float(structure_line[100:112]),
			float(structure_line[112:124]),
			float(structure_line[124:136]),
			structure_line[136:138].strip(),
			structure_line[138:] if len(structure_line[135:]) > 16 else ' ',
		]
	else:
		return None


def parser_electron_number(electron_line):
	"""
	function of parser for electron information

	Args:
		electron_line (str): line

	Returns:
		list: electron information
	"""
	electron_list = parser_split_line_by_length(electron_line.rstrip(), CPF_FORMAT["ELECTRON"]["length"], "int")
	return electron_list


def parser_fragment_bond(bond_line):
	"""
	function of parser for bond information in CPF (the same behavior of electron information)

	Args:
		bond_line (str): bond line

	Returns:
		list: bond information
	"""
	return parser_electron_number(bond_line)


def parser_monomer(monomer_line):
	"""
	function of parser for monomer information

	Args:
		monomer_line (str): monomer line

	Returns:
		list: monomer information
	"""
	monomer_energy = parser_split_line_by_length(monomer_line[0 : 96], 24, "float")
	monomer_orbital = parser_split_line_by_length(monomer_line[96 : 120], 12, "int")
	return monomer_energy + monomer_orbital


def parser_trimer(trimer_line):
	"""
	function of parser for trimer information

	Args:
		trimer_line (str): trimer line

	Returns:
		list: trimer information
	"""
	trimer_ijk = parser_split_line_by_length(line_val[0 : 15], 5, "int")
	trimer_energy = parser_split_line_by_length(line_val[15 : 135], 24, "float")
	return trimer_ijk + trimer_energy


def parser_tetramer(tetramer_line):
	"""
	function of parser for tetramer information

	Args:
		tetramer_line (str): trimer line

	Returns:
		list: tetramer information
	"""
	tetramer_ijk = parser_split_line_by_length(line_val[0 : 20], 5, "int")
	tetramer_energy = parser_split_line_by_length(line_val[20 : 140], 24, "float")
	return tetramer_ijk + tetramer_energy



# =============== class =============== #
class FileCpf:
	""" CPF ファイルクラス """
	def __init__(self, cpf_file = None):
		self._path = None

		self._version = None
		self._n_atom = 0
		self._n_fragment = 0
		self._obj_fragments = []
		self._fragment_number_list = []
		self._basis_set = None
		self._stat = None
		self._method = None
		self._approx = [None, None, None]
		self._energy_total = {
			"repulsion": None,
			"electron": None,
			"whole": None
		}
		self._n_trimer = 0
		self._trimers = []
		self._n_tetramer = 0
		self._tetramers = []
		self._pair_list = []
		self._structure_columns = STRUCTURE_COLUMNS
		self._complete = False
		self.__cache_table = {}

		if cpf_file is not None:
			self._cpf_file = cpf_file
			self.read(cpf_file)

	@property
	def version(self):
		return self._version

	@property
	def n_atom(self):
		return self._n_atom

	@property
	def n_fragment(self):
		return self._n_fragment

	@property
	def structure(self):
		return [v.structure_info for v in self._obj_fragments]

	@property
	def electrons(self):
		return [v.electron for v in self._obj_fragments]

	@property
	def bonds(self):
		return [v.bond for v in self._obj_fragments]

	@property
	def neighbors(self):
		return [v.neighbor for v in self._obj_fragments]

	@property
	def dipole_info(self):
		return [v.dipole for v in self._obj_fragments]

	@property
	def basis_set(self):
		return self._basis_set

	@property
	def stat(self):
		return self._stat

	@property
	def method(self):
		return self._method

	@property
	def approx(self):
		return self._approx

	@property
	def energy_total(self):
		return self._energy_total

	@property
	def monomer_info(self):
		return [v.monomer for v in self._obj_fragments]

	@property
	def n_trimer(self):
		return self._n_trimer

	@property
	def trimers(self):
		return self._trimers

	@property
	def n_tetramer(self):
		return self._n_tetramer

	@property
	def tetramers(self):
		return self._tetramers

	@property
	def structure_columns(self):
		return self._structure_columns

	@property
	def fragments(self):
		return self._obj_fragments

	@property
	def is_completed(self):
		return self._complete


	def read(self, input_file):
		"""
		read .cpf file

		Args:
			input_file (str): .cpf file path

		Returns:
			self
		"""
		max_lines = [float('inf') for _ in range(20)]
		max_lines[0] = 1
		max_lines[1] = 2

		list_idx = 0

		total_electron = 0
		pair_idx = 0
		with open(input_file, "r") as obj_input:
			for line_idx, line_val in enumerate(obj_input, 1):
				if line_val.startswith("END"):
					self._complete = True
					break

				if max_lines[0] == line_idx:
					# バージョン
					self._version = line_val.strip()

				elif max_lines[1] == line_idx:
					# 構造概要
					values = parser_split_line_by_length(line_val.rstrip(), 5, "int")
					self._n_atom = values[0]
					self._n_fragment = values[1]
					self._pair_list = [[i, j] for i in range(1, self._n_fragment + 1) for j in range(1, self._n_fragment + 1) if i > j]
					max_lines[2] = max_lines[1] + self._n_atom
					max_lines[3] = max_lines[2] + np.ceil(self._n_fragment / CPF_FORMAT["ELECTRON"]["number"])
					max_lines[4] = max_lines[3] + np.ceil(self._n_fragment / CPF_FORMAT["ELECTRON"]["number"])

				elif max_lines[1] < line_idx <= max_lines[2]:
					# 原子情報
					structure_info = parser_structure(line_val)
					fragment_number = structure_info[5]
					if fragment_number not in self._fragment_number_list:
						# フラグメントオブジェクトが存在しない場合
						obj_fragment = Fragment(fragment_number)
						obj_fragment.append_atom(structure_info)
						self._fragment_number_list.append(fragment_number)
						self._obj_fragments.append(obj_fragment)
					else:
						# フラグメントオブジェクトが存在する場合
						list_idx_tmp = self._fragment_number_list.index(fragment_number)
						obj_fragment = self._obj_fragments[list_idx_tmp]
						obj_fragment.append_atom(structure_info)

				elif max_lines[2] < line_idx <= max_lines[3]:
					# 電子情報
					electron_info = parser_electron_number(line_val)
					for n_electron in electron_info:
						obj_fragment = self._obj_fragments[list_idx]
						obj_fragment.set_electron(n_electron)
						list_idx += 1

					if list_idx >= self._n_fragment:
						list_idx = 0

				elif max_lines[3] < line_idx <= max_lines[4]:
					# 結合情報
					bond_info = parser_fragment_bond(line_val)
					for n_bond in bond_info:
						obj_fragment = self._obj_fragments[list_idx]
						obj_fragment.set_bond(n_bond)
						list_idx += 1

					if list_idx >= self._n_fragment:
						list_idx = 0

				elif max_lines[4] < line_idx and (max_lines[6] == float('inf') or line_idx <= max_lines[6]):
					values = line_val.strip().split()
					if len(values) == 2:
						# フラグメント間接続
						values = [int(v) for v in values]
						fragment_info = [int(v) - 1 for v in values]
						if fragment_info[0] < 0 or fragment_info[1] < 0:
							sys.stderr.write("ERROR: Unexpected case at fragments `{0[0]}` and `{0[1]}`.\n".format(fragment_info))
							sys.exit(1)

						obj_fragments = [obj_fragment for atom_number in values for obj_fragment in self._obj_fragments if obj_fragment.has_atom(atom_number)]
						obj_fragments[0].append_neighbor(obj_fragments[1].number, values[1])
						obj_fragments[1].append_neighbor(obj_fragments[0].number, values[0])

					elif len(values) == 3:
						# フラグメント間距離
						values = [v for v in line_val.strip().split()]
						values[0] = int(values[0])
						values[1] = int(values[1])
						values[2] = float(values[2])
						self._obj_fragments[values[0] - 1].append_distance(self._obj_fragments[values[1] - 1], values[2])
						self._obj_fragments[values[1] - 1].append_distance(self._obj_fragments[values[0] - 1], values[2])

						if max_lines[5] == float('inf'):
							max_lines[5] = line_idx - 1
							max_lines[6] = max_lines[5] + int(self._n_fragment * (self._n_fragment - 1) / 2)
							max_lines[7] = max_lines[6] + self._n_fragment
							max_lines[7 : 14] = [i for i in range(max_lines[7], max_lines[7] + 8)]
							max_lines[15] = max_lines[14] + self._n_fragment
							max_lines[16] = max_lines[15] + int(self._n_fragment * (self._n_fragment - 1) / 2)
							max_lines[17] = max_lines[16] + 1

				elif max_lines[6] < line_idx <= max_lines[7]:
					# 双極子モーメント
					list_dmoment = []
					try:
						list_dmoment = [float(v) for v in line_val.strip().split()]
					except ValueError:
						break
					self._obj_fragments[list_idx].set_dipole_info(list_dmoment)
					list_idx += 1
					if list_idx >= self._n_fragment:
						list_idx = 0

				elif max_lines[8] == line_idx:
					# 基底関数
					self._basis_set = line_val.strip()

				elif max_lines[9] == line_idx:
					# 電子状態
					self._stat = line_val.strip()

				elif max_lines[10] == line_idx:
					# 手法
					self._method = line_val.strip()

				elif max_lines[11] == line_idx:
					# 近似
					self._approx = line_val.strip().split()

				elif max_lines[12] == line_idx:
					# 核反発エネルギー
					self._energy_total["repulsion"] = float(line_val.strip())

				elif max_lines[13] == line_idx:
					# 全電子エネルギー
					self._energy_total["electron"] = float(line_val.strip())

				elif max_lines[14] == line_idx:
					# 全エネルギー
					self._energy_total["whole"] = float(line_val.strip())

				elif max_lines[14] < line_idx <= max_lines[15]:
					# モノマー
					self._obj_fragments[list_idx].set_monomer_info(parser_monomer(line_val))
					list_idx += 1
					if list_idx >= self._n_fragment:
						list_idx = 0

				elif max_lines[15] < line_idx <= max_lines[16]:
					# IFIE
					info_IFIE = parser_split_line_by_length(line_val.rstrip(), 24, "float")
					obj_fragment1 = self._obj_fragments[self._pair_list[pair_idx][0] - 1]
					obj_fragment2 = self._obj_fragments[self._pair_list[pair_idx][1] - 1]
					obj_fragment1.append_IFIE(obj_fragment2, info_IFIE)
					obj_fragment2.append_IFIE(obj_fragment1, info_IFIE)
					pair_idx += 1

				elif max_lines[17] == line_idx:
					# n_trimer
					self._n_trimer = int(line_val.strip())
					max_lines[18] = max_lines[17] + (self._n_trimer * (self._n_trimer - 1) * (self._n_trimer - 2) / (3 * 2))
					max_lines[19] = max_lines[18] + 1

				elif max_lines[17] < line_idx <= max_lines[18]:
					# trimer data
					self._trimers.append(parser_trimer(line_val))

				elif max_lines[19] == line_idx:
					# n_tetramer
					self._n_tetramer = int(line_val.strip())
					max_lines[20] = max_lines[19] + (self._n_tetramer * (self._n_tetramer - 1) * (self._n_tetramer - 2) / (4 * 3 * 2))

				elif max_lines[19] < line_idx <= max_lines[20]:
					# tetramer
					self._tetramers.append(parser_tetramer(line_val))


	def get_structure_list(self, column_name):
		"""
		method for getting structure information

		Args:
			column_name (str): column name

		Returns:
			list: structure information
		"""
		target_idx = None
		if column_name in self._structure_columns:
			target_idx = self._structure_columns.index(column_name)
		else:
			sys.stderr.write("ERROR: undefined column name.\n")
			sys.exit(1)
		return [v[target_idx] for obj_fragment in self._obj_fragments for v in obj_fragment.structure_info]


	def extract_distance(self, fragment1, fragment2, unit="bohr"):
		"""
		フラグメント間の距離を取得するメソッド

		Args:
			fragment1 (int or obj_Fragment): フラグメント番号か、Fragment オブジェクト
			fragment2 (int or obj_Fragment): フラグメント番号か、Fragment オブジェクト
			unit (str): "bohr" or "angstrom"

		Returns:
			float: フラグメント間距離
		"""
		# フラグメント番号が一致するリストのリストのインデックスを取得する [[list_index, counter_fragment_number], ...]
		if fragment1 == fragment2:
			return 0.0

		obj_fragment1 = fragment1
		if isinstance(fragment1, int):
			obj_fragment1 = self._obj_fragments[fragment1 - 1]
		obj_fragment2 = fragment2
		if isinstance(fragment2, int):
			obj_fragment2 = self._obj_fragments[fragment2 - 1]

		distance = obj_fragment1.get_distance(obj_fragment2)
		if unit == "angstrom":
			return BOHR_RADIUS * distance
		else:
			return distance


	def extract_IFIE_energy(self, fragment, unit="a.u."):
		"""
		method for getting IFIE energy

		Args:
			fragment (int or obj_Fragment): フラグメント番号か、Fragment オブジェクト
			unit (str): "a.u." or "kcal/mol" (Default: "a.u.")

		Returns:
			list: [[fragment_pair_index(int), energy, ...], ...]
		"""
		# フラグメント番号が一致するリストのリストのインデックスを取得する [[list_index, counter_fragment_number], ...]
		obj_fragment = fragment
		if isinstance(fragment, int):
			obj_fragment = self._obj_fragments[fragment - 1]

		values = [[obj_fragment_other.number] + obj_fragment.get_IFIE(obj_fragment_other, no_data="none") for obj_fragment_other in self._obj_fragments if obj_fragment.get_IFIE(obj_fragment_other, no_data="none") is not None]
		if unit == "a.u.":
			return values
		elif unit == "kcal/mol":
			index_dq = IFIE_FORMAT[self._version].index("PIEDA-dq") + 1
			return [[v if i == index_dq else v2 * AU_TO_KCAL for i, v2 in enumerate(v1)] for v1 in values]
		else:
			sys.stderr.write("ERROR: undefined unit.\n")
			sys.exit(1)


	def output_IFIE_format(self, fragment_number, column_list, unit="a.u."):
		"""
		IFIE の結果を指定されたカラムの値で出力するメソッド

		Args:
			fragment_number (int): fragment number
			column_list (list): カラム名のリスト
			unit (str): "a.u." or "kcal/mol" (Default: "a.u.")

		Returns:
			list: [[fragment_pair_index(int), energy, ...], ...]
		"""
		list_energy = self.extract_IFIE_energy(fragment_number, unit)
		list_idx = [IFIE_FORMAT[self._version].index(v) for v in column_list]
		return [[v[0]] + [v[i + 1] for i in list_idx] for v in list_energy]


	def output_IFIE_piedalog_style(self, fragment_number):
		"""
		IFIE の結果をフォーマットして出力するメソッド (piedalog フォーマット)

		Args:
			fragment_number(int): 対象フラグメント

		Returns:
			str
		"""
		print("\t".join(["frag_Num", "Chain", "seq", "RES", "FCHARGE", "MAINSIDE", "DIST", "Total", "ES", "EX", "CT+mix", "DI(MP2)", "q(I=>J)"]))
		for values in self.extract_IFIE_energy(fragment_number, "kcal/mol"):
			fragment_number_pair = values[0]
			fragment_index_pair = self._fragment_number_list.index(fragment_number_pair)
			obj_fragment = self._obj_fragments[fragment_index_pair]

			energy_ES = values[IFIE_FORMAT[self._version].index("HF-ES") + 1]
			energy_EX = values[IFIE_FORMAT[self._version].index("PIEDA-EX") + 1]
			energy_CT = values[IFIE_FORMAT[self._version].index("PIEDA-CT") + 1]
			energy_DI = values[IFIE_FORMAT[self._version].index("MP2-IFIE") + 1]
			dq = values[IFIE_FORMAT[self._version].index("PIEDA-dq") + 1]

			print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}".format(
				fragment_number_pair,
				obj_fragment.chain_name,
				obj_fragment.residue_number,
				obj_fragment.residue_name,
				int(round(obj_fragment.charge, 0)),
				"MainSide",
				round(self.extract_distance(fragment_number, fragment_number_pair, "angstrom"), 3),
				round(energy_ES + energy_EX + energy_CT + energy_DI, 3),
				round(energy_ES, 3),
				round(energy_EX, 3),
				round(energy_CT, 3),
				round(energy_DI, 3),
				round(dq, 3)
			))


	def get_label(self, frag_idx=None):
		"""
		ラベルを返すメソッド (cpf2csv 用メソッド)

		Args:
			frag_idx (int, optional): フラグメントインデックス (Default: None)

		Returns:
			str: ラベル
		"""
		list_fragment = [obj_fragment.number for obj_fragment in self._obj_fragments]
		if frag_idx is not None:
			return list_fragment[frag_idx - 1]
		else:
			return list_fragment


	def get_fragment_atom(self, frag_idx=None):
		"""
		フラグメント構成原子を返すメソッド (cpf2csv 用メソッド)

		Args:
			frag_idx (int, optional): フラグメントインデックス (Default: None)

		Returns:
			list: 構成原子のリスト
		"""
		list_fragment = [obj_fragment.atoms for obj_fragment in self._obj_fragments]
		if frag_idx is not None:
			return list_fragment[frag_idx - 1]
		else:
			return list_fragment


	def get_atom_charge(self, atom_idx=None):
		"""
		原子電荷を返すメソッド (cpf2csv 用メソッド)

		Args:
			atom_idx (int, optional): 原子インデックス (Default: None)

		Returns:
			float or list: 原子電荷
		"""
		list_fragment = [[info[0], info[1].strip(), info[9]] for obj_fragment in self._obj_fragments for info in obj_fragment.structure_info]
		if atom_idx is not None:
			return list_fragment[atom_idx - 1]
		else:
			return list_fragment


	def get_fragment_charge(self, frag_idx=None):
		"""
		フラグメント電荷を返すメソッド (cpf2csv 用メソッド)

		Args:
			frag_idx (int, optional): フラグメントインデックス (Default: None)

		Returns:
			float or list: フラグメント電荷
		"""
		list_charge = [obj_fragment.charge for obj_fragment in self._obj_fragments]
		if frag_idx is not None:
			return list_charge[frag_idx - 1]
		else:
			return list_charge


	def get_energy(self, energy_type="Total", frag_idx=None, unit="kcal/mol"):
		"""
		IFIE エネルギーを返すメソッド (cpf2csv 用メソッド)

		Args:
			energy_type (str, optional): `Total`, `ES`, `EX`, `CT`, `DI` or `Q` (Default: "Total")
			frag_idx (list, optional): [frag_idx_A, frag_idx_B] (Default: None)
			unit (str): "a.u." or "kcal/mol" (Default: "kcal/mol")

		Returns:
			list
		"""
		f = 1
		if unit == "kcal/mol":
			f = AU_TO_KCAL

		energy = None
		if energy_type == "Total":
			list_idx = [IFIE_FORMAT[self._version].index(ENERGY_TYPE[energy_name]) for energy_name in ["ES", "EX", "CT", "DI"]]
			energy = np.array([[sum([obj_fragment1.get_IFIE(obj_fragment2)[idx_energy] for idx_energy in list_idx]) for obj_fragment2 in self._obj_fragments] for obj_fragment1 in self._obj_fragments])

		else:
			idx_energy = IFIE_FORMAT[self._version].index(ENERGY_TYPE[energy_type])
			energy = np.array([[obj_fragment1.get_IFIE(obj_fragment2)[idx_energy] for obj_fragment2 in self._obj_fragments] for obj_fragment1 in self._obj_fragments])
			if energy_type == "Q":
				f = 1

		if frag_idx is None:
			return np.round(energy * f, DIGIT)
		else:
			return np.round((energy * f)[frag_idx[0] - 1][frag_idx[1] - 1], DIGIT)


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
		else:
			result = self.get_energy(energy_type)

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
			output_range (list, optional): 出力するフラグメントラベルリスト (Default: None)

		Returns:
			list
		"""
		result = []
		list_frag = [obj_fragment.number for obj_fragment in self._obj_fragments]
		if output_range is not None:
			list_frag = output_range

		result_frag = [["Fragment index", "Fragment charge"]]
		result_frag += [[obj_fragment.number, obj_fragment.charge] for obj_fragment in self._obj_fragments if obj_fragment.number in list_frag]

		result_atom = [["Fragment index", "Atom index", "Atom", "Atomic charge"]]
		result_atom += [[frag_i, atom_i, atom_name, charge] for frag_i, atom_i, atom_name, charge in zip(
			self.get_structure_list("FragmentNumber"),
			self.get_structure_list("Index"),
			self.get_structure_list("Element"),
			self.get_structure_list("HF_MullikenCharge")
		) if frag_i in list_frag]

		if self.n_fragment > self.n_atom:
			result_atom += [["", "", "", ""] for _ in range(self.n_fragment - self.n_atom)]
		elif self.n_fragment < self.n_atom:
			result_frag += [["", ""] for _ in range(self.n_atom - self.n_fragment)]

		return [frag + [""] + atom for frag, atom in zip(result_frag, result_atom)]


class Fragment:
	def __init__(self, fragment_number:int):
		self._fragment_number = fragment_number
		self._structure_info = []
		self._electron = None
		self._bond_info = []
		self._neighbor_info = {}
		self._dipole_info = []
		self._monomer_info = []

		self._fragment_name = ""
		self._chain_name = ""
		self._residue_number = ""
		self._residue_name = ""
		self._charge = 0.0

		self._distances = {}
		self._info_IFIE = {}

	@property
	def number(self):
		return self._fragment_number

	@property
	def structure_info(self):
		return self._structure_info

	@property
	def atoms(self):
		return [v[0] for v in self._structure_info]

	@property
	def electron(self):
		return self._electron

	@property
	def bond(self):
		return self._bond_info

	@property
	def neighbor(self):
		return self._neighbor_info

	@property
	def dipole(self):
		return self._dipole_info

	@property
	def monomer(self):
		return self._monomer_info

	@property
	def name(self):
		self.__determine_property()
		return self._fragment_name

	@property
	def chain_name(self):
		self.__determine_property()
		return self._chain_name

	@property
	def residue_number(self):
		self.__determine_property()
		return self._residue_number

	@property
	def residue_name(self):
		self.__determine_property()
		return self._residue_name

	@property
	def charge(self):
		self.__determine_property()
		return self._charge


	def __determine_property(self):
		"""
		複数の原子情報から特定の名前を決定するメソッド

		Returns:
			None
		"""
		if self._fragment_name == "":
			list_fragment_name = ["{0}{1}".format(v[3], v[4]) for v in self._structure_info]
			self._fragment_name = sorted(collections.Counter(list_fragment_name).items(), key = lambda x : x[1], reverse = True)[0][0]

			list_chain_name = [v[15] for v in self._structure_info]
			self._chain_name = sorted(collections.Counter(list_chain_name).items(), key = lambda x : x[1], reverse = True)[0][0]

			list_residue_name = [v[3] for v in self._structure_info]
			self._residue_name = sorted(collections.Counter(list_residue_name).items(), key = lambda x : x[1], reverse = True)[0][0]

			list_residue_number = [v[4] for v in self._structure_info]
			self._residue_number = sorted(collections.Counter(list_residue_number).items(), key = lambda x : x[1], reverse = True)[0][0]

			self._charge = sum([v[9] for v in self._structure_info])


	def append_atom(self, structure_info:list):
		"""
		構造の原子情報を追加するメソッド

		Args:
			structure_info (list): 構造情報リスト

		Returns:
			self
		"""
		self._structure_info.append(structure_info)
		return self


	def set_electron(self, electron:int):
		"""
		電子数を設定するメソッド

		Args:
			electron (int): 電子数

		Returns:
			self
		"""
		self._electron = electron
		return self


	def set_bond(self, bonds:int):
		"""
		接続数を設定するメソッド

		Args:
			bonds (int): 接続数

		Returns:
			self
		"""
		self._bond_info = bonds
		return self


	def set_neighbor(self, neighbors:dict):
		"""
		接続原子インデックスリストを設定するメソッド

		Args:
			neighbors (dict): {fragment_number(int): atom_number(int), ...}

		Returns:
			self
		"""
		self._neighbor_info = neighbors
		return self


	def append_neighbor(self, fragment_number:int, atom_number:int):
		"""
		接続原子インデックス情報を追加するメソッド

		Args:
			fragment_number (int): 接続先のフラグメント番号
			atom_number (int): 接続している原子インデックス

		Returns:
			self
		"""
		self._neighbor_info[fragment_number] = atom_number
		return self


	def set_dipole_info(self, dipole:list):
		"""
		双極子モーメント情報を設定するメソッド

		Args:
			dipole (list): 双極子モーメント情報

		Returns:
			self
		"""
		self._dipole = dipole
		return self


	def set_monomer_info(self, monomer_info:list):
		"""
		モノマー情報を設定するメソッド

		Args:
			monomer_info (list): モノマー情報

		Returns:
			self
		"""
		self._monomer_info = monomer_info
		return self


	def has_atom(self, atom_number):
		"""
		フラグメントに指定された原子番号の原子を含むかどうかを返すメソッド

		Args:
			atom_number (int): 原子インデックス

		Returns:
			bool: True: 原子を含む / False: 原子を含まない
		"""
		return atom_number in self.atoms


	def append_distance(self, obj_Fragment, distance):
		"""
		距離情報を追加するメソッド

		Args:
			obj_Fragment (obj_Fragment): 他の Fragment オブジェクト
			distance (float): 距離

		Returns:
			self
		"""
		self._distances[obj_Fragment] = distance
		return self


	def append_IFIE(self, obj_Fragment, values):
		"""
		IFIE 情報を追加するメソッド

		Args:
			obj_Fragment (obj_Fragment): 他の Fragment オブジェクト
			values (list): IFIE 情報

		Returns:
			self
		"""
		self._info_IFIE[obj_Fragment] = values
		return self


	def get_distance(self, obj_Fragment, unit="bohr"):
		"""
		距離を返すメソッド

		Args:
			obj_Fragment (obj_Fragment): 他の Fragment オブジェクト
			unit (str): "bohr" or "angstrom"

		Returns:
			float: 距離
		"""
		if unit == "angstrom":
			return BOHR_RADIUS * self._distances[obj_Fragment]
		else:
			return self._distances[obj_Fragment]


	def get_IFIE(self, obj_Fragment, no_data="zero", raw_data=False):
		"""
		IFIE 情報を返すメソッド

		Args:
			obj_Fragment (obj_Fragment): 他の Fragment オブジェクト
			raw_data (bool): 接続フラグメントの場合、エネルギーの生データにするか (Default: False)
			no_data (str): IFIE データがない場合の値 (`zero` (zero-padding data) or `none` (None)) (Default: `zero`)

		Returns:
			list: IFIE 情報
		"""
		# 接続フラグメントでない場合
		if obj_Fragment in self._info_IFIE.keys():
			# IFIE 情報がある場合
			if self.get_distance(obj_Fragment) == 0:
				# 接続フラグメントの場合
				if raw_data:
					# 生データを返す
					return self._info_IFIE[obj_Fragment]
				else:
					# ゼロ埋めデータを返す
					key = list(self._info_IFIE.keys())[0]
					return [0.0 for _ in range(len(self._info_IFIE[key]))]

			else:
				# 問題のないデータ
				return self._info_IFIE[obj_Fragment]

		if no_data.lower() == "none":
			# IFIE データがない場合で、None を返す
			return None
		elif no_data.lower() == "zero":
			# IFIE データがない場合で、ゼロ埋めデータを返す
			key = list(self._info_IFIE.keys())[0]
			return [0.0 for _ in range(len(self._info_IFIE[key]))]
		else:
			sys.stderr.write("ERROR: undefined `no_data` value.\n")
			sys.exit(1)



# =============== main =============== #
# if __name__ == '__main__':
	# PROGRAM_ROOT = os.path.dirname(os.path.realpath(sys.argv[0]))
	# PROGRAM_NAME = os.path.splitext(os.path.basename(os.path.realpath(sys.argv[0])))[0]
