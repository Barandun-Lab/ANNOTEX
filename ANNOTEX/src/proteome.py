#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 17:35:54 2022

@author: dennis
"""

from glob import glob
from os import path
from chimerax.core.commands import run
import pickle
from Qt.QtWidgets import QTableWidgetItem
from Qt.QtGui import QColor

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, 4)

def load_obj(name):
    if path.exists(name):
        try:
            with open(name, 'rb') as f:
                return pickle.load(f)
        except:
            print('Could not load annotations')

def read_pdb(filename):
    atoms = []
    with open(filename) as file:
        for line in file.readlines():
            line = line.strip()
            line = line.split()
            if len(line) < 12:
                if line[6].count('.') == 3:
                    s = line[6]
                    i1 = s.index('.') + 4
                    i2 = s.index('.', i1) + 4
                    x = s[:i1]
                    y = s[i1:i2]
                    z = s[i2:]
                    del line[6]
                    line[6:6] = [x, y, z]
                elif line[6].count('.') == 2:
                    s = line[6]
                    i1 = s.index('.') + 4
                    x = s[:i1]
                    y = s[i1:]
                    del line[6]
                    line[6:6] = [x, y]
                elif line[7].count('.') == 2:
                    s = line[6]
                    i1 = s.index('.') + 4
                    y = s[:i1]
                    z = s[i1:]
                    del line[7]
                    line[7:7] = [y, z]
            atoms.append(line)
    return atoms

def read_m8(file):
    out = dict()
    if not path.exists(file):
        return out
    with open(file) as f:
        for line in f.readlines():
            line = line.strip().split()
            n1 = line[0][:-4] # TODO add support for mmcif
            
            n2 = line[1]
            if n2.startswith('AF-'):
                i1 = n2.index('-') + 1
                i2 = n2.index('-', i1)
                n2 = n2[i1:i2]
            else:
                n2 = n2[:4] + n2[n2.index('_'):]
            
            e = float(line[-2])
            
            if n1 not in out:
                out[n1] = []
            if len(out[n1]) < 10:
                out[n1].append((n2, e))
    return out

def read_annotations(file):
    annotations = dict()
    with open(file) as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            p = line[0]
            a = line[1]
            annotations[p] = a
            if len(line) > 2:
                print(line)
    return annotations


def size(protein):
    x = [float(atom[6]) for atom in protein]
    # y = [float(atom[7]) for atom in protein]
    # z = [float(atom[8]) for atom in protein]
    return max(x) - min(x)#, max(y) - min(y), max(z) - min(z)

class Proteome:

    def __init__(self, parent):
        # Expects the ScreeningTool as parent
        self.microsporidians = {'Nosema', 'Encephalitozoon', 'Mitosporidium',
                                'Ordospora', 'Nematocida', 'Edhazardia',
                                'Pseudoloma', 'Amphiamblys', 'Vavraia',
                                'Vittaforma', 'Enterocytozoon', 'Trachipleistophora',
                                'Spraguea', 'Anncaliia', 'Hepatospora',
                                'Enterospora', 'Tubulinosema', 'Hamiltosporidium',
                                'Thelohania', 'Microsporidium', 'Papilio xuthus',
                                'Dictyocoela'}
        self.parent = parent
        self.protein_list = []
        self.all = False
        self.current_domains = 0
        self.annotations = dict()
        self.opened = 0
        self.afdb = dict()
        self.pdb = dict()
        self.afsp = dict()
        self.diamond = dict()
        self.afdb_descs = dict()
        self.pdb_descs = dict()
        self.base_opened = 0
        # self.data_folder = ''
    
    def read_proteins(self):
        self._update_inputs()
        self.deeptmhmm = load_obj(path.join(self.data_folder, 'deeptmhmm.pkl'))
        self.protein_list = []
        for file in glob(path.join(self.data_folder, 'pdbs', '*.pdb')): # TODO Add support for mmcif
            name = file.split('/')[-1][:-4]
            self.protein_list.append(name)
        self.protein_list.sort(key=lambda x: int(x[1:x.index('.')]))
        self.parent.result_list.clear()
        for protein in self.protein_list:
            self.parent.result_list.addItem(protein)
        self.open_protein()
    
    def open_protein(self):
        run(self.parent.session, 'close')
        self.opened = 0
        self._update_inputs()
        protein = self.parent.result_list.currentText()
        file_path = path.join(self.data_folder, 'pdbs', f'{protein}.pdb')
        run(self.parent.session, f'open {file_path}')
        self.opened += 1
        
        i = 1
        self.parent.domain_list.addItem('All domains')
        while path.exists(path.join(self.data_folder, 'domains', f'{protein}_domain_{i}.pdb')):
            file_path = path.join(self.data_folder, 'domains', f'{protein}_domain_{i}.pdb')
            run(self.parent.session, f'open {file_path}')
            self.parent.domain_list.addItem(f'Domain {i}')
            i += 1
            self.opened += 1
        self.current_domains = i - 1
        
        self.color()
        self._update_deeptmhmm(protein)
        self._update_eggnog(protein)
        self._read_annotation(protein)
        self._read_foldseek(protein)
        self.update_afdb_table(protein)
        self.update_blast_table(protein)
    
    def color(self):
        run(self.parent.session, 'color bfactor #1 range 50,90 palette ^rainbow')
        # self.parent.msa_label.setText('MSA Coverage: -')
    
    def open_match(self, tableItem):
        file = tableItem.text().split()[-1] # Expects '{e} {desc}' format
        file = file[1:-1]
        if self.opened > self.current_domains + 1:
            run(self.parent.session, f'close #{self.opened}')
            self.opened -= 1
        if '_' in file:
            pdb_code = file[:4]
            chain = file[file.index('_')+1:]
            self.opened += 1
            run(self.parent.session, f'open {pdb_code}; mm #{self.opened}/{chain} to #1')
            run(self.parent.session, f'hide #{self.opened} atoms; show #{self.opened} cartoons')
            run(self.parent.session, f'color #{self.opened}/{chain} green')
        else:
            file = f'AF-{file}-F1-model_v2.pdb'
            file_path = path.join(self.data_folder, 'afdb_needed', file)
            if path.exists(file_path):
                run(self.parent.session, f'open {file_path}')
                self.opened += 1
                run(self.parent.session, f'mm #{self.opened} to #1')
            else:
                print('Tried to open match but file does not exist')
    
    def update_afdb_table(self, protein):
        afdb, afsp, pdb, to_add = [], [], [], []
        afdb_seen, afsp_seen = set(), set()
        if protein in self.afdb:
            afdb = self.afdb[protein]
            for x in afdb:
                afdb_seen.add(x[0])
        if protein in self.afsp:
            afsp = self.afsp[protein]
            for x in afsp:
                afsp_seen.add(x[0])
        if protein in self.pdb:
            pdb = self.pdb[protein]
        to_add = afdb + afsp + pdb
        self.parent.foldseek_table.setRowCount(len(to_add))
        to_add.sort(key=lambda x: x[1])
        for i, x in enumerate(to_add):
            if x[0] in self.afdb_descs:
                s = f'{x[1]} {self.afdb_descs[x[0]]}'
            elif x[0] in self.pdb_descs:
                s = f'{x[1]} {self.pdb_descs[x[0]]}'
            # else:
                # s = f'{x[1]} {x[0]}'
            self.parent.foldseek_table.setItem(i, 0, QTableWidgetItem(s))
            if x[0] in afdb_seen:
                self.parent.foldseek_table.item(i, 0).setBackground(QColor(155, 0, 0))
            elif x[0] in afsp_seen:
                self.parent.foldseek_table.item(i, 0).setBackground(QColor(100, 100, 0))
            else:
                self.parent.foldseek_table.item(i, 0).setBackground(QColor(0, 0, 155))
    
    def update_blast_table(self, protein):
        blast = []
        if protein in self.diamond:
            blast = self.diamond[protein]
        self.parent.blast_table.setRowCount(len(blast))
        for row, x in enumerate(blast):
            name, desc, e = x
            s = f'{e} {desc} {name}'
            self.parent.blast_table.setItem(row, 0, QTableWidgetItem(s))
            microsporidian = False
            if '[' and ']' in desc:
                species = desc[desc.rindex('['):-1]
                for x in self.microsporidians:
                    if x in species:
                        microsporidian = True
            if microsporidian:
                self.parent.blast_table.item(row, 0).setBackground(QColor(155, 0, 0))
    
    def _read_annotation(self, protein):
        if protein in self.annotations:
            self.parent.annotation.setText(self.annotations[protein])
        else:
            self.parent.annotation.setText('N/A')
    
    def _read_foldseek(self, protein):
        pass
    
    def _clear_tables(self):
        self.parent.domain_list.clear()
        self.parent.foldseek_table.setRowCount(0)
        self.parent.blast_table.setRowCount(0)
    
    def save_annotations(self):
        protein = self.parent.result_list.currentText()
        annotation = self.parent.annotation.text()
        self.annotations[protein] = annotation
        annotation_list = list(self.annotations.items())
        annotation_list.sort(key=lambda x: int(x[0][1:x[0].index('.')]))

        with open(path.join(self.data_folder, 'annotations.tsv'), 'w') as f:
            for l in annotation_list:
                p, a = l
                f.write(f'{p}\t{a}\n')
    
    def _update_deeptmhmm(self, protein):
        self.parent.deeptmhmm_label.setText(f'DeepTMHMM: {self.deeptmhmm[protein]}')
    
    def _load_eggnog(self):
        eggnog_file = path.join(self.data_folder, 'eggnog.tsv')
        self.eggnog = dict()
        if path.exists(eggnog_file):
            with open(eggnog_file) as f:
                for line in f.readlines()[5:-3]:
                    line = line.strip().split('\t')
                    gene = line[0]
                    desc = line[7]
                    suggested_name = line[8]
                    e = float(line[2])
                    self.eggnog[gene] = {'desc': desc, 'name': suggested_name, 'e': e}
    
    def _update_eggnog(self, protein):
        if protein in self.eggnog:
            self.parent.eggnog_name_label.setText(f'Name: {self.eggnog[protein]["name"]}')
            self.parent.eggnog_desc_label.setText(f'Desc: {self.eggnog[protein]["desc"]}')
            self.parent.eggnog_e_label.setText(f'E-value: {self.eggnog[protein]["e"]}')
        else:
            self.parent.eggnog_name_label.setText('Name: N/A')
            self.parent.eggnog_desc_label.setText('Desc: N/A')
            self.parent.eggnog_e_label.setText('E-value: N/A')
        
    
    def next_protein(self):
        self.save_annotations()
        current_index = self.parent.result_list.currentIndex()
        items = self.parent.result_list.count()
        new_index = min(current_index + 1, items - 1)
        self.parent.result_list.setCurrentIndex(new_index)
        self._clear_tables()
        self.open_protein()
    
    def previous_protein(self):
        self.save_annotations()
        current_index = self.parent.result_list.currentIndex()
        new_index = max(current_index - 1, 0)
        self.parent.result_list.setCurrentIndex(new_index)
        self._clear_tables()
        self.open_protein()
    
    def _update_inputs(self):
        self.data_folder = self.parent.data_folder.text()
        
        self._load_eggnog()
        
        self.annotations = read_annotations(path.join(self.data_folder, 'annotations.tsv'))
        if not self.afdb:
            self.afdb = read_m8(path.join(self.data_folder, 'foldseek', 'afdb.m8'))
        if not self.afsp:
            self.afsp = read_m8(path.join(self.data_folder, 'foldseek', 'afsp.m8'))
        if not self.pdb:
            self.pdb = read_m8(path.join(self.data_folder, 'foldseek', 'pdb.m8'))
        if not self.diamond:
            self.diamond = load_obj(path.join(self.data_folder, 'diamond.pkl'))
        if not self.afdb_descs:
            self.afdb_descs = load_obj(path.join(self.data_folder, 'afdb_descs.pkl'))
        if not self.pdb_descs:
            self.pdb_descs = load_obj(path.join(self.data_folder, 'pdb_descs.pkl'))
