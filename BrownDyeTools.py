#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# Last modified: 2016-11-18 14:34:49
#
'''BrownDye plugin for Pymol

For documentation see: https://github.io

Author : Robert Konecny
Email: rok@ucsd.edu
Release date: November 2016
License: GNU General Public License v.3

Copyright 2016 Robert Konecny, NBCR

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
from __future__ import print_function
import os, subprocess
import sys, string, filecmp, shutil
import random
import time
#import tkSimpleDialog
#import tkMessageBox
import tkFileDialog
import Tkinter
import Pmw
from threading import Thread
from lxml import etree

DEBUG = 5

__version__ = "0.0.1"

PDB2PQR_PATH = ''
APBS_PATH = ''
BD_PATH = ''
LOGFILE = 'LOGFILE'
DEFAULT_CONTACTS_FILE = 'protein-protein-contacts-default.xml'

if DEBUG > 0:
    PDB2PQR_PATH = '/opt/pdb2pqr-linux-bin64-2.1.0/'
    APBS_PATH = '/export1/Builds/SandBox-2016.10.11-14.05/bin/'
    BD_PATH = '/home/rok/BrownianDynamics/browndye-2016.4.14/bin/'

if "PDB2PQR_PATH" in os.environ: PDB2PQR_PATH = os.environ["PDB2PQR_PATH"]
if "APBS_PATH" in os.environ: APBS_PATH = os.environ["APBS_PATH"]
if "BD_PATH" in os.environ: BD_PATH = os.environ["BD_PATH"]


def __init__(self):
    """BrownDye plugin for PyMol."""
    self.menuBar.addmenuitem('Plugin', 'command',
                             'BrownDye Plugin', label='BrownDye Plugin',
                             command=lambda s=self: BDPlugin(s))


class DummyPymol(object):
    """Dummy pymol class for testing purposes or when running standalone GUI."""
    class Cmd:
        def load(self, name, sel=''):
            pass

        def get_names(self, name):
            return ['mol1', 'mol2', 'map1', 'map2']

        def get_type(self, thing):
            if thing.startswith('mol'):
                return 'object:molecule'
            else:
                return 'object:map'
    cmd = Cmd()

try:
    import pymold
except ImportError:
    print("::: Pymol import failed - Pymol features not available!")
    pymol = DummyPymol()


class BDPlugin(object):
    """ The main BrowDye plugin class."""
    def __init__(self, app):
        self.parent = app.root
        self.createGUI()
        
    def createGUI(self):
        """The main GUI class - sets up all GUI elements."""
        self.dialog = Pmw.Dialog(self.parent, buttons=('Exit',),
                                 title='BrownDye Plugin for PyMOL',
                                 command=self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        self.projectDir = Tkinter.StringVar()
        self.projectDir.set(self.getProjectDir())
        self.pdb2pqr_path = Tkinter.StringVar()
        self.pdb2pqr_path.set(PDB2PQR_PATH)
        self.apbs_path = Tkinter.StringVar()
        self.apbs_path.set(APBS_PATH)
        self.bd_path = Tkinter.StringVar()
        self.bd_path.set(BD_PATH)
                
        # parameters used by pdb2pqr
        self.mol0 = Tkinter.StringVar()
        self.mol1 = Tkinter.StringVar()
        self.mol0.set('mol0.pdb') # for testing only
        self.mol1.set('mol1.pdb')
        self.mol0_object = Tkinter.StringVar()
        self.mol1_object = Tkinter.StringVar()
        self.mol0_object.set(None)
        self.mol1_object.set(None)
        self.pdb2pqr_opt = Tkinter.StringVar()
        self.pdb2pqr_opt.set('--apbs-input')
        self.pqr_assign_only = Tkinter.BooleanVar()
        self.pqr_assign_only.set(True)
        self.pqr_ff = Tkinter.StringVar()
        self.pqr_ff.set('parse')

        # APBS parameters and defaults
        self.dime0 = [Tkinter.IntVar() for _ in range(3)]
        [self.dime0[x].set(0) for x in range(3)]
        self.cglen0 = [Tkinter.DoubleVar() for _ in range(3)]
        [self.cglen0[x].set(0.0) for x in range(3)]
        self.fglen0 = [Tkinter.DoubleVar() for _ in range(3)]
        [self.fglen0[x].set(0.0) for x in range(3)]
        self.dime1 = [Tkinter.IntVar() for _ in range(3)]
        [self.dime1[x].set(0) for x in range(3)]
        self.cglen1 = [Tkinter.DoubleVar() for _ in range(3)]
        [self.cglen1[x].set(0.0) for x in range(3)]
        self.fglen1 = [Tkinter.DoubleVar() for _ in range(3)]
        [self.fglen1[x].set(0.0) for x in range(3)]

        self.apbs_mode = Tkinter.StringVar()
        self.apbs_mode.set('lpbe')
        self.bcfl = Tkinter.StringVar()
        self.bcfl.set('sdh')
        self.ion_charge = [Tkinter.IntVar() for _ in range(2)]
        self.ion_conc = [Tkinter.DoubleVar() for _ in range(2)]
        self.ion_rad = [Tkinter.DoubleVar() for _ in range(2)]
        self.ion_charge[0].set(1)
        self.ion_charge[1].set(-1)
        [self.ion_conc[x].set(0.15) for x in range(2)]
        [self.ion_rad[x].set(1.0) for x in range(2)]

        self.interior_dielectric = Tkinter.DoubleVar()
        self.interior_dielectric.set(4.0)
        self.solvent_dielectric = Tkinter.DoubleVar()
        self.solvent_dielectric.set(78.0)
        self.chgm = Tkinter.StringVar()
        self.chgm.set('spl2')
        self.sdens = Tkinter.DoubleVar()
        self.sdens.set(10.0)
        self.swin = Tkinter.DoubleVar()
        self.swin.set(0.3)
        self.srfm = Tkinter.StringVar()
        self.srfm.set('smol')
        self.srad = Tkinter.DoubleVar()
        self.srad.set(1.4)
        self.system_temp = Tkinter.DoubleVar()
        self.system_temp.set(298.15)

        # reaction criteria
        self.contacts_f = Tkinter.StringVar()
        self.contacts_f.set('protein-protein-contacts.xml') # FIXME
        self.default_contacts_f = Tkinter.BooleanVar()
        self.default_contacts_f.set(False)
        self.rxn_distance = Tkinter.DoubleVar()
        self.rxn_distance.set(5.0)
        self.npairs = Tkinter.IntVar()
        self.npairs.set(3)
        
        # BD parameters and defaults
        self.solvent_eps = Tkinter.DoubleVar()
        self.solvent_eps.set(self.solvent_dielectric.get())
        self.mol0_eps = Tkinter.DoubleVar()
        self.mol0_eps.set(self.interior_dielectric.get())
        self.mol1_eps = Tkinter.DoubleVar()
        self.mol1_eps.set(self.interior_dielectric.get())
        self.debyel = Tkinter.DoubleVar()
        self.debyel.set(7.86) #FIXME
        self.ntraj = Tkinter.IntVar()
        self.ntraj.set(100)
        self.nthreads = Tkinter.IntVar()
        self.nthreads.set(1)
        self.mindx = Tkinter.DoubleVar()
        self.mindx.set(0.2)
        self.ntrajo = Tkinter.IntVar()
        self.ntrajo.set(1)
        self.ncopies = Tkinter.IntVar()
        self.ncopies.set(200)
        self.nbincopies = Tkinter.IntVar()
        self.nbincopies.set(200)
        self.nsteps = Tkinter.IntVar()
        self.nsteps.set(10000)
        self.westeps = Tkinter.IntVar()
        self.westeps.set(10000)
        self.maxnsteps = Tkinter.IntVar()
        self.maxnsteps.set(100000)

        self.run_in_background = Tkinter.BooleanVar()
        self.run_in_background.set(False)
        
        # self.pdb_fn        = Tkinter.StringVar()
        # self.pymol_sel     = Tkinter.StringVar()
        # self.msms_bin      = Tkinter.StringVar()
        # self.pdb2xyzr_bin  = Tkinter.StringVar()
        # self.pdb2xyzrn_bin = Tkinter.StringVar()
        # self.tmp_dir       = Tkinter.StringVar()

        self.cleanup_saved_pymol_sel = Tkinter.BooleanVar()
        self.cleanup_saved_pymol_sel.set(True) # by default, clean up

        # Analysis
        self.traj_f = Tkinter.StringVar()
        self.traj_f.set('traj0.xml')
        self.traj_index_n = Tkinter.IntVar()
        self.traj_index_n.set(1)
        #######################################################################
        # Main code
        w = Tkinter.Label(self.dialog.interior(),
                          text=('\nBrownDye Plugin for PyMOL\n'
                                'Version %s, NBCR 2016\n\n'
                                'Plugin for setting up and running BrownDye '
                                'Browndian dynamics simulations.' % __version__),
                          background='black', foreground='white')
        w.pack(expand=1, fill='both', padx=10, pady=5)

        # make a notebook
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=10, pady=10)

        #####################
        # Tab: Configuration
        #####################
        page = self.notebook.add('Configuration')       
        self.notebook.tab('Configuration').focus_set()
        config = Tkinter.LabelFrame(page, text='Calculation configuration')
        config.pack(fill='both', expand=True, padx=10, pady=10)

        project_path_ent = Pmw.EntryField(config,
                                          label_text='Project directory: ',
                                          labelpos='wn',
                                          entry_textvariable=self.projectDir)
        project_path_b_but = Tkinter.Button(config, text='Browse ...',
                                            command=self.browseProjectDir)
        label = Tkinter.Label(config, text='or')
        project_path_but = Tkinter.Button(config, text='Create',
                                          command=self.createProjectDir)
        pdb2pqr_path_ent = Pmw.EntryField(config,
                                          label_text='Select PDB2PQR location: ',
                                          labelpos='wn',
                                          entry_textvariable=self.pdb2pqr_path)
        pdb2pqr_path_but = Tkinter.Button(config, text='Browse...',
                                          command=self.getPDB2PQRpath)
        apbs_path_ent = Pmw.EntryField(config,
                                       label_text='Select APBS_PATH location:',
                                       labelpos='w',
                                       entry_textvariable=self.apbs_path)
        apbs_path_but = Tkinter.Button(config, text='Browse...',
                                       command=self.getAPBSpath)
        bd_path_ent = Pmw.EntryField(config,
                                     label_text='Select BD_PATH location: ',
                                     labelpos='wn',
                                     entry_textvariable=self.bd_path)
        bd_path_but = Tkinter.Button(config, text='Browse...',
                                     command=self.getBDpath)

        # arrange widgets using grid
        project_path_ent.grid(sticky='we', row=1, column=0, padx=5, pady=5)
        project_path_b_but.grid(sticky='we', row=1, column=1, padx=5, pady=5)
        label.grid(sticky='we', row=1, column=2, padx=5, pady=10)
        project_path_but.grid(sticky='we', row=1, column=3, padx=5, pady=5)
        pdb2pqr_path_ent.grid(sticky='we', row=2, column=0, padx=5, pady=5)
        pdb2pqr_path_but.grid(sticky='we', row=2, column=1, padx=5, pady=5)
        apbs_path_ent.grid(   sticky='we', row=3, column=0, padx=5, pady=5)
        apbs_path_but.grid(   sticky='we', row=3, column=1, padx=5, pady=5)
        bd_path_ent.grid(     sticky='we', row=4, column=0, padx=5, pady=5)
        bd_path_but.grid(     sticky='we', row=4, column=1, padx=5, pady=5)
        config.columnconfigure(0, weight=8)
        config.columnconfigure(1, weight=2)

        #############################
        # Tab: PQR files preparation
        #############################
        page = self.notebook.add('PQR files')
        group_pqr = Tkinter.LabelFrame(page, text='PQR files')
        group_pqr.grid(sticky='eswn',row=0, column=0, columnspan=2, padx=10, pady=5)

        pdb_a_ent = Pmw.EntryField(group_pqr,
                                   label_text='Molecule 0 PDB file:', labelpos='wn',
                                   entry_textvariable=self.mol0)
        pdb_a_but = Tkinter.Button(group_pqr, text='Browse...',
                                   command=self.getPDBMol0)

        label0 = Tkinter.Label(group_pqr, text='or')
        sel_list = []
        self.dialog0 = Pmw.SelectionDialog(page,
                                           title='Molecule 0',
                                           buttons=('OK', 'Cancel'),
                                           defaultbutton='OK',
                                           scrolledlist_labelpos='n',
                                           label_text='Select molecule 0',
                                           scrolledlist_items=sel_list,
                                           command=self.selectMol0)
        self.dialog0.withdraw()
        select0_but = Tkinter.Button(group_pqr, text='Select molecule 0',
                                     command=self.dialog0Call)

        pymol_obj0_opt = Pmw.OptionMenu(group_pqr, labelpos='w',
                                        label_text='Select molecule 0: ',
                                        menubutton_textvariable=self.mol0_object,
                                        menubutton_width=7,
                                        items=(['None'] + pymol.cmd.get_names("all")))
        pdb_b_ent = Pmw.EntryField(group_pqr,
                                   label_text='Molecule 1 PDB file:', labelpos='wn',
                                   entry_textvariable=self.mol1)
        pdb_b_but = Tkinter.Button(group_pqr, text='Browse...',
                                   command=self.getPDBMol1)
        label1 = Tkinter.Label(group_pqr, text='or')
        self.dialog1 = Pmw.SelectionDialog(page,
                                           title='Molecule 1',
                                           buttons=('OK', 'Cancel'),
                                           defaultbutton='OK',
                                           scrolledlist_labelpos='n',
                                           label_text='Select molecule 1',
                                           scrolledlist_items=sel_list,
                                           command=self.selectMol1)
        self.dialog1.withdraw()
        select1_but = Tkinter.Button(group_pqr, text='Select molecule 1',
                                     command=self.dialog1Call)
        pymol_obj1_opt = Pmw.OptionMenu(group_pqr, labelpos='w',
                                        label_text='Select molecule 1: ',
                                        menubutton_textvariable=self.mol1_object,
                                        menubutton_width=7,
                                        items=(['None'] + pymol.cmd.get_names("all")))
        #pqr_opt = Pmw.EntryField(group_pqr,
        #                         label_text='pdb2pqr options:', labelpos='wn',
        #                         value=self.pdb2pqr_opt.get(),
        #                         entry_textvariable=self.pdb2pqr_opt,
        #                         entry_width=30)
        pqr_ff_opt = Pmw.OptionMenu(group_pqr, labelpos='w',
                                    label_text='Force field: ',
                                    menubutton_textvariable=self.pqr_ff,
                                    menubutton_width=7,
                                    items=['parse', 'charmm', 'amber'])
        pqr_an_but = Tkinter.Checkbutton(group_pqr,
                                         text=('Assign charge and radius only '
                                               '(no structure optimization)'),
                                         variable=self.pqr_assign_only,
                                         onvalue=True, offvalue=False)
        pqr_opt_but = Tkinter.Button(page, text='Create PQR files',
                                     command=self.pdb2pqr)

        pdb_a_ent.grid(sticky='we', row=0, column=0, padx=5, pady=1)
        pdb_a_but.grid(sticky='we', row=0, column=1, padx=5, pady=1)
        label0.grid(sticky='we', row=0, column=2, padx=5, pady=1)
        select0_but.grid(sticky='we', row=0, column=3, padx=5, pady=1)
        pymol_obj0_opt.grid(sticky='we', row=0, column=4, padx=5, pady=1)
        pdb_b_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        pdb_b_but.grid(sticky='we', row=1, column=1, padx=5, pady=1)
        label1.grid(sticky='we', row=1, column=2, padx=5, pady=1)
        select1_but.grid(sticky='we', row=1, column=3, padx=5, pady=1)
        pymol_obj1_opt.grid(sticky='we', row=1, column=4, padx=5, pady=1)
        pqr_ff_opt.grid(sticky='we', row=2, column=0, padx=5, pady=1)
        pqr_an_but.grid(sticky='w', row=3, column=0, padx=5, pady=1)
        pqr_opt_but.grid(sticky='we', row=4, column=0, padx=5, pady=1)

        ############
        # Tab: APBS
        ############
        page = self.notebook.add('APBS')
        group_grids = Tkinter.LabelFrame(page, text='Grid size')
        group_apbs = Tkinter.LabelFrame(page, text='APBS options')
        group_grids.grid(sticky='eswn', row=0, column=0, columnspan=3, padx=10, pady=5)
        group_apbs.grid(sticky='eswn', row=0, column=3, columnspan=2, padx=10, pady=5)
        page.columnconfigure(0, weight=2)
        page.columnconfigure(1, weight=1)

        label0 = Tkinter.Label(group_grids, text='Molecule 0')
        dime0_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                     label_text='dime: ',
                                     value=self.dime0[0].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime0[0],
                                     entry_width=5)
        dime1_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                     label_text='',
                                     value=self.dime0[1].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime0[1],
                                     entry_width=5)
        dime2_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                     label_text='',
                                     value=self.dime0[2].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime0[2],
                                     entry_width=5)
        cglen0_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='cglen: ',
                                      value=self.cglen0[0].get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen0[0],
                                      entry_width=8)
        cglen1_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.cglen0[1].get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen0[1],
                                      entry_width=8)
        cglen2_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.cglen0[2].get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen0[2],
                                      entry_width=8)
        fglen0_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='fglen: ',
                                      value=self.fglen0[0].get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen0[0],
                                      entry_width=8)
        fglen1_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.fglen0[1].get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen0[1],
                                      entry_width=8)
        fglen2_0_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.fglen0[2].get(),
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen0[2],
                                      entry_width=8)
        get_size0_but = Tkinter.Button(group_grids,
                                       text="Calculate grid size",
                                       command=self.getSizemol0)

        label1 = Tkinter.Label(group_grids, text='Molecule 1')
        dime0_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                     label_text='dime: ',
                                     value=self.dime1[0].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime1[0],
                                     entry_width=5)
        dime1_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                     label_text='',
                                     value=self.dime1[1].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime1[1],
                                     entry_width=5)
        dime2_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                     label_text='',
                                     value=self.dime1[2].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime1[2],
                                     entry_width=5)

        cglen0_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='cglen: ',
                                      value=self.cglen1[0].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen1[0],
                                      entry_width=8)
        cglen1_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.cglen1[1].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen1[1],
                                      entry_width=8)
        cglen2_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.cglen1[2].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen1[2],
                                      entry_width=8)
        fglen0_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='fglen: ',
                                      value=self.fglen1[0].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen1[0],
                                      entry_width=8)
        fglen1_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.fglen1[1].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen1[1],
                                      entry_width=8)
        fglen2_1_ent = Pmw.EntryField(group_grids, labelpos='w',
                                      label_text='',
                                      value=self.fglen1[2].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen1[2],
                                      entry_width=8)

        get_size1_but = Tkinter.Button(group_grids,
                                       text="Calculate grid size",
                                       command=self.getSizemol1)

        label0.grid(sticky='we', row=0, column=0, padx=5, pady=10)
        dime0_0_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        dime1_0_ent.grid(sticky='we', row=1, column=1, padx=5, pady=1)
        dime2_0_ent.grid(sticky='we', row=1, column=2, padx=5, pady=1)
        cglen0_0_ent.grid(sticky='we', row=2, column=0, padx=5, pady=1)
        cglen1_0_ent.grid(sticky='we', row=2, column=1, padx=5, pady=1)
        cglen2_0_ent.grid(sticky='we', row=2, column=2, padx=5, pady=1)
        fglen0_0_ent.grid(sticky='we', row=3, column=0, padx=5, pady=1)
        fglen1_0_ent.grid(sticky='we', row=3, column=1, padx=5, pady=1)
        fglen2_0_ent.grid(sticky='we', row=3, column=2, padx=5, pady=1)
        get_size0_but.grid(sticky='we', row=4, column=2, padx=5, pady=1)

        label1.grid(sticky='we', row=5, column=0, padx=5, pady=10)
        dime0_1_ent.grid(sticky='we', row=6, column=0, padx=5, pady=1)
        dime1_1_ent.grid(sticky='we', row=6, column=1, padx=5, pady=1)
        dime2_1_ent.grid(sticky='we', row=6, column=2, padx=5, pady=1)
        cglen0_1_ent.grid(sticky='we', row=7, column=0, padx=5, pady=1)
        cglen1_1_ent.grid(sticky='we', row=7, column=1, padx=5, pady=1)
        cglen2_1_ent.grid(sticky='we', row=7, column=2, padx=5, pady=1)
        fglen0_1_ent.grid(sticky='we', row=8, column=0, padx=5, pady=1)
        fglen1_1_ent.grid(sticky='we', row=8, column=1, padx=5, pady=1)
        fglen2_1_ent.grid(sticky='we', row=8, column=2, padx=5, pady=1)

        get_size1_but.grid(sticky='we', row=9, column=2, padx=5, pady=1)

        group_grids.columnconfigure(0, weight=9)
        group_grids.columnconfigure(1, weight=1)
        
        apbs_mode_ent = Pmw.OptionMenu(group_apbs, labelpos='w',
                                       label_text='APBS mode: ',
                                       menubutton_textvariable=self.apbs_mode,
                                       menubutton_width=10,
                                       items=['lpbe', 'npbe'])
        solvent_die_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                         label_text='Solvent eps :',
                                         value=self.solvent_dielectric.get(),
                                         validate={'validator': 'real', 'min': 0.0},
                                         entry_textvariable=self.solvent_dielectric,
                                         entry_width=10)
        solute_die_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                        label_text='Molecule 0/1 eps: ',
                                        value=self.interior_dielectric.get(),
                                        validate={'validator': 'real', 'min': 0.0},
                                        entry_textvariable=self.interior_dielectric,
                                        entry_width=10)
        ion1_charge_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                         label_text='Ion(1) charge: ',
                                         value=self.ion_charge[0].get(),
                                         validate={'validator': 'integer', 'min': -2},
                                         entry_textvariable=self.ion_charge[0],
                                         entry_width=5)
        ion1_conc_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                       label_text='conc.: ',
                                       value=self.ion_conc[0].get(),
                                       validate={'validator': 'real', 'min': 0.0},
                                       entry_textvariable=self.ion_conc[0],
                                       entry_width=5)
        ion1_rad_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                      label_text='radius: ',
                                      value=self.ion_rad[0].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.ion_rad[0],
                                      entry_width=5)
        ion2_charge_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                         label_text='Ion(2) charge: ',
                                         value=self.ion_charge[1].get(),
                                         validate={'validator': 'integer', 'min': -2},
                                         entry_textvariable=self.ion_charge[1],
                                         entry_width=5)
        ion2_conc_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                       label_text='conc.: ',
                                       value=self.ion_conc[1].get(),
                                       validate={'validator': 'real', 'min': 0.0},
                                       entry_textvariable=self.ion_conc[1],
                                       entry_width=5)
        ion2_rad_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                      label_text='radius: ',
                                      value=self.ion_rad[1].get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.ion_rad[1],
                                      entry_width=5)

        sdens_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                   label_text='Surf. sphere density: ',
                                   value=self.sdens.get(),
                                   validate={'validator': 'real', 'min': 0.0},
                                   entry_textvariable=self.sdens, entry_width=5)
        srad_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                  label_text='Solvent radius: ',
                                  value=self.srad.get(),
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.srad, entry_width=5)
        swin_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                  label_text='Spline window: ',
                                  value=self.swin.get(),
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.swin, entry_width=5)
        temp_ent = Pmw.EntryField(group_apbs, labelpos='w',
                                  label_text='Temperature: ',
                                  value=self.system_temp.get(),
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.system_temp, entry_width=5)
        bcfl_ent = Pmw.OptionMenu(group_apbs, labelpos='w',
                                  label_text='Boundary condition: ',
                                  menubutton_textvariable=self.bcfl,
                                  menubutton_width=5,
                                  items=['sdh', 'zero', 'mdh', 'focus'])
        chgm_ent = Pmw.OptionMenu(group_apbs, labelpos='w',
                                  label_text='Charge mapping: ',
                                  menubutton_textvariable=self.chgm,
                                  menubutton_width=5,
                                  items=['spl2', 'spl0', 'spl4'])
        srfm_ent = Pmw.OptionMenu(group_apbs, labelpos='w',
                                  label_text='Diel. surf. calc. method: ',
                                  menubutton_textvariable=self.srfm,
                                  menubutton_width=5,
                                  items=['smol', 'mol', 'spl2', 'spl4'])

        run_apbs_but = Tkinter.Button(group_apbs,
                                      text="Run APBS to generate grids",
                                      command=self.runAPBS)

        solvent_die_ent.grid(sticky='we', row=0, column=0, padx=5, pady=1)
        solute_die_ent.grid( sticky='we', row=1, column=0, padx=5, pady=1)
        sdens_ent.grid(      sticky='we', row=2, column=0, padx=5, pady=1)
        srad_ent.grid(       sticky='we', row=3, column=0, padx=5, pady=1)
        swin_ent.grid(       sticky='we', row=4, column=0, padx=5, pady=1)
        temp_ent.grid(       sticky='we', row=5, column=0, padx=5, pady=1)

        ion1_charge_ent.grid(sticky='we', row=6, column=0, padx=5, pady=1)
        ion1_conc_ent.grid(  sticky='we', row=6, column=1, padx=5, pady=1)
        ion1_rad_ent.grid(   sticky='we', row=6, column=2, padx=5, pady=1)
        ion2_charge_ent.grid(sticky='we', row=7, column=0, padx=5, pady=1)
        ion2_conc_ent.grid(  sticky='we', row=7, column=1, padx=5, pady=1)
        ion2_rad_ent.grid(   sticky='we', row=7, column=2, padx=5, pady=1)

        apbs_mode_ent.grid(sticky='we', row=0, column=1, columnspan=2, padx=5, pady=1)
        bcfl_ent.grid(     sticky='we', row=1, column=1, columnspan=2, padx=5, pady=1)
        chgm_ent.grid(     sticky='we', row=2, column=1, columnspan=2, padx=5, pady=1)
        srfm_ent.grid(     sticky='we', row=3, column=1, columnspan=2, padx=5, pady=1)

        run_apbs_but.grid(sticky='we', row=8, column=0, columnspan=3, padx=5, pady=1)

        ###############################
        # Tab: Reaction criteria setup
        ###############################
        page = self.notebook.add('Reaction citeria')
        group_rxn = Tkinter.LabelFrame(page, text='Setup reaction criteria')
        group_rxn.grid(sticky='eswn', row=0, column=0, columnspan=2, padx=10, pady=5)

        contacts_ent = Pmw.EntryField(group_rxn,
                                      label_text='Contacts file:', labelpos='w',
                                      entry_textvariable=self.contacts_f, entry_width=50)
        contacts_but = Tkinter.Button(group_rxn, text='Browse...',
                                      command=self.getContacts)
        def_contacts_but = Tkinter.Checkbutton(group_rxn,
                                               text='or use default contacts file ',
                                               variable=self.default_contacts_f,
                                               onvalue=True, offvalue=False)
        rxn_distance_ent = Pmw.EntryField(group_rxn, labelpos='wn',
                                          label_text='Reaction distance: ',
                                          value=self.rxn_distance.get(),
                                          validate={'validator': 'real', 'min': 0.01},
                                          entry_textvariable=self.rxn_distance,
                                          entry_width=10)
        npairs_ent = Pmw.EntryField(group_rxn, labelpos='wn',
                                    label_text='Number of reaction pairs: ',
                                    value=self.npairs.get(),
                                    validate={'validator': 'integer', 'min': 1},
                                    entry_textvariable=self.npairs, entry_width=10)
        run_rxn_crit = Tkinter.Button(page, text="Creat rxn files",
                                      command=self.runRxnCrit)

        contacts_ent.grid(    sticky='we', row=0, column=0, padx=5, pady=1)
        contacts_but.grid(    sticky='e',  row=0, column=1, padx=5, pady=1)
        def_contacts_but.grid(sticky='e',  row=0, column=2, padx=5, pady=1)
        rxn_distance_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        npairs_ent.grid(      sticky='we', row=2, column=0, padx=5, pady=1)
        run_rxn_crit.grid(    sticky='we', row=3, column=0, padx=5, pady=1)
        
        ######################
        # Tab: BD input files
        ######################
        page = self.notebook.add('BD setup')
        group_bdinput = Tkinter.LabelFrame(page, text='BD input file')
        group_bdinput.grid(sticky='eswn', row=0, column=0, columnspan=2, padx=10, pady=5)

        solvent_eps_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                         label_text='Solvent eps: ',
                                         value=self.solvent_eps.get(),
                                         validate={'validator': 'real', 'min': 0.0},
                                         entry_textvariable=self.solvent_eps)
        debyel_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                    label_text='Debye length: ',
                                    value=self.debyel.get(),
                                    validate={'validator': 'real', 'min': 0.0},
                                    entry_textvariable=self.debyel)
        mol0_eps_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                      label_text='Molecule 0 eps: ',
                                      value=self.mol0_eps.get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.mol0_eps)
        mol1_eps_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                      label_text='Molecule 1 eps: ',
                                      value=self.mol1_eps.get(),
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.mol1_eps)
        ntraj_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                   label_text='Number of trajectories: ',
                                   value=self.ntraj.get(),
                                   validate={'validator': 'integer', 'min': 1},
                                   entry_textvariable=self.ntraj)
        nthreads_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                      label_text='Number of threads: ',
                                      value=self.nthreads.get(),
                                      validate={'validator': 'integer', 'min': 1},
                                      entry_textvariable=self.nthreads)
        mindx_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                   label_text='Time step tolerance: ',
                                   value=self.mindx.get(),
                                   validate={'validator': 'real', 'min': 0.01},
                                   entry_textvariable=self.mindx)
        ntrajo_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                    label_text='Number of trajectories per output: ',
                                    value=self.ntrajo.get(),
                                    validate={'validator': 'integer', 'min': 1},
                                    entry_textvariable=self.ntrajo)
        ncopies_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                     label_text='Number of copies: ',
                                     value=self.ncopies.get(),
                                     validate={'validator': 'integer', 'min': 1},
                                     entry_textvariable=self.ncopies)
        nbincopies_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                        label_text='Number of bin copies: ',
                                        value=self.nbincopies.get(),
                                        validate={'validator': 'integer', 'min': 1},
                                        entry_textvariable=self.nbincopies)
        nsteps_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                    label_text='Number of steps: ',
                                    value=self.nsteps.get(),
                                    validate={'validator': 'integer', 'min': 1},
                                    entry_textvariable=self.nsteps)
        westeps_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                     label_text='Number of WE steps per output: ',
                                     value=self.westeps.get(),
                                     validate={'validator': 'integer', 'min': 1},
                                     entry_textvariable=self.westeps)
        maxnsteps_ent = Pmw.EntryField(group_bdinput, labelpos='wn',
                                       label_text='Max number of steps: ',
                                       value=self.maxnsteps.get(),
                                       validate={'validator': 'integer', 'min': 1},
                                       entry_textvariable=self.maxnsteps)

        prep_bd_but = Tkinter.Button(page, text="Generate BD input file",
                                     command=self.prepBD)

        solvent_eps_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        debyel_ent.grid(     sticky='we', row=1, column=1, padx=5, pady=1)
        mol0_eps_ent.grid(   sticky='we', row=2, column=0, padx=5, pady=1)
        mol1_eps_ent.grid(   sticky='we', row=2, column=1, padx=5, pady=1)

        ntraj_ent.grid(      sticky='we', row=3, column=0, padx=5, pady=1)
        nthreads_ent.grid(   sticky='we', row=3, column=1, padx=5, pady=1)
        mindx_ent.grid(      sticky='we', row=4, column=0, padx=5, pady=1)
        ntrajo_ent.grid(     sticky='we', row=4, column=1, padx=5, pady=1)
        ncopies_ent.grid(    sticky='we', row=5, column=0, padx=5, pady=1)
        nbincopies_ent.grid( sticky='we', row=5, column=1, padx=5, pady=1)
        nsteps_ent.grid(     sticky='we', row=6, column=0, padx=5, pady=1)
        westeps_ent.grid(    sticky='we', row=6, column=1, padx=5, pady=1)
        maxnsteps_ent.grid(  sticky='we', row=7, column=0, padx=5, pady=1)

        prep_bd_but.grid(    sticky='we', row=8, column=0, padx=5, pady=1)

        ######################
        # Tab: BD Simulation
        ######################
        page = self.notebook.add('BD Simulation')
        group_sim = Tkinter.LabelFrame(page, text='BrownDye Simulation')
        group_sim.grid(sticky='eswn', row=0, column=0, columnspan=2, padx=10, pady=5)

        bkgj_cb = Tkinter.Checkbutton(group_sim,
                                      text='Run job in background', 
                                      variable=self.run_in_background,
                                      onvalue=True, offvalue=False)

        run_bd_but = Tkinter.Button(group_sim,
                                    text="Start BD simulation", command=self.runBD)
        kill_bd_but = Tkinter.Button(group_sim,
                                     text="Stop background job", command=self.killBD)
        self.messagebar1 = Pmw.MessageBar(group_sim,
                                          entry_width=10, entry_relief='sunken',
                                          labelpos='w',
                                          label_text='Trajectories:')
        self.messagebar2 = Pmw.MessageBar(group_sim,
                                          entry_width=20, entry_relief='sunken',
                                          labelpos='w',
                                          label_text='Completed events / escaped / stuck:')

        self.messagebar3 = Pmw.MessageBar(group_sim,
                                          entry_width=20, entry_relief='sunken',
                                          labelpos='w',
                                          label_text=('Calculated reaction rate '
                                                      'and probability:'))

        self.logtxt_ent = Pmw.ScrolledText(group_sim, labelpos='wn',
                                      borderframe=5, 
                                      vscrollmode='dynamic',
                                      hscrollmode='dynamic',
                                      text_width=100, text_height=15,
                                      text_wrap='none',
                                      text_background='#000000',
                                      text_foreground='white')

        bkgj_cb.grid(sticky='w', row=0, column=0, columnspan=2, padx=1, pady=1)
        run_bd_but.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        kill_bd_but.grid(sticky='we', row=1, column=1, padx=5, pady=1)

        self.messagebar1.grid(sticky='we', row=2, column=0, columnspan=1, padx=5, pady=1)
        self.messagebar2.grid(sticky='we', row=2, column=1, columnspan=1, padx=5, pady=1)
        self.messagebar3.grid(sticky='we', row=3, column=0, columnspan=2, padx=5, pady=1)
        self.logtxt_ent.grid(sticky='we', row=4, column=0, columnspan=2, padx=5, pady=1)

        ######################
        # Tab: Analysis
        ######################
        page = self.notebook.add('Analysis')
        group_analysis = Tkinter.LabelFrame(page, text='Analysis and Visualization')
        group_analysis.grid(sticky='eswn', row=0, column=0, columnspan=2, padx=10, pady=5)

        load_traj_ent = Pmw.EntryField(group_analysis,
                                       label_text='Select trajectory file:',
                                       labelpos='w',
                                       entry_textvariable=self.traj_f,
                                       entry_width=30)
        load_traj_but = Tkinter.Button(group_analysis, text='Browse...',
                                       command=self.loadTrajectoryFile)
        analyze_but = Tkinter.Button(group_analysis, text='Analyze',
                                       command=self.analyzeTrajectoryFile)

        self.message_ent = Pmw.ScrolledText(group_analysis, labelpos='wn',
                                            borderframe=5, 
                                            vscrollmode='dynamic',
                                            hscrollmode='dynamic',
                                            text_width=80, text_height=10,
                                            text_wrap='word',
                                            text_background='#000000',
                                            text_foreground='white')

        self.index_list = [None]
        self.dialog_idx = Pmw.SelectionDialog(page,
                                              title='Trajectory index',
                                              buttons=('OK', 'Cancel'),
                                              defaultbutton='OK',
                                              scrolledlist_labelpos='n',
                                              label_text='Select trajectory',
                                              scrolledlist_items=self.index_list,
                                              command=self.selectTrajIndex)
        self.dialog_idx.withdraw()
        select_index_but = Tkinter.Button(group_analysis, text='Select trajectory index',
                                          command=self.dialog_idx.activate)
        self.messagebar_idx = Pmw.MessageBar(group_analysis,
                                             entry_width=20, entry_relief='sunken',
                                             labelpos='w',
                                             label_text='Selected trajectory index:')

        traj_index_n_ent = Pmw.EntryField(group_analysis, labelpos='wn',
                                          label_text='Select trajectory index: ',
                                          value=self.traj_index_n.get(),
                                          validate={'validator': 'integer', 'min': 1},
                                          entry_textvariable=self.traj_index_n)
        convert_but = Tkinter.Button(group_analysis, text='Convert to xyz trajectory',
                                     command=self.convertTrajectoryToXYZ)
        load_xyztraj_but = Tkinter.Button(group_analysis, text='Load xyz trajectory',
                                          command=self.loadTrajectoryFileXYZ)

        load_traj_ent.grid(   sticky='we', row=0, column=0, padx=5, pady=1)
        load_traj_but.grid(   sticky='e',  row=0, column=1, padx=5, pady=1)
        analyze_but.grid(   sticky='e',  row=0, column=2, padx=5, pady=1)
        self.message_ent.grid(sticky='we', row=1, column=0, columnspan=3, padx=5, pady=1)
        select_index_but.grid(sticky='we',  row=2, column=0, padx=5, pady=1)
        self.messagebar_idx.grid(sticky='we', row=2, column=1, padx=5, pady=1)
        traj_index_n_ent.grid(sticky='we',  row=3, column=0, padx=5, pady=1)
        convert_but.grid(sticky='we', row=4, column=0, padx=5, pady=1)
        load_xyztraj_but.grid(sticky='we', row=5, column=0, padx=5, pady=1)
        
        #############
        # Tab: About
        #############
        page = self.notebook.add('About')
        group_about = Tkinter.LabelFrame(page, text='About BrownDye Plugin for PyMOL')
        group_about.grid(sticky='n', row=0, column=0, columnspan=2, padx=10, pady=5)
        about_plugin = (
            'This plugin provides a GUI for setting up and running Brownian '
            'dynamics simulations with BrownDye.\n\n'

            'The plugin requires PDB2PQR, APBS and BrownDye. '
            'To download and install these applications go to:\n\n'
            'http://www.poissonboltzmann.org/ \nand\n'
            'http://browndye.ucsd.edu/\n\n'

            'This software is released under the terms of GNU GPL3 license.\n'
            'For more details please see the accompanying documentation.\n\n'

            '(c) 2016 National Biomedical Computation Resource\n'
            'http://nbcr.ucsd.edu/')

        label_about = Tkinter.Label(group_about, text=about_plugin)
        label_about.grid(sticky='we', row=0, column=2, padx=5, pady=10)

        self.notebook.setnaturalsize()
        return

    def getProjectDir(self):
        """Generate a radnom project directory name."""
        rndID = random.randint(1000, 9999)
        cwd = os.getcwd()
        pDir = '%s/bd-project-%d' % (cwd, rndID)
        return pDir
    
    def browseProjectDir(self):
        """Browse for project directory and chdir there."""
        d = tkFileDialog.askdirectory(
            title='Project directory', initialdir='',
            parent=self.parent)
        self.projectDir.set(d)
        try:
            os.chdir(self.projectDir.get())
        except OSError:
            print("::: %s No such file or directory." % self.projectDir.get())
        return

    def createProjectDir(self):
        """Create project directory. """
        if os.path.exists(self.projectDir.get()):
            print("This directory already exists!")
            return
        os.makedirs(self.projectDir.get())
        os.chdir(self.projectDir.get())
        return

    def runCmd(self, command):
        """Generic wrapper for running arbitrary shell commands."""
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=True)
        p.wait()
        stdout, stderr = p.communicate()
        if DEBUG > 2: print("returncode = %d" % p.returncode)
        if DEBUG > 2: print("stdout:\n%s\n" % stdout)
        if DEBUG > 2: print("stderr:\n%s\n" % stderr)
        if p.returncode > 0:
            sys.stderr.write("::: Non-zero return code!\n") 
            sys.stderr.write("::: Failed command: \n\n")
            sys.stderr.write("%s \n\n" % command)
            sys.stderr.write(stderr)
            # sys.exit(1)
        return(p.returncode)

    def getPDB2PQRpath(self):
        """Get PDB2PQR binary path."""
        d = tkFileDialog.askdirectory(
            title='PDB2PQR binary directory', initialdir='',
            parent=self.parent)
        self.pdb2pqr_path.set(d)
        return
        
    def getAPBSpath(self):
        """Get APBS binary path."""
        d = tkFileDialog.askdirectory(
            title='APBS binary directory', initialdir='',
            parent=self.parent)
        self.apbs_path.set(d)
        return
        
    def getBDpath(self):
        """Get BrownDye binary path."""
        d = tkFileDialog.askdirectory(
            title='BrownDye binary directory', initialdir='',
            parent=self.parent)
        self.bd_path.set(d)
        return

    def getPDBMol0(self):
        """Get molecule 0 filename."""
        file_name = tkFileDialog.askopenfilename(
            title='PDB File', initialdir='',
            filetypes=[('pdb files', '*.pdb *.ent'), ('all files', '*')],
            parent=self.parent)
        self.mol0.set(file_name)
        return
        
    def getPDBMol1(self):
        """Get molecule 1 filename."""
        file_name = tkFileDialog.askopenfilename(
            title='PDB File', initialdir='',
            filetypes=[('pdb files', '*.pdb *.ent'), ('all files', '*')],
            parent=self.parent)
        self.mol1.set(file_name)
        return

    def selectMol0(self, result):
        sel = self.dialog0.getcurselection()
        if len(sel) > 0: print("::: Selection: %s" % sel)
        self.dialog0.deactivate(result)
        return

    def selectMol1(self, result):
        sel = self.dialog1.getcurselection()
        if len(sel) > 0: print("::: Selection: %s" % sel)
        self.dialog1.deactivate(result)
        return

    def dialog0Call(self):
        """Populate the selection list with pymol objects."""
        self.dialog0.delete(0, 'end')
        for i in pymol.cmd.get_names():
            self.dialog0.insert('end', i)
        self.dialog0.activate()
        return
    
    def dialog1Call(self):
        """Populate the selection list with pymol objects."""
        self.dialog1.delete(0, 'end')
        for i in pymol.cmd.get_names():
            self.dialog1.insert('end', i)
        self.dialog1.activate()
        return
    
    def getSizemol0(self):
        """Calculate APBS grid dimensions for molecule 0."""
        pqr_filename = 'mol0.pqr'
        if not os.path.isfile(pqr_filename):
            print("::: %s does not exist!" % pqr_filename)
            return
        psize = Psize()
        psize.runPsize(pqr_filename)
        #print(psize.getCharge())
        grid_points = psize.getFineGridPoints()
        cglen = psize.getCoarseGridDims()
        fglen = psize.getFineGridDims()
        [self.dime0[x].set(grid_points[x]) for x in range(3)]
        [self.cglen0[x].set(cglen[x]) for x in range(3)]
        [self.fglen0[x].set(fglen[x]) for x in range(3)]
        return

    def getSizemol1(self):
        """Calculate APBS grid dimensions for molecule 1."""
        pqr_filename = 'mol1.pqr'
        if not os.path.isfile(pqr_filename):
            print("::: %s does not exist!" % pqr_filename)
            return
        psize = Psize()
        psize.runPsize(pqr_filename)
        #print(psize.getCharge())
        grid_points = psize.getFineGridPoints()
        cglen = psize.getCoarseGridDims()
        fglen = psize.getFineGridDims()
        [self.dime1[x].set(grid_points[x]) for x in range(3)]
        [self.cglen1[x].set(cglen[x]) for x in range(3)]
        [self.fglen1[x].set(fglen[x]) for x in range(3)]
        return
    
    def getContacts(self):
        """Get contacts file."""
        file_name = tkFileDialog.askopenfilename(
            title='Contacts File', initialdir='',
            filetypes=[('xml files', '*.xml'), ('all files', '*')],
            parent=self.parent)
        self.contacts_f.set(file_name)
        return
    
    def pdb2pqr(self):
        """Convert PDB to PQR."""
        if self.mol0_object.get() == 'None':
            if not filecmp.cmp(self.mol0.get(), 'mol0.pdb'):
                try:
                    shutil.copyfile(self.mol0.get(), 'mol0.pdb')
                except:
                    print("::: Creating of mol0.pdb failed!")
                    return
        else:
            # save pymol mol0_object to pdb file
            self.mol0.set('mol0.pdb')
        if self.mol1_object.get() == 'None':
            if not filecmp.cmp(self.mol1.get(), 'mol1.pdb'):
                try:
                    shutil.copyfile(self.mol1.get(), 'mol1.pdb')
                except:
                    print("::: Creating of mol1.pdb failed!")
                    return
        else:
            # save pymol mol1_object to pdb file
            self.mol1.set('mol1.pdb')

        assign_only = ''
        if self.pqr_assign_only.get(): assign_only = '--assign-only'
        pqr_options = ('%s %s --ff=%s' %
                       (assign_only, self.pdb2pqr_opt.get(), self.pqr_ff.get()))
        for i in ['mol0', 'mol1']:
            command = ('%s/pdb2pqr %s %s.pdb %s.pqr' %
                       (self.pdb2pqr_path.get(), pqr_options, i, i))
            if DEBUG > 2: print(command)
            print("::: Running pdb2pqr on %s ..." % i)
            rc = self.runCmd(command)
            if rc == 0:
                print("::: Done.")
            else:
                print("::: Command failed: %s" % command)

        return

    def runAPBS(self):
        """Run APBS calculations on molecule 0 and molecule 1"""
        apbs_template = """# APBS template for BD grids
read
    mol pqr %s
end
elec 
    mg-auto
    dime %d %d %d
    cglen %f %f %f
    fglen %f %f %f
    cgcent mol 1
    fgcent mol 1
    mol 1
    %s # lpbe/npbe
    bcfl %s # sdh
    ion charge 1 conc %f radius %f
    ion charge -1 conc %f radius %f
    pdie %f
    sdie %f
    chgm %s # spl2
    sdens %f
    srfm %s # smol
    srad %f
    swin %f
    temp %f
    calcenergy total
    calcforce no
    write pot dx %s
end
print elecEnergy 1 end
quit
"""

        for i in ['mol0', 'mol1']:
            pqr_filename = '%s.pqr' % i
            # psize = Psize()
            # psize.runPsize(pqr_filename)
            # print(psize.getCharge())
            # grid_points = psize.getFineGridPoints()
            # cglen = psize.getCoarseGridDims()
            # fglen = psize.getFineGridDims()

            if i == 'mol0':
                grid_points = [self.dime0[x].get() for x in range(3)]
                cglen = [self.cglen0[x].get() for x in range(3)]
                fglen = [self.fglen0[x].get() for x in range(3)]
            if i == 'mol1':
                grid_points = [self.dime1[x].get() for x in range(3)]
                cglen = [self.cglen1[x].get() for x in range(3)]
                fglen = [self.fglen1[x].get() for x in range(3)]
            if grid_points[0] == 0 or grid_points[1] == 0 or grid_points[2] == 0:
                print("::: %s - no grid points defined!" % i)
                return
                
            dx_filename = i
            fnout = '%s.in' % i
            with open(fnout, "w") as fout:
                fout.write((apbs_template % 
                            (pqr_filename,
                             grid_points[0], grid_points[1], grid_points[2],
                             cglen[0], cglen[1], cglen[2],
                             fglen[0], fglen[1], fglen[2], 
                             self.apbs_mode.get(), self.bcfl.get(),
                             self.ion_conc[0].get(), self.ion_rad[0].get(),
                             self.ion_conc[1].get(), self.ion_rad[1].get(),
                             self.interior_dielectric.get(),
                             self.solvent_dielectric.get(),
                             self.chgm.get(), self.sdens.get(),
                             self.srfm.get(), self.srad.get(),
                             self.swin.get(), self.system_temp.get(),
                             dx_filename)))

            command = '%s/apbs %s.in' % (self.apbs_path.get(), i)
            if DEBUG > 2:
                print(command)
                print(grid_points)
                print(cglen)
                print(fglen)
            print("::: Running apbs on %s ... " % i)
            gmem = 200.0 * grid_points[0] * grid_points[1] * grid_points[2] / 1024 / 1024
            print("::: Estimated memory requirements: %.3f MB" % gmem)
            rc = self.runCmd(command)
            if rc == 0:
                print("::: Done.")
            else:
                print("::: Failed: %s" % command)

        return

    def runPqr2xml(self):
        """Run pqr2xml on mol0 and mol1 PQR files. """
        for i in ['mol0', 'mol1']:
            command = ('%s/pqr2xml < %s.pqr > %s-atoms.pqrxml'
                       % (self.bd_path.get(), i, i))
            if DEBUG > 2: print(command)
            print("::: Running pqr2xml on %s ..." % i)
            rc = self.runCmd(command)
            if rc == 0:
                print("::: Done.")
            else:
                print("::: Failed: %s" % command)
        return

    def createDefaultContactsFile(self):
        contacts_template = """
<!-- Default protein/protein contacts file -->
<contacts>
  <combinations>
    <molecule0>
      <contact>
        <atom> SG </atom> <residue> CYS </residue>
      </contact>
      <contact>
        <atom> NE2 </atom> <residue> HIS </residue>
      </contact>
      <contact>
        <atom> ND1 </atom> <residue> HIS </residue>
      </contact>
      <contact>
        <atom> NZ </atom> <residue> LYS </residue>
      </contact>
      <contact>
        <atom> ND2 </atom> <residue> ASN </residue>
      </contact>
      <contact>
        <atom> NE2 </atom> <residue> GLN </residue>
      </contact>
      <contact>
        <atom> NH1 </atom> <residue> ARG </residue>
      </contact>
      <contact>
        <atom> NZ </atom> <residue> ARG </residue>
      </contact>
      <contact>
        <atom> NH2 </atom> <residue> ARG </residue>
      </contact>
      <contact>
        <atom> OG </atom> <residue> SER </residue>
      </contact>
      <contact>
        <atom> OH </atom> <residue> TYR </residue>
      </contact>
      <contact>
        <atom> NE1 </atom> <residue> TRP </residue>
      </contact>
      <contact> <atom> N </atom> <residue> ALA </residue> </contact>
      <contact> <atom> N </atom> <residue> ARG </residue> </contact>
      <contact> <atom> N </atom> <residue> ASN </residue> </contact>
      <contact> <atom> N </atom> <residue> ASP </residue> </contact>
      <contact> <atom> N </atom> <residue> CYS </residue> </contact>
      <contact> <atom> N </atom> <residue> GLU </residue> </contact>
      <contact> <atom> N </atom> <residue> GLN </residue> </contact>
      <contact> <atom> N </atom> <residue> GLY </residue> </contact>
      <contact> <atom> N </atom> <residue> HIS </residue> </contact>
      <contact> <atom> N </atom> <residue> ILE </residue> </contact>
      <contact> <atom> N </atom> <residue> LEU </residue> </contact>
      <contact> <atom> N </atom> <residue> LYS </residue> </contact>
      <contact> <atom> N </atom> <residue> MET </residue> </contact>
      <contact> <atom> N </atom> <residue> PHE </residue> </contact>
      <contact> <atom> N </atom> <residue> PRO </residue> </contact>
      <contact> <atom> N </atom> <residue> SER </residue> </contact>
      <contact> <atom> N </atom> <residue> THR </residue> </contact>
      <contact> <atom> N </atom> <residue> TRP </residue> </contact>
      <contact> <atom> N </atom> <residue> TYR </residue> </contact>
      <contact> <atom> N </atom> <residue> VAL </residue> </contact>
    </molecule0>
    <molecule1>
      <contact>
        <atom> SG </atom> <residue> CYS </residue>
      </contact>
      <contact>
        <atom> OD1 </atom> <residue> ASP </residue>
      </contact>
      <contact>
        <atom> OD2 </atom> <residue> ASP </residue>
      </contact>
      <contact>
        <atom> OE1 </atom> <residue> GLU </residue>
      </contact>
      <contact>
        <atom> OE2 </atom> <residue> GLU </residue>
      </contact>
      <contact>
        <atom> ND1 </atom> <residue> HIS </residue>
      </contact>
      <contact>
        <atom> SD </atom> <residue> MET </residue>
      </contact>
      <contact>
        <atom> OD1 </atom> <residue> ASN </residue>
      </contact>
      <contact>
        <atom> OE1 </atom> <residue> GLN</residue>
      </contact>
      <contact>
        <atom> OG </atom> <residue> SER </residue>
      </contact>
      <contact>
        <atom> OG1 </atom> <residue> THR </residue>
      </contact>
      <contact>
        <atom> OH </atom> <residue> TYR </residue>
      </contact>
      <contact> <atom> O </atom> <residue> ALA </residue> </contact>
      <contact> <atom> O </atom> <residue> ARG </residue> </contact>
      <contact> <atom> O </atom> <residue> ASN </residue> </contact>
      <contact> <atom> O </atom> <residue> ASP </residue> </contact>
      <contact> <atom> O </atom> <residue> CYS </residue> </contact>
      <contact> <atom> O </atom> <residue> GLU </residue> </contact>
      <contact> <atom> O </atom> <residue> GLN </residue> </contact>
      <contact> <atom> O </atom> <residue> GLY </residue> </contact>
      <contact> <atom> O </atom> <residue> HIS </residue> </contact>
      <contact> <atom> O </atom> <residue> ILE </residue> </contact>
      <contact> <atom> O </atom> <residue> LEU </residue> </contact>
      <contact> <atom> O </atom> <residue> LYS </residue> </contact>
      <contact> <atom> O </atom> <residue> MET </residue> </contact>
      <contact> <atom> O </atom> <residue> PHE </residue> </contact>
      <contact> <atom> O </atom> <residue> PRO </residue> </contact>
      <contact> <atom> O </atom> <residue> SER </residue> </contact>
      <contact> <atom> O </atom> <residue> THR </residue> </contact>
      <contact> <atom> O </atom> <residue> TRP </residue> </contact>
      <contact> <atom> O </atom> <residue> TYR </residue> </contact>
      <contact> <atom> O </atom> <residue> VAL </residue> </contact>
    </molecule1>
  </combinations>
  <combinations>
    <molecule1>
      <contact>
        <atom> SG </atom> <residue> CYS </residue>
      </contact>
      <contact>
        <atom> NE2 </atom> <residue> HIS </residue>
      </contact>
      <contact>
        <atom> ND1 </atom> <residue> HIS </residue>
      </contact>
      <contact>
        <atom> NZ </atom> <residue> LYS </residue>
      </contact>
      <contact>
        <atom> ND2 </atom> <residue> ASN </residue>
      </contact>
      <contact>
        <atom> NE2 </atom> <residue> GLN </residue>
      </contact>
      <contact>
        <atom> NH1 </atom> <residue> ARG </residue>
      </contact>
      <contact>
        <atom> NZ </atom> <residue> ARG </residue>
      </contact>
      <contact>
        <atom> NH2 </atom> <residue> ARG </residue>
      </contact>
      <contact>
        <atom> OG </atom> <residue> SER </residue>
      </contact>
      <contact>
        <atom> OH </atom> <residue> TYR </residue>
      </contact>
      <contact>
        <atom> NE1 </atom> <residue> TRP </residue>
      </contact>
      <contact> <atom> N </atom> <residue> ALA </residue> </contact>   
      <contact> <atom> N </atom> <residue> ARG </residue> </contact> 
      <contact> <atom> N </atom> <residue> ASN </residue> </contact>  
      <contact> <atom> N </atom> <residue> ASP </residue> </contact> 
      <contact> <atom> N </atom> <residue> CYS </residue> </contact>
      <contact> <atom> N </atom> <residue> GLU </residue> </contact>
      <contact> <atom> N </atom> <residue> GLN </residue> </contact>
      <contact> <atom> N </atom> <residue> GLY </residue> </contact>
      <contact> <atom> N </atom> <residue> HIS </residue> </contact>
      <contact> <atom> N </atom> <residue> ILE </residue> </contact>
      <contact> <atom> N </atom> <residue> LEU </residue> </contact>
      <contact> <atom> N </atom> <residue> LYS </residue> </contact>
      <contact> <atom> N </atom> <residue> MET </residue> </contact>
      <contact> <atom> N </atom> <residue> PHE </residue> </contact>
      <contact> <atom> N </atom> <residue> PRO </residue> </contact>
      <contact> <atom> N </atom> <residue> SER </residue> </contact>
      <contact> <atom> N </atom> <residue> THR </residue> </contact>
      <contact> <atom> N </atom> <residue> TRP </residue> </contact>
      <contact> <atom> N </atom> <residue> TYR </residue> </contact>
      <contact> <atom> N </atom> <residue> VAL </residue> </contact>
    </molecule1>
    <molecule0>
      <contact>
        <atom> SG </atom> <residue> CYS </residue>
      </contact>
      <contact>
        <atom> OD1 </atom> <residue> ASP </residue>
      </contact>
      <contact>
        <atom> OD2 </atom> <residue> ASP </residue>
      </contact>
      <contact>
        <atom> OE1 </atom> <residue> GLU </residue>
      </contact>
      <contact>
        <atom> OE2 </atom> <residue> GLU </residue>
      </contact>
      <contact>
        <atom> ND1 </atom> <residue> HIS </residue>
      </contact>
      <contact>
        <atom> SD </atom> <residue> MET </residue>
      </contact>
      <contact>
        <atom> OD1 </atom> <residue> ASN </residue>
      </contact>
      <contact>
        <atom> OE1 </atom> <residue> GLN</residue>
      </contact>
      <contact>
        <atom> OG </atom> <residue> SER </residue>
      </contact>
      <contact>
        <atom> OG1 </atom> <residue> THR </residue>
      </contact>
      <contact>
        <atom> OH </atom> <residue> TYR </residue>
      </contact>
      <contact> <atom> O </atom> <residue> ALA </residue> </contact>   
      <contact> <atom> O </atom> <residue> ARG </residue> </contact> 
      <contact> <atom> O </atom> <residue> ASN </residue> </contact>  
      <contact> <atom> O </atom> <residue> ASP </residue> </contact> 
      <contact> <atom> O </atom> <residue> CYS </residue> </contact>
      <contact> <atom> O </atom> <residue> GLU </residue> </contact>
      <contact> <atom> O </atom> <residue> GLN </residue> </contact>
      <contact> <atom> O </atom> <residue> GLY </residue> </contact>
      <contact> <atom> O </atom> <residue> HIS </residue> </contact>
      <contact> <atom> O </atom> <residue> ILE </residue> </contact>
      <contact> <atom> O </atom> <residue> LEU </residue> </contact>
      <contact> <atom> O </atom> <residue> LYS </residue> </contact>
      <contact> <atom> O </atom> <residue> MET </residue> </contact>
      <contact> <atom> O </atom> <residue> PHE </residue> </contact>
      <contact> <atom> O </atom> <residue> PRO </residue> </contact>
      <contact> <atom> O </atom> <residue> SER </residue> </contact>
      <contact> <atom> O </atom> <residue> THR </residue> </contact>
      <contact> <atom> O </atom> <residue> TRP </residue> </contact>
      <contact> <atom> O </atom> <residue> TYR </residue> </contact>
      <contact> <atom> O </atom> <residue> VAL </residue> </contact>
    </molecule0>
  </combinations>
</contacts>
        """
        fout = DEFAULT_CONTACTS_FILE
        with open(fout, 'w') as f:
            f.write(contacts_template)
        return

    def makeRxnCriteria(self):
        """Create rxn criteria files."""
        if self.default_contacts_f.get():
            self.createDefaultContactsFile()
            self.contacts_f.set(DEFAULT_CONTACTS_FILE)
        if not os.path.isfile(self.contacts_f.get()):
            print("::: File not found: %s" % self.contacts_f.get())
            return
        command = ('%s/make_rxn_pairs '
                   '-nonred -mol0 mol0-atoms.pqrxml -mol1 mol1-atoms.pqrxml '
                   '-ctypes %s  -dist %f > mol0-mol1-rxn-pairs.xml'
                   % (self.bd_path.get(), self.contacts_f.get(),
                      self.rxn_distance.get()))
        if DEBUG > 2: print(command)
        print("::: Running make_rxn_pairs ...")
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
        command =('%s/make_rxn_file '
                  '-pairs mol0-mol1-rxn-pairs.xml -distance %f '
                  ' -nneeded %d > mol0-mol1-rxns.xml'
                  % (self.bd_path.get(), self.rxn_distance.get(),
                     self.npairs.get()))
        if DEBUG > 2: print(command)
        print("::: Running make_rxn_file ...")
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
        return
            
    def prepareInputFile(self):
        """Create BD input file."""
        nam_simulation_template = """
<root>
 <protein> true </protein>
 <solvent>
    <dielectric> %f </dielectric>
    <debye-length>  %f </debye-length>
  </solvent>
  <output> results.xml </output>
  <n-trajectories> %d </n-trajectories>
  <n-threads> %d </n-threads>
  <molecule0>
    <prefix> mol0 </prefix>
    <atoms>  mol0-atoms.pqrxml </atoms>
    <apbs-grids>
       <grid> mol0.dx </grid>
    </apbs-grids>
    <solute-dielectric> %f </solute-dielectric>
  </molecule0>
  <molecule1>
    <prefix> mol1 </prefix>
    <atoms>  mol1-atoms.pqrxml </atoms>
    <all-in-surface> false </all-in-surface>
    <apbs-grids>
       <grid> mol1.dx </grid>
    </apbs-grids>
    <solute-dielectric> %f </solute-dielectric>
  </molecule1>
  <include-desolvation-forces> true </include-desolvation-forces>
  <time-step-tolerances>
    <minimum-dx> %f </minimum-dx>
  </time-step-tolerances>
  <reactions> mol0-mol1-rxns.xml </reactions>
  <seed> 11111113 </seed>
  <n-trajectories-per-output> %d </n-trajectories-per-output>
  <n-copies> %d </n-copies>
  <n-bin-copies> %d </n-bin-copies>
  <n-steps> %d </n-steps>
  <n-we-steps-per-output> %d </n-we-steps-per-output>
  <max-n-steps> %d </max-n-steps>
  <trajectory-file> traj </trajectory-file>
  <n-steps-per-output> 10000 </n-steps-per-output>
</root>
"""
        with open('input.xml', "w") as fout:
            fout.write((nam_simulation_template % 
                       (self.solvent_eps.get(), self.debyel.get(),
                        self.ntraj.get(), self.nthreads.get(),
                        self.mol0_eps.get(), self.mol1_eps.get(),
                        self.mindx.get(), self.ntrajo.get(),
                        self.ncopies.get(), self.nbincopies.get(),
                        self.nsteps.get(), self.westeps.get(),
                        self.maxnsteps.get())))

        # FIXME check for .dx files
        command = ('PATH=%s:${PATH} bd_top input.xml'
                   % self.bd_path.get())
        if DEBUG > 2: print(command)
        print("::: Running bd_top (this will take a couple of minutes) ...")
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
        return

    def runRxnCrit(self):
        self.runPqr2xml()
        self.makeRxnCriteria()
        return
    
    def prepBD(self):
        self.prepareInputFile()
        return

    def runBD(self):
        """Start BrownDye simulation either in foreground or background."""
        print("::: Starting BrownDye simulation ...")
        command = ('%s/nam_simulation mol0-mol1-simulation.xml >& %s'
                   % (self.bd_path.get(), LOGFILE))
        if self.run_in_background.get():
            p = subprocess.Popen(" nohup %s" % command, shell=True)
            self.jobPID = p.pid
            self.notebook.selectpage('BD Simulation')
            self.logtxt_ent.insert('end', "::: Starting BrownDye simulation in background.\n")
            self.logtxt_ent.insert('end', "::: Job PID: %d \n" % self.jobPID)
        else:
            p = RunThread(self.projectDir.get(), command, self.logtxt_ent)
            p.start()
            self.notebook.selectpage('BD Simulation')
            self.logtxt_ent.insert('end', "::: Starting BrownDye simulation ...\n")
            time.sleep(1)
            tl = MonitorThread(LOGFILE, p, 3600, self.ntraj.get(), self.logtxt_ent,
                               self.messagebar1, self.messagebar2,
                               self.messagebar3, self.bd_path.get())
            tl.start()
        return

    def killBD(self):
        """Terminate background BD job, if exists."""
        try:
            os.system('pkill -9 -P %d' % self.jobPID)
        except AttributeError:
            print("No background job running!")
        return

    def loadTrajectoryFile(self):
        """Load trajectory file."""
        file_name = tkFileDialog.askopenfilename(title='Trajectory File',
                                                 initialdir='',
                                                 filetypes=[('xml files', '*.xml'),
                                                            ('all files', '*')],
                                                 parent=self.parent)
        self.traj_f.set(file_name)
        return
    
    def analyzeTrajectoryFile(self):
        """Process trajectory file and print various stats."""
        if not os.path.isfile(self.traj_f.get()):
            print("::: %s does not exist, check the path."
                  % self.traj_f.get())
            return
        traj_f_base = self.traj_f.get().strip('.xml')
        command = 'grep escaped %s |wc -l' % self.traj_f.get()
        escaped = subprocess.check_output(command, shell=True)
        command = 'grep stuck %s |wc -l' % self.traj_f.get()
        stuck = subprocess.check_output(command, shell=True)
        command = 'grep reacted %s |wc -l' % self.traj_f.get()
        reacted = subprocess.check_output(command, shell=True)
        escaped = int(escaped) / 2
        stuck = int(stuck) / 2
        reacted = int(reacted) / 2
        logline = ('Trajectory file %s contains %d stuck, '
                   '%d escaped and %d reacted events.\n' 
                   % (self.traj_f.get(), stuck, escaped, reacted))
        self.message_ent.insert('end', "%s" % logline)
        if reacted == 0:
            print("::: No association events found in %s trajectory file."
                  % self.traj_f.get())
            return
        command = ('%s/process_trajectories -traj %s '
                   '-index %s.index.xml -srxn association'
                   % (self.bd_path.get(), self.traj_f.get(),
                       traj_f_base))
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=True)
        p.wait()
        stdout, stderr = p.communicate()
        if p.returncode > 0:
            print("::: Error processing %s trajectory file." % self.traj_f.get())
            return
        t = etree.fromstring(stdout)
        tr = t.xpath('//trajectory/number')
        traj_index = [int(tr[x].text.strip()) for x in range(len(tr))]
        logline = ('%d association event trajectories found (index numbers: %s)\n'
                   % (len(traj_index), str(traj_index)))
        self.message_ent.insert('end', "%s" % logline)
        # number of frames
        for i in traj_index:
            with open(self.traj_f.get(), 'r') as f:
                d = etree.parse(f)
                xpath_str = 'trajectory[n-traj=" %d "]//s/n' % i
                j = d.xpath(xpath_str)
                logline = ('trajectory %d: %s frames\n' 
                           % (i, j[0].text.strip()))
                self.message_ent.insert('end', "%s" % logline)
                self.dialog_idx.insert('end', i)
        return

    def selectTrajIndex(self, result):
        """Select trajectory index number from Dialog window"""
        sel = self.dialog_idx.getcurselection()
        if len(sel) > 0:
            if DEBUG > 1: print("::: Selection: %s" % sel)
            self.traj_index_n.set(sel)
            self.messagebar_idx.insert(sel)
        self.dialog_idx.deactivate(result)
        return
    
    def convertTrajectoryToXYZ(self):
        """Convert XML trajectory to XYZ format."""
        traj_f_base = self.traj_f.get().strip('.xml')
        ofile = 'trajectory-%d.xml' % self.traj_index_n.get()
        command = ('%s/process_trajectories -traj %s '
                   '-index %s.index.xml -n %d > %s'
                   % (self.bd_path.get(), self.traj_f.get(),
                      traj_f_base, self.traj_index_n.get(), ofile))
        print("::: Processing trajectory %s ..." % self.traj_f.get())
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
            return
        self.xyz_ofile = 'trajectory-%d.xyz' % self.traj_index_n.get()
        command = ('%s/xyz_trajectory -mol0 mol0-atoms.pqrxml '
                   '-mol1 mol1-atoms.pqrxml -trajf %s '
                   '> %s'
                   % (self.bd_path.get(), ofile, self.xyz_ofile))
        print("::: Converting %s to XYZ trajectory ..." % ofile)
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
            return
        return

    def convertTrajectoryToPQR(self):
        """Convert XML trajectory to PQR format."""
        return
    
    def loadTrajectoryFileXYZ(self):
        """Load XYZ trajectory to pymol."""
        xyz_trajectory_object = 'trajectory-%d' % self.traj_index_n.get()
        try:
            pymol.cmd.load(self.xyz_ofile, xyz_trajectory_object)
        except:
            e = sys.exc_info()[0]
            print("::: Loading xyz trajectory failed!")
            print("::: Error: %s" % e)
        return
        
    def execute(self, result):
        """Quid BD plugin."""
        print("Exiting BrownDye Plugin ...")
        if __name__ == '__main__':
            self.parent.destroy()
        else:
            self.dialog.withdraw()
        print("Done.")
        return


class RunThread(Thread):
    """Thread management class."""
    def __init__(self, work_dir, command, page):
        Thread.__init__(self)
        self.page = page
        self.command = command
        if DEBUG > 2: print("%s" % (command))
        self.work_dir = work_dir
        if DEBUG > 2: print("Work directory: %s" % (work_dir))
        self.pid = 0
        return

    def run(self):
        print("::: Project directory: %s" % (self.work_dir))
        # current_dir = os.getcwd()
        os.chdir(self.work_dir)
        self.page.yview('moveto', 1.0)
        p = subprocess.Popen(self.command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=True)
        p.wait()
        stdout, stderr = p.communicate()
        self.outlog = stdout
        self.status = p.returncode
	# try:
	#  self.page.insert('end',"%s %s" % (out, err))
        #  self.page.yview('moveto', 1.0)
        # except:
	#  self.screwup("Fatal --> No results, so something is wrong")
   	#  return
        # os.chdir(current_dir)
        return


class MonitorThread(Thread):
    """Monitor runing BD job and print out progress information."""
    def __init__(self, logfile, mythread, timeout, totntraj, page,
                 messagebar1, messagebar2, messagebar3, bd_path):
        Thread.__init__(self)
        self.logfile = logfile.replace("\\","/")
        self.totntraj = totntraj
        self.page = page
        self.messagebar1 = messagebar1
        self.messagebar2 = messagebar2
        self.messagebar3 = messagebar3
        self.bd_path = bd_path
        if DEBUG > 2: print("::: logfile: %s" % logfile)
        if DEBUG > 2: print("::: self.logfile: %s" % self.logfile)
        self.mythread = mythread
        self.timeout = timeout
    def run(self):
        seconds = 0
        readcount = 0 
        while self.mythread.is_alive() and seconds < self.timeout:
            time.sleep(5)
            if os.path.getsize(self.logfile) > readcount+20:
                with open(self.logfile,'r') as f:
                    f.seek(readcount,0)
                    while readcount < os.path.getsize(self.logfile):
                        logline = f.readline()
                        if not logline: break
                        print(logline, end='')
                        self.page.insert('end', "%s" % logline)
                        readcount = f.tell()

                seconds = seconds+5
                #transfer_status['log'] = 'Running for %d seconds'%seconds
            with open('results.xml', 'r') as f:
                results = etree.parse(f)
                ntraj = results.xpath('//reactions/n-trajectories')[0].text.strip()
                stuck = results.xpath('//reactions/stuck')[0].text.strip()
                escaped = results.xpath('//reactions/escaped')[0].text.strip()
                ncompleted = results.xpath('//reactions/completed/name')[0].text.strip()
                completed = results.xpath('//reactions/completed/n')[0].text.strip()
                mymessage = '%s out of %d' % (ntraj, self.totntraj)
                self.messagebar1.message('state', mymessage)
                mymessage= '%s / %s / %s' % (completed, escaped, stuck) 
                self.messagebar2.message('state', mymessage)
                
            command = ('cat results.xml | %s/compute_rate_constant'
                       % (self.bd_path))
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, shell=True)
            p.wait()
            stdout, stderr = p.communicate()
            rates = etree.fromstring(stdout)
            rate_constant = rates.xpath('//rate-constant/mean')[0].text.strip()
            rxn_probability = rates.xpath('//reaction-probability/mean')[0].text.strip()
            mymessage = ('%s / %s' % (rate_constant, rxn_probability))
            self.messagebar3.message('state', mymessage)
            
        time.sleep(2)
        print("::: BrownDye simulation finished.")
        self.page.insert('end', "::: BrownDye simulation finished\n")
        return
    
class StopThread(Thread):
    """Kill running thread."""
    def __init__(self, mythread):
        self.pid = mythread.pid
        os.system('kill -9 %d' % self.pid)
        return
    
class Psize:
    """Calculate grid size dimensions for a pqr molecule.
    
    This is based on pdb2pqr version of psize. All licensing info applies.

    Note: CFAC and FADD defaults changed to accomodate BrownDye requirements.
    
    CFAC 1.7 -> 3.0
    FADD 20 -> 50

    This significantly increases the grid size and thus memory requirements.
    """
    def __init__(self):
        self.constants = {"CFAC": 3.0, "FADD": 50, "SPACE": 0.50,
                          "GMEMFAC": 200, "GMEMCEIL": 400, "OFAC": 0.1,
                          "REDFAC": 0.25, "TFAC_ALPHA": 9e-5,
                          "TFAC_XEON": 3e-4, "TFAC_SPARC": 5e-4}
        self.minlen = [360.0, 360.0, 360.0]
        self.maxlen = [0.0, 0.0, 0.0]
        self.q = 0.0
        self.gotatom = 0
        self.gothet = 0
        self.olen = [0.0, 0.0, 0.0]
        self.cen = [0.0, 0.0, 0.0]
        self.clen = [0.0, 0.0, 0.0]
        self.flen = [0.0, 0.0, 0.0]
        self.n = [0, 0, 0]
        self.np = [0.0, 0.0, 0.0]
        self.nsmall = [0, 0, 0]
        self.nfocus = 0

    def parseInput(self, filename):
        """Parse input structure file in PDB or PQR format."""
        file = open(filename, "r")
        self.parseLines(file.readlines())

    def parseLines(self, lines):
        """Parse the lines."""
        for line in lines:
            if string.find(line, "ATOM") == 0:
                subline = string.replace(line[30:], "-", " -")
                ## words = string.split(subline) ## this is a hack
                ## adhering to lovely PDB format definition (fixed space)
                words = (line[30:38], line[38:46], line[46:54], line[54:63],
                         line[63:69], line[72:76], line[76:78])
                if len(filter(string.strip, words)) < 4:    
                    sys.stderr.write("Can't parse following line:\n")
                    sys.stderr.write("%s\n" % line)
                    sys.exit(2)
                    continue
                self.gotatom = self.gotatom + 1
                try:
                    self.q = self.q + float(words[3])
                except ValueError as ve:
                    print("Error parsing line", line)
                    raise ve
                rad = float(words[4])
                center = []
                for word in words[0:3]:
                    center.append(float(word))
                for i in range(3):
                    self.minlen[i] = min(center[i] - rad, self.minlen[i])
                    self.maxlen[i] = max(center[i] + rad, self.maxlen[i])
            elif string.find(line, "HETATM") == 0:
                self.gothet = self.gothet + 1

    def setConstant(self, name, value):
        """ Set a constant to a value; returns 0 if constant not found """
        try:
            self.constants[name] = value
            return 1
        except KeyError:
            return 0

    def getConstant(self, name):
        """ Get a constant value; raises KeyError if constant not found """
        return self.constants[name]

    def setLength(self, maxlen, minlen):
        """ Compute molecule dimensions """
        for i in range(3):
            self.olen[i] = maxlen[i] - minlen[i]
            if self.olen[i] < 0.1:
                self.olen[i] = 0.1
        return self.olen

    def setCoarseGridDims(self, olen):
        """ Compute coarse mesh dimensions """
        for i in range(3):
            self.clen[i] = self.constants["CFAC"] * olen[i]
        return self.clen

    def setFineGridDims(self, olen, clen):
        """ Compute fine mesh dimensions """
        for i in range(3):
            self.flen[i] = olen[i] + self.constants["FADD"]
            if self.flen[i] > clen[i]:
                #str = "WARNING: Fine length (%.2f) cannot be larger than coarse length (%.2f)\n" % (self.flen[i], clen[i])
                #str = str + "         Setting fine grid length equal to coarse grid length\n"
                #stdout.write(str)
                self.flen[i] = clen[i]
        return self.flen

    def setCenter(self, maxlen, minlen):
        """ Compute molecule center """
        for i in range(3):
            self.cen[i] = (maxlen[i] + minlen[i]) / 2
        return self.cen
    
    def setFineGridPoints(self, flen):
        """ Compute mesh grid points, assuming 4 levels in MG hierarchy """
        tn = [0, 0, 0]
        for i in range(3):
            tn[i] = int(flen[i]/self.constants["SPACE"] + 0.5)
            self.n[i] = 32*(int((tn[i] - 1) / 32.0 + 0.5)) + 1
            if self.n[i] < 33:
                self.n[i] = 33
        return self.n

    def setAll(self):
        """ Set up all of the things calculated individually above. """
        maxlen = self.getMax()
        minlen = self.getMin()
        self.setLength(maxlen, minlen)
        olen = self.getLength()
        
        self.setCoarseGridDims(olen)
        clen = self.getCoarseGridDims()        
        
        self.setFineGridDims(olen, clen)
        flen = self.getFineGridDims()
        
        self.setCenter(maxlen, minlen)
        cen = self.getCenter()
        
        self.setFineGridPoints(flen)
        n = self.getFineGridPoints()
        
    def getMax(self): return self.maxlen
    def getMin(self): return self.minlen
    def getCharge(self): return self.q
    def getLength(self): return self.olen
    def getCoarseGridDims(self): return self.clen
    def getFineGridDims(self): return self.flen
    def getCenter(self): return self.cen
    def getFineGridPoints(self): return self.n

    def runPsize(self, filename):
        """ Parse input PQR file and set parameters. """
        self.parseInput(filename)
        self.setAll()

#############################################
#
# Create root window for testing.
#
##############################################
if __name__ == '__main__':
    
    class App:
        def my_show(self, *args, **kwargs):
            pass

    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Root window')

    widget = BDPlugin(app)
    app.root.mainloop()
