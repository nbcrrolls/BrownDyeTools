#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Last modified: 2019-02-26 19:41:38
#
# pylint: disable=no-member,import-error,missing-docstring,invalid-name,multiple-statements
#
"""BrownDye Tools plugin for Pymol

For documentation see: https://github.com/rokdev/BrownDyeTools

Author : Robert Konecny <rok@ucsd.edu>
Release date: November 2016
License: GNU General Public License version 3

Copyright 2016-19 Robert Konecny, NBCR

This is free software, licensed under the GNU General Public License
version 3. You should have received a copy of the GNU General Public
License along with this program.
If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
import os
import subprocess
if sys.version_info[0] < 3:
    #import urllib2
    #from urllib2 import URLError, HTTPError
    #from urllib import quote
    #from idlelib.TreeWidget import TreeItem, TreeNode
    #from Tkinter import *
    import Tkinter
    #import tkSimpleDialog
    import tkFileDialog
    #import tkMessageBox
    #import tkColorChooser
    #import Queue
    #import ttk
    #from StringIO import StringIO
else:
    #import urllib.request as urllib2
    #from urllib.error import URLError, HTTPError
    #from urllib.parse import quote
    #from idlelib.tree import TreeItem, TreeNode
    #from tkinter import *
    import tkinter as Tkinter
    #from tkinter import simpledialog as tkSimpleDialog
    from tkinter import filedialog as tkFileDialog
    #import tkinter.messagebox as tkMessageBox
    #import tkinter.colorchooser as tkColorChooser
    #import queue as Queue
    #import tkinter.ttk as ttk
    #from io import StringIO

#import string
import re
##import filecmp
import json
import shutil
import random
import time
import datetime
##import tkSimpleDialog
##import tkMessageBox
#import Tkinter
#import tkFileDialog
from threading import Thread
##from lxml import etree
from xml.etree import ElementTree as etree
##import importlib
import Pmw

DEBUG = 0

__version__ = '1.0.0'
__author__ = 'Robert Konecny <rok@ucsd.edu>'

PDB2PQR_PATH = None
APBS_PATH = None
BD_PATH = None
CONFIG_FILE = 'bd-config.json'
DEFAULT_CONTACTS_FILE = 'protein-protein-contacts-default.xml'
MOL = ['mol0', 'mol1']
APBS_EXE = 'apbs'
PDB2PQR_EXE = 'pdb2pqr'
BDTOP_EXE = 'bd_top'
NAM_SIMULATION_EXE = 'nam_simulation'
JOBPID = None

PQR_DEFAULTS = {
    'pqr_ff': 'parse',
    'pqr_opts': '--apbs-input --chain',
    'pqr_ph': 7.0,
    'pqr_assign_only': True,
    'pqr_use_propka': False,
}
PSIZE_DEFAULTS = {
    'cfac': 3.0,
    'fadd': 50.0,
    'gspace': 0.5,
}
APBS_DEFAULTS = {
    'apbs_mode': 'lpbe',
    'apbs_pdie': 4.0,
    'apbs_sdie': 78.0,
    'apbs_srfm': 'smol',
    'apbs_chgm': 'spl2',
    'apbs_bcfl': 'sdh',
    'apbs_sdens': 10.,
    'apbs_swin': 0.3,
    'apbs_srad': 1.4,
    'apbs_temp': 298.15,
    'apbs_ion_charge': [1, -1],
    'apbs_ion_conc': [0.15, 0.15],
    'apbs_ion_radius': [1.0, 1.0],
}
BD_DEFAULTS = {
    'contacts_f': DEFAULT_CONTACTS_FILE,
    'default_contacts_f': True,
    'ntraj': 100,
    'nthreads': 1,
    'mindx': 0.2,
    'sdie': APBS_DEFAULTS['apbs_sdie'],
    'pdie': [APBS_DEFAULTS['apbs_pdie'], APBS_DEFAULTS['apbs_pdie']],
    'debyel': [0.0, 0.0],
    'ntrajo': 1,
    'ncopies': 200,
    'nbincopies': 200,
    'nsteps': 10000,
    'westeps': 10000,
    'maxnsteps': 10000,
    'nsteps_per_output': 2,
    'rxn_distance': 5.0,
    'npairs': 3,
}

if DEBUG > 0:
    PDB2PQR_PATH = '/opt/pdb2pqr-linux-bin64-2.1.0/'
    APBS_PATH = '/export1/Builds/SandBox-2018-03-08T15.28.17/bin/'
    BD_PATH = '/home/rok/BrownianDynamics/browndye-2016.4.14/bin/'

if "PDB2PQR_PATH" in os.environ: PDB2PQR_PATH = os.environ["PDB2PQR_PATH"]
if "APBS_PATH" in os.environ: APBS_PATH = os.environ["APBS_PATH"]
if "BD_PATH" in os.environ: BD_PATH = os.environ["BD_PATH"]

def __init__(self):
    """BrownDyeTools plugin for PyMol."""
    self.menuBar.addmenuitem('Plugin', 'command',
                             'BrownDye Tools', label='BrownDye Tools',
                             command=lambda s=self: BDPlugin(s))

class DummyPymol(object):
    """Dummy pymol class when running standalone GUI."""
    class Cmd:
        def load(self, name, sel=''):
            pass

        def get_names(self, name):
            return ['mol0', 'mol1', 'map0', 'map1']

        def get_type(self, thing):
            if thing.startswith('mol'):
                return 'object:molecule'
            else:
                return 'object:map'
    cmd = Cmd()

if 'pymol.gui' in sys.modules:
    try:
        import pymol
    except ImportError:
        print("::: Pymol import failed - Pymol features not available!")
        pymol = DummyPymol()
else:
    print("::: Pymol import failed - Pymol features will not available!")
    pymol = DummyPymol()

class BDPlugin(object):
    """ The main BrowDye plugin class."""
    def __init__(self, app):
        self.parent = app.root
        Tkinter.Grid.rowconfigure(self.parent, 0, weight=1)
        Tkinter.Grid.columnconfigure(self.parent, 0, weight=1)
        self.createGUI()

    def createGUI(self):
        """The main GUI class - sets up all GUI elements."""
        self.dialog = Pmw.Dialog(self.parent, buttons=('Exit',),
                                 title='BrownDyeTools Plugin for PyMOL',
                                 command=self.exitBDPlugin)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        self.projectDir = Tkinter.StringVar()
        self.projectDir.set(self.getProjectDir())
        self.pdb2pqr_path = Tkinter.StringVar()
        self.pdb2pqr_path.set(PDB2PQR_PATH)
        self.apbs_path = Tkinter.StringVar()
        self.apbs_path.set(APBS_PATH)
        self.bd_path = Tkinter.StringVar()
        self.bd_path.set(BD_PATH)
        self.config_file = Tkinter.StringVar()
        self.config_file.set(CONFIG_FILE)

        # parameters used by pdb2pqr
        self.mol = [Tkinter.StringVar(), Tkinter.StringVar()]
        self.mol[0].set(None)
        self.mol[1].set(None)
        self.pqr = [Tkinter.StringVar(), Tkinter.StringVar()]
        self.pqr[0].set(None)
        self.pqr[1].set(None)
        self.mol_object = [Tkinter.StringVar(), Tkinter.StringVar()]
        self.mol_object[0].set(None)
        self.mol_object[1].set(None)
        self.pqr_opts = Tkinter.StringVar()
        self.pqr_opts.set(PQR_DEFAULTS['pqr_opts'])
        self.pqr_assign_only = Tkinter.BooleanVar()
        self.pqr_assign_only.set(PQR_DEFAULTS['pqr_assign_only'])
        self.pqr_use_propka = Tkinter.BooleanVar()
        self.pqr_use_propka.set(PQR_DEFAULTS['pqr_use_propka'])
        self.pqr_ph = Tkinter.DoubleVar()
        self.pqr_ph.set(PQR_DEFAULTS['pqr_ph'])
        self.pqr_ff = Tkinter.StringVar()
        self.pqr_ff.set(PQR_DEFAULTS['pqr_ff'])

        # psize/APBS parameters and defaults
        self.gspace = Tkinter.DoubleVar()
        self.gspace.set(PSIZE_DEFAULTS['gspace'])
        self.cfac = Tkinter.DoubleVar()
        self.cfac.set(PSIZE_DEFAULTS['cfac'])
        self.fadd = Tkinter.DoubleVar()
        self.fadd.set(PSIZE_DEFAULTS['fadd'])
        self.dime = [[Tkinter.IntVar() for _ in range(3)],
                     [Tkinter.IntVar() for _ in range(3)]]
        self.cglen = [[Tkinter.DoubleVar() for _ in range(3)],
                      [Tkinter.DoubleVar() for _ in range(3)]]
        self.fglen = [[Tkinter.DoubleVar() for _ in range(3)],
                      [Tkinter.DoubleVar() for _ in range(3)]]
        self.apbs_mode = Tkinter.StringVar()
        self.apbs_mode.set(APBS_DEFAULTS['apbs_mode'])
        self.apbs_bcfl = Tkinter.StringVar()
        self.apbs_bcfl.set(APBS_DEFAULTS['apbs_bcfl'])
        self.apbs_ion_charge = [Tkinter.IntVar() for _ in range(2)]
        self.apbs_ion_conc = [Tkinter.DoubleVar() for _ in range(2)]
        self.apbs_ion_radius = [Tkinter.DoubleVar() for _ in range(2)]
        for i in range(2):
            self.apbs_ion_charge[i].set(APBS_DEFAULTS['apbs_ion_charge'][i])
            self.apbs_ion_conc[i].set(APBS_DEFAULTS['apbs_ion_conc'][i])
            self.apbs_ion_radius[i].set(APBS_DEFAULTS['apbs_ion_radius'][i])
        self.apbs_pdie = Tkinter.DoubleVar()
        self.apbs_pdie.set(APBS_DEFAULTS['apbs_pdie'])
        self.apbs_sdie = Tkinter.DoubleVar()
        self.apbs_sdie.set(APBS_DEFAULTS['apbs_sdie'])
        self.apbs_chgm = Tkinter.StringVar()
        self.apbs_chgm.set(APBS_DEFAULTS['apbs_chgm'])
        self.apbs_sdens = Tkinter.DoubleVar()
        self.apbs_sdens.set(APBS_DEFAULTS['apbs_sdens'])
        self.apbs_swin = Tkinter.DoubleVar()
        self.apbs_swin.set(APBS_DEFAULTS['apbs_swin'])
        self.apbs_srfm = Tkinter.StringVar()
        self.apbs_srfm.set(APBS_DEFAULTS['apbs_srfm'])
        self.apbs_srad = Tkinter.DoubleVar()
        self.apbs_srad.set(APBS_DEFAULTS['apbs_srad'])
        self.apbs_temp = Tkinter.DoubleVar()
        self.apbs_temp.set(APBS_DEFAULTS['apbs_temp'])

        # reaction criteria
        self.contacts_f = Tkinter.StringVar()
        self.contacts_f.set(BD_DEFAULTS['contacts_f'])
        self.default_contacts_f = Tkinter.BooleanVar()
        self.default_contacts_f.set(BD_DEFAULTS['default_contacts_f'])
        self.rxn_distance = Tkinter.DoubleVar()
        self.rxn_distance.set(BD_DEFAULTS['rxn_distance'])
        self.npairs = Tkinter.IntVar()
        self.npairs.set(BD_DEFAULTS['npairs'])

        # BD parameters and defaults
        self.sdie = Tkinter.DoubleVar()
        self.sdie.set(self.apbs_sdie.get())
        self.pdie = [Tkinter.DoubleVar() for _ in range(2)]
        self.debyel = [Tkinter.DoubleVar() for _ in range(2)]
        for i in range(2):
            self.pdie[i].set(self.apbs_pdie.get())
            self.debyel[i].set(BD_DEFAULTS['debyel'][i])
        self.ntraj = Tkinter.IntVar()
        self.ntraj.set(BD_DEFAULTS['ntraj'])
        self.nthreads = Tkinter.IntVar()
        self.nthreads.set(BD_DEFAULTS['nthreads'])
        self.mindx = Tkinter.DoubleVar()
        self.mindx.set(BD_DEFAULTS['mindx'])
        self.ntrajo = Tkinter.IntVar()
        self.ntrajo.set(BD_DEFAULTS['ntrajo'])
        self.ncopies = Tkinter.IntVar()
        self.ncopies.set(BD_DEFAULTS['ncopies'])
        self.nbincopies = Tkinter.IntVar()
        self.nbincopies.set(BD_DEFAULTS['nbincopies'])
        self.nsteps = Tkinter.IntVar()
        self.nsteps.set(BD_DEFAULTS['nsteps'])
        self.westeps = Tkinter.IntVar()
        self.westeps.set(BD_DEFAULTS['westeps'])
        self.maxnsteps = Tkinter.IntVar()
        self.maxnsteps.set(BD_DEFAULTS['maxnsteps'])
        self.nsteps_per_output = Tkinter.IntVar()
        self.nsteps_per_output.set(BD_DEFAULTS['nsteps_per_output'])

        self.run_in_background = Tkinter.BooleanVar()
        self.run_in_background.set(False)

        # Analysis
        self.traj_f = Tkinter.StringVar()
        self.traj_f.set('traj0.xml')
        self.traj_index_n = Tkinter.IntVar()
        self.traj_index_n.set(1)
        #######################################################################
        # Main code
        pref = dict(padx=5, pady=5)
        w = Tkinter.Label(self.dialog.interior(),
                          text=('\nBrownDye Tools for PyMOL\n'
                                'Version %s, NBCR 2016-19\n\n'
                                'Plugin for setting up and running '
                                'Brownian dynamics simulations with BrownDye.'
                                % __version__),
                          background='black', foreground='white')
        w.pack(expand=1, fill='both', **pref)
        # create a notebook
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, **pref)
        self.status_bar = Pmw.MessageBar(self.dialog.interior(), entry_width=40,
                                         entry_relief='sunken',
                                         labelpos='w', label_text='Status: ')

        self.status_bar.pack(fill='both', expand=1, **pref)
        self.status_bar.message('state', 'Idle')
        self.balloon = Pmw.Balloon(self.dialog.interior())

        #####################
        # Tab: Configuration
        #####################
        page = self.notebook.add('Configuration')
        self.notebook.tab('Configuration').focus_set()
        config = Tkinter.LabelFrame(page, text='Calculation configuration')
        config.pack(fill='both', expand=True, **pref)

        project_path_ent = Pmw.EntryField(config,
                                          label_text='Create project directory:',
                                          labelpos='wn',
                                          entry_textvariable=self.projectDir)
        project_path_but = Tkinter.Button(config, text='Create',
                                          command=self.createProjectDir)
        self.balloon.bind(project_path_but, 'Create a new project directory')
        label = Tkinter.Label(config, text='or')
        project_path_b_but = Tkinter.Button(config, text='Browse ...',
                                            command=self.browseProjectDir)
        pdb2pqr_path_ent = Pmw.EntryField(config,
                                          label_text='Select PDB2PQR_PATH location:',
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
                                     label_text='Select BD_PATH location:',
                                     labelpos='wn',
                                     entry_textvariable=self.bd_path)
        bd_path_but = Tkinter.Button(config, text='Browse...',
                                     command=self.getBDpath)
        label1 = Tkinter.Label(config, text='')
        config_ent = Pmw.EntryField(config,
                                    label_text='Save or load calculation configuration: ',
                                    labelpos='wn',
                                    entry_textvariable=self.config_file)
        config_but = Tkinter.Button(config, text='Browse...',
                                    command=self.getConfigPath)
        load_config_but = Tkinter.Button(config, text='Load configuration',
                                         command=self.loadConfig)
        self.balloon.bind(load_config_but, 'Load calculation configuration from a file')
        save_config_but = Tkinter.Button(config, text='Save configuration',
                                         command=self.saveConfig)
        self.balloon.bind(save_config_but, 'Save calculation configuration to a file')

        # arrange widgets using grid
        project_path_ent.grid(sticky='we', row=1, column=0, **pref)
        project_path_but.grid(sticky='we', row=1, column=3, **pref)
        label.grid(sticky='we', row=1, column=2, **pref)
        project_path_b_but.grid(sticky='we', row=1, column=1, **pref)
        pdb2pqr_path_ent.grid(sticky='we', row=2, column=0, **pref)
        pdb2pqr_path_but.grid(sticky='we', row=2, column=1, **pref)
        apbs_path_ent.grid(sticky='we', row=3, column=0, **pref)
        apbs_path_but.grid(sticky='we', row=3, column=1, **pref)
        bd_path_ent.grid(sticky='we', row=4, column=0, **pref)
        bd_path_but.grid(sticky='we', row=4, column=1, **pref)
        label1.grid(sticky='we', row=5, column=2, **pref)
        config_ent.grid(sticky='we', row=6, column=0, **pref)
        config_but.grid(sticky='we', row=6, column=1, **pref)
        load_config_but.grid(sticky='e', row=6, column=2, **pref)
        save_config_but.grid(sticky='w', row=6, column=3, **pref)
        config.columnconfigure(0, weight=8)
        config.columnconfigure(1, weight=2)

        #############################
        # Tab: PQR files preparation
        #############################
        page = self.notebook.add('PQR files')
        grp_pqr = Tkinter.LabelFrame(page, text='PQR files')
        grp_pqr.grid(sticky='eswn', row=0, column=0, columnspan=3, **pref)

        pdb0_ent = Pmw.EntryField(grp_pqr,
                                  label_text='Molecule 0 PDB file:', labelpos='wn',
                                  entry_textvariable=self.mol[0])
        pdb0_but = Tkinter.Button(grp_pqr, text='Browse...',
                                  command=lambda: self.getPDBMol(0))
        label0 = Tkinter.Label(grp_pqr, text='or')
        pymol_obj0_opt = Pmw.OptionMenu(grp_pqr, labelpos='w',
                                        label_text='Select molecule 0: ',
                                        menubutton_textvariable=self.mol_object[0],
                                        menubutton_width=7,
                                        items=(['None'] + pymol.cmd.get_names("all")))
        pdb1_ent = Pmw.EntryField(grp_pqr,
                                  label_text='Molecule 1 PDB file:', labelpos='wn',
                                  entry_textvariable=self.mol[1])
        pdb1_but = Tkinter.Button(grp_pqr, text='Browse...',
                                  command=lambda: self.getPDBMol(1))
        label1 = Tkinter.Label(grp_pqr, text='or')
        pymol_obj1_opt = Pmw.OptionMenu(grp_pqr, labelpos='w',
                                        label_text='Select molecule 1: ',
                                        menubutton_textvariable=self.mol_object[1],
                                        menubutton_width=7,
                                        items=(['None'] + pymol.cmd.get_names("all")))
        pqr_ff_opt = Pmw.OptionMenu(grp_pqr, labelpos='w',
                                    label_text='Force field: ',
                                    menubutton_textvariable=self.pqr_ff,
                                    menubutton_width=7,
                                    items=['parse', 'charmm', 'amber', 'swanson'])
        pqr_an_but = Tkinter.Checkbutton(grp_pqr,
                                         text=('Assign charge and radius only '
                                               '(no structure optimization)'),
                                         variable=self.pqr_assign_only,
                                         onvalue=True, offvalue=False,
                                         command=lambda: propka_but.toggle())
        propka_but = Tkinter.Checkbutton(grp_pqr,
                                         text=('Use PROPKA to assign protonation states'),
                                         variable=self.pqr_use_propka,
                                         onvalue=True, offvalue=False,
                                         command=lambda: pqr_an_but.toggle())
        ph_ent = Pmw.EntryField(grp_pqr, labelpos='w',
                                label_text='pH: ',
                                validate={'validator': 'real', 'min': 0.00},
                                entry_textvariable=self.pqr_ph,
                                entry_width=4)
        pqr_opt_but = Tkinter.Button(page, text='Create PQR files',
                                     command=self.pdb2pqr)
        label2 = Tkinter.Label(page, text='or load your PQR files:')
        pqr_0_ent = Pmw.EntryField(page,
                                   label_text='Molecule 0 PQR file:', labelpos='wn',
                                   entry_textvariable=self.pqr[0])
        pqr_0_but = Tkinter.Button(page, text='Browse...',
                                   command=lambda: self.getPQRMol(0))
        pqr_1_ent = Pmw.EntryField(page,
                                   label_text='Molecule 1 PQR file:', labelpos='wn',
                                   entry_textvariable=self.pqr[1])
        pqr_1_but = Tkinter.Button(page, text='Browse...',
                                   command=lambda: self.getPQRMol(1))

        pdb0_ent.grid(sticky='we', row=0, column=0, **pref)
        pdb0_but.grid(sticky='we', row=0, column=1, **pref)
        label0.grid(sticky='we', row=0, column=2, **pref)
        pymol_obj0_opt.grid(sticky='we', row=0, column=3, **pref)
        pdb1_ent.grid(sticky='we', row=1, column=0, **pref)
        pdb1_but.grid(sticky='we', row=1, column=1, **pref)
        label1.grid(sticky='we', row=1, column=2, **pref)
        pymol_obj1_opt.grid(sticky='we', row=1, column=3, **pref)
        pqr_ff_opt.grid(sticky='we', row=2, column=0, **pref)
        pqr_an_but.grid(sticky='we', row=3, column=0, **pref)
        propka_but.grid(sticky='we', row=3, column=1, columnspan=2, **pref)
        ph_ent.grid(sticky='we', row=3, column=3, **pref)
        pqr_opt_but.grid(sticky='we', row=4, column=0, **pref)
        label2.grid(sticky='we', row=5, column=0, **pref)
        pqr_0_ent.grid(sticky='we', row=6, column=0, **pref)
        pqr_0_but.grid(sticky='we', row=6, column=1, **pref)
        pqr_1_ent.grid(sticky='we', row=7, column=0, **pref)
        pqr_1_but.grid(sticky='we', row=7, column=1, **pref)
        page.columnconfigure(0, weight=1)
        page.columnconfigure(1, weight=1)

        ############
        # Tab: APBS
        ############
        page = self.notebook.add('APBS')
        grp_grids = Tkinter.LabelFrame(page, text='Grid size')
        grp_apbs = Tkinter.LabelFrame(page, text='APBS options')
        grp_grids.grid(sticky='eswn', row=0, column=0, columnspan=3, **pref)
        grp_apbs.grid(sticky='eswn', row=0, column=3, columnspan=2, **pref)
        #page.columnconfigure(0, weight=2)
        #page.columnconfigure(1, weight=1)


        gspace_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                    label_text='Grid spacing (A): ',
                                    validate={'validator': 'real', 'min': 0.00},
                                    entry_textvariable=self.gspace,
                                    entry_width=4)
        fadd_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                  label_text='Grid span beyond the molecule (in A) : ',
                                  validate={'validator': 'real', 'min': 0.00},
                                  entry_textvariable=self.fadd,
                                  entry_width=4)
        label0 = Tkinter.Label(grp_grids, text='Molecule 0')
        dime0_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='dime: ',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime[0][0],
                                     entry_width=5)
        dime0_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime[0][1],
                                     entry_width=5)
        dime0_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime[0][2],
                                     entry_width=5)
        cglen0_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='cglen: ',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen[0][0],
                                      entry_width=8)
        cglen0_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen[0][1],
                                      entry_width=8)
        cglen0_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen[0][2],
                                      entry_width=8)
        fglen0_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='fglen: ',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen[0][0],
                                      entry_width=8)
        fglen0_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen[0][1],
                                      entry_width=8)
        fglen0_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen[0][2],
                                      entry_width=8)
        get_size0_but = Tkinter.Button(grp_grids,
                                       text="Set grid",
                                       command=lambda: self.getSizemol(0))
        label1 = Tkinter.Label(grp_grids, text='Molecule 1')
        dime1_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='dime: ',
                                     # value=self.dime1[0].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime[1][0],
                                     entry_width=5)
        dime1_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime[1][1],
                                     entry_width=5)
        dime1_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime[1][2],
                                     entry_width=5)
        cglen1_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='cglen: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen[1][0],
                                      entry_width=8)
        cglen1_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen[1][1],
                                      entry_width=8)
        cglen1_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen[1][2],
                                      entry_width=8)
        fglen1_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='fglen: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen[1][0],
                                      entry_width=8)
        fglen1_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen[1][1],
                                      entry_width=8)
        fglen1_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen[1][2],
                                      entry_width=8)
        get_size1_but = Tkinter.Button(grp_grids,
                                       text="Set grid",
                                       command=lambda: self.getSizemol(1))

        gspace_ent.grid(sticky='we', row=0, column=0, **pref)
        fadd_ent.grid(sticky='we', row=0, column=1, columnspan=2, **pref)
        label0.grid(sticky='we', row=1, column=0, **pref)
        dime0_0_ent.grid(sticky='we', row=2, column=0, **pref)
        dime0_1_ent.grid(sticky='we', row=2, column=1, **pref)
        dime0_2_ent.grid(sticky='we', row=2, column=2, **pref)
        cglen0_0_ent.grid(sticky='we', row=3, column=0, **pref)
        cglen0_1_ent.grid(sticky='we', row=3, column=1, **pref)
        cglen0_2_ent.grid(sticky='we', row=3, column=2, **pref)
        fglen0_0_ent.grid(sticky='we', row=4, column=0, **pref)
        fglen0_1_ent.grid(sticky='we', row=4, column=1, **pref)
        fglen0_2_ent.grid(sticky='we', row=4, column=2, **pref)
        get_size0_but.grid(sticky='we', row=5, column=2, **pref)

        label1.grid(sticky='we', row=6, column=0, **pref)
        dime1_0_ent.grid(sticky='we', row=7, column=0, **pref)
        dime1_1_ent.grid(sticky='we', row=7, column=1, **pref)
        dime1_2_ent.grid(sticky='we', row=7, column=2, **pref)
        cglen1_0_ent.grid(sticky='we', row=8, column=0, **pref)
        cglen1_1_ent.grid(sticky='we', row=8, column=1, **pref)
        cglen1_2_ent.grid(sticky='we', row=8, column=2, **pref)
        fglen1_0_ent.grid(sticky='we', row=9, column=0, **pref)
        fglen1_1_ent.grid(sticky='we', row=9, column=1, **pref)
        fglen1_2_ent.grid(sticky='we', row=9, column=2, **pref)

        get_size1_but.grid(sticky='we', row=10, column=2, **pref)

        #grp_grids.columnconfigure(0, weight=9)
        #grp_grids.columnconfigure(1, weight=1)

        apbs_mode_ent = Pmw.OptionMenu(grp_apbs, labelpos='w',
                                       label_text='APBS mode: ',
                                       menubutton_textvariable=self.apbs_mode,
                                       menubutton_width=10,
                                       items=['lpbe', 'npbe'])
        solvent_die_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                         label_text='Solvent eps :',
                                         validate={'validator': 'real', 'min': 0.0},
                                         entry_textvariable=self.apbs_sdie,
                                         entry_width=10)
        solute_die_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                        label_text='Molecule 0/1 eps: ',
                                        validate={'validator': 'real', 'min': 0.0},
                                        entry_textvariable=self.apbs_pdie,
                                        entry_width=10)
        ion1_charge_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                         label_text='Ion(1) charge: ',
                                         validate={'validator': 'integer', 'min': -2},
                                         entry_textvariable=self.apbs_ion_charge[0],
                                         entry_width=5)
        ion1_conc_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                       label_text='conc.: ',
                                       validate={'validator': 'real', 'min': 0.0},
                                       entry_textvariable=self.apbs_ion_conc[0],
                                       entry_width=5)
        ion1_rad_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                      label_text='radius: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.apbs_ion_radius[0],
                                      entry_width=5)
        ion2_charge_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                         label_text='Ion(2) charge: ',
                                         validate={'validator': 'integer', 'min': -2},
                                         entry_textvariable=self.apbs_ion_charge[1],
                                         entry_width=5)
        ion2_conc_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                       label_text='conc.: ',
                                       validate={'validator': 'real', 'min': 0.0},
                                       entry_textvariable=self.apbs_ion_conc[1],
                                       entry_width=5)
        ion2_rad_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                      label_text='radius: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.apbs_ion_radius[1],
                                      entry_width=5)
        sdens_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                   label_text='Surf. sphere density: ',
                                   validate={'validator': 'real', 'min': 0.0},
                                   entry_textvariable=self.apbs_sdens, entry_width=5)
        srad_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                  label_text='Solvent radius: ',
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.apbs_srad, entry_width=5)
        swin_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                  label_text='Spline window: ',
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.apbs_swin, entry_width=5)
        temp_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                  label_text='Temperature: ',
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.apbs_temp, entry_width=5)
        bcfl_ent = Pmw.OptionMenu(grp_apbs, labelpos='w',
                                  label_text='Boundary condition: ',
                                  menubutton_textvariable=self.apbs_bcfl,
                                  menubutton_width=5,
                                  items=['sdh', 'zero', 'mdh', 'focus'])
        chgm_ent = Pmw.OptionMenu(grp_apbs, labelpos='w',
                                  label_text='Charge mapping: ',
                                  menubutton_textvariable=self.apbs_chgm,
                                  menubutton_width=5,
                                  items=['spl2', 'spl0', 'spl4'])
        srfm_ent = Pmw.OptionMenu(grp_apbs, labelpos='w',
                                  label_text='Diel. surf. calc. method: ',
                                  menubutton_textvariable=self.apbs_srfm,
                                  menubutton_width=5,
                                  items=['smol', 'mol', 'spl2', 'spl4'])
        run_apbs_but = Tkinter.Button(page,
                                      text="Run APBS to generate grids",
                                      command=self.runAPBS)

        solvent_die_ent.grid(sticky='we', row=0, column=0, **pref)
        solute_die_ent.grid(sticky='we', row=1, column=0, **pref)
        sdens_ent.grid(sticky='we', row=2, column=0, **pref)
        srad_ent.grid(sticky='we', row=3, column=0, **pref)
        swin_ent.grid(sticky='we', row=4, column=0, **pref)
        temp_ent.grid(sticky='we', row=5, column=0, **pref)

        ion1_charge_ent.grid(sticky='we', row=6, column=0, **pref)
        ion1_conc_ent.grid(sticky='we', row=6, column=1, **pref)
        ion1_rad_ent.grid(sticky='we', row=6, column=2, **pref)
        ion2_charge_ent.grid(sticky='we', row=7, column=0, **pref)
        ion2_conc_ent.grid(sticky='we', row=7, column=1, **pref)
        ion2_rad_ent.grid(sticky='we', row=7, column=2, **pref)

        apbs_mode_ent.grid(sticky='we', row=0, column=1, columnspan=2, **pref)
        bcfl_ent.grid(sticky='we', row=1, column=1, columnspan=2, **pref)
        chgm_ent.grid(sticky='we', row=2, column=1, columnspan=2, **pref)
        srfm_ent.grid(sticky='we', row=3, column=1, columnspan=2, **pref)

        run_apbs_but.grid(sticky='e', row=8, column=4, **pref)
        page.columnconfigure(0, weight=1)
        page.columnconfigure(1, weight=1)

        ###############################
        # Tab: Reaction criteria setup
        ###############################
        page = self.notebook.add('Reaction citeria')
        grp_rxn = Tkinter.LabelFrame(page, text='Setup reaction criteria')
        grp_rxn.grid(sticky='eswn', row=0, column=0, columnspan=2, **pref)

        contacts_ent = Pmw.EntryField(grp_rxn,
                                      label_text='Contacts file:', labelpos='w',
                                      entry_textvariable=self.contacts_f, entry_width=50)
        contacts_but = Tkinter.Button(grp_rxn, text='Browse...',
                                      command=self.getContacts)
        def_contacts_but = Tkinter.Checkbutton(grp_rxn,
                                               text='or use default contacts file ',
                                               variable=self.default_contacts_f,
                                               onvalue=True, offvalue=False)
        rxn_distance_ent = Pmw.EntryField(grp_rxn, labelpos='wn',
                                          label_text='Reaction distance: ',
                                          # value=self.rxn_distance.get(),
                                          validate={'validator': 'real', 'min': 0.01},
                                          entry_textvariable=self.rxn_distance,
                                          entry_width=10)
        npairs_ent = Pmw.EntryField(grp_rxn, labelpos='wn',
                                    label_text='Number of reaction pairs: ',
                                    # value=self.npairs.get(),
                                    validate={'validator': 'integer', 'min': 1},
                                    entry_textvariable=self.npairs, entry_width=10)
        run_rxn_crit = Tkinter.Button(page, text="Create reaction files",
                                      command=self.runRxnCrit)

        contacts_ent.grid(sticky='we', row=0, column=0, **pref)
        contacts_but.grid(sticky='e', row=0, column=1, **pref)
        def_contacts_but.grid(sticky='e', row=0, column=2, **pref)
        rxn_distance_ent.grid(sticky='we', row=1, column=0, **pref)
        npairs_ent.grid(sticky='we', row=2, column=0, **pref)
        run_rxn_crit.grid(sticky='e', row=3, column=0, **pref)
        grp_rxn.columnconfigure(0, weight=1)
        page.columnconfigure(0, weight=1)

        ######################
        # Tab: BD input files
        ######################
        page = self.notebook.add('BD setup')
        grp_bdinput = Tkinter.LabelFrame(page, text='BrownDye input file')
        grp_bdinput.grid(sticky='eswn', row=0, column=0, columnspan=2, **pref)

        sdie_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                  label_text='Solvent eps: ',
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.sdie, entry_width=5)
        debyel_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                    label_text='Solvent Debye length: ',
                                    validate={'validator': 'real', 'min': 0.0},
                                    entry_textvariable=self.debyel[0], entry_width=5)
        get_debyel_but = Tkinter.Button(grp_bdinput, text='Get Debye length',
                                        command=self.getDebyeLength)
        pdie0_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                   label_text='Molecule 0 eps: ',
                                   validate={'validator': 'real', 'min': 0.0},
                                   entry_textvariable=self.pdie[0], entry_width=5)
        pdie1_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                   label_text='Molecule 1 eps: ',
                                   validate={'validator': 'real', 'min': 0.0},
                                   entry_textvariable=self.pdie[1], entry_width=5)
        ntraj_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                   label_text='Number of trajectories: ',
                                   validate={'validator': 'integer', 'min': 1},
                                   entry_textvariable=self.ntraj, entry_width=5)
        nthreads_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                      label_text='Number of threads: ',
                                      validate={'validator': 'integer', 'min': 1},
                                      entry_textvariable=self.nthreads, entry_width=5)
        mindx_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                   label_text='Time step tolerance: ',
                                   validate={'validator': 'real', 'min': 0.01},
                                   entry_textvariable=self.mindx, entry_width=5)
        ntrajo_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                    label_text='Number of trajectories per output: ',
                                    validate={'validator': 'integer', 'min': 1},
                                    entry_textvariable=self.ntrajo, entry_width=5)
        ncopies_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                     label_text='Number of copies: ',
                                     validate={'validator': 'integer', 'min': 1},
                                     entry_textvariable=self.ncopies, entry_width=5)
        nbincopies_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                        label_text='Number of bin copies: ',
                                        validate={'validator': 'integer', 'min': 1},
                                        entry_textvariable=self.nbincopies,
                                        entry_width=5)
        nsteps_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                    label_text='Number of steps: ',
                                    validate={'validator': 'integer', 'min': 1},
                                    entry_textvariable=self.nsteps, entry_width=5)
        westeps_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                     label_text='Number of WE steps per output: ',
                                     validate={'validator': 'integer', 'min': 1},
                                     entry_textvariable=self.westeps, entry_width=5)
        maxnsteps_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                       label_text='Max number of steps: ',
                                       validate={'validator': 'integer', 'min': 1},
                                       entry_textvariable=self.maxnsteps, entry_width=5)
        nsteps_per_output_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                               label_text='Number of steps per output: ',
                                               validate={'validator': 'integer', 'min': 1},
                                               entry_textvariable=self.nsteps_per_output,
                                               entry_width=5)
        prep_bd_but = Tkinter.Button(page, text="Generate BrownDye input files",
                                     command=self.prepBD)

        #self.status_bar = Pmw.MessageBar(page, entry_width=20, entry_relief='sunken',
        #                                 labelpos='w', label_text='Status:')
        #self.status_bar.message('state', '')

        sdie_ent.grid(sticky='we', row=1, column=0, **pref)
        debyel_ent.grid(sticky='we', row=2, column=0, **pref)
        get_debyel_but.grid(sticky='w', row=2, column=1, **pref)
        pdie0_ent.grid(sticky='we', row=3, column=0, **pref)
        pdie1_ent.grid(sticky='we', row=3, column=1, **pref)

        ntraj_ent.grid(sticky='we', row=4, column=0, **pref)
        nthreads_ent.grid(sticky='we', row=4, column=1, **pref)
        mindx_ent.grid(sticky='we', row=5, column=0, **pref)
        ntrajo_ent.grid(sticky='we', row=5, column=1, **pref)
        ncopies_ent.grid(sticky='we', row=6, column=0, **pref)
        nbincopies_ent.grid(sticky='we', row=6, column=1, **pref)
        nsteps_ent.grid(sticky='we', row=7, column=0, **pref)
        westeps_ent.grid(sticky='we', row=7, column=1, **pref)
        maxnsteps_ent.grid(sticky='we', row=8, column=0, **pref)
        nsteps_per_output_ent.grid(sticky='we', row=8, column=1, **pref)
        prep_bd_but.grid(sticky='e', row=9, column=0, **pref)
        #self.status_bar.grid(sticky='e', row=10, column=0, **pref)
        grp_bdinput.columnconfigure(0, weight=1)
        grp_bdinput.columnconfigure(1, weight=1)
        page.columnconfigure(0, weight=1)

        ######################
        # Tab: BD Simulation
        ######################
        page = self.notebook.add('BD simulation')
        grp_sim = Tkinter.LabelFrame(page, text='BrownDye simulation')
        grp_sim.grid(sticky='eswn', row=0, column=0, columnspan=2, **pref)

        bkgj_cb = Tkinter.Checkbutton(grp_sim,
                                      text='Run job in background',
                                      variable=self.run_in_background,
                                      onvalue=True, offvalue=False)
        run_bd_but = Tkinter.Button(grp_sim,
                                    text="Start BrownDye simulation", command=self.runBD)
        kill_bd_but = Tkinter.Button(grp_sim,
                                     text="Stop background job", command=self.killBD)
        self.msgbar1 = Pmw.MessageBar(grp_sim,
                                      entry_width=10, entry_relief='sunken',
                                      labelpos='w',
                                      label_text='Trajectories: ')
        self.msgbar2 = Pmw.MessageBar(grp_sim,
                                      entry_width=10, entry_relief='sunken',
                                      labelpos='w',
                                      label_text='Completed / escaped / stuck events: ')
        self.msgbar3 = Pmw.MessageBar(grp_sim,
                                      entry_width=10, entry_relief='sunken',
                                      labelpos='w',
                                      label_text=('Calculated reaction rate '
                                                  'and probability: '))
        self.logtxt_ent = Pmw.ScrolledText(grp_sim, labelpos='wn',
                                           borderframe=5,
                                           vscrollmode='dynamic',
                                           hscrollmode='dynamic',
                                           text_width=100, text_height=15,
                                           text_wrap='none',
                                           text_background='#000000',
                                           text_foreground='white')

        bkgj_cb.grid(sticky='w', row=0, column=0, columnspan=2, **pref)
        run_bd_but.grid(sticky='we', row=1, column=0, **pref)
        kill_bd_but.grid(sticky='we', row=1, column=1, **pref)

        self.msgbar1.grid(sticky='we', row=2, column=0, columnspan=1, **pref)
        self.msgbar2.grid(sticky='we', row=2, column=1, columnspan=1, **pref)
        self.msgbar3.grid(sticky='we', row=3, column=0, columnspan=2, **pref)
        self.logtxt_ent.grid(sticky='we', row=4, column=0, columnspan=2, **pref)
        grp_sim.columnconfigure(0, weight=1)
        grp_sim.columnconfigure(1, weight=1)

        page.columnconfigure(0, weight=1)

        ######################
        # Tab: Analysis
        ######################
        page = self.notebook.add('Analysis')
        grp_analysis = Tkinter.LabelFrame(page, text='Analysis and Visualization')
        grp_analysis.grid(sticky='eswn', row=0, column=0, **pref)

        load_traj_ent = Pmw.EntryField(grp_analysis,
                                       label_text='Select trajectory file: ',
                                       labelpos='w',
                                       entry_textvariable=self.traj_f,
                                       entry_width=20)
        load_traj_but = Tkinter.Button(grp_analysis, text='Browse...',
                                       command=self.loadTrajectoryFile)
        analyze_but = Tkinter.Button(grp_analysis, text='Analyze',
                                     command=self.analyzeTrajectoryFile)
        self.msg_ent = Pmw.ScrolledText(grp_analysis, labelpos='wn',
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
        select_index_but = Tkinter.Button(grp_analysis, text='Select trajectory index',
                                          command=self.dialog_idx.activate)
        self.msgbar_idx = Pmw.MessageBar(grp_analysis,
                                         entry_width=20, entry_relief='sunken',
                                         labelpos='w',
                                         label_text='Selected trajectory: ')
        convert_but = Tkinter.Button(grp_analysis, text='Convert to xyz trajectory',
                                     command=self.convertTrajectoryToXYZ)
        load_xyztraj_but = Tkinter.Button(grp_analysis, text='Load xyz trajectory',
                                          command=self.loadTrajectoryFileXYZ)

        load_traj_ent.grid(sticky='we', row=0, column=0, **pref)
        load_traj_but.grid(sticky='w', row=0, column=1, **pref)
        analyze_but.grid(sticky='w', row=0, column=2, **pref)
        self.msg_ent.grid(sticky='we', row=1, column=0, columnspan=3, **pref)
        select_index_but.grid(sticky='we', row=2, column=0, **pref)
        self.msgbar_idx.grid(sticky='we', row=2, column=1, **pref)
        #traj_index_n_ent.grid(sticky='we', row=3, column=0, **pref)
        convert_but.grid(sticky='we', row=4, column=0, **pref)
        load_xyztraj_but.grid(sticky='we', row=5, column=0, **pref)
        grp_analysis.columnconfigure(0, weight=1)
        page.columnconfigure(0, weight=1)

        #############
        # Tab: About
        #############
        page = self.notebook.add('About')
        grp_about = Tkinter.LabelFrame(page, text='About BrownDyeTools Plugin for PyMOL')
        # grp_about.grid(sticky='n', row=0, column=0, columnspan=2, **pref)
        grp_about.pack(fill='both', expand=True, **pref)
        about_plugin = ("""
        This plugin provides a GUI for setting up and running Brownian 
        dynamics simulations with BrownDye.

        The plugin requires PDB2PQR, APBS and BrownDye. 
        To download and install these applications go to:

        http://www.poissonboltzmann.org/
        and
        http://browndye.ucsd.edu/


        This software is released under the terms of GNU GPL3 license.
        For more details please see the accompanying documentation.

        (c) 2016-19 National Biomedical Computation Resource
        http://nbcr.ucsd.edu/""")

        label_about = Tkinter.Label(grp_about, text=about_plugin)
        # label_about.grid(sticky='we', row=0, column=2, **pref)
        label_about.pack(fill='both', expand=True, **pref)
        self.notebook.setnaturalsize()
        return

    def getProjectDir(self):
        """Generate a random project directory name."""
        rnd_id = random.randint(10000, 99999)
        cwd = os.getcwd()
        p_dir = '%s/bd-project-%d' % (cwd, rnd_id)
        return p_dir

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
        print("::: Created project directory: %s" % self.projectDir.get())
        return

    def runCmd(self, command):
        """Generic wrapper for running shell commands."""
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             universal_newlines=True,
                             stderr=subprocess.PIPE, shell=True)
        #p.wait()
        stdout, stderr = p.communicate()
        if DEBUG > 2: print("returncode = %d" % p.returncode)
        if DEBUG > 2: print("stdout:\n%s\n" % stdout)
        if DEBUG > 2: print("stderr:\n%s\n" % stderr)
        if p.returncode > 0:
            sys.stderr.write("::: Non-zero return code!\n")
            sys.stderr.write("::: Failed command: \n\n")
            sys.stderr.write("%s \n\n" % command)
            sys.stderr.write(stderr)
        return p.returncode

    def getPDB2PQRpath(self):
        """Get PDB2PQR binary path."""
        d = tkFileDialog.askdirectory(title='PDB2PQR binary directory',
                                      initialdir='',
                                      parent=self.parent)
        self.pdb2pqr_path.set(d)
        return

    def getAPBSpath(self):
        """Get APBS binary path."""
        d = tkFileDialog.askdirectory(title='APBS binary directory',
                                      initialdir='',
                                      parent=self.parent)
        self.apbs_path.set(d)
        return

    def getBDpath(self):
        """Get BrownDye binary path."""
        d = tkFileDialog.askdirectory(title='BrownDye binary directory',
                                      initialdir='',
                                      parent=self.parent)
        self.bd_path.set(d)
        return

    def getConfigPath(self):
        """Ask for path to the configuration json file."""
        f = tkFileDialog.askopenfilename(title='Configuration file',
                                         initialdir='',
                                         parent=self.parent)
        self.config_file.set(f)
        return

    def loadConfig(self):
        """Load configuration from a file, update all affected variables."""
        configf = self.config_file.get()
        with open(configf) as f:
            data = json.load(f)

        self.pdb2pqr_path.set(data['paths']['PDB2PQR_PATH'])
        self.apbs_path.set(data['paths']['APBS_PATH'])
        self.bd_path.set(data['paths']['BD_PATH'])
        self.projectDir.set(data['paths']['ProjectDir'])

        pqr_config = data['pdb2pqr']
        for k in PQR_DEFAULTS:
            getattr(self, k).set(pqr_config[k])
        psize_config = data['psize']
        for k in PSIZE_DEFAULTS:
            getattr(self, k).set(psize_config[k])
        apbs_config = data['apbs']
        for k in APBS_DEFAULTS:
            try:
                getattr(self, k).set(apbs_config[k])
            except AttributeError:
                getattr(self, k)[0].set(apbs_config[k][0])
                getattr(self, k)[1].set(apbs_config[k][1])
        bd_config = data["browndye"]
        for k in BD_DEFAULTS:
            try:
                getattr(self, k).set(bd_config[k])
            except AttributeError:
                getattr(self, k)[0].set(bd_config[k][0])
                getattr(self, k)[1].set(bd_config[k][1])
        return

    def saveConfig(self):
        """Save configuration information to a json file."""
        configf = self.config_file.get()
        data = {}
        paths_config = {
            "PDB2PQR_PATH": self.pdb2pqr_path.get(),
            "APBS_PATH": self.apbs_path.get(),
            "BD_PATH": self.bd_path.get(),
            "ProjectDir": self.projectDir.get()
        }
        pqr_config = {}
        for k in PQR_DEFAULTS:
            pqr_config[k] = getattr(self, k).get()
        psize_config = {}
        for k in PSIZE_DEFAULTS:
            psize_config[k] = getattr(self, k).get()
        apbs_config = {}
        for k in APBS_DEFAULTS:
            try:
                apbs_config[k] = getattr(self, k).get()
            except AttributeError:
                apbs_config[k] = [getattr(self, k)[0].get(),
                                  getattr(self, k)[1].get()]
        bd_config = {}
        for k in BD_DEFAULTS:
            try:
                bd_config[k] = getattr(self, k).get()
            except AttributeError:
                bd_config[k] = [getattr(self, k)[0].get(),
                                getattr(self, k)[1].get()]
        bdtools_config = {
            "version": __version__,
            "date": datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        }
        data['paths'] = paths_config
        data['pdb2pqr'] = pqr_config
        data['psize'] = psize_config
        data['apbs'] = apbs_config
        data['browndye'] = bd_config
        data['_bdtools'] = bdtools_config
        with open(configf, 'w') as f:
            json.dump(data, f, indent=4, sort_keys=True)
        return

    def getPDBMol(self, n):
        """Get molecule 0/1 PDB filename."""
        fname = tkFileDialog.askopenfilename(title='Select PDB File',
                                             initialdir='',
                                             filetypes=[
                                                 ('pdb files', '*.pdb *.ent'),
                                                 ('all files', '*')],
                                             parent=self.parent)
        self.mol[n].set(fname)
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

    def getPQRMol(self, n):
        """Get molecule 0/1 PQR filename."""
        fname = tkFileDialog.askopenfilename(title='Select PQR File',
                                             initialdir='',
                                             filetypes=[
                                                 ('pqr files', '*.pqr'),
                                                 ('all files', '*')],
                                             parent=self.parent)
        if len(fname) > 0:
            self.pqr[n].set(fname)
            target_f = '%s/%s.pqr' % (self.projectDir.get(), MOL[n])
            if fname != target_f:
                if os.path.isfile(target_f): os.remove(target_f)
                shutil.copyfile(self.pqr[n].get(), target_f)
        return

    def getSizemol(self, n):
        """Calculate APBS grid dimensions for molecule 0/1."""
        pqr_fname = '%s.pqr' % MOL[n]
        if not os.path.isfile(pqr_fname):
            print("::: %s does not exist!" % pqr_fname)
            return
        psize = Psize(self)
        psize.runPsize(pqr_fname)
        #print(psize.getCharge())
        grid_points = psize.getFineGridPoints()
        cglen = psize.getCoarseGridDims()
        fglen = psize.getFineGridDims()
        for i in range(3):
            self.dime[n][i].set(grid_points[i])
            self.cglen[n][i].set(cglen[i])
            self.fglen[n][i].set(fglen[i])
        return

    def getContacts(self):
        """Get contacts file."""
        fname = tkFileDialog.askopenfilename(
            title='Contacts File', initialdir='',
            filetypes=[('xml files', '*.xml'), ('all files', '*')],
            parent=self.parent)
        self.contacts_f.set(fname)
        return

    def pdb2pqr(self):
        """Convert PDB to PQR."""
        if not os.path.exists(self.projectDir.get()):
            print("::: Project directory does not exist!")
            print("::: You need to set it first.")
            return
        for i in range(2):
            target_f = '%s.pdb' % MOL[i]
            if self.mol_object[i].get() == 'None':
                # if not filecmp.cmp(self.mol0.get(), target_f):
                try:
                    shutil.copyfile(self.mol[i].get(), target_f)
                except:
                    e = sys.exc_info()[0]
                    print(e)
                    print("::: Creating of %s failed!" % target_f)
                    print("::: The pqr file already exists.")
                    return
            else:
                pymol.cmd.save(filename=target_f,
                               selection=self.mol_object[i].get())
        if self.pqr_assign_only.get():
            assign_only = '--assign-only'
            use_propka = ''
        if self.pqr_use_propka.get():
            assign_only = ''
            use_propka = ('--ph-calc-method=propka --with-ph=%s --drop-water'
                          % self.pqr_ph.get())
            self.pqr_ff.set('parse')
            print("::: Using PROPKA to assign protonation states and "
                  "optimizing the structure.\n"
                  "The force field is set to PARSE.")

        pqr_options = ('%s %s %s --ff=%s' %
                       (assign_only, use_propka, self.pqr_opts.get(),
                        self.pqr_ff.get()))
        pdb2pqr_exe = ('%s/%s' % (self.pdb2pqr_path.get(), PDB2PQR_EXE))
        if not self.checkExe(pdb2pqr_exe): return
        for i in MOL:
            command = ('%s %s %s.pdb %s.pqr' %
                       (pdb2pqr_exe, pqr_options, i, i))
            if DEBUG > 2: print(command)
            print("::: Running pdb2pqr on %s ..." % i)
            self.status_bar.message('state', 'Busy: Running pdb2pqr. Please wait ...')
            self.dialog.update()
            rc = self.runCmd(command)
            if rc == 0:
                print("::: Done.")
            else:
                print("::: Command failed: %s" % command)
        self.status_bar.message('state', 'Idle')
        return

    def runAPBS2(self):
        """Run APBS calculations on molecule 0 and molecule 1.
        Non-threaded version.
        """
        apbs_template = """
# APBS template for BrownDye grids
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
    ion charge %d conc %f radius %f
    ion charge %d conc %f radius %f
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
        apbs_exe = '%s/%s' % (self.apbs_path.get(), APBS_EXE)
        if not self.checkExe(apbs_exe): return
        for i in range(2):
            pqr_filename = '%s.pqr' % MOL[i]
            grid_points = [self.dime[i][x].get() for x in range(3)]
            cglen = [self.cglen[i][x].get() for x in range(3)]
            fglen = [self.fglen[i][x].get() for x in range(3)]
            if grid_points[0] == 0 or grid_points[1] == 0 or grid_points[2] == 0:
                print("::: %s - no grid points defined!" % MOL[i])
                return
            dx_filename = MOL[i]
            fout = '%s.in' % MOL[i]
            with open(fout, "w") as f:
                f.write((apbs_template %
                         (pqr_filename,
                          grid_points[0], grid_points[1], grid_points[2],
                          cglen[0], cglen[1], cglen[2],
                          fglen[0], fglen[1], fglen[2],
                          self.apbs_mode.get(), self.apbs_bcfl.get(),
                          self.apbs_ion_charge[0].get(),
                          self.apbs_ion_conc[0].get(), self.apbs_ion_radius[0].get(),
                          self.apbs_ion_charge[1].get(),
                          self.apbs_ion_conc[1].get(), self.apbs_ion_radius[1].get(),
                          self.apbs_pdie.get(),
                          self.apbs_sdie.get(),
                          self.apbs_chgm.get(), self.apbs_sdens.get(),
                          self.apbs_srfm.get(), self.apbs_srad.get(),
                          self.apbs_swin.get(), self.apbs_temp.get(),
                          dx_filename)))

            command = 'MCSH_HOME=. %s %s.in' % (apbs_exe, MOL[i])
            if DEBUG > 2:
                print(command)
                print(grid_points)
                print(cglen, fglen)
            print("::: Running apbs on %s ... " % MOL[i])
            gmem = 200.0 * grid_points[0] * grid_points[1] * grid_points[2] / 1024 / 1024
            print("::: Estimated memory requirements: %.3f MB" % gmem)
            self.status_bar.message('state', 'Busy: Running apbs. Please wait ...')
            self.dialog.update()
            rc = self.runCmd(command)
            if rc == 0:
                iomc = '%s-io.mc' % MOL[i]
                if os.path.isfile(iomc): os.remove(iomc)
                shutil.copyfile('io.mc', iomc)
                print("::: Done.")
            else:
                print("::: Failed: %s" % command)
        self.status_bar.message('state', 'Idle')
        return

    def runAPBS(self):
        """Run APBS calculations on molecule 0 and molecule 1.
        Threaded version.
        """
        apbs_template = """
# APBS template for BrownDye grids
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
    ion charge %d conc %f radius %f
    ion charge %d conc %f radius %f
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
        apbs_exe = '%s/%s' % (self.apbs_path.get(), APBS_EXE)
        if not self.checkExe(apbs_exe): return
        for i in range(2):
            pqr_filename = '%s.pqr' % MOL[i]
            grid_points = [self.dime[i][x].get() for x in range(3)]
            cglen = [self.cglen[i][x].get() for x in range(3)]
            fglen = [self.fglen[i][x].get() for x in range(3)]
            if grid_points[0] == 0 or grid_points[1] == 0 or grid_points[2] == 0:
                print("::: %s - no grid points defined!" % MOL[i])
                return
            dx_filename = MOL[i]
            fout = '%s.in' % MOL[i]
            with open(fout, "w") as f:
                f.write((apbs_template %
                         (pqr_filename,
                          grid_points[0], grid_points[1], grid_points[2],
                          cglen[0], cglen[1], cglen[2],
                          fglen[0], fglen[1], fglen[2],
                          self.apbs_mode.get(), self.apbs_bcfl.get(),
                          self.apbs_ion_charge[0].get(),
                          self.apbs_ion_conc[0].get(), self.apbs_ion_radius[0].get(),
                          self.apbs_ion_charge[1].get(),
                          self.apbs_ion_conc[1].get(), self.apbs_ion_radius[1].get(),
                          self.apbs_pdie.get(),
                          self.apbs_sdie.get(),
                          self.apbs_chgm.get(), self.apbs_sdens.get(),
                          self.apbs_srfm.get(), self.apbs_srad.get(),
                          self.apbs_swin.get(), self.apbs_temp.get(),
                          dx_filename)))

            command = 'MCSH_HOME=. %s %s.in' % (apbs_exe, MOL[i])
            if DEBUG > 2:
                print(command)
                print(grid_points)
                print(cglen, fglen)
            print("::: Prepping apbs run for %s ... " % MOL[i])
            gmem = 200.0 * grid_points[0] * grid_points[1] * grid_points[2] / 1024 / 1024
            print("::: Estimated memory requirements: %.3f MB" % gmem)
        thread = APBSRunner(self, apbs_exe)
        thread.start()
        return

    def getDebyeLength(self):
        dl = [[], []]
        for i in range(2):
            fname = '%s-io.mc' % MOL[i]
            if DEBUG > 2: print("Parsing %s for Debye length ..." % fname)
            if not os.path.isfile(fname):
                print("::: File %s does not exist!" % fname)
                print("::: You need to run APBS first.")
                return
            with open(fname) as f:
                for line in f:
                    l = re.findall(r'Debye length = .*', line)
                    if l:
                        dl[i].append(l[0].split(' ')[3])
            if len(dl[i]) == 0:
                print("::: No Debye length found in %s for %s!" %
                      (fname, MOL[i]))
                return
        if DEBUG > 2: print("Debye lengths: %s %s" %
                            (dl[0][-1], dl[1][-1]))
        if dl[0][-1] != dl[1][-1]:
            print("::: Debye lengths don't match!")
            print("%s / %s" % (dl[0][-1], dl[1][-1]))
            return
        self.debyel[0].set(float(dl[0][-1]))
        return

    def runPqr2xml(self):
        """Run pqr2xml on mol0 and mol1 PQR files. """
        pqr2xml_exe = '%s/pqr2xml' % self.bd_path.get()
        if not self.checkExe(pqr2xml_exe): return
        for i in MOL:
            command = ('%s < %s.pqr > %s-atoms.pqrxml'
                       % (pqr2xml_exe, i, i))
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
        <atom> OE1 </atom> <residue> GLN </residue>
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
        <atom> OE1 </atom> <residue> GLN </residue>
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
                   '-nonred -mol0 %s-atoms.pqrxml -mol1 %s-atoms.pqrxml '
                   '-ctypes %s  -dist %f > %s-%s-rxn-pairs.xml'
                   % (self.bd_path.get(), MOL[0], MOL[1],
                      self.contacts_f.get(),
                      self.rxn_distance.get(), MOL[0], MOL[1]))
        if DEBUG > 2: print(command)
        print("::: Running make_rxn_pairs ...")
        self.status_bar.message('state', 'Busy: Processing rxn criteria. Please wait ...')
        self.dialog.update()
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
        command = ('%s/make_rxn_file '
                   '-pairs %s-%s-rxn-pairs.xml -distance %f '
                   '-nneeded %d > %s-%s-rxns.xml'
                   % (self.bd_path.get(), MOL[0], MOL[1],
                      self.rxn_distance.get(),
                      self.npairs.get(), MOL[0], MOL[1]))
        if DEBUG > 2: print(command)
        print("::: Running make_rxn_file ...")
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
        self.status_bar.message('state', 'Idle')
        return

    def prepBD(self):
        """Create BrownDye input file and run bd_top."""
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
    <prefix> %s </prefix>
    <atoms>  %s-atoms.pqrxml </atoms>
    <apbs-grids>
       <grid> %s.dx </grid>
    </apbs-grids>
    <solute-dielectric> %f </solute-dielectric>
  </molecule0>
  <molecule1>
    <prefix> %s </prefix>
    <atoms>  %s-atoms.pqrxml </atoms>
    <all-in-surface> false </all-in-surface>
    <apbs-grids>
       <grid> %s.dx </grid>
    </apbs-grids>
    <solute-dielectric> %f </solute-dielectric>
  </molecule1>
  <include-desolvation-forces> true </include-desolvation-forces>
  <time-step-tolerances>
    <minimum-dx> %f </minimum-dx>
  </time-step-tolerances>
  <reactions> %s-%s-rxns.xml </reactions>
  <seed> 11111113 </seed>
  <n-trajectories-per-output> %d </n-trajectories-per-output>
  <n-copies> %d </n-copies>
  <n-bin-copies> %d </n-bin-copies>
  <n-steps> %d </n-steps>
  <n-we-steps-per-output> %d </n-we-steps-per-output>
  <max-n-steps> %d </max-n-steps>
  <trajectory-file> traj </trajectory-file>
  <n-steps-per-output> %d </n-steps-per-output>
</root>
"""
        bdtop_input = 'input.xml'
        if self.debyel[0].get() == 0.0:
            print("::: Invalid Debye length (%s)!" %
                  self.debyel[0].get())
            print("::: Acquire Debye length first")
            return
        with open(bdtop_input, "w") as f:
            f.write((nam_simulation_template %
                     (self.sdie.get(), self.debyel[0].get(),
                      self.ntraj.get(), self.nthreads.get(),
                      MOL[0], MOL[0], MOL[0], self.pdie[0].get(),
                      MOL[1], MOL[1], MOL[1], self.pdie[1].get(),
                      self.mindx.get(), MOL[0], MOL[1],
                      self.ntrajo.get(),
                      self.ncopies.get(), self.nbincopies.get(),
                      self.nsteps.get(), self.westeps.get(),
                      self.maxnsteps.get(),
                      self.nsteps_per_output.get())))
        for i in MOL:
            dxfile = '%s.dx' % i
            if not os.path.isfile(dxfile):
                print("::: %s DX file does not exist!" % dxfile)
                print("::: Run APBS first to create it.")
                return
        command = ('PATH=%s:${PATH} %s %s'
                   % (self.bd_path.get(), BDTOP_EXE, bdtop_input))
        if DEBUG > 2: print(command)
        print("::: Running bd_top (this will take a couple of minutes) ...")
        thread = BDTopRunner(self, command)
        thread.start()
        return

    def runRxnCrit(self):
        self.runPqr2xml()
        self.makeRxnCriteria()
        return

    def runBD(self):
        """Start BrownDye simulation either in foreground or background."""
        print("::: Starting BrownDye simulation ...")
        nam_simulation_exe = ('%s/%s' %
                              (self.bd_path.get(), NAM_SIMULATION_EXE))
        if not self.checkExe(nam_simulation_exe): return
        command = ('%s %s-%s-simulation.xml'
                   % (nam_simulation_exe, MOL[0], MOL[1]))
        if self.run_in_background.get():
            p = subprocess.Popen(" nohup %s" % command,
                                 universal_newlines=True, shell=True)
            global JOBPID
            JOBPID = p.pid
            self.notebook.selectpage('BD simulation')
            self.logtxt_ent.insert(
                'end', "::: Starting BrownDye simulation in background.\n")
            self.logtxt_ent.insert('end', "::: Job PID: %d \n" % JOBPID)
        else:
            thread = BDRunner(self, command)
            thread.start()
            self.notebook.selectpage('BD simulation')
            self.logtxt_ent.insert('end', "::: Starting BrownDye simulation ...\n")
            time.sleep(1)
            tm = MonitorThread(self, thread, 3600)
            tm.start()
        return

    def killBD(self):
        """Terminate background BrownDye job, if exists."""
        try:
            command = 'pkill -9 -P %d' % JOBPID
            s = subprocess.Popen(command, shell=True,
                                 universal_newlines=True,
                                 stdout=subprocess.PIPE).stdout.read()
            if DEBUG > 1: print(command, s)
            command = 'kill -9 %d' % JOBPID
            s = subprocess.Popen(command, shell=True,
                                 universal_newlines=True,
                                 stdout=subprocess.PIPE).stdout.read()
            if DEBUG > 1: print(command, s)
            print("::: Background job (PID %d) terminated." % JOBPID)
        except AttributeError:
            print("::: No background job running!")
        return

    def loadTrajectoryFile(self):
        """Load trajectory file."""
        fname = tkFileDialog.askopenfilename(title='Trajectory File',
                                             initialdir='',
                                             filetypes=[('xml files', '*.xml'),
                                                        ('all files', '*')],
                                             parent=self.parent)
        self.traj_f.set(fname)
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
        self.msg_ent.insert('end', "%s" % logline)
        if reacted == 0:
            print("::: No association events found in %s trajectory file."
                  % self.traj_f.get())
            return
        command = ('%s/process_trajectories -traj %s '
                   '-index %s.index.xml -srxn association'
                   % (self.bd_path.get(), self.traj_f.get(),
                      traj_f_base))
        self.status_bar.message('state', 'Busy: Processing trajectory. Please wait ...')
        self.dialog.update()
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             universal_newlines=True,
                             stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        if p.returncode > 0:
            print("::: Error processing %s trajectory file." % self.traj_f.get())
            return
        t = etree.fromstring(stdout)
        tr = t.findall('.//trajectory/number')
        traj_index = [int(tr[x].text.strip()) for x in range(len(tr))]
        logline = (
            '%d association event trajectories found (index numbers: %s)\n'
            % (len(traj_index), str(traj_index)))
        self.msg_ent.insert('end', "%s" % logline)
        # number of frames
        self.dialog_idx.clear()
        for i in traj_index:
            with open(self.traj_f.get(), 'r') as f:
                d = etree.parse(f)
                # xpath_str = 'trajectory[n-traj=" %d "]//s/n' % i
                xpath_str = 'trajectory[n-traj=" %d "]//s' % i
                j = d.findall(xpath_str)
                logline = ('trajectory %d: approx. %s frames\n'
                           % (i, len(j) * 20))
                # % (i, j[0].text.strip()))
                self.msg_ent.insert('end', "%s" % logline)
                self.dialog_idx.insert('end', i)
        self.status_bar.message('state', 'Idle')
        return

    def selectTrajIndex(self, result):
        """Select trajectory index number from Dialog window"""
        sel = self.dialog_idx.getcurselection()
        if len(sel) > 0:
            if DEBUG > 1: print("::: Selection: %s" % sel)
            self.traj_index_n.set(sel[0])
            # self.msgbar_idx.insert(sel)
            self.msgbar_idx.message('state', str(sel[0]))
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
        command = ('%s/xyz_trajectory -mol0 %s-atoms.pqrxml '
                   '-mol1 %s-atoms.pqrxml -trajf %s > %s'
                   % (self.bd_path.get(), MOL[0], MOL[1], ofile,
                      self.xyz_ofile))
        print("::: Converting %s to XYZ trajectory ..." % ofile)
        self.status_bar.message('state', 'Busy: Processing trajectory. Please wait ...')
        self.dialog.update()
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
            return
        self.status_bar.message('state', 'Idle')
        return

    def convertTrajectoryToPQR(self):
        """Convert XML trajectory to PQR format."""
        return

    def loadTrajectoryFileXYZ(self):
        """Load XYZ trajectory to pymol."""
        xyz_trajectory_object = 'trajectory-%d' % self.traj_index_n.get()
        print("::: Loading XYZ trajectory ...")
        self.status_bar.message('state', 'Busy: Loading trajectory. Please wait ...')
        self.dialog.update()
        try:
            pymol.cmd.load(self.xyz_ofile, xyz_trajectory_object)
            print("::: Done.")
        except:
            e = sys.exc_info()[0]
            print("::: Loading xyz trajectory failed!")
            print("::: Error: %s" % e)
        self.status_bar.message('state', 'Idle')
        return

    def checkExe(self, command):
        """Check if command exists and is executable."""
        if os.path.isfile(command) and os.access(command, os.X_OK):
            is_exe = True
        else:
            is_exe = False
            print("::: Can't find %s!" % command)
        return is_exe

    def exitBDPlugin(self, result):
        """Quit BrownDye Tools plugin."""
        #if tkMessageBox.askokcancel("Really quit?",
        #                            "Exit BrownDye Tools Plugin now?") == 0:
        #    return 0
        print("Exiting BrownDye Tools Plugin ...")
        if __name__ == '__main__':
            self.parent.destroy()
        else:
            self.dialog.withdraw()
        print("Done.")
        return

class APBSRunner(Thread):
    """Run APBS in a thread."""
    def __init__(self, my_inst, apbs_exe):
        Thread.__init__(self)
        self.apbs_exe = apbs_exe
        self.status_bar = my_inst.status_bar
        if DEBUG > 2: print("%s" % (self.apbs_exe))
        self.work_dir = my_inst.projectDir.get()
        if DEBUG > 2: print("Work directory: %s" % (self.work_dir))
        return

    def run(self):
        print("::: Project directory: %s" % (self.work_dir))
        os.chdir(self.work_dir)
        for i in range(2):
            command = 'MCSH_HOME=. %s %s.in' % (self.apbs_exe, MOL[i])
            if DEBUG > 2: print(command)
            self.status_bar.message('state',
                                    ('Busy: Running apbs on %s. Please wait ...' % MOL[i]))
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True,
                                 shell=True, bufsize=1)
            global JOBPID
            JOBPID = p.pid
            self.pid = p.pid
            if DEBUG > 1: print('JOBPID: %d' % self.pid)
            for line in iter(p.stdout.readline, ''):
                print(line)
                ##self.page.insert('end', "%s" % line)
            stdout, stderr = p.communicate()
            self.outlog = stdout
            self.status = p.returncode
            self.status_bar.message('state', ('Finished apbs on %s.' % MOL[i]))
            if self.status == 0:
                iomc = '%s-io.mc' % MOL[i]
                if os.path.isfile(iomc): os.remove(iomc)
                shutil.copyfile('io.mc', iomc)
                print("::: Done.")
            else:
                print("::: Failed: %s" % command)
        self.status_bar.message('state', 'Idle')
        return

    def kill(self):
        os.system('kill -9 %d' % self.pid)
        return


class BDTopRunner(Thread):
    """bd_top thread management class."""
    def __init__(self, my_inst, command):
        Thread.__init__(self)
        self.status_bar = my_inst.status_bar
        self.command = command
        if DEBUG > 2: print("%s" % (command))
        return

    def run(self):
        self.status_bar.message('state', 'Busy: Running bd_top. Please wait ...')
        p = subprocess.Popen(self.command, stdout=subprocess.PIPE,
                             universal_newlines=True,
                             stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        self.outlog = stdout
        self.status = p.returncode
        try:
            self.status_bar.message('state', 'Finished bd_top ...')
            time.sleep(3)
            self.status_bar.message('state', 'Idle')
            print(stdout, stderr, p.returncode)
        except:
            print("::: Thread error!")
        return


class BDRunner(Thread):
    """BD thread management class."""
    def __init__(self, my_inst, command):
        Thread.__init__(self)
        self.page = my_inst.logtxt_ent
        self.status_bar = my_inst.status_bar
        self.command = command
        if DEBUG > 2: print("%s" % (command))
        self.work_dir = my_inst.projectDir.get()
        if DEBUG > 2: print("Work directory: %s" % (self.work_dir))
        return

    def run(self):
        print("::: Project directory: %s" % (self.work_dir))
        # current_dir = os.getcwd()
        os.chdir(self.work_dir)
        #self.page.yview('moveto', 1.0)
        self.status_bar.message('state', 'Busy: Running simulation. Please wait ...')
        p = subprocess.Popen(self.command, stdout=subprocess.PIPE,
                             universal_newlines=True,
                             stderr=subprocess.PIPE, shell=True, bufsize=1)
        global JOBPID
        JOBPID = p.pid
        self.pid = p.pid
        if DEBUG > 1: print("JOBPID: %d" % self.pid)
        for line in iter(p.stdout.readline, ''):
            print(line,)
            self.page.insert('end', "%s" % line)
        stdout, stderr = p.communicate()
        self.outlog = stdout
        self.status = p.returncode
        self.status_bar.message('state', 'Idle')
        try:
            self.page.insert('end', "%s %s" % (stdout, stderr))
            # self.page.yview('moveto', 1.0)
        except:
            print("::: Thread error!")
        return

    def kill(self):
        os.system('kill -9 %d' % self.pid)
        return

class MonitorThread(Thread):
    """Monitor runing BrownDye job and print out progress information."""
    def __init__(self, bd_inst, thrd, timeout):
        Thread.__init__(self)
        self.totntraj = bd_inst.ntraj.get()
        self.page = bd_inst.logtxt_ent
        self.msgbar1 = bd_inst.msgbar1
        self.msgbar2 = bd_inst.msgbar2
        self.msgbar3 = bd_inst.msgbar3
        self.bd_path = bd_inst.bd_path.get()
        self.mythread = thrd
        self.timeout = timeout

    def run(self):
        seconds = 0
        #readcount = 0
        #log_fileok = False
        results_file = 'results.xml'
        results_fileok = False
        while self.mythread.is_alive() and seconds < self.timeout:
            time.sleep(3)
            try:
                if os.stat(results_file).st_size > 0:
                    results_fileok = True
            except:
                if DEBUG > 2: print("::: Waiting for %s ..." % results_file)
                results_fileok = False
            if results_fileok:
                with open(results_file, 'r') as f:
                    results = etree.parse(f)
                    ntraj = results.findall('.//reactions/n-trajectories')[0].text.strip()
                    stuck = results.findall('.//reactions/stuck')[0].text.strip()
                    escaped = results.findall('.//reactions/escaped')[0].text.strip()
                    ncompleted = results.findall('.//reactions/completed/name')[0].text.strip()
                    completed = results.findall('.//reactions/completed/n')[0].text.strip()
                    mymsg = '%s out of %d' % (ntraj, self.totntraj)
                    self.msgbar1.message('state', mymsg)
                    mymsg = '%s / %s / %s' % (completed, escaped, stuck)
                    self.msgbar2.message('state', mymsg)
                command = ('cat %s | %s/compute_rate_constant'
                           % (results_file, self.bd_path))
                p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                     universal_newlines=True,
                                     stderr=subprocess.PIPE, shell=True)
                # p.wait()
                stdout, stderr = p.communicate()
                try:
                    rates = etree.fromstring(stdout)
                    rate_constant = rates.findall('.//rate-constant/mean')[0].text.strip()
                    rxn_probability = rates.findall('.//reaction-probability/mean')[0].text.strip()
                    mymsg = ('%s / %s' % (rate_constant, rxn_probability))
                    self.msgbar3.message('state', mymsg)
                except:
                    e = sys.exc_info()[0]
                    print("::: rates etree parse XMLSyntaxError")
                    print("::: Error: %s" % e)
                    raise

        time.sleep(2)
        if DEBUG > 2: print(self.mythread.is_alive())
        print("::: BrownDye simulation finished.")
        self.page.insert('end', "::: BrownDye simulation finished\n")
        return


class Psize(object):
    """Calculate grid size dimensions for a pqr molecule.

    This is based on pdb2pqr version of psize. All licensing info applies.

    Note: CFAC and FADD defaults changed to accomodate BrownDye requirements.
    CFAC 1.7 -> 3.0
    FADD 20 -> 50

    This significantly increases the grid size and thus memory requirements, however.
    """
    def __init__(self, bd_instance):
        self.constants = {"GMEMFAC": 200, "GMEMCEIL": 400, "OFAC": 0.1,
                          "REDFAC": 0.25, "TFAC_ALPHA": 9e-5,
                          "TFAC_XEON": 3e-4, "TFAC_SPARC": 5e-4}

        self.constants['CFAC'] = PSIZE_DEFAULTS['cfac']
        self.constants['FADD'] = bd_instance.fadd.get()
        self.constants['SPACE'] = bd_instance.gspace.get()
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
        with open(filename, "r") as f:
            self.parseLines(f.readlines())

    def parseLines(self, lines):
        """Parse the lines."""
        for line in lines:
            if line.find("ATOM") == 0:
                ## subline = line[30:].replace("-", " -")
                ## words = subline.split ## this is a hack
                ## adhering to lovely PDB format definition (fixed space)
                words = (line[30:38], line[38:46], line[46:54], line[54:63],
                         line[63:69], line[72:76], line[76:78])
                # rok: shouldn't the following be a map?
                #if len(list(filter(str.strip, words))) < 4:
                # like this?
                #list(map(lambda x: x.strip(), words))
                # anyway, the following should work just as well
                if len(words) < 4:
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
            elif line.find("HETATM") == 0:
                self.gothet = self.gothet + 1

    def setConstant(self, name, value):
        """Set a constant to a value; returns 0 if constant not found."""
        try:
            self.constants[name] = value
            return 1
        except KeyError:
            return 0

    def getConstant(self, name):
        """Get a constant value; raises KeyError if constant not found."""
        return self.constants[name]

    def setLength(self, maxlen, minlen):
        """Compute molecule dimensions."""
        for i in range(3):
            self.olen[i] = maxlen[i] - minlen[i]
            if self.olen[i] < 0.1:
                self.olen[i] = 0.1
        return self.olen

    def setCoarseGridDims(self, olen):
        """Compute coarse mesh dimensions."""
        for i in range(3):
            self.clen[i] = self.constants["CFAC"] * olen[i]
        return self.clen

    def setFineGridDims(self, olen, clen):
        """Compute fine mesh dimensions."""
        for i in range(3):
            self.flen[i] = olen[i] + self.constants["FADD"]
            if self.flen[i] > clen[i]:
                self.flen[i] = clen[i]
        return self.flen

    def setCenter(self, maxlen, minlen):
        """Compute molecule center."""
        for i in range(3):
            self.cen[i] = (maxlen[i] + minlen[i]) / 2
        return self.cen

    def setFineGridPoints(self, flen):
        """Compute mesh grid points, assuming 4 levels in MG hierarchy."""
        tn = [0, 0, 0]
        for i in range(3):
            tn[i] = int(flen[i]/self.constants["SPACE"] + 0.5)
            self.n[i] = 32*(int((tn[i] - 1) / 32.0 + 0.5)) + 1
            if self.n[i] < 33:
                self.n[i] = 33
        return self.n

    def setAll(self):
        """Set up all of the things calculated individually above."""
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
        """Parse input PQR file and set parameters."""
        self.parseInput(filename)
        self.setAll()

#############################################
#
# Create root window for standalone GUI.
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
