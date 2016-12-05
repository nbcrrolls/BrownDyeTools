#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# Last modified: 2016-12-05 13:49:39
#
'''BrownDye Tools plugin for Pymol

For documentation see: https://github.com/rokdev/BrownDyeTools

Author : Robert Konecny <rok@ucsd.edu>
Release date: November 2016
License: GNU General Public License version 3

Copyright 2016 Robert Konecny, NBCR

This is free software, licensed under the GNU General Public License
version 3. You should have received a copy of the GNU General Public
License along with this program.
If not, see <http://www.gnu.org/licenses/>.
'''

from __future__ import print_function
import os, subprocess
import sys, string, filecmp, shutil, re
import random
import time
#import tkSimpleDialog
import tkMessageBox
import tkFileDialog
import Tkinter
import Pmw
from threading import Thread
from lxml import etree

DEBUG = 5

__version__ = '0.1.0'
__author__ = 'Robert Konecny <rok@ucsd.edu>'

PDB2PQR_PATH = None
APBS_PATH = None
BD_PATH = None
DEFAULT_CONTACTS_FILE = 'protein-protein-contacts-default.xml'
MOL0 = 'mol0'
MOL1 = 'mol1'
APBS_EXE = 'apbs'
PDB2PQR_EXE = 'pdb2pqr'

pqr_defaults = {
    'force_field': 'parse',
    'options': '--apbs-input',
}

psize_defaults = {
    'cfac': 3.0,
    'fadd': 50.0,
    'space': 0.5,
}
apbs_defaults = {
    'mode': 'lpbe',
    'pdie': 4.0,
    'sdie': 78.0,
    'srfm': 'smol',
    'chgm': 'spl2',
    'bcfl': 'sdh',
    'sdens': 10.,
    'swin': 0.3,
    'srad': 1.4,
    'temp': 298.15,
    'ion_charge': [1, -1],
    'ion_conc': [0.15, 0.15],
    'ion_radius': [1.0, 1.0],
}

bd_defaults = {
    'ntraj': 100,
    'nthreads': 1,
    'mindx': 0.2,
    'sdie': apbs_defaults['sdie'],
    'pdie0': apbs_defaults['pdie'],
    'pdie1': apbs_defaults['pdie'],
    'debyel': 0.0,
    'ntrajo': 1,
    'ncopies': 200,
    'nbincopies': 200,
    'nsteps': 10000,
    'westeps': 10000,
    'maxnsteps': 10000,
    'rxn_distance': 5.0,
    'npairs': 3,
}

if DEBUG > 0:
    PDB2PQR_PATH = '/opt/pdb2pqr-linux-bin64-2.1.0/'
    APBS_PATH = '/export1/Builds/SandBox-2016.10.11-14.05/bin/'
    BD_PATH = '/home/rok/BrownianDynamics/browndye-2016.4.14/bin/'

if "PDB2PQR_PATH" in os.environ: PDB2PQR_PATH = os.environ["PDB2PQR_PATH"]
if "APBS_PATH" in os.environ: APBS_PATH = os.environ["APBS_PATH"]
if "BD_PATH" in os.environ: BD_PATH = os.environ["BD_PATH"]

JOBPID = None

def __init__(self):
    """BrownDye plugin for PyMol."""
    self.menuBar.addmenuitem('Plugin', 'command',
                             'BrownDye Tools', label='BrownDye Tools',
                             command=lambda s=self: BDPlugin(s))


class DummyPymol(object):
    """Dummy pymol class when running standalone GUI."""
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
        Tkinter.Grid.rowconfigure(self.parent, 0, weight=1)
        Tkinter.Grid.columnconfigure(self.parent, 0, weight=1)
        self.createGUI()
        
    def createGUI(self):
        """The main GUI class - sets up all GUI elements."""
        self.dialog = Pmw.Dialog(self.parent, buttons=('Exit',),
                                 title='BrownDye Plugin for PyMOL',
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

        # parameters used by pdb2pqr
        self.mol0 = Tkinter.StringVar()
        self.mol1 = Tkinter.StringVar()
        self.mol0.set(None)
        self.mol1.set(None)
        self.pqr0 = Tkinter.StringVar()
        self.pqr1 = Tkinter.StringVar()
        self.pqr0.set(None)
        self.pqr1.set(None)
        self.mol0_object = Tkinter.StringVar()
        self.mol1_object = Tkinter.StringVar()
        self.mol0_object.set(None)
        self.mol1_object.set(None)
        self.pdb2pqr_opt = Tkinter.StringVar()
        self.pdb2pqr_opt.set(pqr_defaults['options'])
        self.pqr_assign_only = Tkinter.BooleanVar()
        self.pqr_assign_only.set(True)
        self.pqr_ff = Tkinter.StringVar()
        self.pqr_ff.set(pqr_defaults['force_field'])

        # APBS parameters and defaults
        self.gspace = Tkinter.DoubleVar()
        self.gspace.set(psize_defaults['space'])
        self.fadd = Tkinter.DoubleVar()
        self.fadd.set(psize_defaults['fadd'])
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
        self.apbs_mode.set(apbs_defaults['mode'])
        self.bcfl = Tkinter.StringVar()
        self.bcfl.set(apbs_defaults['bcfl'])
        self.ion_charge = [Tkinter.IntVar() for _ in range(2)]
        self.ion_conc = [Tkinter.DoubleVar() for _ in range(2)]
        self.ion_rad = [Tkinter.DoubleVar() for _ in range(2)]
        self.ion_charge[0].set(apbs_defaults['ion_charge'][0])
        self.ion_charge[1].set(apbs_defaults['ion_charge'][1])
        [self.ion_conc[x].set(apbs_defaults['ion_conc'][x]) for x in range(2)]
        [self.ion_rad[x].set(apbs_defaults['ion_radius'][x]) for x in range(2)]

        self.interior_dielectric = Tkinter.DoubleVar()
        self.interior_dielectric.set(apbs_defaults['pdie'])
        self.solvent_dielectric = Tkinter.DoubleVar()
        self.solvent_dielectric.set(apbs_defaults['sdie'])
        self.chgm = Tkinter.StringVar()
        self.chgm.set(apbs_defaults['chgm'])
        self.sdens = Tkinter.DoubleVar()
        self.sdens.set(apbs_defaults['sdens'])
        self.swin = Tkinter.DoubleVar()
        self.swin.set(apbs_defaults['swin'])
        self.srfm = Tkinter.StringVar()
        self.srfm.set(apbs_defaults['srfm'])
        self.srad = Tkinter.DoubleVar()
        self.srad.set(apbs_defaults['srad'])
        self.system_temp = Tkinter.DoubleVar()
        self.system_temp.set(apbs_defaults['temp'])

        # reaction criteria
        self.contacts_f = Tkinter.StringVar()
        self.contacts_f.set('protein-protein-contacts.xml')
        self.default_contacts_f = Tkinter.BooleanVar()
        self.default_contacts_f.set(False)
        self.rxn_distance = Tkinter.DoubleVar()
        self.rxn_distance.set(bd_defaults['rxn_distance'])
        self.npairs = Tkinter.IntVar()
        self.npairs.set(bd_defaults['npairs'])
        
        # BD parameters and defaults
        self.solvent_eps = Tkinter.DoubleVar()
        self.solvent_eps.set(self.solvent_dielectric.get())
        self.mol_eps = [Tkinter.DoubleVar() for _ in range(2)]
        [self.mol_eps[x].set(self.interior_dielectric.get())  for x in range(2)]
        self.debyel = [Tkinter.DoubleVar() for _ in range(2)]
        [self.debyel[x].set(bd_defaults['debyel']) for x in range(2)]
        self.ntraj = Tkinter.IntVar()
        self.ntraj.set(bd_defaults['ntraj'])
        self.nthreads = Tkinter.IntVar()
        self.nthreads.set(bd_defaults['nthreads'])
        self.mindx = Tkinter.DoubleVar()
        self.mindx.set(bd_defaults['mindx'])
        self.ntrajo = Tkinter.IntVar()
        self.ntrajo.set(bd_defaults['ntrajo'])
        self.ncopies = Tkinter.IntVar()
        self.ncopies.set(bd_defaults['ncopies'])
        self.nbincopies = Tkinter.IntVar()
        self.nbincopies.set(bd_defaults['nbincopies'])
        self.nsteps = Tkinter.IntVar()
        self.nsteps.set(bd_defaults['nsteps'])
        self.westeps = Tkinter.IntVar()
        self.westeps.set(bd_defaults['westeps'])
        self.maxnsteps = Tkinter.IntVar()
        self.maxnsteps.set(bd_defaults['maxnsteps'])

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
                                'Version %s, NBCR 2016\n\n'
                                'Plugin for setting up and running BrownDye '
                                'Brownian dynamics simulations.' % __version__),
                          background='black', foreground='white')
        w.pack(expand=1, fill='both', **pref)

        # create a notebook
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, **pref)

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

        # arrange widgets using grid
        project_path_ent.grid(sticky='we', row=1, column=0, **pref)
        project_path_but.grid(sticky='we', row=1, column=1, **pref)
        label.grid(           sticky='we', row=1, column=2, **pref)
        project_path_b_but.grid(sticky='we', row=1, column=3, **pref)
        pdb2pqr_path_ent.grid(sticky='we', row=2, column=0, **pref)
        pdb2pqr_path_but.grid(sticky='we', row=2, column=1, **pref)
        apbs_path_ent.grid(   sticky='we', row=3, column=0, **pref)
        apbs_path_but.grid(   sticky='we', row=3, column=1, **pref)
        bd_path_ent.grid(     sticky='we', row=4, column=0, **pref)
        bd_path_but.grid(     sticky='we', row=4, column=1, **pref)
        config.columnconfigure(0, weight=8)
        config.columnconfigure(1, weight=2)

        #############################
        # Tab: PQR files preparation
        #############################
        page = self.notebook.add('PQR files')
        grp_pqr = Tkinter.LabelFrame(page, text='PQR files')
        grp_pqr.grid(sticky='eswn', row=0, column=0, columnspan=3, **pref)

        pdb_0_ent = Pmw.EntryField(grp_pqr,
                                   label_text='Molecule 0 PDB file:', labelpos='wn',
                                   entry_textvariable=self.mol0)
        pdb_0_but = Tkinter.Button(grp_pqr, text='Browse...',
                                   command=lambda: self.getPDBMol(0))
        label0 = Tkinter.Label(grp_pqr, text='or')
        pymol_obj0_opt = Pmw.OptionMenu(grp_pqr, labelpos='w',
                                        label_text='Select molecule 0: ',
                                        menubutton_textvariable=self.mol0_object,
                                        menubutton_width=7,
                                        items=(['None'] + pymol.cmd.get_names("all")))
        pdb_1_ent = Pmw.EntryField(grp_pqr,
                                   label_text='Molecule 1 PDB file:', labelpos='wn',
                                   entry_textvariable=self.mol1)
        pdb_1_but = Tkinter.Button(grp_pqr, text='Browse...',
                                   command=lambda: self.getPDBMol(1))
        label1 = Tkinter.Label(grp_pqr, text='or')
        pymol_obj1_opt = Pmw.OptionMenu(grp_pqr, labelpos='w',
                                        label_text='Select molecule 1: ',
                                        menubutton_textvariable=self.mol1_object,
                                        menubutton_width=7,
                                        items=(['None'] + pymol.cmd.get_names("all")))
        pqr_ff_opt = Pmw.OptionMenu(grp_pqr, labelpos='w',
                                    label_text='Force field: ',
                                    menubutton_textvariable=self.pqr_ff,
                                    menubutton_width=7,
                                    items=['parse', 'charmm', 'amber'])
        pqr_an_but = Tkinter.Checkbutton(grp_pqr,
                                         text=('Assign charge and radius only '
                                               '(no structure optimization)'),
                                         variable=self.pqr_assign_only,
                                         onvalue=True, offvalue=False)
        pqr_opt_but = Tkinter.Button(page, text='Create PQR files',
                                     command=self.pdb2pqr)
        label2 = Tkinter.Label(page, text='or load your PQR files:')
        pqr_0_ent = Pmw.EntryField(page,
                                   label_text='Molecule 0 PQR file:', labelpos='wn',
                                   entry_textvariable=self.pqr0)
        pqr_0_but = Tkinter.Button(page, text='Browse...',
                                   command=self.getPQRMol0)
        pqr_1_ent = Pmw.EntryField(page,
                                   label_text='Molecule 1 PQR file:', labelpos='wn',
                                   entry_textvariable=self.pqr1)
        pqr_1_but = Tkinter.Button(page, text='Browse...',
                                   command=self.getPQRMol1)

        pdb_0_ent.grid(     sticky='we', row=0, column=0, **pref)
        pdb_0_but.grid(     sticky='we', row=0, column=1, **pref)
        label0.grid(        sticky='we', row=0, column=2, **pref)
        pymol_obj0_opt.grid(sticky='we', row=0, column=3, **pref)
        pdb_1_ent.grid(     sticky='we', row=1, column=0, **pref)
        pdb_1_but.grid(     sticky='we', row=1, column=1, **pref)
        label1.grid(        sticky='we', row=1, column=2, **pref)
        pymol_obj1_opt.grid(sticky='we', row=1, column=3, **pref)
        pqr_ff_opt.grid(    sticky='we', row=2, column=0, **pref)
        pqr_an_but.grid(    sticky='we', row=3, column=0, **pref)
        pqr_opt_but.grid(   sticky='we', row=4, column=0, **pref)
        label2.grid(        sticky='we', row=5, column=0, **pref)
        pqr_0_ent.grid(     sticky='we', row=6, column=0, **pref)
        pqr_0_but.grid(     sticky='we', row=6, column=1, **pref)
        pqr_1_ent.grid(     sticky='we', row=7, column=0, **pref)
        pqr_1_but.grid(     sticky='we', row=7, column=1, **pref)
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
                                     entry_textvariable=self.dime0[0],
                                     entry_width=5)
        dime0_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime0[1],
                                     entry_width=5)
        dime0_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime0[2],
                                     entry_width=5)
        cglen0_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='cglen: ',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen0[0],
                                      entry_width=8)
        cglen0_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen0[1],
                                      entry_width=8)
        cglen0_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.cglen0[2],
                                      entry_width=8)
        fglen0_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='fglen: ',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen0[0],
                                      entry_width=8)
        fglen0_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen0[1],
                                      entry_width=8)
        fglen0_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.00},
                                      entry_textvariable=self.fglen0[2],
                                      entry_width=8)
        get_size0_but = Tkinter.Button(grp_grids,
                                       text="Calculate grid size",
                                       command=self.getSizemol0)

        label1 = Tkinter.Label(grp_grids, text='Molecule 1')
        dime1_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='dime: ',
                                     # value=self.dime1[0].get(),
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime1[0],
                                     entry_width=5)
        dime1_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime1[1],
                                     entry_width=5)
        dime1_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                     label_text='',
                                     validate={'validator': 'integer', 'min': 0},
                                     entry_textvariable=self.dime1[2],
                                     entry_width=5)

        cglen1_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='cglen: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen1[0],
                                      entry_width=8)
        cglen1_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen1[1],
                                      entry_width=8)
        cglen1_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.cglen1[2],
                                      entry_width=8)
        fglen1_0_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='fglen: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen1[0],
                                      entry_width=8)
        fglen1_1_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen1[1],
                                      entry_width=8)
        fglen1_2_ent = Pmw.EntryField(grp_grids, labelpos='w',
                                      label_text='',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.fglen1[2],
                                      entry_width=8)

        get_size1_but = Tkinter.Button(grp_grids,
                                       text="Calculate grid size",
                                       command=self.getSizemol1)

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
                                         entry_textvariable=self.solvent_dielectric,
                                         entry_width=10)
        solute_die_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                        label_text='Molecule 0/1 eps: ',
                                        validate={'validator': 'real', 'min': 0.0},
                                        entry_textvariable=self.interior_dielectric,
                                        entry_width=10)
        ion1_charge_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                         label_text='Ion(1) charge: ',
                                         validate={'validator': 'integer', 'min': -2},
                                         entry_textvariable=self.ion_charge[0],
                                         entry_width=5)
        ion1_conc_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                       label_text='conc.: ',
                                       validate={'validator': 'real', 'min': 0.0},
                                       entry_textvariable=self.ion_conc[0],
                                       entry_width=5)
        ion1_rad_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                      label_text='radius: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.ion_rad[0],
                                      entry_width=5)
        ion2_charge_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                         label_text='Ion(2) charge: ',
                                         validate={'validator': 'integer', 'min': -2},
                                         entry_textvariable=self.ion_charge[1],
                                         entry_width=5)
        ion2_conc_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                       label_text='conc.: ',
                                       validate={'validator': 'real', 'min': 0.0},
                                       entry_textvariable=self.ion_conc[1],
                                       entry_width=5)
        ion2_rad_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                      label_text='radius: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.ion_rad[1],
                                      entry_width=5)

        sdens_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                   label_text='Surf. sphere density: ',
                                   validate={'validator': 'real', 'min': 0.0},
                                   entry_textvariable=self.sdens, entry_width=5)
        srad_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                  label_text='Solvent radius: ',
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.srad, entry_width=5)
        swin_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                  label_text='Spline window: ',
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.swin, entry_width=5)
        temp_ent = Pmw.EntryField(grp_apbs, labelpos='w',
                                  label_text='Temperature: ',
                                  validate={'validator': 'real', 'min': 0.0},
                                  entry_textvariable=self.system_temp, entry_width=5)
        bcfl_ent = Pmw.OptionMenu(grp_apbs, labelpos='w',
                                  label_text='Boundary condition: ',
                                  menubutton_textvariable=self.bcfl,
                                  menubutton_width=5,
                                  items=['sdh', 'zero', 'mdh', 'focus'])
        chgm_ent = Pmw.OptionMenu(grp_apbs, labelpos='w',
                                  label_text='Charge mapping: ',
                                  menubutton_textvariable=self.chgm,
                                  menubutton_width=5,
                                  items=['spl2', 'spl0', 'spl4'])
        srfm_ent = Pmw.OptionMenu(grp_apbs, labelpos='w',
                                  label_text='Diel. surf. calc. method: ',
                                  menubutton_textvariable=self.srfm,
                                  menubutton_width=5,
                                  items=['smol', 'mol', 'spl2', 'spl4'])

        run_apbs_but = Tkinter.Button(page,
                                      text="Run APBS to generate grids",
                                      command=self.runAPBS)

        solvent_die_ent.grid(sticky='we', row=0, column=0, **pref)
        solute_die_ent.grid( sticky='we', row=1, column=0, **pref)
        sdens_ent.grid(      sticky='we', row=2, column=0, **pref)
        srad_ent.grid(       sticky='we', row=3, column=0, **pref)
        swin_ent.grid(       sticky='we', row=4, column=0, **pref)
        temp_ent.grid(       sticky='we', row=5, column=0, **pref)

        ion1_charge_ent.grid(sticky='we', row=6, column=0, **pref)
        ion1_conc_ent.grid(  sticky='we', row=6, column=1, **pref)
        ion1_rad_ent.grid(   sticky='we', row=6, column=2, **pref)
        ion2_charge_ent.grid(sticky='we', row=7, column=0, **pref)
        ion2_conc_ent.grid(  sticky='we', row=7, column=1, **pref)
        ion2_rad_ent.grid(   sticky='we', row=7, column=2, **pref)

        apbs_mode_ent.grid(sticky='we', row=0, column=1, columnspan=2, **pref)
        bcfl_ent.grid(     sticky='we', row=1, column=1, columnspan=2, **pref)
        chgm_ent.grid(     sticky='we', row=2, column=1, columnspan=2, **pref)
        srfm_ent.grid(     sticky='we', row=3, column=1, columnspan=2, **pref)

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

        contacts_ent.grid(    sticky='we', row=0, column=0, **pref)
        contacts_but.grid(    sticky='e',  row=0, column=1, **pref)
        def_contacts_but.grid(sticky='e',  row=0, column=2, **pref)
        rxn_distance_ent.grid(sticky='we', row=1, column=0, **pref)
        npairs_ent.grid(      sticky='we', row=2, column=0, **pref)
        run_rxn_crit.grid(    sticky='e', row=3, column=0, **pref)
        grp_rxn.columnconfigure(0, weight=1)
        page.columnconfigure(0, weight=1)
        
        ######################
        # Tab: BD input files
        ######################
        page = self.notebook.add('BD setup')
        grp_bdinput = Tkinter.LabelFrame(page, text='BrownDye input file')
        grp_bdinput.grid(sticky='eswn', row=0, column=0, columnspan=2, **pref)

        solvent_eps_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                         label_text='Solvent eps: ',
                                         validate={'validator': 'real', 'min': 0.0},
                                         entry_textvariable=self.solvent_eps, entry_width=5)
        debyel_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                    label_text='Solvent Debye length: ',
                                    validate={'validator': 'real', 'min': 0.0},
                                    entry_textvariable=self.debyel[0], entry_width=5)
        get_debyel_but = Tkinter.Button(grp_bdinput, text='Get Debye length',
                                        command=self.getDebyeLength)
        mol0_eps_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                      label_text='Molecule 0 eps: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.mol_eps[0], entry_width=5)
        mol1_eps_ent = Pmw.EntryField(grp_bdinput, labelpos='wn',
                                      label_text='Molecule 1 eps: ',
                                      validate={'validator': 'real', 'min': 0.0},
                                      entry_textvariable=self.mol_eps[1], entry_width=5)
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
                                        entry_textvariable=self.nbincopies, entry_width=5)
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

        prep_bd_but = Tkinter.Button(page, text="Generate BrownDye input files",
                                     command=self.prepBD)

        self.status_bar = Pmw.MessageBar(page, entry_width=20, entry_relief='sunken',
                                         labelpos='w', label_text='Status:')
        self.status_bar.message('state', '')

        solvent_eps_ent.grid(sticky='we', row=1, column=0, **pref)
        debyel_ent.grid(     sticky='we', row=2, column=0, **pref)
        get_debyel_but.grid( sticky='w',  row=2, column=1, **pref)
        mol0_eps_ent.grid(   sticky='we', row=3, column=0, **pref)
        mol1_eps_ent.grid(   sticky='we', row=3, column=1, **pref)

        ntraj_ent.grid(      sticky='we', row=4, column=0, **pref)
        nthreads_ent.grid(   sticky='we', row=4, column=1, **pref)
        mindx_ent.grid(      sticky='we', row=5, column=0, **pref)
        ntrajo_ent.grid(     sticky='we', row=5, column=1, **pref)
        ncopies_ent.grid(    sticky='we', row=6, column=0, **pref)
        nbincopies_ent.grid( sticky='we', row=6, column=1, **pref)
        nsteps_ent.grid(     sticky='we', row=7, column=0, **pref)
        westeps_ent.grid(    sticky='we', row=7, column=1, **pref)
        maxnsteps_ent.grid(  sticky='we', row=8, column=0, **pref)

        prep_bd_but.grid(    sticky='e', row=9, column=0, **pref)
        self.status_bar.grid(sticky='e', row=10, column=0, **pref)
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
        grp_about = Tkinter.LabelFrame(page, text='About BrownDye Plugin for PyMOL')
        # grp_about.grid(sticky='n', row=0, column=0, columnspan=2, **pref)
        grp_about.pack(fill='both', expand=True, **pref)
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

        label_about = Tkinter.Label(grp_about, text=about_plugin)
        # label_about.grid(sticky='we', row=0, column=2, **pref)
        label_about.pack(fill='both', expand=True, **pref)
        self.notebook.setnaturalsize()
        return

    def getProjectDir(self):
        """Generate a random project directory name."""
        rndID = random.randint(10000, 99999)
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
        print("::: Created project directory: %s" % self.projectDir.get())
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
        return(p.returncode)

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

    def getPDBMol(self, n):
        """Get molecule 0/1 PDB filename."""
        fname = tkFileDialog.askopenfilename(title='PDB File',
                                             initialdir='',
                                             filetypes=[
                                                 ('pdb files', '*.pdb *.ent'),
                                                 ('all files', '*')],
                                             parent=self.parent)
        if n == 0:
            self.mol0.set(fname)
        else:
            self.mol1.set(fname)
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

    def getPQRMol0(self):
        """Get molecule 0 PQR filename."""
        fname = tkFileDialog.askopenfilename(title='PQR File',
                                             initialdir='',
                                             filetypes=[
                                                 ('pqr files', '*.pqr'),
                                                 ('all files', '*')],
                                             parent=self.parent)
        if len(fname) > 0:
            self.pqr0.set(fname)
            target_f = '%s/%s.pqr' % (self.projectDir.get(), MOL0)
            if fname != target_f:
                if os.path.isfile(target_f): os.remove(target_f)
                shutil.copyfile(self.pqr0.get(), target_f)
        return

    def getPQRMol1(self):
        """Get molecule 1 PQR filename."""
        fname = tkFileDialog.askopenfilename(title='PQR File',
                                             initialdir='',
                                             filetypes=[
                                                 ('pqr files', '*.pqr'),
                                                 ('all files', '*')],
                                             parent=self.parent)
        if len(fname) > 0:
            self.pqr1.set(fname)
            target_f = '%s/%s.pqr' % (self.projectDir.get(), MOL1)
            if fname != target_f:
                if os.path.isfile(target_f): os.remove(target_f)
                shutil.copyfile(self.pqr1.get(), target_f)
        return

    def getSizemol0(self):
        """Calculate APBS grid dimensions for molecule 0."""
        pqr_fname = '%s.pqr' % MOL0
        if not os.path.isfile(pqr_fname):
            print("::: %s does not exist!" % pqr_fname)
            return
        psize = Psize(self)
        psize.runPsize(pqr_fname)
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
        pqr_fname = '%s.pqr' % MOL1
        if not os.path.isfile(pqr_fname):
            print("::: %s does not exist!" % pqr_fname)
            return
        psize = Psize(self)
        psize.runPsize(pqr_fname)
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
        fname = tkFileDialog.askopenfilename(
            title='Contacts File', initialdir='',
            filetypes=[('xml files', '*.xml'), ('all files', '*')],
            parent=self.parent)
        self.contacts_f.set(fname)
        return

    def pdb2pqr(self):
        """Convert PDB to PQR."""
        target_f = '%s.pdb' % MOL0
        if self.mol0_object.get() == 'None':
            # if not filecmp.cmp(self.mol0.get(), target_f):
            try:
                shutil.copyfile(self.mol0.get(), target_f)
            except:
                e = sys.exc_info()[0]
                print(e)
                print("::: Creating of %s failed!" % target_f)
                return
        else:
            pymol.cmd.save(filename=target_f,
                           selection=self.mol0_object.get())
        target_f = '%s.pdb' % MOL1
        if self.mol1_object.get() == 'None':
            # if not filecmp.cmp(self.mol1.get(), target_f):
            try:
                shutil.copyfile(self.mol1.get(), target_f)
            except:
                e = sys.exc_info()[0]
                print(e)
                print("::: Creating of %s failed!" % target_f)
                return
        else:
            pymol.cmd.save(filename=target_f,
                           selection=self.mol1_object.get())

        assign_only = ''
        if self.pqr_assign_only.get(): assign_only = '--assign-only'
        pqr_options = ('%s %s --ff=%s' %
                       (assign_only, self.pdb2pqr_opt.get(), self.pqr_ff.get()))
        pdb2pqr_exe = ('%s/%s' % (self.pdb2pqr_path.get(), PDB2PQR_EXE))
        if not self.check_exe(pdb2pqr_exe): return
        for i in [MOL0, MOL1]:
            command = ('%s %s %s.pdb %s.pqr' %
                       (pdb2pqr_exe, pqr_options, i, i))
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
        apbs_template = """# APBS template for BrownDye grids
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
        if not self.check_exe(apbs_exe): return
        for i in [MOL0, MOL1]:
            pqr_filename = '%s.pqr' % i
            if i == MOL0:
                grid_points = [self.dime0[x].get() for x in range(3)]
                cglen = [self.cglen0[x].get() for x in range(3)]
                fglen = [self.fglen0[x].get() for x in range(3)]
            if i == MOL1:
                grid_points = [self.dime1[x].get() for x in range(3)]
                cglen = [self.cglen1[x].get() for x in range(3)]
                fglen = [self.fglen1[x].get() for x in range(3)]
            if grid_points[0] == 0 or grid_points[1] == 0 or grid_points[2] == 0:
                print("::: %s - no grid points defined!" % i)
                return
            dx_filename = i
            fout = '%s.in' % i
            with open(fout, "w") as f:
                f.write((apbs_template %
                         (pqr_filename,
                          grid_points[0], grid_points[1], grid_points[2],
                          cglen[0], cglen[1], cglen[2],
                          fglen[0], fglen[1], fglen[2],
                          self.apbs_mode.get(), self.bcfl.get(),
                          self.ion_charge[0].get(),
                          self.ion_conc[0].get(), self.ion_rad[0].get(),
                          self.ion_charge[1].get(),
                          self.ion_conc[1].get(), self.ion_rad[1].get(),
                          self.interior_dielectric.get(),
                          self.solvent_dielectric.get(),
                          self.chgm.get(), self.sdens.get(),
                          self.srfm.get(), self.srad.get(),
                          self.swin.get(), self.system_temp.get(),
                          dx_filename)))

            command = 'MCSH_HOME=. %s %s.in' % (apbs_exe, i)
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
                iomc = '%s-io.mc' % i
                if os.path.isfile(iomc): os.remove(iomc)
                shutil.copyfile('io.mc', iomc)
                print("::: Done.")
            else:
                print("::: Failed: %s" % command)
        return

    def getDebyeLength(self):
        dl = [[], []]
        for i in range(2):
            if i == 0: mol = MOL0
            if i == 1: mol = MOL1
            fname = '%s-io.mc' % mol
            if DEBUG > 2: print("Parsing %s for Debye length ..." % fname)
            if not os.path.isfile(fname):
                print("::: File %s does not exist!" % fname )
                print("::: Run APBS first.")
                return
            with open(fname) as f:
                for line in f:
                    l = re.findall(r'Debye length = .*', line)
                    if l:
                        dl[i].append(l[0].split(' ')[3])
            if len(dl[i]) == 0:
                print("::: No Debye length found in %s for %s!" %
                      (fname, mol))
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
        if not self.check_exe(pqr2xml_exe): return
        for i in [MOL0, MOL1]:
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
                   '-nonred -mol0 %s-atoms.pqrxml -mol1 %s-atoms.pqrxml '
                   '-ctypes %s  -dist %f > %s-%s-rxn-pairs.xml'
                   % (self.bd_path.get(), MOL0, MOL1, self.contacts_f.get(),
                      self.rxn_distance.get(), MOL0, MOL1))
        if DEBUG > 2: print(command)
        print("::: Running make_rxn_pairs ...")
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
        command =('%s/make_rxn_file '
                  '-pairs %s-%s-rxn-pairs.xml -distance %f '
                  ' -nneeded %d > %s-%s-rxns.xml'
                  % (self.bd_path.get(), MOL0, MOL1, self.rxn_distance.get(),
                     self.npairs.get(), MOL0, MOL1))
        if DEBUG > 2: print(command)
        print("::: Running make_rxn_file ...")
        rc = self.runCmd(command)
        if rc == 0:
            print("::: Done.")
        else:
            print("::: Failed: %s" % command)
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
  <n-steps-per-output> 10000 </n-steps-per-output>
</root>
"""
        bdtop_input = 'input.xml'
        BDTOP_EXE = 'bd_top'
        with open(bdtop_input, "w") as f:
            f.write((nam_simulation_template %
                     (self.solvent_eps.get(), self.debyel[0].get(),
                      self.ntraj.get(), self.nthreads.get(),
                      MOL0, MOL0, MOL0, self.mol_eps[0].get(),
                      MOL1, MOL1, MOL1, self.mol_eps[1].get(),
                      self.mindx.get(), MOL0, MOL1,
                      self.ntrajo.get(),
                      self.ncopies.get(), self.nbincopies.get(),
                      self.nsteps.get(), self.westeps.get(),
                      self.maxnsteps.get())))
        for i in [MOL0, MOL1]:
            dxfile = '%s.dx' % i
            if not os.path.isfile(dxfile):
                print("::: %s DX file does not exist!" % dxfile)
                print("::: Run APBS first to create it.")
                return
        command = ('PATH=%s:${PATH} %s %s'
                   % (self.bd_path.get(), BDTOP_EXE, bdtop_input))
        if DEBUG > 2: print(command)
        print("::: Running bd_top (this will take a couple of minutes) ...")
        self.status_bar.message('state', ' Starting bd_top ...')
        thread = RunBDTopThread(self, command)
        thread.start()
        return

    def runRxnCrit(self):
        self.runPqr2xml()
        self.makeRxnCriteria()
        return
    
    def runBD(self):
        """Start BrownDye simulation either in foreground or background."""
        print("::: Starting BrownDye simulation ...")
        nam_simulation_exe = '%s/nam_simulation' % self.bd_path.get()
        if not self.check_exe(nam_simulation_exe): return
        command = ('%s %s-%s-simulation.xml'
                   % (nam_simulation_exe, MOL0, MOL1))
        if self.run_in_background.get():
            p = subprocess.Popen(" nohup %s" % command, shell=True)
            global JOBPID
            JOBPID = p.pid
            self.notebook.selectpage('BD simulation')
            self.logtxt_ent.insert('end', "::: Starting BrownDye simulation in background.\n")
            self.logtxt_ent.insert('end', "::: Job PID: %d \n" % JOBPID)
        else:
            thread = RunThread(self, command)
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
                                 stdout=subprocess.PIPE).stdout.read()
            if DEBUG > 1: print(command, s)
            command = 'kill -9 %d' % JOBPID
            s = subprocess.Popen(command, shell=True,
                                 stdout=subprocess.PIPE).stdout.read()
            if DEBUG > 1: print(command, s)
            print("::: Background job (PID %d) terminated." % JOPBPID)
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
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=True)
        # p.wait()
        stdout, stderr = p.communicate()
        if p.returncode > 0:
            print("::: Error processing %s trajectory file." % self.traj_f.get())
            return
        t = etree.fromstring(stdout)
        tr = t.xpath('//trajectory/number')
        traj_index = [int(tr[x].text.strip()) for x in range(len(tr))]
        logline = ('%d association event trajectories found (index numbers: %s)\n'
                   % (len(traj_index), str(traj_index)))
        self.msg_ent.insert('end', "%s" % logline)
        # number of frames
        self.dialog_idx.clear()
        for i in traj_index:
            with open(self.traj_f.get(), 'r') as f:
                d = etree.parse(f)
                xpath_str = 'trajectory[n-traj=" %d "]//s/n' % i
                j = d.xpath(xpath_str)
                logline = ('trajectory %d: %s frames\n' 
                           % (i, j[0].text.strip()))
                self.msg_ent.insert('end', "%s" % logline)
                self.dialog_idx.insert('end', i)
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
                   % (self.bd_path.get(), MOL0, MOL1, ofile, self.xyz_ofile))
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

    def check_exe(self, command):
        """Check if command exists and is executable."""
        if os.path.isfile(command) and os.access(command, os.X_OK):
            is_exe = True
        else:
            is_exe = False
            print("::: Can't find %s!" % command)
        return is_exe
        
    def exitBDPlugin(self, result):
        """Quit BrownDye plugin."""
        if tkMessageBox.askokcancel("Really quit?",
                                    "Exit BrownDye Plugin now?") == 0:
            return 0
        print("Exiting BrownDye Plugin ...")
        if __name__ == '__main__':
            self.parent.destroy()
        else:
            self.dialog.withdraw()
        print("Done.")
        return

class RunBDTopThread(Thread):
    """bd_top thread management class."""
    def __init__(self, my_inst, command):
        Thread.__init__(self)
        self.status_bar = my_inst.status_bar
        self.command = command
        if DEBUG > 2: print("%s" % (command))
        return

    def run(self):
        self.status_bar.message('state', ' Running bd_top ...')
        p = subprocess.Popen(self.command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        self.outlog = stdout
        self.status = p.returncode
	try:
            self.status_bar.message('state', ' Finished bd_top ...')
            time.sleep(3)
            self.status_bar.message('state', ' Done.')
            print(stdout, stderr, p.returncode)
        except:
            print("::: Thread error!")
        return


class RunThread(Thread):
    """Thread management class."""
    def __init__(self, my_inst, command):
        Thread.__init__(self)
        self.page = my_inst.logtxt_ent
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
        p = subprocess.Popen(self.command, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, shell=True, bufsize=1)
        global JOBPID
        JOBPID = p.pid
        if DEBUG > 1: print('JOBPID: %d' % JOBPID)
        for line in iter(p.stdout.readline, ''):
            print(line,)
            self.page.insert('end',"%s" % line)
        stdout, stderr = p.communicate()
        self.outlog = stdout
        self.status = p.returncode
	try:
            self.page.insert('end',"%s %s" % (stdout, stderr))
            # self.page.yview('moveto', 1.0)
        except:
            print("::: Thread error!")
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
        readcount = 0 
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
                    ntraj = results.xpath('//reactions/n-trajectories')[0].text.strip()
                    stuck = results.xpath('//reactions/stuck')[0].text.strip()
                    escaped = results.xpath('//reactions/escaped')[0].text.strip()
                    ncompleted = results.xpath('//reactions/completed/name')[0].text.strip()
                    completed = results.xpath('//reactions/completed/n')[0].text.strip()
                    mymsg = '%s out of %d' % (ntraj, self.totntraj)
                    self.msgbar1.message('state', mymsg)
                    mymsg= '%s / %s / %s' % (completed, escaped, stuck) 
                    self.msgbar2.message('state', mymsg)
                command = ('cat %s | %s/compute_rate_constant'
                           % (results_file, self.bd_path))
                p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True)
                # p.wait()
                stdout, stderr = p.communicate()
                rates = etree.fromstring(stdout)
                rate_constant = rates.xpath('//rate-constant/mean')[0].text.strip()
                rxn_probability = rates.xpath('//reaction-probability/mean')[0].text.strip()
                mymsg = ('%s / %s' % (rate_constant, rxn_probability))
                self.msgbar3.message('state', mymsg)

        time.sleep(2)
        if DEBUG > 2: print(self.mythread.is_alive())
        print("::: BrownDye simulation finished.")
        self.page.insert('end', "::: BrownDye simulation finished\n")
        return

class StopThread(Thread):
    """Kill running thread."""
    def __init__(self, mythread):
        self.pid = mythread.pid
        os.system('kill -9 %d' % self.pid)
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

        self.constants['CFAC'] = psize_defaults['cfac']
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
        with open(filename, "r") as file:
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
