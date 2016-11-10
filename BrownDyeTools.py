#!/usr/bin/env python2
#-*- coding: utf-8 -*-
#
# Last modified: 2016-11-10 12:36:22

import os, subprocess
import sys, string
#from sys import stdout
#from math import log
import time
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
import tkColorChooser
import Tkinter
import Pmw

DEBUG = 5

PDB2PQR_PATH = ''
APBS_PATH = ''
BD_PATH = ''

if DEBUG > 0:
    PDB2PQR_PATH = '/opt/pdb2pqr-linux-bin64-2.1.0/'
    APBS_PATH = '/export1/Builds/SandBox-2016.10.11-14.05/bin/'
    BD_PATH = '/home/rok/BrownianDynamics/browndye-2016.4.14/bin/'

if "PDB2PQR_PATH" in os.environ: PDB2PQR_PATH = os.environ["PDB2PQR_PATH"]
if "APBS_PATH" in os.environ: APBS_PATH = os.environ["APBS_PATH"]
if "BD_PATH" in os.environ: BD_PATH = os.environ["BD_PATH"]


def __init__(self):
    """ BD plugin for PyMol
    """
    self.menuBar.addmenuitem('Plugin', 'command',
                             'BrownDye Plugin', label='BrownDye Plugin',
                             command=lambda s=self: BDPlugin(s))
    
class BDPlugin:
    def __init__(self, app):
        self.parent = app.root
        self.createGUI()
        
    def createGUI(self):
        self.dialog = Pmw.Dialog(self.parent, buttons=('Exit',),
                                 title='BrownDye Plugin for PyMOL',
                                command=self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))


        self.pdb2pqr_path = Tkinter.StringVar()
        self.pdb2pqr_path.set(PDB2PQR_PATH)
        self.apbs_path = Tkinter.StringVar()
        self.apbs_path.set(APBS_PATH)
        self.bd_path = Tkinter.StringVar()
        self.bd_path.set(BD_PATH)
                
        # parameters used by BD
        self.mol0 = Tkinter.StringVar()
        self.mol1 = Tkinter.StringVar()
        self.mol0.set('mol0.pdb') # for testing only
        self.mol1.set('mol1.pdb')
        self.pdb2pqr_opt = Tkinter.StringVar()
        self.pdb2pqr_opt.set('--apbs-input')
        self.pqr_assign_only = Tkinter.BooleanVar()
        self.pqr_assign_only.set(True)
        self.pqr_ff = Tkinter.StringVar()
        self.pqr_ff.set('parse')

        # APBS parameters and defaults
        self.dime0 = [Tkinter.IntVar(), Tkinter.IntVar(), Tkinter.IntVar()]
        self.dime0[0].set(0)
        self.dime0[1].set(0)
        self.dime0[2].set(0)
        self.cglen0 = [Tkinter.DoubleVar(), Tkinter.DoubleVar(), Tkinter.DoubleVar()]
        self.cglen0[0].set(0.0)
        self.cglen0[1].set(0.0)
        self.cglen0[2].set(0.0)
        self.fglen0 = [Tkinter.DoubleVar(), Tkinter.DoubleVar(), Tkinter.DoubleVar()]
        self.fglen0[0].set(0.0)
        self.fglen0[1].set(0.0)
        self.fglen0[2].set(0.0)
        # APBS parameters and defaults
        self.dime1 = [Tkinter.IntVar(), Tkinter.IntVar(), Tkinter.IntVar()]
        self.dime1[0].set(0)
        self.dime1[1].set(0)
        self.dime1[2].set(0)
        self.cglen1 = [Tkinter.DoubleVar(), Tkinter.DoubleVar(), Tkinter.DoubleVar()]
        self.cglen1[0].set(0.0)
        self.cglen1[1].set(0.0)
        self.cglen1[2].set(0.0)
        self.fglen1 = [Tkinter.DoubleVar(), Tkinter.DoubleVar(), Tkinter.DoubleVar()]
        self.fglen1[0].set(0.0)
        self.fglen1[1].set(0.0)
        self.fglen1[2].set(0.0)

        self.apbs_mode = Tkinter.StringVar()
        self.apbs_mode.set('lpbe')
        self.bcfl = Tkinter.StringVar()
        self.bcfl.set('sdh')
        self.ion_plus_one_conc = Tkinter.DoubleVar()
        self.ion_plus_one_conc.set(0.15)
        self.ion_plus_one_rad = Tkinter.DoubleVar()
        self.ion_plus_one_rad.set(0.2)
        self.ion_minus_one_conc = Tkinter.DoubleVar()
        self.ion_minus_one_conc.set(0.15)
        self.ion_minus_one_rad = Tkinter.DoubleVar()
        self.ion_minus_one_rad.set(0.2)
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
        self.contacts = Tkinter.StringVar()
        self.contacts.set('protein-protein-contacts.xml')
        self.rxn_distance = Tkinter.DoubleVar()
        self.rxn_distance.set(5.0)
        self.npairs= Tkinter.IntVar()
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
        self.ntraj.set(1000)
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
        
        self.pdb_fn        = Tkinter.StringVar()
        self.pymol_sel     = Tkinter.StringVar()
        self.msms_bin      = Tkinter.StringVar()
##         self.pdb2xyzr_bin  = Tkinter.StringVar()
        self.pdb2xyzrn_bin = Tkinter.StringVar()
        self.tmp_dir       = Tkinter.StringVar()

        self.cleanup_saved_pymol_sel = Tkinter.BooleanVar()
        self.cleanup_saved_pymol_sel.set(True) # by default, clean up

        w = Tkinter.Label(self.dialog.interior(),
                          text = '\nBrownDye Plugin for PyMOL\nrok@NBCR, 2016\n\n\
   Plugin for setting up and running BrownDye Browndian dynamics simulations.',
                          background = 'black', foreground = 'green'
                          )
        w.pack(expand=1, fill='both', padx=10, pady=5)

        # make a few tabs within the dialog
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=10, pady=10)

        #####################
        # Tab: Configuration
        #####################
        page = self.notebook.add('Configuration')       
        self.notebook.tab('Configuration').focus_set()
        config = Tkinter.LabelFrame(page, text='Calculation configuration')
        config.pack(fill='both', expand=True, padx=10, pady=10)

        pymol_sel_ent = Pmw.EntryField(config,
                                       label_text='PyMOL selection:',
                                       labelpos='wn',
                                       entry_textvariable=self.pymol_sel
                                       )
        clean_cb = Tkinter.Checkbutton(config,
                                       text='Clean up tmp pdb (saved PyMOL selection) in the temp dir.', 
                                       variable=self.cleanup_saved_pymol_sel,
                                       onvalue=True, offvalue=False)
        label = Tkinter.Label(config, text='or')


        pdb2pqr_path_ent = Pmw.EntryField(config,
                                          label_text='Select PDB2PQR location: ', labelpos='wn',
                                          entry_textvariable=self.pdb2pqr_path)
        pdb2pqr_path_but = Tkinter.Button(config, text = 'Browse...', command=self.getPDB2PQRpath)
        apbs_path_ent = Pmw.EntryField(config,
                                       label_text = 'Select APBS_PATH location:', labelpos='w',
                                       entry_textvariable=self.apbs_path)
        apbs_path_but = Tkinter.Button(config, text = 'Browse...', command = self.getAPBSpath)
        bd_path_ent = Pmw.EntryField(config,
                                    label_text='Select BD_PATH location: ', labelpos='wn',
                                     entry_textvariable=self.bd_path)
        bd_path_but = Tkinter.Button(config, text = 'Browse...', command=self.getBDpath)

        # arrange widgets using grid
        pymol_sel_ent.grid(sticky='we', row=0, column=0, columnspan=2, padx=5, pady=5)
        clean_cb.grid(sticky='w', row=1, column=0, columnspan=2, padx=1, pady=1)
        label.grid(sticky='we', row=2, column=0, columnspan=2, padx=5, pady=10)
        pdb2pqr_path_ent.grid(sticky='we', row=3, column=0, padx=5, pady=5)
        pdb2pqr_path_but.grid(sticky='we', row=3, column=1, padx=5, pady=5)
        apbs_path_ent.grid(   sticky='we', row=4, column=0, padx=5, pady=5)
        apbs_path_but.grid(   sticky='we', row=4, column=1, padx=5, pady=5)
        bd_path_ent.grid(     sticky='we', row=5, column=0, padx=5, pady=5)
        bd_path_but.grid(     sticky='we', row=5, column=1, padx=5, pady=5)
        config.columnconfigure(0, weight=9)
        config.columnconfigure(1, weight=1)

        #############################
        # Tab: PQR files preparation
        #############################
        page = self.notebook.add('PQR files')
        group_pqr = Tkinter.LabelFrame(page, text='PQR files')
        group_pqr.grid(sticky='eswn',row=0, column=0, columnspan=2, padx=10, pady=5)

        pdb_a_ent = Pmw.EntryField(group_pqr,
                                   label_text = 'Molecule 0 PDB file:', labelpos='wn',
                                   entry_textvariable=self.mol0)
        pdb_a_but = Tkinter.Button(group_pqr, text = 'Browse...',
                                   command = self.getPDBmol0)
        pdb_b_ent = Pmw.EntryField(group_pqr,
                                   label_text = 'Molecule 1 PDB file:', labelpos='wn',
                                   entry_textvariable=self.mol1)
        pdb_b_but = Tkinter.Button(group_pqr, text = 'Browse...',
                                   command = self.getPDBmol1)
        #pqr_opt = Pmw.EntryField(group_pqr,
        #                         label_text='pdb2pqr options:', labelpos='wn',
        #                         value=self.pdb2pqr_opt.get(),
        #                         entry_textvariable=self.pdb2pqr_opt,
        #                         entry_width=30)
        pqr_ff_opt = Pmw.OptionMenu(group_pqr, labelpos='w',
                                    label_text='Force field: ',
                                    menubutton_textvariable=self.pqr_ff,
                                    menubutton_width=7,
                                    items=['parse', 'charmm'])


        pqr_an_but = Tkinter.Checkbutton(group_pqr,
                                         text='Assign charge and radius only (no structure optimization)', 
                                         variable=self.pqr_assign_only,
                                         onvalue=True, offvalue=False)


        pqr_opt_but = Tkinter.Button(page, text='Create PQR files',
                                      command=self.pdb2pqr)

        pdb_a_ent.grid(sticky='we', row=0, column=0, padx=5, pady=1)
        pdb_a_but.grid(sticky='we', row=0, column=1, padx=5, pady=1)
        pdb_b_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        pdb_b_but.grid(sticky='we', row=1, column=1, padx=5, pady=1)
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
        dime0_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                     label_text='dime: ',
                                     value=self.dime0[0].get(),
                                     validate = {'validator':'integer', 'min':0},
                                     entry_textvariable=self.dime0[0],
                                     entry_width=5)
        dime1_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                     label_text='',
                                     value=self.dime0[1].get(),
                                     validate = {'validator':'integer', 'min':0},
                                     entry_textvariable=self.dime0[1],
                                     entry_width=5)
        dime2_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                     label_text='',
                                     value=self.dime0[2].get(),
                                     validate = {'validator':'integer', 'min':0},
                                     entry_textvariable=self.dime0[2],
                                     entry_width=5)

        cglen0_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='cglen: ',
                                    value=self.cglen0[0].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                      entry_textvariable=self.cglen0[0],
                                    entry_width=8)
        cglen1_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.cglen0[1].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                      entry_textvariable=self.cglen0[1],
                                    entry_width=8)
        cglen2_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.cglen0[2].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                      entry_textvariable=self.cglen0[2],
                                    entry_width=8)
        fglen0_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='fglen: ',
                                    value=self.fglen0[0].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.fglen0[0],
                                    entry_width=8)
        fglen1_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.fglen0[1].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.fglen0[1],
                                    entry_width=8)
        fglen2_0_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.fglen0[2].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.fglen0[2],
                                    entry_width=8)
        get_size0_but = Tkinter.Button(group_grids,
                                      text="Get grid size for molecule 0",
                                       command=self.getSizemol0)

        label1 = Tkinter.Label(group_grids, text='Molecule 1')
        dime0_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                     label_text='dime: ',
                                     value=self.dime1[0].get(),
                                     validate = {'validator':'integer', 'min':0},
                                     entry_textvariable=self.dime1[0],
                                     entry_width=5)
        dime1_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                     label_text='',
                                     value=self.dime1[1].get(),
                                     validate = {'validator':'integer', 'min':0},
                                     entry_textvariable=self.dime1[1],
                                     entry_width=5)
        dime2_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                     label_text='',
                                     value=self.dime1[2].get(),
                                     validate = {'validator':'integer', 'min':0},
                                     entry_textvariable=self.dime1[2],
                                     entry_width=5)

        cglen0_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='cglen: ',
                                    value=self.cglen1[0].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.cglen1[0],
                                    entry_width=8)
        cglen1_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.cglen1[1].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.cglen1[1],
                                    entry_width=8)
        cglen2_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.cglen1[2].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.cglen1[2],
                                    entry_width=8)
        fglen0_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='fglen: ',
                                    value=self.fglen1[0].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.fglen1[0],
                                    entry_width=8)
        fglen1_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.fglen1[1].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.fglen1[1],
                                    entry_width=8)
        fglen2_1_ent = Pmw.EntryField(group_grids, labelpos = 'w',
                                    label_text='',
                                    value=self.fglen1[2].get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.fglen1[2],
                                    entry_width=8)

        get_size1_but = Tkinter.Button(group_grids,
                                      text="Get grid size for molecule",
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
        solvent_die_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                         label_text='Solvent eps :',
                                         value=self.solvent_dielectric.get(),
                                         validate = {'validator':'real', 'min':0.00},
                                         entry_textvariable=self.solvent_dielectric,
                                         entry_width=10)
        solute_die_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                        label_text='Molecule 0/1 eps: ',
                                        value=self.interior_dielectric.get(),
                                        validate = {'validator':'real', 'min':0.00},
                                        entry_textvariable=self.interior_dielectric,
                                        entry_width=10)

        ion_plus_one_conc_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                               label_text='Ion +1 conc.: ',
                                               value=self.ion_plus_one_conc.get(),
                                               validate = {'validator':'real', 'min':0.00},
                                               entry_textvariable=self.ion_plus_one_conc,
                                               entry_width=5)
        ion_plus_one_rad_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                              label_text='radius: ',
                                              value=self.ion_plus_one_rad.get(),
                                              validate = {'validator':'real', 'min':0.00},
                                              entry_textvariable=self.ion_plus_one_rad,
                                              entry_width=5)
        ion_minus_one_conc_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                               label_text='Ion -1 conc.: ',
                                               value=self.ion_minus_one_conc.get(),
                                               validate = {'validator':'real', 'min':0.00},
                                               entry_textvariable=self.ion_minus_one_conc,
                                                entry_width=5)
        ion_minus_one_rad_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                               label_text='radius: ',
                                               value=self.ion_minus_one_rad.get(),
                                               validate = {'validator':'real', 'min':0.00},
                                               entry_textvariable=self.ion_minus_one_rad,
                                               entry_width=5)

        sdens_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                   label_text='Surf. sphere density: ',
                                   value=self.sdens.get(),
                                   validate = {'validator':'real', 'min':0.00},
                                   entry_textvariable=self.sdens, entry_width=5)
        srad_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                   label_text='Solvent radius: ',
                                   value=self.srad.get(),
                                   validate = {'validator':'real', 'min':0.00},
                                   entry_textvariable=self.srad, entry_width=5)
        swin_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                   label_text='Spline window: ',
                                   value=self.swin.get(),
                                   validate = {'validator':'real', 'min':0.00},
                                   entry_textvariable=self.swin, entry_width=5)
        temp_ent = Pmw.EntryField(group_apbs, labelpos = 'w',
                                   label_text='Temperature: ',
                                   value=self.system_temp.get(),
                                   validate = {'validator':'real', 'min':0.00},
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

        ion_plus_one_conc_ent.grid( sticky='we', row=6, column=0, padx=5, pady=1)
        ion_plus_one_rad_ent.grid(  sticky='we', row=6, column=1, padx=5, pady=1)
        ion_minus_one_conc_ent.grid(sticky='we', row=7, column=0, padx=5, pady=1)
        ion_minus_one_rad_ent.grid( sticky='we', row=7, column=1, padx=5, pady=1)

        apbs_mode_ent.grid(  sticky='we', row=0, column=1, padx=5, pady=1)
        bcfl_ent.grid(              sticky='we', row=1, column=1, padx=5, pady=1)
        chgm_ent.grid(              sticky='we', row=2, column=1, padx=5, pady=1)
        srfm_ent.grid(              sticky='we', row=3, column=1, padx=5, pady=1)

        run_apbs_but.grid(sticky='we', row=8, column=0, columnspan=2, padx=5, pady=1)

        ###############################
        # Tab: Reaction criteria setup
        ###############################
        page = self.notebook.add('Reaction citeria')
        group_rxn = Tkinter.LabelFrame(page, text='Setup reaction criteria')
        group_rxn.grid(sticky='eswn', row=0, column=0, columnspan=2, padx=10, pady=5)

        contacts_ent = Pmw.EntryField(group_rxn,
                                   label_text = 'Contacts file:', labelpos='w',
                                      entry_textvariable=self.contacts, entry_width=50)
        contacts_but = Tkinter.Button(group_rxn, text = 'Browse...',
                                   command = self.getContacts)
        rxn_distance_ent = Pmw.EntryField(group_rxn, labelpos = 'wn',
                                      label_text='Reaction distance: ',
                                      value=self.rxn_distance.get(),
                                      validate = {'validator':'real', 'min':0.01},
                                      entry_textvariable=self.rxn_distance, entry_width=10)
        npairs_ent = Pmw.EntryField(group_rxn, labelpos='wn',
                                    label_text='Number of reaction pairs: ',
                                    value=self.npairs.get(),
                                    validate={'validator': 'integer', 'min': 1},
                                    entry_textvariable=self.npairs, entry_width=10)
        run_rxn_crit = Tkinter.Button(page, text="Creat rxn files", command=self.run_rxn_crit)

        contacts_ent.grid(    sticky='we', row=0, column=0, padx=5, pady=1)
        contacts_but.grid(    sticky='e',  row=0, column=1, padx=5, pady=1)
        rxn_distance_ent.grid(sticky='we', row=1, column=0, padx=5, pady=1)
        npairs_ent.grid(      sticky='we', row=2, column=0, padx=5, pady=1)
        run_rxn_crit.grid(    sticky='we', row=3, column=0, padx=5, pady=1)
        
        ######################
        # Tab: BD input files
        ######################
        page = self.notebook.add('BD setup')
        group_bdinput = Tkinter.LabelFrame(page, text='BD input file')
        group_bdinput.grid(sticky='eswn', row=0, column=0, columnspan=2, padx=10, pady=5)

        solvent_eps_ent = Pmw.EntryField(group_bdinput, labelpos = 'wn',
                                      label_text='Solvent eps: ',
                                      value=self.solvent_eps.get(),
                                      validate = {'validator':'real', 'min':0.00},
                                      entry_textvariable=self.solvent_eps
                                      )
        debyel_ent = Pmw.EntryField(group_bdinput, labelpos = 'wn',
                                    label_text='Debye length: ',
                                    value=self.debyel.get(),
                                    validate = {'validator':'real', 'min':0.00},
                                    entry_textvariable=self.debyel)
        mol0_eps_ent = Pmw.EntryField(group_bdinput, labelpos = 'wn',
                                      label_text='Molecule 0 eps: ',
                                      value=self.mol0_eps.get(),
                                      validate = {'validator':'real', 'min':0.00},
                                      entry_textvariable=self.mol0_eps)
        mol1_eps_ent = Pmw.EntryField(group_bdinput, labelpos = 'wn',
                                      label_text='Molecule 1 eps: ',
                                      value=self.mol1_eps.get(),
                                      validate = {'validator':'real', 'min':0.00},
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
        mindx_ent = Pmw.EntryField(group_bdinput, labelpos = 'wn',
                                   label_text='Time step tolerance: ',
                                   value=self.mindx.get(),
                                   validate = {'validator':'real', 'min':0.01},
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

        prep_bd_but = Tkinter.Button(page, text="Generate BD input file", command=self.prepBD)
        run_bd_but = Tkinter.Button(page, text="Start BD simulation", command=self.runBD)

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
        run_bd_but.grid(     sticky='we', row=8, column=1, padx=5, pady=1)
        
        #############
        # Tab: About
        #############
        page = self.notebook.add('About')
        group_about = Tkinter.LabelFrame(page, text='About BrownDye Plugin for PyMOL')
        group_about.grid(sticky='n', row=0, column=0, columnspan=2, padx=10, pady=5)
        about_plugin = ''' This plugin provides a GUI for setting up and running Brownian \
dynamics simulations with BrownDye.

The plugin requires PDB2PQR, APBS and BrownDye. To download and \
install these applications go to:

http://www.poissonboltzmann.org/
 and 
http://browndye.ucsd.edu/
'''

        label_about = Tkinter.Label(group_about, text=about_plugin)
        label_about.grid(sticky='we', row=0, column=2, padx=5, pady=10)

        self.notebook.setnaturalsize()
        return

    def runcmd(self, command):
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        p.wait()
        (out,err) = p.communicate();

        if (DEBUG > 2): print("returncode = %d" % p.returncode)
        if (DEBUG > 2): print("stdout:\n%s\n" % out)
        if (DEBUG > 2): print("stderr:\n%s\n" % err)
        if p.returncode > 0:
            sys.stderr.write("::: Non-zero return code!\n") 
            sys.stderr.write("::: Failed command: \n\n")
            sys.stderr.write(command + "\n\n")
            sys.stderr.write(err)
            #sys.exit(1)
        return(p.returncode)

    def getPDB2PQRpath(self):
        d = tkFileDialog.askdirectory(
            title='PDB2PQR binary directory', initialdir='',
            parent=self.parent)
        self.pdb2pqr_path.set(d)
        return
        
    def getAPBSpath(self):
        d = tkFileDialog.askdirectory(
            title='APBS binary directory', initialdir='',
            parent=self.parent)
        self.apbs_path.set(d)
        return
        
    def getBDpath(self):
        d = tkFileDialog.askdirectory(
            title='BrownDye binary directory', initialdir='',
            parent=self.parent)
        self.bd_path.set(d)
        return

    def getPDBmol0(self):
        file_name = tkFileDialog.askopenfilename(
            title='PDB File', initialdir='',
            filetypes=[('pdb files', '*.pdb *.ent'), ('all files', '*')],
            parent=self.parent)
        self.mol0.set(file_name)
        return
        
    def getPDBmol1(self):
        file_name = tkFileDialog.askopenfilename(
            title='PDB File', initialdir='',
            filetypes=[('pdb files', '*.pdb *.ent'), ('all files', '*')],
            parent=self.parent)
        self.mol1.set(file_name)
        return

    def getSizemol0(self):
        pqr_filename = 'mol0.pqr'
        psize = Psize()
        psize.runPsize(pqr_filename)
        #print(psize.getCharge())
        grid_points = psize.getFineGridPoints()
        cglen = psize.getCoarseGridDims()
        fglen = psize.getFineGridDims()
        self.dime0[0].set(grid_points[0])
        self.dime0[1].set(grid_points[1])
        self.dime0[2].set(grid_points[2])
        self.cglen0[0].set(cglen[0])
        self.cglen0[1].set(cglen[1])
        self.cglen0[2].set(cglen[2])
        self.fglen0[0].set(fglen[0])
        self.fglen0[1].set(fglen[0])
        self.fglen0[2].set(fglen[0])
        return

    def getSizemol1(self):
        pqr_filename = 'mol1.pqr'
        psize = Psize()
        psize.runPsize(pqr_filename)
        #print(psize.getCharge())
        grid_points = psize.getFineGridPoints()
        cglen = psize.getCoarseGridDims()
        fglen = psize.getFineGridDims()
        self.dime1[0].set(grid_points[0])
        self.dime1[1].set(grid_points[1])
        self.dime1[2].set(grid_points[2])
        self.cglen1[0].set(cglen[0])
        self.cglen1[1].set(cglen[1])
        self.cglen1[2].set(cglen[2])
        self.fglen1[0].set(fglen[0])
        self.fglen1[1].set(fglen[0])
        self.fglen1[2].set(fglen[0])
        return
    
    def getContacts(self):
        file_name = tkFileDialog.askopenfilename(
            title='Contacts File', initialdir='',
            filetypes=[('xml files', '*.xml'), ('all files', '*')],
            parent=self.parent)
        self.contacts.set(file_name)
        return
    
    def pdb2pqr(self):
        os.system("cp %s ./mol0.pdb" % self.mol0.get())
        os.system("cp %s ./mol1.pdb" % self.mol1.get())
        an = ''
        if self.pqr_assign_only.get(): an = '--assign-only'
        
        pqr_options = an + ' ' + self.pdb2pqr_opt.get() + ' --ff=' + self.pqr_ff.get()
        
        for i in ['mol0', 'mol1']:
            command = self.pdb2pqr_path.get() + '/pdb2pqr ' + \
                      pqr_options + ' ' + i + '.pdb ' + i + '.pqr'
            if (DEBUG > 0): print(command)
            print("::: Running pdb2pqr on " + i + " ...")
            rc = self.runcmd(command)
            if rc != 0:
                print("::: Failed: " + command)
        return

    def runAPBS(self):
        apbs_template = '''# APBS template for BD grids
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
'''

        for i in ['mol0', 'mol1']:
            pqr_filename = i + '.pqr'
            #psize = Psize()
            #psize.runPsize(pqr_filename)
            #print(psize.getCharge())
            #grid_points = psize.getFineGridPoints()
            #cglen = psize.getCoarseGridDims()
            #fglen = psize.getFineGridDims()

            if i == 'mol0':
                grid_points = [self.dime0[x].get() for x in range(3)]
                cglen = [self.cglen0[x].get() for x in range(3)]
                fglen = [self.fglen0[x].get() for x in range(3)]
            if i == 'mol1':
                grid_points = [self.dime1[x].get() for x in range(3)]
                cglen = [self.cglen1[x].get() for x in range(3)]
                fglen = [self.fglen1[x].get() for x in range(3)]
            if grid_points[0] == 0 or grid_points[1] == 0 or grid_points[2] == 0:
                print("::: " + i + " - no grid points defined!")
                return
                
            dx_filename = i
            fnout = i + '.in'
            fout = open(fnout, "w")
            fout.write(apbs_template % \
                       (pqr_filename, grid_points[0], grid_points[1], grid_points[2],
                        cglen[0], cglen[1], cglen[2],
                        fglen[0], fglen[1], fglen[2], 
                        self.apbs_mode.get(), self.bcfl.get(),
                        self.ion_plus_one_conc.get(), self.ion_plus_one_rad.get(),
                        self.ion_minus_one_conc.get(), self.ion_minus_one_rad.get(),
                        self.interior_dielectric.get(), self.solvent_dielectric.get(),
                        self.chgm.get(), self.sdens.get(), self.srfm.get(),
                        self.srad.get(), self.swin.get(), self.system_temp.get(),
                        dx_filename))
            fout.close()
            command = self.apbs_path.get() + '/apbs ' + i + '.in'
            if (DEBUG > 0):
                print(command)
                print(grid_points)
                print(cglen)
                print(fglen)
            print("::: Running apbs on " + i + " ...")
            rc = self.runcmd(command)
            if rc != 0:
                print("::: Failed: " + command)

        return

    def run_pqr2xml(self):
        for i in ['mol0', 'mol1']:
            command = self.bd_path.get() + '/pqr2xml < ' + i + '.pqr > ' + i + '-atoms.pqrxml'
            if (DEBUG > 0): print(command)
            print("::: Running pqr2xml on " + i + " ...")
            rc = self.runcmd(command)
            if rc != 0:
                print("::: Failed: " + command)
        return
                
    def make_rxn_criteria(self):
        command = self.bd_path.get() + '/make_rxn_pairs ' + \
                  '-nonred -mol0 mol0-atoms.pqrxml -mol1 mol1-atoms.pqrxml \
                  -ctypes ' + self.contacts.get() + ' -dist ' \
                  + str(self.rxn_distance.get()) + ' > mol0-mol1-rxn-pairs.xml'
        if (DEBUG > 0): print(command)
        print("::: Running make_rxn_pairs ...")
        rc = self.runcmd(command)
        if rc != 0:
            print("::: Failed: " + command)
        command = self.bd_path.get() + '/make_rxn_file ' + \
                  '-pairs mol0-mol1-rxn-pairs.xml -distance ' + \
                  str(self.rxn_distance.get()) + ' -nneeded ' + str(self.npairs.get()) + ' > mol0-mol1-rxns.xml'
        if (DEBUG > 0): print(command)
        print("::: Running make_rxn_file ...")
        rc = self.runcmd(command)
        if rc != 0:
            print("::: Failed: " + command)
        return
            
    def prepareInputFile(self):
        nam_simulation_template = '''
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

</root>

'''
        fout = open('input.xml', "w")
        fout.write(nam_simulation_template % \
                   (self.solvent_eps.get(), self.debyel.get(), self.ntraj.get(), self.nthreads.get(),
                    self.mol0_eps.get(), self.mol1_eps.get(), self.mindx.get(), self.ntrajo.get(), self.ncopies.get(),
                    self.nbincopies.get(), self.nsteps.get(), self.westeps.get(), self.maxnsteps.get()))

        fout.close()

        # FIXME check for .dx files
        command = 'PATH=' + self.bd_path.get() + ':${PATH}' + ' bd_top' + ' input.xml'
        if (DEBUG > 0): print(command)
        print("::: Running bd_top ...")
        rc = self.runcmd(command)
        if rc != 0:
            print("::: Failed: " + command)
        return

    def run_rxn_crit(self):
        self.run_pqr2xml()
        self.make_rxn_criteria()
        return
    
    def prepBD(self):
        self.prepareInputFile()
        return

    def runBD(self):
        return

    def execute(self, result):
        print 'Exiting BD Plugin ...'
        if __name__ == '__main__':
            self.parent.destroy()
        else:
            self.dialog.withdraw()
        print 'Done.'
        return

class Psize:
    def __init__(self):
        self.constants = {"CFAC":1.7, "FADD":20, "SPACE":0.50, "GMEMFAC":200, "GMEMCEIL":400,
                          "OFAC":0.1, "REDFAC":0.25, "TFAC_ALPHA":9e-5,
                          "TFAC_XEON":3e-4, "TFAC_SPARC": 5e-4}
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
        self.nsmall = [0,0,0]
        self.nfocus = 0

    def parseInput(self, filename):
        """ Parse input structure file in PDB or PQR format """
        file = open(filename, "r")
        self.parseLines(file.readlines())

    def parseLines(self, lines):
        """ Parse the lines """
        for line in lines:
            if string.find(line,"ATOM") == 0:
                subline = string.replace(line[30:], "-", " -")
                ## words = string.split(subline) ## this is a hack
                ## adhering to lovely PDB format definition (fixed space)
                #words = line[30:38], line[38:46], line[46:54], line[54:60], line[60:66], line[72:76], line[76:78] 
                words = line[30:38], line[38:46], line[46:54], line[54:63], line[63:69], line[72:76], line[76:78] 
                if len(filter(string.strip, words)) < 4:    
                    sys.stderr.write("Can't parse following line:\n")
                    sys.stderr.write("%s\n" % line)
                    sys.exit(2)
                    continue
                self.gotatom = self.gotatom + 1
                try:
                    self.q = self.q + float(words[3])
                except ValueError, ve:                    
                    print("Error parsing line", line)
                    #ATOM     12  CG  HIS A   0      41.299 139.172 108.417  1.00100.00           C
                    raise ve
                rad = float(words[4])
                center = []
                for word in words[0:3]:
                    center.append(float(word))
                for i in range(3):
                    self.minlen[i] = min(center[i]-rad, self.minlen[i])
                    self.maxlen[i] = max(center[i]+rad, self.maxlen[i])
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
        tn = [0,0,0]
        for i in range(3):
            tn[i] = int(flen[i]/self.constants["SPACE"] + 0.5)
            self.n[i] = 32*(int((tn[i] - 1) / 32.0 + 0.5)) + 1
            if self.n[i] < 33:
                self.n[i] = 33
        return self.n

    def setAll(self):
        """ Set up all of the things calculated individually above """
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
        """ Parse input PQR file and set parameters """
        self.parseInput(filename)
        self.setAll()

#############################################
#
#
# Create root window for testing.
#
#
##############################################
if __name__ == '__main__':
    
    class App:
        def my_show(self,*args,**kwargs):
            pass

    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)   
    app.root.title('Root window')

    widget = BDPlugin(app)
    app.root.mainloop()
