# vim: set expandtab shiftwidth=4 softtabstop=4:

# === UCSF ChimeraX Copyright ===
# Copyright 2016 Regents of the University of California.
# All rights reserved.  This software provided pursuant to a
# license agreement containing restrictions on its disclosure,
# duplication and use.  For details see:
# http://www.rbvi.ucsf.edu/chimerax/docs/licensing.html
# This notice must be embedded in or attached to all copies,
# including partial copies, of the software or any revisions
# or derivations thereof.
# === UCSF ChimeraX Copyright ===

from chimerax.core.tools import ToolInstance
from Qt.QtWidgets import QLabel, QLineEdit, QGridLayout, QPushButton
from Qt.QtWidgets import QHBoxLayout, QComboBox, QCheckBox, QTableWidget, QHeaderView, QRadioButton
from Qt.QtGui import QDesktopServices
from Qt.QtCore import QUrl

from .proteome import Proteome

class ScreeningTool(ToolInstance):

    # Inheriting from ToolInstance makes us known to the ChimeraX tool mangager,
    # so we can be notified and take appropriate action when sessions are closed,
    # saved, or restored, and we will be listed among running tools and so on.
    #
    # If cleaning up is needed on finish, override the 'delete' method
    # but be sure to call 'delete' from the superclass at the end.

    SESSION_ENDURING = False    # Does this instance persist when session closes
    SESSION_SAVE = True         # We do save/restore in sessions
    # help = "help:user/tools/tutorial.html"
                                # Let ChimeraX know about our help page

    def __init__(self, session, tool_name):
        # 'session'   - chimerax.core.session.Session instance
        # 'tool_name' - string

        # Initialize base class.
        super().__init__(session, tool_name)

        # Set name displayed on title bar (defaults to tool_name)
        # Must be after the superclass init, which would override it.
        self.display_name = "ANNOTEX - Annotating AF proteomes"

        # Create the main window for our tool.  The window object will have
        # a 'ui_area' where we place the widgets composing our interface.
        # The window isn't shown until we call its 'manage' method.
        #
        # Note that by default, tool windows are only hidden rather than
        # destroyed when the user clicks the window's close button.  To change
        # this behavior, specify 'close_destroys=True' in the MainToolWindow
        # constructor.
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)

        # We will be adding an item to the tool's context menu, so override
        # the default MainToolWindow fill_context_menu method
        self.tool_window.fill_context_menu = self.fill_context_menu

        self.proteome = Proteome(self)
        self._build_ui()
        

    def _build_ui(self):

        layout = QGridLayout()
        
        layout_text1 = QHBoxLayout()
        layout_text1.addWidget(QLabel('Data folder:'))
        self.data_folder = QLineEdit('/PATH/TO/ANNOTATION_DATA/')
        layout_text1.addWidget(self.data_folder)
        layout.addLayout(layout_text1, 0, 0, 1, 3)
        
        
        submit = QPushButton('Search')
        submit.clicked.connect(self.proteome.read_proteins)
        layout.addWidget(submit, 3, 1, 1, 1)
        
        
        self.result_list = QComboBox()
        self.result_list.activated.connect(self.proteome.open_protein)
        layout.addWidget(self.result_list, 5, 0, 1, 3)
        
        previous = QPushButton('Previous')
        previous.clicked.connect(self.proteome.previous_protein)
        layout.addWidget(previous, 6, 0, 1, 1)
        next_protein = QPushButton('Next')
        next_protein.clicked.connect(self.proteome.next_protein)
        layout.addWidget(next_protein, 6, 2, 1, 1)
        
        self.domain_list = QComboBox()
        layout.addWidget(self.domain_list, 7, 0, 1, 3)
        
        layout_annotation = QHBoxLayout()
        layout_annotation.addWidget(QLabel('Annotation:'))
        self.annotation = QLineEdit('N/A')
        layout_annotation.addWidget(self.annotation)
        layout.addLayout(layout_annotation, 9, 0, 1, 3)
        
        
        self.eggnog_name_label = QLabel('Eggnog: N/A')
        self.eggnog_name_label.setWordWrap(True)
        layout.addWidget(self.eggnog_name_label, 10, 0, 1, 3)
        
        self.eggnog_desc_label = QLabel('Eggnog: N/A')
        self.eggnog_desc_label.setWordWrap(True)
        layout.addWidget(self.eggnog_desc_label, 11, 0, 1, 3)
        
        self.eggnog_e_label = QLabel('Eggnog: N/A')
        self.eggnog_e_label.setWordWrap(True)
        layout.addWidget(self.eggnog_e_label, 12, 0, 1, 3)
        
        self.msa_label = QLabel('MSA Coverage: N/A')
        layout.addWidget(self.msa_label, 13, 0, 1, 1)
        
        self.deeptmhmm_label = QLabel('DeepTMHMM: N/A')
        layout.addWidget(self.deeptmhmm_label, 14, 0, 1, 1)
        
        
        self.foldseek_table = QTableWidget()
        self.foldseek_table.setRowCount(1)
        self.foldseek_table.setColumnCount(1)
        self.foldseek_table.horizontalHeader().setStretchLastSection(True)
        self.foldseek_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.foldseek_table.setHorizontalHeaderLabels(['Structural matches'])
        layout.addWidget(self.foldseek_table, 15, 0, 1, 3)
        self.foldseek_table.itemClicked.connect(self.proteome.open_match)

        
        self.blast_table = QTableWidget()
        self.blast_table.setRowCount(1)
        self.blast_table.setColumnCount(1)
        self.blast_table.horizontalHeader().setStretchLastSection(True)
        self.blast_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch)
        self.blast_table.setHorizontalHeaderLabels(['Diamond'])
        layout.addWidget(self.blast_table, 16, 0, 1, 3)
        self.blast_table.itemClicked.connect(self.test)
        
        
        self.tool_window.ui_area.setLayout(layout)
        self.tool_window.manage('side')
    
    def test(self, tableItem):
        print(tableItem.text())
        
    def link(self, linkStr):
        QDesktopServices.openUrl(QUrl(linkStr))

    def fill_context_menu(self, menu, x, y):
        # Add any tool-specific items to the given context menu (a QMenu instance).
        # The menu will then be automatically filled out with generic tool-related actions
        # (e.g. Hide Tool, Help, Dockable Tool, etc.) 
        #
        # The x,y args are the x() and y() values of QContextMenuEvent, in the rare case
        # where the items put in the menu depends on where in the tool interface the menu
        # was raised.
        from Qt.QtWidgets import QAction
        clear_action = QAction("Clear", menu)
        clear_action.triggered.connect(lambda *args: self.data_folder.clear())
        menu.addAction(clear_action)

    def take_snapshot(self, session, flags):
        return {
            'version': 1,
            'data_text': self.data_folder.text(),
            'protein_text': self.protein_search.text()
        }

    @classmethod
    def restore_snapshot(class_obj, session, data):
        # Instead of using a fixed string when calling the constructor below, we could
        # have saved the tool name during take_snapshot() (from self.tool_name, inherited
        # from ToolInstance) and used that saved tool name.  There are pros and cons to
        # both approaches.
        inst = class_obj(session, "annotex")
        inst.line_edit.setText(data['data_text'])
        inst.protein_search.setText(data['protein_text'])
        return inst
