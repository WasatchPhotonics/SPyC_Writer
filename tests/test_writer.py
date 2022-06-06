import os
import sys
import datetime

import spc_spectra as spc
import numpy as np
from hypothesis import given, strategies as st

print(os.getcwd())

from src.SPyC_Writer.SPCFileWriter import SPCFileWriter
from src.SPyC_Writer.SPCEnums import SPCFileType, SPCTechType

sys.path.append(os.path.join(os.getcwd(),"src")) # Adds higher directory to python modules path.
class TestWritingParse:
    OUTPUT_DIR = os.path.join(os.getcwd(), "SPC_Output")

    @given(st.lists(st.floats(width=16), min_size=2, max_size=30))
    def test_y_format(self, y_values: list[float]) -> bool:
        writer = SPCFileWriter(SPCTechType.SPCTechGen,
                               )
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, "my_output.spc")
        if not os.path.isdir(self.OUTPUT_DIR):
            os.makedirs(self.OUTPUT_DIR)
        writer.write_spc_file(file_output, y_values=y_values)

        data = spc.File(file_output)

        assert np.array_equal(data.sub[0].y, y_values, equal_nan = True)

