import os
import sys
import datetime

import pytest
import numpy as np
import spc_spectra as spc
from hypothesis import given, strategies as st

print(os.getcwd())

from src.SPyC_Writer.SPCFileWriter import SPCFileWriter
from src.SPyC_Writer.SPCEnums import SPCFileType

sys.path.append(os.path.join(os.getcwd(),"src")) # Adds higher directory to python modules path.

@pytest.fixture
def setup():
    if not os.path.isdir(TestWritingParse.OUTPUT_DIR):
        os.makedirs(TestWritingParse.OUTPUT_DIR)
    yield
    TestWritingParse.clean_up_files()

@pytest.mark.usefixtures("setup")
class TestWritingParse:
    OUTPUT_DIR = os.path.join(os.getcwd(), "SPC_Output")

    @classmethod
    def clean_up_files(cls):
        for file_name in os.listdir(cls.OUTPUT_DIR):
            os.remove(os.path.join(cls.OUTPUT_DIR, file_name))

    def comapre_spc_file_array(self, file_name: str, original_arr: list):
        data = spc.File(file_name)
        match len(data.sub):
            case 0:
                file_ys = np.empty(shape=(0,))
            case 1:
                file_ys = data.sub[0].y
            case _:
                file_ys = np.asarray([subf.y for subf in data.sub])

        if np.array_equal(file_ys, original_arr, equal_nan = True):
            return True
        else:
            print(f"arrays {file_ys} and {original_arr}")
            return False

    @given(st.lists(st.floats(width=16), min_size=2, max_size=30))
    def test_y_format(self, y_values: list[float]) -> None:
        writer = SPCFileWriter(SPCFileType.DEFAULT)
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "y format spc output array doesn't match parse"


    @given(
        st.lists(
            st.lists(
                st.floats(width=16), min_size=10, max_size=10),
                min_size=2))
    def test_yyy_format(self, y_values: list[float]) -> None:
        writer = SPCFileWriter(SPCFileType.TMULTI)
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "yyy format spc output array doesn't match parse"
