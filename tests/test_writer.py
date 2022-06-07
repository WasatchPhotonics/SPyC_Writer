import os
import sys
import logging
import datetime

import pytest
import numpy as np
import spc_spectra as spc
from hypothesis import settings, given, strategies as st

sys.path.append(os.path.join(os.getcwd(),"src")) # Adds higher directory to python modules path.

from SPyC_Writer.SPCFileWriter import SPCFileWriter
from SPyC_Writer.SPCEnums import SPCFileType


# from hyptohesis docs. Allows for multidim lists with same length
# https://hypothesis.readthedocs.io/en/latest/data.html
rectangle_lists = st.floats(width=16).flatmap(
    lambda n: st.lists(st.lists(st.floats(width=16), min_size=2, max_size=8), min_size=2,max_size=2)
).filter(lambda x: len(x[0]) == len(x[1]))

log = logging.getLogger(__name__)
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

    def comapre_spc_file_array(self, file_name: str, original_arr: list, axis: str = "y"):
        data = spc.File(file_name)
        match len(data.sub):
            case 0:
                file_arr = np.empty(shape=(0,))
            case 1:
                if axis == "x-xy":
                    file_arr = data.x
                elif axis == "xy":
                    file_arr = data.sub[0].x
                else:
                    file_arr = data.sub[0].y
            case _:
                if axis == "x-xyxy":
                    file_arr = data.x # should be one common x array
                elif axis == "xyxy":
                    file_arr = np.asarray([subf.x for subf in data.sub])
                else:
                    file_arr = np.asarray([subf.y for subf in data.sub])

        if np.array_equal(file_arr, original_arr, equal_nan = True):
            return True
        else:
            log.debug(f"failed values arrays  file parse: {file_arr} and  generated value: {original_arr}")
            log.debug(f"file info is {data.__dict__}\n\n")
            log.debug(f"sub is {data.sub[0].__dict__}")
            return False

    @given(st.lists(st.floats(width=16, allow_infinity=False, allow_nan=False), min_size=2, max_size=30))
    @settings(deadline=None, max_examples=30)
    def test_y_format(self, y_values: list[float]) -> None:
        writer = SPCFileWriter(SPCFileType.DEFAULT)
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "y format spc output array doesn't match parse"


    @given(rectangle_lists)
    @settings(deadline=None, max_examples=30)
    def test_xy_format(self, arr_values: list[float]) -> None:
        log.debug("testing xy format")
        y_values = arr_values[0]
        x_values = arr_values[1]
        writer = SPCFileWriter(SPCFileType.TXVALS)
        y_values = np.asarray(y_values)
        x_values = np.asarray(x_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values, x_values=x_values)
        assert self.comapre_spc_file_array(file_output, x_values, axis="x-xy") and self.comapre_spc_file_array(file_output, y_values), "xy format spc output array doesn't match parse"

    @given(
        st.lists(
            st.lists(
                st.floats(width=16), min_size=10, max_size=10),
                min_size=2))
    @settings(deadline=None, max_examples=30)
    def test_yyy_format(self, y_values: list[float]) -> None:
        """
        File type is tested even though KIA lists invalid while others parse successfully
        """
        writer = SPCFileWriter(SPCFileType.TMULTI)
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "yyy format spc output array doesn't match parse"
