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

from SPyC_Writer.common import RES_DESC_LIMIT, SRC_INSTRUMENT_LIMIT


# from hyptohesis docs. Allows for multidim lists with same length
# https://hypothesis.readthedocs.io/en/latest/data.html
two_rectangle_lists = st.integers(min_value=2, max_value=10).flatmap(
    lambda n: st.lists(st.lists(st.floats(width=16), min_size=2, max_size=n), min_size=2,max_size=2)
).filter(lambda x: len(x[0]) == len(x[1]))

# flatmap on a tuple, first element is number of columns, second element is max number of rows
rectangle_lists = st.tuples(st.integers(min_value=2, max_value=10), st.integers(min_value=2, max_value=10)).flatmap(
    lambda n : st.lists(st.lists(st.floats(width=16), min_size=n[0], max_size=n[0]), min_size=2,max_size=n[1])
)

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
        if len(data.sub) == 0:
            file_arr = np.empty(shape=(0,))
        elif len(data.sub) == 1:
            if axis == "x-xy":
                file_arr = data.x
            elif axis == "xy":
                file_arr = data.sub[0].x
            else:
                file_arr = data.sub[0].y
        else:
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

    @given(st.lists(st.floats(width=16, allow_nan=False), min_size=2, max_size=30))
    @settings(deadline=None, max_examples=30)
    def test_y_format(self, y_values: list[float]) -> None:
        writer = SPCFileWriter(SPCFileType.DEFAULT)
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "y format spc output array doesn't match parse"

    @given(two_rectangle_lists)
    @settings(deadline=None, max_examples=30)
    def test_xy_format(self, arr_values: list[float]) -> None:
        log.debug("testing xy format")
        y_values = arr_values[0]
        x_values = arr_values[1]
        writer = SPCFileWriter(SPCFileType.TXVALS)
        y_values = np.asarray(y_values)
        x_values = np.asarray(x_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_xy_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values, x_values=x_values)
        assert self.comapre_spc_file_array(file_output, x_values, axis="x-xy") and self.comapre_spc_file_array(file_output, y_values), "xy format spc output array doesn't match parse"

    @given(rectangle_lists)
    @settings(deadline=None, max_examples=30)
    def test_yyy_format(self, y_values: list[float]) -> None:
        """
        File type is tested even though KIA reports invalid while others parse successfully
        """
        writer = SPCFileWriter(SPCFileType.TMULTI)
        y_values = np.array(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_yyy_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "yyy format spc output array doesn't match parse"

    @given(st.lists(st.floats(width=16, allow_nan=False), min_size=2, max_size=30),
           st.integers(min_value=0, max_value=31),
           st.integers(min_value=0, max_value=30),
           st.integers(min_value=0, max_value=31),
           # hypothesis generates all sorts of unicode, which we expect to be invalid
           # check cases that should always be valid by generating ASCII
           st.text(alphabet=st.characters(min_codepoint = 32, max_codepoint = 126), max_size = RES_DESC_LIMIT),
           st.text(alphabet=st.characters(min_codepoint = 32, max_codepoint = 126), max_size = SRC_INSTRUMENT_LIMIT),
           )
    @settings(deadline=None, max_examples=30)
    def test_header_data(self,
                         y_values: list[float],
                         x_units: int,
                         y_units: int,
                         z_units: int,
                         res_desc: str,
                         src_instrument_desc: str,
                         ) -> None:
        time = datetime.datetime.now()
        writer = SPCFileWriter(
            SPCFileType.DEFAULT,
            compress_date = time,
            x_units = x_units,
            y_units = y_units,
            z_units = z_units,
            res_desc = res_desc,
            src_instrument_desc = src_instrument_desc,
            )
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        data = spc.File(file_output)
        comparison_inputs = [time.year, time.month, time.day, x_units, y_units, z_units, res_desc, src_instrument_desc]
        comparison_outputs = [data.year, data.month, data.day, data.fxtype, data.fytype, data.fztype, data.fres.decode("utf-8").replace("\x00",""), data.fsource.decode("utf-8").replace("\x00","")]
        comparison = [inputs == outputs for inputs, outputs in zip(comparison_inputs, comparison_outputs)]
        assert all(comparison), f"input data {comparison_inputs}\n and output data {comparison_outputs}\nheaders didn't match"

