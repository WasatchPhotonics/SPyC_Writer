import os
import sys
import ast
import shutil
import logging
import datetime
import platform
import subprocess

import pytest
import pprint
import numpy as np
import spc_spectra as spc
from hypothesis import assume, settings, given, strategies as st

sys.path.append(os.path.join(os.getcwd(),"src")) # Adds higher directory to python modules path.

from SPyC_Writer.SPCFileWriter import SPCFileWriter
from SPyC_Writer.SPCEnums import *

from SPyC_Writer.common import RES_DESC_LIMIT, SRC_INSTRUMENT_LIMIT

HALF_PRECISION = 16
CWD = os.getcwd()
DATA_DIR = os.path.join(CWD,"tests","data")
REPO_URL = "https://github.com/lincosamides/spc-sdk.git"
TEMP_CLONE_DIR = "temp_repo"

# from hyptohesis docs. Allows for multidim lists with same length
# https://hypothesis.readthedocs.io/en/latest/data.html
@st.composite
def two_rectangle_lists(draw):
    list_len = draw(st.integers(min_value=2, max_value=10))
    number = st.floats(width=HALF_PRECISION)
    list = st.lists(number, min_size=2, max_size=list_len)
    xs = draw(list)
    ys = draw(list)
    assume(len(xs) == len(ys))
    return (xs,ys)

# flatmap on a tuple, first element is number of columns, second element is max number of rows
@st.composite
def rectangle_lists(draw):
    cols = draw(st.integers(min_value=2, max_value=10))
    rows = draw(st.integers(min_value=2, max_value=10))
    number = st.floats(width=HALF_PRECISION)
    row = st.lists(number, min_size=cols, max_size=cols)
    rect = draw(st.lists(row, min_size=2,max_size=rows))
    return rect

def get_spc_sdk_data():
    try:
        temp_dir = os.path.join(CWD,"tests", TEMP_CLONE_DIR)
        subprocess.run(["git", "clone", REPO_URL, temp_dir], check=True)
        shutil.move(os.path.join(temp_dir, "data"), DATA_DIR)
        if platform.system() == "Windows":
            subprocess.run(["rmdir", "/S", "/Q", temp_dir], shell=True, check=True)
        else:
            subprocess.run(["rm", "-rf", temp_dir], check=True)
    except Exception as e:
        print(f"Failed to clone repo or extract 'data': {e}")

# some sdk data use fields that can't be derived and aren't exposed, for example the subheader flags, see specific sdk test for explanation for thos cases
def compare_file_output(file_name: str, ignore_bytes: list[int] = None):
    if(ignore_bytes is None):
        ignore_bytes = []
    spc_xy = spc.File(os.path.join(DATA_DIR, file_name))
    bin_data = bytes()
    log_str = ""
    if(hasattr(spc_xy,"logbins")):
        bin_data = spc_xy.logbins
    if(hasattr(spc_xy,"logtxto")):
        log_str = spc_xy.logtxto
    log.debug(spc_xy.ftflg)
    log.debug(spc_xy.fversn)
    log.debug(f"fcmnt type is {type(spc_xy.fcmnt)}")
    cmnt = ast.literal_eval(spc_xy.fcmnt).decode("utf-8")
    custom = spc_xy.fcatxt.decode("utf-8").split("\x00")[:3]
    log.debug(f"{spc_xy.fdate=}")
    if(spc_xy.fdate > 0): # s_evenx shows date can be null and python doesn't support year 0 and errors
        date = datetime.datetime(minute=spc_xy.minute, hour=spc_xy.hour, day=spc_xy.day, month=spc_xy.month, year=spc_xy.year)
    else:
        date = None
    writer = SPCFileWriter(
        file_type=SPCFileType(int.from_bytes(spc_xy.ftflg)),
        num_pts=len(spc_xy.x),
        compress_date=date,
        file_version=int.from_bytes(spc_xy.fversn),
        experiment_type=spc_xy.fexper,
        exponent=spc_xy.fexp,
        first_x=spc_xy.ffirst,
        last_x=spc_xy.flast,
        x_units=spc_xy.fxtype,
        y_units=spc_xy.fytype,
        z_units=spc_xy.fztype,
        res_desc=spc_xy.fres.decode("utf-8"),
        src_instrument_desc=spc_xy.fsource.decode("utf-8"),
        custom_units=custom,
        memo=cmnt, # spc doesn't decode right so you get "b'your msg\x00\x00'", need to strip byte type and interior ""
        sample_inject=spc_xy.fsampin,
        method_file=spc_xy.fmethod.decode("utf-8"),
        spectra_mod_flag=spc_xy.fmods,
        process_code=int.from_bytes(spc_xy.fprocs),
        z_subfile_inc=spc_xy.fzinc,
        num_w_planes=spc_xy.fwplanes,
        w_plane_inc=spc_xy.fwinc,
        w_units=int.from_bytes(spc_xy.fwtype),
        log_data=bin_data,
        log_text=log_str
    )
    if spc_xy.dat_fmt.endswith('-xy'):
        writer.write_spc_file("testOutput.spc", np.asarray([s.y for s in spc_xy.sub]), np.asarray([s.x for s in spc_xy.sub]))
    else:
        writer.write_spc_file("testOutput.spc", spc_xy.sub[0].y, spc_xy.x)

    with open(os.path.join(DATA_DIR,file_name), "rb") as ref, open("testOutput.spc", "rb") as output:
        refB = (b for b in ref.read())
        outB = (b for b in output.read())
        diffs = [(idx, pair[0], pair[1]) for idx, pair in enumerate(zip(refB, outB)) if pair[0] != pair[1] and idx not in ignore_bytes]
        if(len(diffs) != 0):
            pprint.pprint(diffs)
            assert False, f"File contents did not match for written output based on {file_name}"

log = logging.getLogger(__name__)
@pytest.fixture
def setup():
    if not os.path.isdir(TestWritingParse.OUTPUT_DIR):
        os.makedirs(TestWritingParse.OUTPUT_DIR)
    if not os.path.isdir(DATA_DIR):
        get_spc_sdk_data()
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

    @given(st.lists(st.floats(width=HALF_PRECISION, allow_nan=False), min_size=2, max_size=30))
    @settings(deadline=None, max_examples=30)
    def test_y_format(self, y_values: list[float]) -> None:
        writer = SPCFileWriter(SPCFileType.DEFAULT)
        y_values = np.asarray(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_y_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "y format spc output array doesn't match parse"

    @given(two_rectangle_lists())
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

    @given(rectangle_lists())
    def test_yyy_format(self, y_values: list[float]) -> None:
        """
        File type is tested even though KIA reports invalid while others parse successfully
        """
        writer = SPCFileWriter(SPCFileType.TMULTI)
        y_values = np.array(y_values)
        file_output = os.path.join(self.OUTPUT_DIR, f"spc_yyy_output{datetime.datetime.now().strftime('%H_%M_%S_%f')}.spc")
        writer.write_spc_file(file_output, y_values=y_values)
        assert self.comapre_spc_file_array(file_output, y_values), "yyy format spc output array doesn't match parse"

    @given(st.lists(st.floats(width=HALF_PRECISION, allow_nan=False), min_size=2, max_size=30),
           st.integers(min_value=SPCXType.SPCXArb, max_value=SPCXType.SPCXAngst),
           st.integers(min_value=SPCYType.SPCYArb, max_value=SPCYType.SPCYSRot),
           st.integers(min_value=SPCXType.SPCXArb, max_value=SPCXType.SPCXAngst), # metoffice docs explicitly say x and z use same enum
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

    def test_evenx_sdk(self):
        # ignored bytes
        # 512 - This is the subheader flag bytes, here it indicates modified by arithmetic
        # None of the specs elaborates what this means, I'm assuming maybe preprocesssing was done before saving
        # since this can't be derived ignore for now, will likely support passing end user made subheaders in the future
        compare_file_output("s_evenx.spc", [512])

    def test_xy_sdk(self):
        # ignored bytes
        # 2576 - This is the subheader subnpts for sub 0, spc_spectra indicates this value is 4
        # no clue what's going on here, this isn't TXYXYS so spec says this should be null
        # also the number of points is 512 and since subnpts is usually x AND y this would be 1024
        # 1024 would be 00 00 04 00 though not 04 00 00 00 so since only off by one byte that isn't clear from the spec I call this good
        compare_file_output("s_xy.spc", [2576])