"""
MIT License

Copyright (c) 2022 Wasatch Photonics

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
# See the following refernces for details on the .spc file format implementation
# https://github.com/bz-dev/spc-sdk
# https://www.yumpu.com/en/document/read/40416248/a-brief-guide-to-spc-file-format-and-using-gspcio
# https://ensembles-eu.metoffice.gov.uk/met-res/aries/technical/GSPC_UDF.PDF

# This file is meant to write the SPC file format.
# It has been tested with file formats of Y, XY, XYYY, and XYXYXY
# It was verified by parsing the outputs using the following tools
# https://github.com/rohanisaac/spc
# https://www.effemm2.de/spectragryph/
# https://sciencesolutions.wiley.com/knowitall-spectroscopy-software/
# XYXYXY files are successfully parsed by spectragryph and rohanisaac
# KIA states they are invalid though.

import math
import logging
from struct import pack
from datetime import datetime
from dataclasses import field

import numpy as np

from .SPCLog import SPCLog
from .SPCDate import SPCDate
from .SPCHeader import SPCHeader
from .SPCSubheader import SPCSubheader
from .SPCEnums import SPCFileType, SPCModFlags, \
   SPCTechType, SPCXType, SPCYType

log = logging.getLogger(__name__)

class SPCFileWriter:
    """
    Created based on the spc file format.
    the primary method for writing files is via the write_spc_file method
    x data points and y data points are passed in via 2 different 2d arrays
    Only new format is currently supported.
    """
    def __init__(self,
                 file_type: SPCFileType,
                 num_pts: int = 0, # according to the docs if given num_pts, first_x, and last_x then it will calculate an evenly spaced x axis
                 compress_date: datetime = datetime.now(),
                 file_version: int = 0x4B,
                 experiment_type: SPCTechType = SPCTechType.SPCTechGen,
                 exponent: int = 0, # available but not supported
                 first_x: int = 0,
                 last_x: int = 0,
                 x_units: SPCXType = SPCXType.SPCXArb,
                 y_units: SPCYType = SPCYType.SPCYArb, 
                 z_units: SPCXType = SPCXType.SPCXArb,
                 res_desc: str = "",
                 src_instrument_desc: str = "",
                 custom_units: list[str] = field(default_factory=list),
                 memo: str = "",
                 custom_axis_str: str = "",
                 spectra_mod_flag: SPCModFlags = SPCModFlags.UNMOD,
                 z_subfile_inc: float = 1.0,
                 num_w_planes: float = 0,
                 w_plane_inc: float = 1.0,
                 w_units: SPCXType = SPCXType.SPCXArb,
                 log_data: bytes = bytes(),
                 log_text: str = "",
                 )-> None:
        """
        According to the formatting document the following file types are most common.
        For the respective file type, the corresponding initialization parameters 
        and arrays to the write function should be passed at a minimum
        ------------------
        Single File Even X
        ------------------
        Pass num_pts, first_x, last_x, and a 1D numpy array for the Y values

        ------------------
        Multifile Even X
        ------------------
        Pass num_pts, first_x, last_x, and a 2D numpy array for the Y values, each row represents a new subfile

        ------------------
        Single File Uneven X
        ------------------
        Pass a 1D nump array for the X values and a 1D numpy array for the Y values

        ------------------
        Multifile Uneven X common X
        ------------------
        Pass a 1D numpy array for the X values and a 2D numpy array for the Y values

        ------------------
        Multifile Uneven X common X
        ------------------
        Pass a 2D numpy array for the X values and a 2D numpy array for the Y values
        """
        self.file_type = file_type
        self.num_pts = num_pts
        self.compress_date = compress_date
        self.file_version = file_version
        self.experiment_type = experiment_type
        self.exponent = exponent
        self.first_x = first_x
        self.last_x = last_x
        self.x_units = x_units
        self.y_units = y_units
        self.z_units = z_units
        self.res_desc = res_desc
        self.src_instrument_desc = src_instrument_desc
        self.custom_units = custom_units
        self.memo = memo
        self.custom_axis_str = custom_axis_str
        self.spectra_mod_flag = spectra_mod_flag
        self.z_subfile_inc = z_subfile_inc
        self.num_w_planes = num_w_planes
        self.w_plane_inc = w_plane_inc
        self.w_units = w_units
        self.log_data = log_data
        self.log_text = log_text

    def validate_inputs(self, x_values, y_values, z_values, w_values) -> bool:
        if x_values.size != 0 and y_values.size != 0 and not (x_values.shape[-1] == y_values.shape[-1]):
            log.error(f"got x and y values of different size. Arrays must be so same length.")
            return False
        if x_values.size == 0:
            if self.file_type == SPCFileType.TXVALS:
                log.error(f"no x values received but file type is a shared x values type")
                return False
            self.first_x = 0
            self.last_x = len(y_values)
        else:
            self.first_x = np.amin(x_values)
            self.last_x = np.amax(x_values)
        return True


    def write_spc_file(self,
                       file_name: str, 
                       y_values: np.ndarray,
                       x_values: np.ndarray = np.empty(shape=(0)),
                       z_values: np.ndarray = np.empty(shape=(0)),
                       w_values: np.ndarray = np.empty(shape=(0)),
                       ) -> bool:
        file_output = b""
        generate_log = False
        if not self.validate_inputs(x_values, y_values, z_values, w_values):
            log.error(f"invalid inputs, returning false")
            return False
        if not (self.file_type & SPCFileType.TMULTI):
            points_count = len(y_values)
        elif self.file_type & SPCFileType.TMULTI and not (self.file_type & SPCFileType.TXYXYS):
            points_count = len(y_values[0]) # since x values are evenly spaced y values shouldn't be jagged array
        else:
            # num_points for XYXYXY is instead supposed to be the byte offset to the directory
            # or null and there is no directory
            points_count = 0 
            self.exponent = self.calculate_exponent(x_values, y_values)
        if len(y_values.shape) == 1:
            num_traces = 1
        else:
            num_traces = y_values.shape[0]

        if w_values.size != 0:
            if num_traces % len(w_values) != 0:
                log.error(f"w_values should divide evenly into the number of sub files")
                return False

        if len(self.log_data) > 0 or len(self.log_text) > 0:
            generate_log = True

        By_values = self.convert_points(y_values, np.single)
        Bx_values = self.convert_points(x_values, np.single)

        header = SPCHeader(
            file_type = self.file_type,
            num_points = points_count,
            compress_date = SPCDate(self.compress_date),
            x_values = x_values,
            y_values = y_values,
            experiment_type = self.experiment_type,
            first_x = self.first_x,
            last_x = self.last_x,
            num_subfiles = num_traces,
            x_units = self.x_units,
            y_units = self.y_units,
            z_units = self.z_units,
            res_desc = self.res_desc,
            src_instrument_desc = self.src_instrument_desc,
            memo = self.memo,
            custom_axes = self.custom_units,
            spectra_mod_flag = self.spectra_mod_flag,
            z_subfile_inc = self.z_subfile_inc,
            num_w_planes = self.num_w_planes,
            w_plane_inc = self.w_plane_inc,
            w_units = self.w_units,
            generate_log = generate_log
            )
        file_header = header.generate_header()
        file_output = b"".join([file_output, file_header])

        if (self.file_type & SPCFileType.TXVALS) and not (self.file_type & SPCFileType.TXYXYS):
            file_output = b"".join([file_output, Bx_values]) # x values should be a flat array so shouldn't be any issues with this

        dir_pointers = []
        for i in range(num_traces):
            subfile = b""
            w_val = 0
            if w_values.size != 0:
                w_val = w_values[math.floor(i/w_values(len))]
            if SPCFileType.TXYXYS & self.file_type:
                points_count = len(y_values[i]) # inverse from header, header it is 0, here it's the length of a specific y input
            sub_header = b""
            if len(z_values) == 0:
                z_val = 0
            elif len(z_values) == 1:
                z_val = z_values[0]
            else:
                try:
                    z_val = z_values[i]
                except:
                    z_val = 0

            subheader = SPCSubheader(start_z = z_val,
                                   sub_index = i, 
                                   num_points = points_count,
                                   w_axis_value = w_val)
            if self.file_type & SPCFileType.TXYXYS:
                bx = self.convert_points(x_values[i], "<f4") #self.convert_points(np.ones(shape=(1952,)), "<f4")#self.convert_points(x_values[i], "<f4")
                by = self.convert_points(y_values[i], "<f4")
                sub_head = subheader.generate_subheader()
                subfile = b"".join([sub_head, bx, by])
            elif self.file_type & SPCFileType.TMULTI and not (self.file_type & SPCFileType.TXYXYS):
                sub_head = subheader.generate_subheader()
                subfile = b"".join([sub_head, self.convert_points(y_values[i], "<f4")])
            else:
                sub_head = subheader.generate_subheader()
                subfile = b"".join([sub_head, By_values])

            pointer = self.generate_dir_pointer(len(file_output), len(subfile), z_val)
            dir_pointers.append(pointer)
            file_output = b"".join([file_output, subfile])

        if self.file_type & SPCFileType.TXVALS and self.file_type & SPCFileType.TXYXYS:
            file_output = b"".join([file_output, b"".join(dir_pointers)])

        if generate_log:
            log.debug(f"generating spc log")
            log_head = SPCLog(self.log_data, self.log_text)
            log_header = log_head.generate_log_header()
            file_output = b"".join([file_output, log_header, self.log_data, self.log_text.encode()])

        try:
            with open(file_name, 'wb') as f:
                f.write(file_output)
                return True
        except Exception as e:
            log.error(f"error in creating spc file of {e}")
            return False

    def generate_dir_pointer(self, offset: int, sub_size: int, z_val: float) -> bytes:
        Boffset = offset.to_bytes(4, byteorder="little")
        Bsub_size = sub_size.to_bytes(4, byteorder="little")
        Bz_value = pack("f", z_val)
        pointer = b"".join([Boffset, Bsub_size, Bz_value])
        return pointer

    def calculate_exponent(self, x_values: np.ndarray, y_values) -> int:
        """
        Exploits the fact that we are on a 64 bit architecture.
        any value greater than 1 results in a number greater than what a 32 bit int can hold.
        We can hold that value and keep dividing by 2 until we get smaller than 32 bit.
        This will then be our exponent. Since that final integer over 2^32 results in some decimal.
        The left shift of that decimal should be our original number.
        A max_x that is only a very small decimal should be very rare so shouldn't need to worry about
        shifting the decimal right.
        """
        max_x = abs(np.amax(x_values))
        max_y = abs(np.amax(y_values))
        max_num = max([max_x, max_y])
        product = max_num * (2**32)
        exponent = 0
        while product > 2**32:
            product /= 2
            exponent += 1
            if exponent > 127:
                log.error(f"exponent is only a signed byte. Cannot store greater than 127")
                raise
        return exponent

    def convert_points(self, data_points: np.ndarray, conversion: np.dtype) -> bytes:
        """
        Takes a numpy array of data points and converts them to single precision floats.
        Then converts them to a string of bytes. Currently only supports the single precision.
        Does not support the spc specific exponent representation of floating point numbers.
        """
        data_points = data_points.astype(conversion)
        return data_points.tobytes()

            


