import logging
from struct import pack
from dataclasses import dataclass, field

import numpy as np

from .SPCDate import SPCDate
from .common import fit_byte_block, RES_DESC_LIMIT, SPARE_LIMIT, \
    MEMO_LIMIT, AXES_LIMIT, SRC_INSTRUMENT_LIMIT, RESERVE_LIMIT, METHOD_FILE_LIMIT
from .SPCEnums import SPCFileType, SPCModFlags, SPCPostDisposition, \
   SPCProcessCode, SPCTechType, SPCXType, SPCYType

log = logging.getLogger(__name__)

@dataclass
class SPCHeader:
    """
    Corresponds to the spc file header.
    minimum fields are file type, num points and compression date.
    Source code comments detail corresponding spc.h name.
    """
    file_type: SPCFileType # (ftflags)
    num_points: int # (fnpts) or directory position for XYXYXY files
    compress_date: SPCDate # (fdate)
    x_values: np.ndarray = np.empty(shape=(0))
    y_values: np.ndarray = np.empty(shape=(0))
    file_version: int = 0x4B # (fversn)
    experiment_type: SPCTechType = SPCTechType.SPCTechGen # (fexper)
    exponent: int = -128 # (fexp)
    first_x: float = 0 # (ffirst)
    last_x: float = 0 # (flast)
    num_subfiles: int = 0 # (fnsub)
    x_units: SPCXType = SPCXType.SPCXArb # (fxtype)
    y_units: SPCYType = SPCYType.SPCYArb # (fytype)
    z_units: SPCXType = SPCXType.SPCXArb # (fztype)
    post_disposition: SPCPostDisposition = SPCPostDisposition.PSTDEFT # (fpost) should normally be null according to old format doc
    res_desc: str = "" # (fres)
    src_instrument_desc: str = "" # (fsource)
    peak_point: float = 0 # (fpeakpt) Interferogram peak point, associated with y_units = 2 
    memo: str = "" # (fcmnt)
    custom_axes: list[str] = field(default_factory=list) # (fcatxt)
    spectra_mod_flag: SPCModFlags = SPCModFlags.UNMOD # (fmods)
    # For proc codes see https://github.com/bz-dev/spc-sdk/blob/master/GRAMSDDE.H#L104
    # There are two defines for the value 1. This is explained more in their comments.
    # Rather than repeat this, I just use PPCOMP for a default of 1.
    process_code: SPCProcessCode = SPCProcessCode.PPCOMP 
    calib_plus_one: int = b"\x00" # (flevel) old format doc says galactic internal use and should be null
    sample_inject: int = b"\x00\x00" # (fsampin) spc.h lists 1 as valid, old format doc says only for galactic internal use and should be null
    data_mul: float = b"\x00\x00\x00\x00" # (ffactor) old format doc says galactic internal use only and should be null
    method_file: str = b"\x00" # (fmethod) according to pdf it seems to just be the string rep of a file name for program data. Although old doc also says this should be null
    z_subfile_inc: float = 1.0 # (fzinc)
    num_w_planes: float = 0 # (fwplanes)
    w_plane_inc: float = 0.0 # (fwinc)
    w_units: SPCXType = SPCXType.SPCXArb # (freserv)
    generate_log: bool = False

    def generate_header(self) -> bytes:
        Bfile_type = self.file_type.to_bytes(1, byteorder="little")
        Bfile_version = self.file_version.to_bytes(1, byteorder = "little")
        Bexperiment_type = self.experiment_type.to_bytes(1, byteorder = "little")
        Bexponent = self.exponent.to_bytes(1, byteorder = "little", signed=True)
        if not (self.file_type & SPCFileType.TXYXYS):
            log.debug(f"header is even spaced or not XYXYXY, setting num points to count {self.num_points}")
            Bnum_points = self.num_points.to_bytes(4, byteorder="little")
        else:
            log.debug(f"uneven, x will be specified, setting points count to 0")
            Bnum_points = b"\x00\x00\x00\x00"
        Bfirst_x = pack("d", self.first_x)
        Blast_x = pack("d", self.last_x)
        Bnum_subfiles = pack("l", self.num_subfiles)
        Bx_units = self.x_units.to_bytes(1, "little")
        By_units = self.y_units.to_bytes(1, "little")
        Bz_units = self.z_units.to_bytes(1, "little")
        Bpost_disposition = self.post_disposition.to_bytes(1, "little")
        # compress date is already formated to bytes by the object
        Bres_desc = bytearray(self.res_desc, encoding="utf-8")
        Bres_desc = fit_byte_block(Bres_desc, RES_DESC_LIMIT)
        Bsrc_instrument_desc = bytearray(self.src_instrument_desc, encoding="utf-8")
        Bsrc_instrument_desc = fit_byte_block(Bsrc_instrument_desc, SRC_INSTRUMENT_LIMIT)
        Bpeak_point = pack("e", self.peak_point)
        spare = bytes(bytearray(b"\x00\x00\x00\x00"*SPARE_LIMIT))
        Bmemo = bytearray(self.memo, encoding="utf-8")
        Bmemo = fit_byte_block(Bmemo, MEMO_LIMIT)
        Bcustom_axes = b"\x00".join([bytes(ax, encoding="utf-8") for ax in self.custom_axes.default_factory()])
        Bcustom_axes = fit_byte_block(bytearray(Bcustom_axes), AXES_LIMIT)
        log_offset = self.calc_log_offset(self.file_type, self.num_subfiles, self.x_values, self.y_values) # (flogoff)
        if self.file_type & SPCFileType.TMULTI and self.file_type & SPCFileType.TXVALS:
            # dir offset is log offset without dir since dir comes before log
            # per the spec num_points is the dir offset when there is a dir
            log.debug(f"setting dir offset to {log_offset}, num subfiles is {self.num_subfiles}")
            if self.file_type & SPCFileType.TXYXYS:
                Bnum_points = log_offset.to_bytes(4, byteorder="little") # only TXYXYS changes num_points to dir offset
            log_offset += self.num_subfiles * 12 # spc.h defines each dir entry as 12 bytes, one entry per subfile
        if self.generate_log:
            Blog_offset = pack("l", log_offset)
        else:
            Blog_offset = b"\x00\x00\x00\x00"
        Blog_offset = b"\x00\x00\x00\x00"
        Bspectra_mod_flag = pack("l", self.spectra_mod_flag)
        Bprocess_code = self.process_code.to_bytes(1, byteorder = "little")
        Bmethod_file = fit_byte_block(bytearray(self.method_file), METHOD_FILE_LIMIT)
        Bz_subfile_inc = pack("f", self.z_subfile_inc)
        Bnum_w_planes = pack("l", self.num_w_planes)
        Bw_plane_inc = pack("f", self.w_plane_inc)
        Bw_units = self.w_units.to_bytes(1, byteorder="little")
        Breserved = bytes(bytearray(b"\x00"*RESERVE_LIMIT))
        field_bytes = [
            Bfile_type, Bfile_version, Bexperiment_type, Bexponent, Bnum_points,
            Bfirst_x, Blast_x, Bnum_subfiles, Bx_units, By_units, Bz_units,
            Bpost_disposition, self.compress_date.get(), Bres_desc, Bsrc_instrument_desc, 
            Bpeak_point, spare, Bmemo, Bcustom_axes, Blog_offset, Bspectra_mod_flag,
            Bprocess_code, self.calib_plus_one, self.sample_inject, self.data_mul, Bmethod_file,
            Bz_subfile_inc, Bnum_w_planes, Bw_plane_inc, Bw_units, Breserved
            ]
        file_header = b"".join(field_bytes)

        if len(file_header) != 512:
            log.critical(f"file_header length is not 512 length was {len(file_header)}") # this shouldn't happen
            raise RuntimeError("Header incorrect length")
        return file_header

    def calc_log_offset(self, 
                        file_type: SPCFileType,
                        num_subfiles: int, 
                        x_values: np.ndarray,
                        y_values: np.ndarray,
                        ) -> int:
        """
        Calculates the log offset in bytes,
        assumes no directory present. If there is a directory
        that value is added later in the generate_header function
        """
        log_offset = 0
        log_offset += 512 # length of header
        log_offset += num_subfiles * 32 # add 32 bytes for each subheader
        if file_type & SPCFileType.TXVALS and not (file_type & SPCFileType.TXYXYS):
            log_offset += (len(x_values) * 4) # number of points for common x axis
        if file_type & SPCFileType.TMULTI and file_type & SPCFileType.TXYXYS:
            # multiply by 4 due to 32bit float
            for idx in range(len(x_values)):
                log_offset += (len(x_values[idx]) * 4) # chose clarity over conciseness here
                log_offset += (len(y_values[idx]) * 4) # first line represnets each x axis in subfile, second is each y
        elif file_type & SPCFileType.TMULTI:
            for idx in range(len(y_values)):
                log_offset += (len(y_values[idx]) * 4) * num_subfiles
        else:
            log_offset += (len(y_values) * 4) * num_subfiles # only a y axis in each sub file
        return int(log_offset) # should be int but just to be safe