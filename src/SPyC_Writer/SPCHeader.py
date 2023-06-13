import logging
from struct import pack
from dataclasses import dataclass, field
import numpy as np
from .SPCDate import SPCDate
from .common import fit_byte_block, RES_DESC_LIMIT, SPARE_LIMIT, MEMO_LIMIT, AXES_LIMIT, SRC_INSTRUMENT_LIMIT, RESERVE_LIMIT, METHOD_FILE_LIMIT
from .SPCEnums import SPCFileType, SPCModFlags, SPCPostDisposition, SPCProcessCode, SPCTechType, SPCXType, SPCYType
log = logging.getLogger(__name__)

class SPCHeader:
    """

    Corresponds to the spc file header.

    minimum fields are file type, num points and compression date.

    Source code comments detail corresponding spc.h name.

    """
    x_values = np.empty(shape=0)
    y_values = np.empty(shape=0)
    file_version = 75
    experiment_type = SPCTechType.SPCTechGen
    exponent = -128
    first_x = 0
    last_x = 0
    num_subfiles = 0
    x_units = SPCXType.SPCXArb
    y_units = SPCYType.SPCYArb
    z_units = SPCXType.SPCXArb
    post_disposition = SPCPostDisposition.PSTDEFT
    res_desc = ''
    src_instrument_desc = ''
    peak_point = 0
    memo = ''
    custom_axes = field(default_factory=list)
    spectra_mod_flag = SPCModFlags.UNMOD
    process_code = SPCProcessCode.PPCOMP
    calib_plus_one = b'\x00'
    sample_inject = b'\x00\x00'
    data_mul = b'\x00\x00\x00\x00'
    method_file = b'\x00'
    z_subfile_inc = 1.0
    num_w_planes = 0
    w_plane_inc = 0.0
    w_units = SPCXType.SPCXArb
    generate_log = False

    def generate_header(self):
        Bfile_type = self.file_type.to_bytes(1, byteorder='little')
        Bfile_version = self.file_version.to_bytes(1, byteorder='little')
        Bexperiment_type = self.experiment_type.to_bytes(1, byteorder='little')
        Bexponent = self.exponent.to_bytes(1, byteorder='little', signed=True)
        if not self.file_type & SPCFileType.TXYXYS:
            log.debug(
                f'header is even spaced or not XYXYXY, setting num points to count {self.num_points}'
                )
            Bnum_points = self.num_points.to_bytes(4, byteorder='little')
        else:
            log.debug(f'uneven, x will be specified, setting points count to 0'
                )
            Bnum_points = b'\x00\x00\x00\x00'
        Bfirst_x = pack('d', self.first_x)
        Blast_x = pack('d', self.last_x)
        Bnum_subfiles = pack('l', self.num_subfiles)
        Bx_units = self.x_units.to_bytes(1, 'little')
        By_units = self.y_units.to_bytes(1, 'little')
        Bz_units = self.z_units.to_bytes(1, 'little')
        Bpost_disposition = self.post_disposition.to_bytes(1, 'little')
        Bres_desc = bytearray(self.res_desc, encoding='utf-8')
        Bres_desc = fit_byte_block(Bres_desc, RES_DESC_LIMIT)
        Bsrc_instrument_desc = bytearray(self.src_instrument_desc, encoding
            ='utf-8')
        Bsrc_instrument_desc = fit_byte_block(Bsrc_instrument_desc,
            SRC_INSTRUMENT_LIMIT)
        Bpeak_point = pack('e', self.peak_point)
        spare = bytes(bytearray(b'\x00\x00\x00\x00' * SPARE_LIMIT))
        Bmemo = bytearray(self.memo, encoding='utf-8')
        Bmemo = fit_byte_block(Bmemo, MEMO_LIMIT)
        Bcustom_axes = b'\x00'.join([bytes(ax, encoding='utf-8') for ax in
            self.custom_axes.default_factory()])
        Bcustom_axes = fit_byte_block(bytearray(Bcustom_axes), AXES_LIMIT)
        log_offset = self.calc_log_offset(self.file_type, self.num_subfiles,
            self.x_values, self.y_values)
        if (self.file_type & SPCFileType.TMULTI and self.file_type &
            SPCFileType.TXVALS):
            log.debug(
                f'setting dir offset to {log_offset}, num subfiles is {self.num_subfiles}'
                )
            if self.file_type & SPCFileType.TXYXYS:
                Bnum_points = log_offset.to_bytes(4, byteorder='little')
            log_offset += self.num_subfiles * 12
        if self.generate_log:
            Blog_offset = pack('l', log_offset)
        else:
            Blog_offset = b'\x00\x00\x00\x00'
        Blog_offset = b'\x00\x00\x00\x00'
        Bspectra_mod_flag = pack('l', self.spectra_mod_flag)
        Bprocess_code = self.process_code.to_bytes(1, byteorder='little')
        Bmethod_file = fit_byte_block(bytearray(self.method_file),
            METHOD_FILE_LIMIT)
        Bz_subfile_inc = pack('f', self.z_subfile_inc)
        Bnum_w_planes = pack('l', self.num_w_planes)
        Bw_plane_inc = pack('f', self.w_plane_inc)
        Bw_units = self.w_units.to_bytes(1, byteorder='little')
        Breserved = bytes(bytearray(b'\x00' * RESERVE_LIMIT))
        field_bytes = [Bfile_type, Bfile_version, Bexperiment_type,
            Bexponent, Bnum_points, Bfirst_x, Blast_x, Bnum_subfiles,
            Bx_units, By_units, Bz_units, Bpost_disposition, self.
            compress_date.get(), Bres_desc, Bsrc_instrument_desc,
            Bpeak_point, spare, Bmemo, Bcustom_axes, Blog_offset,
            Bspectra_mod_flag, Bprocess_code, self.calib_plus_one, self.
            sample_inject, self.data_mul, Bmethod_file, Bz_subfile_inc,
            Bnum_w_planes, Bw_plane_inc, Bw_units, Breserved]
        file_header = b''.join(field_bytes)
        if len(file_header) != 512:
            log.critical(
                f'file_header length is not 512 length was {len(file_header)}')
            raise RuntimeError('Header incorrect length')
        return file_header

    def calc_log_offset(self, file_type, num_subfiles, x_values, y_values):
        """

        Calculates the log offset in bytes,

        assumes no directory present. If there is a directory

        that value is added later in the generate_header function

        """
        log_offset = 0
        log_offset += 512
        log_offset += num_subfiles * 32
        if (file_type & SPCFileType.TXVALS and not file_type & SPCFileType.
            TXYXYS):
            log_offset += len(x_values) * 4
        if file_type & SPCFileType.TMULTI and file_type & SPCFileType.TXYXYS:
            for idx in range(len(x_values)):
                log_offset += len(x_values[idx]) * 4
                log_offset += len(y_values[idx]) * 4
        elif file_type & SPCFileType.TMULTI:
            for idx in range(len(y_values)):
                log_offset += len(y_values[idx]) * 4 * num_subfiles
        else:
            log_offset += len(y_values) * 4 * num_subfiles
        return int(log_offset)
