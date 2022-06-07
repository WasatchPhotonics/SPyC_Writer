import logging
from struct import pack
from dataclasses import dataclass

from .SPCEnums import SPCSubfileFlags

log = logging.getLogger(__name__)

@dataclass
class SPCSubheader:
    subfile_flags: SPCSubfileFlags = SPCSubfileFlags.SUBNONE # (subflgs)
    exponent: int = -128 # (subexp) -128 so it will be an IEEE 32bit float
    sub_index: int = 0 # (subindx)
    start_z: float = 0.0 # (subtime) looking at spc.h, z appears to be a time value so won't pass array of data points like x or y
    end_z: float = None # (subnext)
    noise_value: float = None # (subnois) should be null according to old format
    num_points: int = None # (subnpts) only needed for xyxy multifile
    num_coadded: int = None  # (subscan) should be null according to old format
    w_axis_value: float = 0.0 # (subwlevel)

    def generate_subheader(self) -> bytes:
        Bsubfile_flags = self.subfile_flags.to_bytes(1, byteorder="little")
        Bexponent = pack("<b", self.exponent)
        Bsub_index = pack("<H", self.sub_index)
        Bstart_z = pack("<f", self.start_z)
        if self.end_z is None:
            Bend_z = pack("<f", self.start_z)
        else:
            Bend_z = pack("<f", self.end_z)
        if self.noise_value is None:
            Bnoise_value = b"\x00\x00\x00\x00"
        else:
            Bnoise_value = pack("f", self.noise_value)
        if self.num_points is None:
            Bnum_points = b"\x00\x00\x00\x00"
        else:
            Bnum_points = pack("l", self.num_points)
        if self.num_coadded is None:
            Bnum_coadded = b"\x00\x00\x00\x00"
        else:
            Bnum_coadded = pack("l", self.num_coadded)
        Bw_axis_value = pack("f", self.w_axis_value)
        extra = b"\x00\x00\x00\x00" # (subresv) 4 null bytes for the reserved portion
        subheader = b''.join([Bsubfile_flags, 
                              Bexponent, 
                              Bsub_index, 
                              Bstart_z, Bend_z, 
                              Bnoise_value, 
                              Bnum_points,
                              Bnum_coadded,
                              Bw_axis_value,
                              extra
                              ])
        if len(subheader) != 32:
            log.critical(f"This shouldn't happen. Subheader length wasn't 32. Was {len(subheader)}")
            raise RuntimeError("Subheader invalid length")
        return subheader