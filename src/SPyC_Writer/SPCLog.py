from struct import pack
from dataclasses import dataclass

from .common import LOG_RESERVE_LIMIT

@dataclass
class SPCLog:
    log_data: bytes
    log_text: str

    def generate_log_header(self) -> bytes:
        text_len = len(self.log_text.encode())
        data_len = len(self.log_data)
        block_size = 64 + text_len + data_len
        mem_block = 4096 * round(block_size/4096)
        text_offset = 64 + data_len
        Bblock_size = pack("l", block_size)
        Bmem_size = pack("l", mem_block)
        Btext_offset = pack("l", text_offset)
        Bdata_len = pack("l", data_len)
        Bdisk_len = b"\x00\x00\x00\x00"
        Breserved = bytes(bytearray(b"\x00"*LOG_RESERVE_LIMIT))
        return b"".join([Bblock_size, Bmem_size, Btext_offset, Bdata_len, Bdisk_len, Breserved])