def fit_byte_block(field: bytearray, limit: int) -> bytes:
    while len(field) < limit:
        field.extend(bytearray(b"\x00"))
    if len(field) > limit:
        field = field[:limit]
        field[-1] = 0 # most fields are null terminated so just make that the default
    return bytes(field)

RES_DESC_LIMIT = 9
SRC_INSTRUMENT_LIMIT = 9
MEMO_LIMIT = 130
AXES_LIMIT = 30
METHOD_FILE_LIMIT = 48
SPARE_LIMIT = 8
RESERVE_LIMIT = 187
LOG_RESERVE_LIMIT = 44