from datetime import datetime
from typing import Optional

class SPCDate:
    def __init__(self, time: Optional[datetime] = None) -> None:
        if time is not None: # based on s_evenx.spc date can be 0 which python doesn't support for year
            minutes = int(time.minute) 
            hour = int(time.hour) << 6
            day = int(time.day) << 11
            month = int(time.month) << 16
            year = int(time.year) << 20
        else:
            minutes = 0
            hour = 0
            day = 0
            month = 0
            year = 0
        formatted_date = (minutes | hour | day | month | year)
        
        self.compressed_date = formatted_date.to_bytes(4, byteorder = "little")

    def get(self):
        return self.compressed_date