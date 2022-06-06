import datetime

class SPCDate:
    def __init__(self, time: datetime = None) -> None:
        if time is None:
            time = datetime.datetime.now()

        minutes = int(time.minute) 
        hour = int(time.hour) << 6
        day = int(time.day) << 11
        month = int(time.month) << 16
        year = int(time.year) << 20
        formatted_date = (minutes | hour | day | month | year)
        
        self.compressed_date = formatted_date.to_bytes(4, byteorder = "little")

    def get(self):
        return self.compressed_date