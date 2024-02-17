
class MH_Exception(Exception):
    pass


class FileNotFoundMH(MH_Exception):
    print('Invalid directory!')
    pass


class InvalidCardId(MH_Exception):
    pass