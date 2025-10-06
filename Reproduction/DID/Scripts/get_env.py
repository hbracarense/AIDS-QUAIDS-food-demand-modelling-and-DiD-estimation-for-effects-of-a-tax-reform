import os
from dotenv import load_dotenv

load_dotenv()

class path_class:
    path_abs = os.getenv('PATH_ABS')
    path_data = os.getenv('PATH_DATA')
    path_exports = os.getenv('PATH_EXPORTS')

class files_class:
    file_main = os.getenv('FILE_MAIN')
    file_elas = os.getenv('FILE_ELAS')
    file_pre = os.getenv('FILE_PRE')
    file_1 = os.getenv('FILE_1')
    file_2 = os.getenv('FILE_2')