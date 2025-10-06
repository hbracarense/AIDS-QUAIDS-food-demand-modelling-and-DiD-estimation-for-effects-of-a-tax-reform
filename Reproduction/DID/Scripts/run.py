
import os
from main import main
from get_env import path_class, files_class



main(base_path = os.path.join(path_class.path_abs, path_class.path_data, files_class.file_main),
     elas_path = os.path.join(path_class.path_abs, path_class.path_data, files_class.file_elas),
     pre_path = os.path.join(path_class.path_abs, path_class.path_data, files_class.file_pre),
     s1_path  = os.path.join(path_class.path_abs, path_class.path_data, files_class.file_1),
     s2_path  = os.path.join(path_class.path_abs, path_class.path_data, files_class.file_2),
     outdir  = os.path.join(path_class.path_abs, path_class.path_exports),
     boot_B=2000,
     placebo_runs=2000,
     seed=42)