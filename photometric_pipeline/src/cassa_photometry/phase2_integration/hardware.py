# phase2_hardware.py
import os
import sys
import math

try:
    import psutil
except ImportError:
    print("[ERROR] 'psutil' is required. Run: pip install psutil")
    sys.exit(1)

class HardwareManager:
    def __init__(self):
        self.total_cores = os.cpu_count() or 1
        self.sys_ram_gb = psutil.virtual_memory().total / (1024**3)
        self.allocated_cores = 1

    def allocate_resources(self):
        print("\n📡 SCANNING HPC HARDWARE...")
        print(f" [✔] CPU Detected: {self.total_cores} Cores")
        print(f" [✔] System RAM  : {self.sys_ram_gb:.1f} GB")
        
        choice = input("\nSelect CPU Resource Level (1: 25%, 2: 50%, 3: 75%, 4: 100%) [Default: 2]: ").strip()
        multiplier = { '1': 0.25, '2': 0.50, '3': 0.75, '4': 1.0 }.get(choice, 0.50)
        self.allocated_cores = max(1, int(self.total_cores * multiplier))
        
        os.environ["OMP_NUM_THREADS"] = str(self.allocated_cores)
        os.environ["OPENBLAS_NUM_THREADS"] = str(self.allocated_cores)
        os.environ["MKL_NUM_THREADS"] = str(self.allocated_cores)
        
        return self.allocated_cores

    def print_estimations(self, num_files, max_stack, total_size_mb):
        """Restored from phase2_main.py: Warns user of potential RAM limits."""
        avg_file_mb = total_size_mb / num_files if num_files > 0 else 0
        peak_ram_gb = ((avg_file_mb * 2) * max_stack * 5) / 1024 + 1.5
        
        base_time = 3.5 
        thread_eff = math.sqrt(self.allocated_cores / 4) if self.allocated_cores > 4 else (self.allocated_cores / 4)
        est_time_s = (num_files * base_time) / thread_eff
        
        print("\n=======================================================")
        print(" PIPELINE RESOURCE ESTIMATION")
        print("=======================================================")
        print(f" - Total Valid Frames   : {num_files}")
        print(f" - Max Concurrent Stack : {max_stack} frames")
        print(f" - Est. Peak RAM Req.   : ~{peak_ram_gb:.1f} GB")
        print(f" - Est. Processing Time : ~{int(est_time_s//60)}m {int(est_time_s%60)}s")
        print("=======================================================")
        
        if peak_ram_gb > self.sys_ram_gb:
            print("\n[WARNING] Estimated RAM exceeds available system memory! The process may crash or swap heavily.")