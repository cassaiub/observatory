# phase2_models.py

class TargetGroup:
    """Holds the state and filepaths for a specific object/filter/hardware combo."""
    def __init__(self, group_key, meta):
        self.group_key = group_key
        self.meta = meta  
        self.raw_files = []
        
        # Pipeline State
        self.anchor_filepath = None
        self.master_filepath = None
        self.successful_frames = 0
        self.wcs_solved = False