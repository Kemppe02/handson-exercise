import md
import os.path

# Init run
md.run_md()

# Check if run properly
assert os.path.exists('cu.traj')

assert os.path.getsize('cu.traj') > 0
