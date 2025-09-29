"""Demonstrates molecular dynamics with constant energy."""

from asap3 import Trajectory

from ase import units
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet

def calcenergy(a):
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    inst_temp = ekin / (1.5 * units.kB)
    etot = epot + ekin
    return epot, ekin, inst_temp, etot
    

def run_md():
    # Use Asap for a huge performance increase if it is installed
    use_asap = True
    
    if use_asap:
        from asap3 import EMT
        from asap3 import LennardJones
    
        size = 6
    else:
        from ase.calculators.emt import EMT
    
        size = 3
    
    # Set up a crystal
    atoms = FaceCenteredCubic(
        directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        symbol='Ar',
        size=(size, size, size),
        pbc=True,
    )
    
    # Describe the interatomic interactions with the Effective Medium Theory
    atoms.calc = LennardJones([18], [0.010323], [3.4], rCut=6.625, modified=True)
    
    # Set the momenta corresponding to T=40K
    MaxwellBoltzmannDistribution(atoms, temperature_K=40)
    
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 1 * units.fs)  # 5 fs time step.
    
    traj = Trajectory("argon.traj", "w", atoms)
    dyn.attach(traj.write, interval=10)
    
    
    def printenergy(a=atoms):  # store a reference to atoms in the definition.
        """Function to print the potential, kinetic and total energy."""
        epot, ekin, inst_temp, etot = calcenergy(a)
        print(
            f'Energy per atom: Epot ={epot:6.3f}eV  Ekin = {ekin:.3f}eV '
            f'(T={inst_temp:3.0f}K) Etot = {etot:.3f}eV'
        )
    
    
    # Now run the dynamics
    dyn.attach(printenergy, interval=100)
    printenergy()
    dyn.run(20000)

if __name__ == "__main__":
    run_md()
