import itertools
from clease.settings import Concentration
import numpy as np
from ase.visualize import view
import os
from clease.settings import CECrystal
import json
import sys

from ase.optimize import BFGS
from ase.filters import ExpCellFilter
from ase.filters import StrainFilter
import ase.io.cif
from ase.formula import Formula
import pandas as pd
from ase.io.trajectory import Trajectory

from mace.calculators import MACECalculator
from ase.db import connect
from clease.tools import update_db

from ase.db import connect
from clease.tools import update_db

def all_combinations(input_list, r):
    return [list(combo) for combo in itertools.combinations(input_list, r)]

from ase.geometry.analysis import Analysis
def find_first_min_pos_after_max(arr,noise):
    if len(arr) < 2:
        return None  # Array too short to have both a max and a min
    
    # Initialize variables
    first_max_found = False
    first_max_value = 0.0
    first_max_pos = 0
    first_min_after_value = 0.0
    min_pos_after_max = 0
    
    for i in range(len(arr)):
        if not first_max_found:
            # Find the first maximum
            if arr[i] > first_max_value:
                first_max_value = arr[i]
                first_max_pos = i
            elif arr[i] < (first_max_value - noise):
                first_max_found = True
                continue
        else:
            if arr[i] < arr[i-1]:
                first_min_after_value = arr[i]
                min_pos_after_max = i
            elif arr[i] > (first_min_after_value + noise):
                break
              
    print("First max           = (",first_max_pos,",",first_max_value,")")
    print("First min after max = (",min_pos_after_max,",",first_min_after_value,")")   
    return min_pos_after_max

def get_sigma_bond(atoms,rmax,noise):
    atomrdf=Analysis(atoms)
    #Analyze to find first minimum with larger bins
    nbins = round(rmax * 10.0)
    rdf_listpeak = atomrdf.get_rdf(rmax=rmax,nbins=nbins)
    rdf_distpeak = np.arange(0.0,rmax,rmax/nbins,dtype=float)
    current_peak = 0.0
    rdf_valpeak = rdf_listpeak[0].copy()
    minpos = find_first_min_pos_after_max(rdf_valpeak,noise)
    if minpos == 0:
        minpos = 25
    #Now do integral for distribution with smaller bins
    nbins = round(rmax * 100.0)
    minpos=minpos*10
    rdf_list = atomrdf.get_rdf(rmax=rmax,nbins=nbins)
    rdf_dist = np.arange(0.0,rmax,rmax/nbins,dtype=float)
    current_peak = 0.0
    rdf_val = rdf_list[0].copy()
    
    dr = rmax/nbins
    rtot = 0.0
    num = 1E-10
    
    for i in range(minpos):
        rtot = rtot + rdf_dist[i]*(rdf_val[i]* rdf_dist[i]**2 * dr)
        num =  num + (rdf_val[i]* rdf_dist[i]**2 * dr)
    
    ravg = rtot/num
    print(ravg)
    bondsum = 0.0
    for i in range(minpos):
        bondsum = bondsum + ((rdf_dist[i]-ravg)**2)*(rdf_val[i]* rdf_dist[i]**2 * dr)
    
    sigma_bond = np.sqrt(bondsum/num)
    print(sigma_bond)
    return ravg, sigma_bond

# print(sys.argv[1],sys.argv[2])
startnum = int(sys.argv[1])
finnum = int(sys.argv[2])
print(startnum,finnum)

# Example usage:
ncat = 4
input_list1 = ['Ge', 'Ti', 'Sn', 'Hf', 'Zr', 'Pb', 'Ce','Rh','Ru','Ir','Pt','Nb','V','Mn'] 

list_comp = all_combinations(input_list1,ncat) # Generate all set of combinations

print(list_comp)
strucname = "aPbO2"
nsample = 1

first = True

directory = "./trajectories/"

if not os.path.exists(directory):
    os.makedirs(directory)

directory_Ei = "./Ei_dicts/"

if not os.path.exists(directory_Ei):
    os.makedirs(directory_Ei)

directory_csv = "./csv_files/"

if not os.path.exists(directory_csv):
    os.makedirs(directory_csv)

for i in range(startnum,finnum):
    basis_elements=[list_comp[i],['O']]
    
    delim=''
    res = delim.join([str(ele) for ele in basis_elements[0]])
    db_name = strucname+"-"+str(ncat)+"comp_"+str(res)+".db"
    print(db_name)
    
    lst = ([1] * ncat) + [2.0]
    A_eq = [lst]
    b_eq = [1.0]
    
    concrange=[ncat*[(0,(1.0/ncat))],[(1.0,1.0)]]
    #print(concrange)
    
    conc = Concentration(basis_elements=basis_elements, A_eq=A_eq, b_eq=b_eq)
    conc.set_conc_ranges(ranges=concrange)

    #alpha-PbO2 structure
    settings = CECrystal(concentration=conc,
        spacegroup=60,
        basis=[(0.00000, 0.67652, 0.25000), (0.22868, 0.61936, -0.08125)],
        cell=[4.71, 5.61, 5.06, 90, 90, 90],
        size=[5,4,4],
        db_name=db_name,
        max_cluster_dia=(3.0, 3.0, 3.0))
    
    from clease.structgen import NewStructures
    ns = NewStructures(settings, generation_number=0, struct_per_gen=nsample)
    ns.generate_random_structures()
    
    db = connect(db_name)
    for j in range(nsample):   
        istring = str(i).zfill(4)
        jstring = str(j).zfill(4)
        for row in db.select(id=(j+2)):
            atoms = row.toatoms()
        
        sample_description = istring+"-"+jstring+str(res)+strucname+str(ncat)
        outfilestr = directory+sample_description+"comp.traj"
        reltraj = Trajectory(outfilestr,'w',atoms)

        print(f"Iniating calculator for {istring} {jstring}")
        calculator = MACECalculator( \
          model_paths='../../universal_mace_repo/mace-medium-density-agnesi-stress.model', \
          device='cuda')
        
        maxstep=200
        
        print("Attaching calculator")
        atoms.calc = calculator
        print(Formula(atoms.get_chemical_formula()).count())
        print('Volume before = ',atoms.get_volume())
        print('Energy = ',atoms.get_total_energy())
        ecf = ExpCellFilter(atoms)
        dyn = BFGS(ecf)
        try:
            dyn.run(fmax=0.05,steps=maxstep)
        except:
            print(f"Skipping iteration {istring} due to exception")
            continue
        
        reltraj.write(atoms)
        print('Volume after =',atoms.get_volume())
        energy = atoms.get_potential_energy()
        print('Energy =', energy)
        steps=dyn.get_number_of_steps()
        print("steps =",steps)
        if steps >= 200:
            continue
        bond_length, sigma_bond = get_sigma_bond(atoms,5.0,0.1)
    
        #Save individual atom energies Ei as a json file
        atoms.calc.results
    
        node_energies = atoms.calc.results['node_energy']
    
        Eioutfilestr = directory_Ei + sample_description+"comp_Ei.json"
    
        with open(Eioutfilestr, "w") as file:
            json.dump(node_energies.tolist(), file)
        
        atoms.cell.cellpar()
        
        summary_dict={}
        summary_dict['Formula'] = basis_elements[0]
        summary_dict['Energy'] = energy
        summary_dict['sigma bond length'] = sigma_bond
        summary_dict['Bond length'] = bond_length
        csvfilestr = directory_csv+istring+"-"+jstring+str(res)+str(ncat)+"comp_hbond.csv"
        if first == True:
            df_summary = pd.DataFrame([summary_dict])
            df_summary.to_csv(csvfilestr, index=False)
            first=False
        else:
            df_summarynew = pd.DataFrame([summary_dict])
            df_summarynew.to_csv(csvfilestr, index=False)
            df_summary = pd.concat([df_summary,df_summarynew])
        
        df_summary.to_csv("mixed_energy_summary.csv", index=False)

