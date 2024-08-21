import pandas as pd
import numpy as np
import cobra
import re
import csv
from datetime import datetime
from scipy.stats import qmc
print(cobra.show_versions())

import swiglpk
print(swiglpk.glp_version())
from swiglpk import glp_adv_basis

import logging
FORMAT = '%(asctime)s :: %(levelname)s :: %(message)s'
logging.basicConfig(format=FORMAT,level=logging.INFO)
log = logging.getLogger("calc-FBA-Sobol-logger")
log.info("Hello, world")

import os
chunk = int(os.environ['SLURM_ARRAY_TASK_ID'])
chunkSize = int(os.environ['chunkSize'])
log.info(f"chunk={chunk}, chunkSize={chunkSize}")

odir = os.environ['odir']
fname = os.environ['model']
file_name, file_extension = os.path.splitext(fname)
#objname = os.environ['objective']
log.info(f'model:{fname}')

seed = int(os.environ['seed'])
log.info(f'seed:{seed}')

if file_extension == '.xml':
    pcModel=cobra.io.read_sbml_model(fname)
elif file_extension == '.mat':
    pcModel=cobra.io.load_matlab_model(fname)
elif file_extension == '.yml':
    pcModel=cobra.io.load_yaml_model(fname)
else :
    log.error(f'model:{fname} has unknown extension.\n')
    exit

#pcModel.objective=objname

#lreR = [r for r in pcModel.reactions if r.boundary]
ex_re=re.compile('EX_.+')
bm_re=re.compile('.*(B|b)iomass.*')
ureR = [r for r in pcModel.reactions if ex_re.search(r.id) is not None and bm_re.search(r.id) is None and r.objective_coefficient==0]
if len(ureR)==0:
    ureR = [r for r in pcModel.exchanges if bm_re.search(r.id) is None]
rnames=[r.id for r in ureR]
lreR = [r for r in ureR if r.reversibility]
flnames=[(r.id+'_fl') for r in pcModel.reactions]+[(r.id+'_rc') for r in pcModel.reactions]
lre = [(r.id+'_l') for r in lreR]
lreC = [-1000 for r in lreR]

ure = [(r.id+'_u') for r in ureR]
ureC = [1000 for r in ureR]
multiplier = lreC + ureC
columns=lre+ure
solcols=['']+columns+['Solution']
#ofs = open(odir+'/fba_sobol_solution_'+str(chunk)+'_'+str(chunkSize)+'.csv', mode='a')
ofs = odir+'/fba_sobol_solution_'+str(chunk)+'_'+str(chunkSize)+'.csv'
if os.path.exists(ofs):
  prev_table = pd.read_csv(ofs)
  optFlag = False
else:
  prev_table = pd.DataFrame([],columns=solcols)
  optFlag = True
  with open(ofs, mode='a') as file:
    csvwriter = csv.writer(file)
    csvwriter.writerow(solcols)

#off = open(odir+'/fba_sobol_fluxes_'+str(chunk)+'_'+str(chunkSize)+'.csv', mode='a')
off = odir+'/fba_sobol_fluxes_'+str(chunk)+'_'+str(chunkSize)+'.csv'
if not os.path.exists(off):
  with open(off, mode='a') as file:
    csvwriter = csv.writer(file)
    csvwriter.writerow(solcols+flnames)
sampler = qmc.Sobol(d=len(columns),seed=seed)
skip=(chunk-1)*chunkSize
sampler = sampler.fast_forward(skip)
dt = sampler.random(chunkSize)

solutions = []
fluxes = []
reduced_costs = []
with pcModel:
  for i in range(chunkSize):
    if not optFlag:
      if i < prev_table.shape[0]:
        log.info(f'i={i}, solution skipped')
      elif i == prev_table.shape[0]:
        optFlag=True
        prms=[i]+[dt[i,j] for j in range(dt.shape[1])]+[float('inf')]
        fls=[float('inf') for r in pcModel.reactions]
        rcs=[float('inf') for r in pcModel.reactions]
        output_table = pd.DataFrame([prms])#,columns=[columns)
        output_table.to_csv(ofs, mode='a', index=False, header=False)
        output_table1 = pd.DataFrame([prms+fls+rcs])#,columns=[columns)
        output_table1.to_csv(off, mode='a', index=False, header=False)
        log.info(f'i={i}, solution infeasible')
    else:
      j = 0
      for lr in lreR:
        r = pcModel.reactions.get_by_id(lr.id)
        r.lower_bound=multiplier[j]*dt[i,j]
        j = j+1
      for ur in ureR:
        r = pcModel.reactions.get_by_id(ur.id)
        r.upper_bound=max(r.lower_bound,multiplier[j]*dt[i,j])
        j = j+1
      glp_adv_basis(pcModel.solver.problem, 0)
      log.info(f'Dataset {i}, Optimisation starts.')
      s = pcModel.optimize()
      log.info(f'Dataset {i}, Optimisation done.')
      solutions.append(s.objective_value)
      prms=[i]+[dt[i,j] for j in range(dt.shape[1])]+[s.objective_value]
      output_table = pd.DataFrame([prms])#,columns=[columns)
      output_table.to_csv(ofs, mode='a', index=False, header=False)
      fls=[pcModel.reactions.get_by_id(r.id).flux for r in pcModel.reactions]
      fluxes.append(fls)
      rcs=[pcModel.reactions.get_by_id(r.id).reduced_cost for r in pcModel.reactions]
      output_table1 = pd.DataFrame([prms+fls+rcs])#,columns=[columns)
      output_table1.to_csv(off, mode='a', index=False, header=False)
      reduced_costs.append(rcs)
      now = datetime.now()
      log.info(f'i={i}, solution={s}')

log.info('Calculation finished.')
#output_table = pd.DataFrame(dt,columns=columns)
#output_table['Solution'] = solutions
#dt=pd.DataFrame(fluxes,columns=flnames)
#dt1=pd.DataFrame(reduced_costs,columns=[(n+'_rc') for n in rnames])
#output_table.to_csv(of)
#ot=output_table.merge(dt,left_index=True, right_index=True)
#ot=ot.merge(dt1,left_index=True, right_index=True)
#ot.to_csv(of1)

#log.info(f'Calculation saved to {of}.')
