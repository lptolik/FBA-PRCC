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
flnames=[(n+'_fl') for n in rnames]+[(n+'_rc') for n in rnames]
lre = [(r.id+'_l') for r in lreR]
lreC = [-1000 for r in lreR]

ure = [(r.id+'_u') for r in ureR]
ureC = [1000 for r in ureR]
multiplier = lreC + ureC
columns=lre+ure
solcols=['']+columns+['Solution']
#ofs = open(odir+'/fba_sobol_solution_'+str(chunk)+'_'+str(chunkSize)+'.csv', mode='a')
ofs = odir+'/fba_sobol_isactive_'+str(chunk)+'_'+str(chunkSize)+'.csv'
ifs = odir+'/fba_sobol_solution_'+str(chunk)+'_'+str(chunkSize)+'.csv.gz'
if os.path.exists(ifs):
  dt = pd.read_csv(ifs)
  optFlag = True
  with open(ofs, mode='w') as file:
    csvwriter = csv.writer(file)
    csvwriter.writerow(solcols)
else:
  log.info(f'file:{ifs} does not exist.')
  exit

solutions = []
reduced_costs = []
isactive = []
with pcModel:
    for i in range(chunkSize):
        s0=dt.loc[i,'Solution']
        if(float('inf') == s0):
          isactive = [0 for c in columns]
          s = s0
        else:
          j = 0
          for lr in lreR:
            r = pcModel.reactions.get_by_id(lr.id)
            r.lower_bound=multiplier[j]*dt.loc[i,(r.id+'_l')]
            j = j+1
          for ur in ureR:
            r = pcModel.reactions.get_by_id(ur.id)
            r.upper_bound=max(r.lower_bound,multiplier[j]*dt.loc[i,(r.id+'_u')])
            j = j+1
          glp_adv_basis(pcModel.solver.problem, 0)
          log.info(f'Dataset {i}, Optimisation starts.')
          s = pcModel.slim_optimize()
          if abs(dt.loc[i,'Solution']-s) > 1e-7 :
              log.info(f'Dateset {i}, Solutions do not match {s}, {s0} ({abs(s0-s)}).')
          else:
              isactive = []
              for lr in lreR:
                  with pcModel:
                      r = pcModel.reactions.get_by_id(lr.id)
                      r.lower_bound=r.lower_bound*1.1
                      s1 = pcModel.slim_optimize()
                      if (s1-s)/s > 1e-5 :
                          log.info(f'Dateset {i}, Bound {r.id}_l, active: {s1} > {s} ({(s1-s)/s}).')
                          isactive.append(1)
                      else:
                          #log.info(f'Dateset {i}, Bound {r.id}_l, not active: {s1} <= {s} ({(s1-s)/s}).')
                          isactive.append(0)
                          if abs((s1-s)/s) > 1e-5 :
                              log.info(f'Dateset {i}, Bound {r.id}_l, changes: {s1} , {s} ({(s1-s)/s}).')
              for ur in ureR:
                  with pcModel:
                      r = pcModel.reactions.get_by_id(ur.id)
                      r.upper_bound=r.upper_bound*1.1
                      s1 = pcModel.slim_optimize()
                      if (s1-s)/s > 1e-5 :
                          log.info(f'Dateset {i}, Bound {r.id}_u, active: {s1} > {s} ({(s1-s)/s}).')
                          isactive.append(1)
                      else:
                          #log.info(f'Dateset {i}, Bound {r.id}_u, not active: {s1} <= {s} ({(s1-s)/s}).')
                          isactive.append(0)
                          if abs((s1-s)/s) > 1e-5 :
                              log.info(f'Dateset {i}, Bound {r.id}_u, changes: {s1} , {s} ({(s1-s)/s}).')
        log.info(f'Dataset {i}, Optimisation done.({sum(isactive)}), {s0}, {s}')
        prms=[i]+isactive+[s]
        output_table = pd.DataFrame([prms])#,columns=[columns)
        output_table.to_csv(ofs, mode='a', index=False, header=False)
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
