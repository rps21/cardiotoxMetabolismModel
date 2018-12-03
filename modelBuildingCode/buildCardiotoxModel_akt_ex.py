import pickle
import tempfile
from indra.tools.small_model_tools.modelScope.buildPriorModel import build_prior
from indra.statements import *
from indra.tools import assemble_corpus as ac
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildExperimentalStatements
from indra.tools.small_model_tools.modelScope.buildDrugTreatmentStatements import buildDrugTargetStmts
from indra.tools.small_model_tools.modelScope.buildUseCaseModel import expandModel
from indra.tools.small_model_tools.modelScope.addBNGLParameters import addObservables
from indra.tools.small_model_tools.modelScope.addBNGLParameters import addSimParamters
from indra.tools.small_model_tools.modelContext import extraModelReductionTools as ex
from indra.assemblers import PysbAssembler
import pysb


#Build Prior


model_types = [Phosphorylation,Dephosphorylation,ActiveForm,IncreaseAmount,DecreaseAmount,Complex] #Gef,GtpActivation]#Check sos and raf stmts 
drugTargets = ['PDGFRA']#,'KDR','FLT3']
modifiedNodes = ['RPS6','PKM','HIF1A']# ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC']
#modifiedNodes = ['PKM']# ['JUN','STAT1','PKM','RPS6','AURKA','HIF1A','MYC']
otherNodes =[]# ['MYC','JUN']#ligands = ['PDGFA','PDGF','VEGF','VEGFA','FLT3LG']
model_genes = drugTargets + modifiedNodes + otherNodes


reach_stmts = './indraReading/raw_stmts/reach_output.pkl'
trips_stmts = './indraReading/raw_stmts/trips_output.pkl'
prior_stmts = build_prior(model_genes,model_types,dbs=True,additional_stmts_files=[])

#Prior reduction

#Map to model_types
mods_to_keep = list(map(lambda obj: obj.__name__.lower(), model_types))

stmts_to_remove = []
for st in prior_stmts:
    full_mods = list(map(lambda obj: obj.mods, st.agent_list())) #gives mods for a stmt, in list of lists. Actually need mod types. [(phosphorylation, Y, 589)]
    flattened_mods = list(itertools.chain.from_iterable(full_mods))
    #mod_types = list(map(lambda obj: obj.mod_type, mods) #mod types in lol 
    mod_types = list(map(lambda obj: obj.mod_type, flattened_mods))
    if any([el for el in mod_types if el not in mods_to_keep]):
        stmts_to_remove.append(st)
prior_model_stmts = [st for st in prior_stmts if st not in stmts_to_remove]


stmts_to_remove = []
for st in prior_model_stmts:
    for ag in st.agent_list():
        for mod in ag.mods:
            if mod.mod_type not in mods_to_keep:
                stmts_to_remove.append(st)

prior_model_stmts = [st for st in prior_model_stmts if st not in stmts_to_remove]


prior_model_stmts = ex.removeMutations(prior_model_stmts)
prior_model_stmts = ex.removeDimers(prior_model_stmts)


#Add drug-target interactions 
drugSentences = 'SORAFENIB dephosphorylates PDGFRA. SORAFENIB dephosphorylates KDR. SORAFENIB dephosphorylates FLT3'
drug_stmts = buildDrugTargetStmts(drugSentences)
prior_model_stmts = prior_stmts + drug_stmts

#save prior model 
ac.dump_statements(prior_model_stmts,'../final_testing_models/cardiotoxPrior_akt.pkl')



###########################################################################

#Expand model to match experimental observations

prior_model_stmts = ac.load_statements('../final_testing_models/cardiotoxPrior_akt.pkl')
#otherNodes = ['MYC','JUN']#ligands = ['PDGFA','PDGF','VEGF','VEGFA','FLT3LG']

expSentences = 'SORAFENIB dephosphorylates RPS6. SORAFENIB phosphorylates PKM. SORAFENIB transcribes HIF1A.'
exp_stmts = buildExperimentalStatements(expSentences)


expObservations = {'RPS6':['phosphorylation',exp_stmts[0]],'PKM':['dephosphorylation',exp_stmts[1]],'HIF1A':['decreaseamount',exp_stmts[2]]}
#expObservations = {'PKM':['phosphorylation',exp_stmts[1]]}
#expObservations = {'PKM':['dephosphorylation',exp_stmts[1]]}
drug = 'SORAFENIB'
#drugTargets = ['FLT3','PDGFRA','KDR']
drugTargets = ['PDGFRA']
initialStmts = prior_model_stmts
initialNodes = []
for st in prior_model_stmts:
    for ag in st.agent_list():
        if ag.name not in initialNodes:
            initialNodes.append(ag.name)


modelStmts = expandModel(expObservations,drug,drugTargets,initialStmts,initialNodes,otherNodes)


###########################################################################

#Build and write bngl model 

pa = PysbAssembler()
pa.add_statements(modelStmts)
pysbModel = pa.make_model()
#WARNING: [2018-08-29 16:12:59] indra/pysb_assembler - HIF1A transcribes itself, skipping???

model_obs = addObservables(pysbModel,bound=True)


bngl_model_filter = pysb.export.export(model_obs,'bngl')
bngl_file_filter = open('../final_testing_models/cardiotoxPrior_akt.bngl','w')
bngl_file_filter.write(bngl_model_filter+'\n')
bngl_file_filter.close()
model_actions = addSimParamters(method='ode',equil=True,equilSpecies=['SORAFENIB()'],viz=True)
bngl_file_filter = open('../final_testing_models/cardiotoxPrior_akt.bngl','w')
bngl_file_filter.write(model_actions)
bngl_file_filter.close()



#from indra.util import _require_python3
#import os
#import json
#import time
#import pickle

## CREATE A JSON FILE WITH THIS INFORMATION, E.G., a file consisting of:
## {"basename": "fallahi_eval", "basedir": "output"}
#with open('config.json', 'rt') as f:
#    config = json.load(f)
## This is the base name used for all files created/saved
#basen = config['basename']
## This is the base folder to read/write (potentially large) files from/to
## MODIFY ACCORDING TO YOUR OWN SETUP
#based = config['basedir']

## This makes it easier to make standardized pickle file paths
#prefixed_pkl = lambda suffix: os.path.join(based, basen + '_' + suffix + '.pkl')

#def pkldump(suffix, content):
#    fname = prefixed_pkl(suffix)
#    with open(fname, 'wb') as fh:
#        pickle.dump(content, fh)

#def pklload(suffix):
#    fname = prefixed_pkl(suffix)
#    print('Loading %s' % fname)
#    ts = time.time()
#    with open(fname, 'rb') as fh:
#        content = pickle.load(fh)
#    te = time.time()
#    print('Loaded %s in %.1f seconds' % (fname, te-ts))
#    return content




