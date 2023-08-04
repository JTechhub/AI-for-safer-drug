import os
import pandas as pd

wdir = r'/directory/for/ChEMBLNTD_curation'
csv = 'set1_raw.csv'
df = pd.read_csv(os.path.join(wdir, csv))

df['relation'] = '='
df['type']='in vitro'
df['model type']='HepG2'
df['dataset']='set1'
df['BOX']='Malaria'

chemical_info = df[['SMILES','active ingredient','salt','relation','type','model type','dataset','IFI','MW_FREEBASE','ALOGP','HBA','HBD','PSA','SOURCES','BOX']]
chemical_info.columns = ['smiles code','active ingredient','salt','relation','type','model type','dataset','IFI','MW','clogp parent','hba parent','hbd parent','tpsa parent','ChEMBL_CMPD_NUMBER','BOX']

#treatment concetnration was obtained from the research article

PCT_INB_3D7 = chemical_info.copy(deep=True)
PCT_INB_3D7['assay value'] = df['PCT_IHB_3D7']
PCT_INB_3D7['target disease']='Marlaria plasmodium falciparum 3D7'
PCT_INB_3D7['treatment concentration']=2
PCT_INB_3D7['treatment unit']='uM'
PCT_INB_3D7['description']='growth inhibition'
PCT_INB_3D7['measurement']='(%) inhibition 3D7'

PCT_INHB_DD2 = chemical_info.copy(deep=True)
PCT_INHB_DD2['assay value'] = df['PCT_INHB_DD2']
PCT_INHB_DD2['target disease']='Malaria plasmodium falciparum DD2'
PCT_INHB_DD2['treatment concentration']=2
PCT_INHB_DD2['treatment unit']='uM'
PCT_INHB_DD2['description']='growth inhibition'
PCT_INHB_DD2['measurement']='(%) inhibition DD2'

PCT_INHIB_3D7_PFLDH = chemical_info.copy(deep=True)
PCT_INHIB_3D7_PFLDH['assay value']= df['PCT_INHIB_3D7_PFLDH']
PCT_INHIB_3D7_PFLDH['target disease']='Marlaria plasmodium falciparum 3D7'
PCT_INHIB_3D7_PFLDH['protein target']='PFLDH(Lactate dehydrogenase)'
PCT_INHIB_3D7_PFLDH['treatment concentration']=2
PCT_INHIB_3D7_PFLDH['treatment unit']='uM'
PCT_INHIB_3D7_PFLDH['description']='enzyme inhibition'
PCT_INHIB_3D7_PFLDH['measurement']='(%) enzyme inhibition 3D7 PFLDH'

pXC50_3D7 = chemical_info.copy(deep=True)
pXC50_3D7['assay value']=df['pXC50_3D7']
pXC50_3D7['target disease']='Marlaria plasmodium falciparum 3D7'
pXC50_3D7['description']='growth inhibition'
pXC50_3D7['measurement']='PXC50 3D7'

PCT_INHIB_HEPG2 = chemical_info.copy(deep=True)
PCT_INHIB_HEPG2['assay value']=df['PCT_INHIB_HEPG2']
PCT_INHIB_HEPG2['target disease']='Human cytotoxicity (HepG2)'
PCT_INHIB_HEPG2['treatment concentration']=2
PCT_INHIB_HEPG2['treatment unit']='uM'
PCT_INHIB_HEPG2['description']='growth inhibition'
PCT_INHIB_HEPG2['measurement']='(%) inhibition HpeG2'

df_final = pd.concat([PCT_INB_3D7, PCT_INHB_DD2, PCT_INHIB_3D7_PFLDH, pXC50_3D7, PCT_INHIB_HEPG2], ignore_index=True)
df_final = df_final[['smiles code','active ingredient','salt','relation','assay value','target disease','protein target','treatment concentration','treatment unit','measurement','type','model type','description','dataset','IFI','MW','clogp parent','hba parent','hbd parent','tpsa parent','ChEMBL_CMPD_NUMBER','BOX']]
df_final.to_csv(os.path.join(wdir,'set1.csv'), index=False)
