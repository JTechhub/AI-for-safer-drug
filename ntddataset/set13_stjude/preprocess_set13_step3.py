import pandas as pd

df = pd.read_excel('set13_result.xlsx')

df_dose_resp = df[df['DataSet']=='Dose Response']

dose_resp_list = []
for index, row in df_dose_resp.iterrows():
#test_row = df.loc[0,:]
    conc = row['conc_uM'].split()
    resp = row['Growth Inhibition(%)'].split()

    for i in range(len(conc)):
        each_row = {}

        each_row['smiles code']=row['smiles code']
        each_row['active ingredient'] = row['active ingredient']
        each_row['salt'] = row['salt']
        each_row['relation'] = '='
        each_row['assay value'] = float(resp[i])
        each_row['target disease'] = row['target disease']
        each_row['protein target'] = row['protein target']
        each_row['treatment concentration'] = float(conc[i])
        each_row['treatment unit'] = 'uM'
        each_row['exposure time (h)'] = row['exposure time (h)']
        each_row['measurement'] = 'Growth Inhibition(%)'
        each_row['type'] = row['type']
        each_row['model type'] = row['model type']
        each_row['description'] = row['description']
        each_row['dataset'] = row['dataset']
        each_row['Assay Type'] = row['Assay Type']
        each_row['Assay Src Description'] = row['Assay Src Description']
        each_row['Assay Organism'] =  row['Assay Organism']
        each_row['Target Type'] = row['Target Type']
        each_row['Target Name'] = row['Target Name']
        each_row['Reference'] = row['Reference']
        each_row['CMPD_ID'] = row['CMPD_ID']
        each_row['STJUDE_ID'] = row['STJUDE_ID']
        each_row['conc_uM'] = ''
        each_row['Growth Inhibition(%)'] = ''
        each_row['DataSet'] = ''
        each_row['ID'] = ''
        
        dose_resp_list.append(each_row)

    df.loc[index, 'conc_uM'] = ''
    df.loc[index, 'Growth Inhibition(%)'] = ''

df_final = pd.DataFrame(dose_resp_list)

df_check = pd.concat([df,df_final],axis=0)
df_check['Assay Type']='in vitro'
df_check.to_csv('set13_result2.csv', index=False)