import pandas as pd

#엑셀 파일(첫번째 시트만 사용)
datafile = 'TCMDC_TB_NTD.xls'
df = pd.read_excel(datafile,sheet_name=0,header=0,usecols=[0,1,2,3,4,5,6,9])

# 'list2'에서 'TBset' 값이 비어있는 행만 선택
# 'TBset'가 있는 행은 TB set sheet에 중복됨.
df_BCG = df[df['list2'].isna()].copy()

df_BCG['smiles code'] = df['SMILES version']
df_BCG['active ingredient'] = df['SMILES parent']

def process_salt_names(df):
    df['Salt name'].fillna('.', inplace=True)
    salt_to_smiles = {
        "1 Hydrochloride": "Cl"
        # ... 추가적인 매핑은 이곳에 추가
        }
    df['salt'] = df['Salt name'].apply(lambda x: salt_to_smiles.get(x, x))
    df['salt'] = df['salt'].apply(lambda x: x if x in salt_to_smiles.values() else '.')
    return df

df_BCG = process_salt_names(df_BCG)

df_BCG['relation'] = '='
df_BCG['assay value'] = ''
df_BCG['target disease'] = 'Tuberculosis'
df_BCG['protein target'] = ''
df_BCG['treatment concentration'] = df_BCG['MIC90 BCG (µM)']
df_BCG['treatment unit'] = 'µM'
df_BCG['exposure time (h)'] = ''
df_BCG['measurement'] = 'MIC90'
df_BCG['type'] = ''
df_BCG['model type'] = ''
df_BCG['description'] = ''
df_BCG['dataset'] = 'set9'
df_BCG['Assay Type'] = ''
df_BCG['Assay Src Description'] = df_BCG['GSKnumber']
df_BCG['Assay Organism'] = ''
df_BCG['Target Type'] = 'ORGANISM'
df_BCG['Target Name'] = 'Mycobacterium bovis Bacillus Calmette-Guérin'
df_BCG['Reference'] = '.'
df_BCG['CMPD_ID'] = df_BCG['DATABase number']

#엑셀 파일(두번째 시트만 사용)
df_tb = pd.read_excel(datafile,sheet_name=1,header=0,usecols=[0,1,2,3,4,5,9,10,12,13,14,15])

df_tb = process_salt_names(df_tb)

df_tb['smiles code'] = df_tb['SMILES version']
df_tb['active ingredient'] = df_tb['SMILES parent']
df_tb['relation'] = '='
df_tb['assay value'] = ''
df_tb['target disease'] = 'Tuberculosis'
df_tb['protein target'] = ''
df_tb['treatment concentration'] = df_tb['MIC90 BCG (µM)']
df_tb['treatment unit'] = 'µM'
df_tb['exposure time (h)'] = ''
df_tb['measurement'] = 'MIC90'
df_tb['type'] = ''
df_tb['model type'] = 'BCG'
df_tb['description'] = 'growth inhibition'
df_tb['dataset'] = 'set9'
df_tb['Assay Type'] = ''
df_tb['Assay Src Description'] = df_tb['GSKnumber']
df_tb['Assay Organism'] = ''
df_tb['Target Type'] = 'ORGANISM'
df_tb['Target Name'] = 'Mycobacterium bovis Bacillus Calmette-Guérin'
df_tb['Reference'] = '.'
df_tb['CMPD_ID'] = df_tb['DATABase number']

df_tb_ATP = df_tb.copy()
df_tb_ATP['treatment concentration'] = df_tb['MIC90 H37Rv uM (ATP)']
df_tb_ATP['treatment unit'] = 'µM'
df_tb_ATP['measurement'] = 'MIC90'
df_tb_ATP['model type'] = 'H37Rv'
df_tb_ATP['description'] = 'growth inhibition'
df_tb_ATP['Target Type'] = 'ORGANISM'
df_tb_ATP['Target Name'] = 'Mycobacterium tuberculosis'

df_tb_MABA = df_tb.copy()
df_tb_MABA['treatment concentration'] = df_tb['MIC H37Rv uM (MABA)']
df_tb_MABA['treatment unit'] = 'µM'
df_tb_MABA['measurement'] = 'MIC'
df_tb_MABA['model type'] = 'H37Rv'
df_tb_MABA['description'] = 'growth inhibition'
df_tb_MABA['Target Type'] = 'ORGANISM'
df_tb_MABA['Target Name'] = 'Mycobacterium tuberculosis'

df_tb_HepG2 = df_tb.copy()
df_tb_MABA['relation'] = df_tb['HepG2 PIC50 MOD']
df_tb_HepG2['treatment concentration'] = df_tb['HepG2 PIC50']
df_tb_HepG2['measurement'] = 'PIC50'
df_tb_HepG2['model type'] = 'HepG2'
df_tb_HepG2['description'] = 'growth inhibition'
df_tb_HepG2['Target Type'] = 'ORGANISM'
df_tb_HepG2['Target Name'] = 'Mycobacterium tuberculosis'

df_tb_Clogp = df_tb.copy()
df_tb_Clogp['treatment concentration'] = df_tb['Clogp']


df_final = pd.concat([df_BCG, df_tb, df_tb_ATP, df_tb_MABA, df_tb_HepG2, df_tb_Clogp], ignore_index=True)

final_columns = ['smiles code', 'active ingredient', 'salt', 'relation', 'assay value', 'target disease',
                 'protein target', 'treatment concentration', 'treatment unit', 'measurement', 'type',
                 'model type', 'description', 'dataset', 'Assay Type', 'Assay Src Description',
                 'Assay Organism', 'Target Type', 'Target Name', 'Reference', 'CMPD_ID']

df_final = df_final[final_columns]
df_final.to_csv('set9_result.csv', index=False)