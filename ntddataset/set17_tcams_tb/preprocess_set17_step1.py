import os
import string
import pandas as pd

#엑셀 파일(두번째 시트만 사용-숨어있는 보호된 첫번째 시트가 존재)
datafile = 'GSK_TC_TBset-final.xlsx'

#'TOX50 HepG2 MOD' 추가했습니다. Tox50 HepG2 값에 대한 'relation'에 해당되서요.
#엑셀 데이터에 ACT in interacellular pXC50, ACT in resistant strains MIC (uM) 값이 있긴 한데, 아무리 찾아도 ACT가 무슨 뜻인지 찾을 수가 없네요. 그래서 데이터에서 제외되는게 맞는 것 같습니다.
selected_columns = ['DATABase number', 'Parent SMILES', 'H37Rv MIC  (uM)', 'BCG pIC50', 'TOX50 HepG2 MOD', 'TOX50 HepG2 (uM)', 'Salt Name']
df = pd.read_excel(datafile,sheet_name=1,header=0,usecols=selected_columns)

list_to_remove = []

def process_salt_names(df):
    df['Salt Name'].fillna('.', inplace=True)
    salt_to_smiles = {
        "Hydrochloride": "Cl",
        "Sodium salt": "Na",
        "Trifluoroacetic acid salt": "OC(=O)C(F)(F)F"
        }

    def convert_to_smiles(salt_name):
        # Check if the salt name starts with "1" or "2"
        if salt_name.startswith("1 "):
            count = 1
            salt_name = salt_name[2:]
        elif salt_name.startswith("2 "):
            count = 2
            salt_name = salt_name[2:]
        else:
            count = 1  # default count

        smiles_str = salt_to_smiles.get(salt_name, salt_name)
        return '.'.join([smiles_str] * count)

    df['salt'] = df['Salt Name'].apply(convert_to_smiles)

    return df

df = process_salt_names(df)

df = df.drop(list_to_remove, axis=0)

df['smiles code'] = df['Parent SMILES'] + "." + df['salt']
df['smiles code'] = df['smiles code'].apply(lambda x: x.rstrip('.'))
df['active ingredient'] = df['Parent SMILES']
df['salt'] = df['salt']
#BCG pIC50는 아무리 논문도 찾아보고 했는데 unit 언급이 없네요.
df['relation'] = '='
df['assay value'] = df['BCG pIC50']
df['target disease'] = 'Tuberculosis'
df['protein target'] = ''
df['treatment concentration'] = ''
df['treatment unit'] = ''
df['exposure time (h)'] = ''
df['measurement'] = 'pIC50'
df['type'] = ''
df['model type'] = ''
#description에는 컬럼 이름을 추가하는게 제일 좋을 것 같아요!
df['description'] = 'BCG pIC50'
df['dataset'] = 'set17'
df['CMPD_ID'] = df['DATABase number']

df_Assay2 = df.copy()
df_Assay2['assay value'] = df['H37Rv MIC  (uM)']
df_Assay2['treatment unit'] = 'uM'
df_Assay2['measurement'] = 'MIC'
df_Assay2['description'] = 'H37Rv MIC (uM)'

df_Assay3 = df_Assay2.copy()
df_Assay3['assay value'] = df['TOX50 HepG2 (uM)']
#'TOX50 HepG2 MOD' 추가했습니다. Tox50 HepG2 값에 대한 'relation'에 해당되서요.
df_Assay3['relation'] = df['TOX50 HepG2 MOD']
df_Assay3['measurement'] = 'TOX50'
df_Assay3['description'] = 'TOX50 HepG2 (uM)'
df_Assay3['target disease'] = 'Human toxicity'

df_final = pd.concat([df, df_Assay2, df_Assay3], ignore_index=True)

final_columns = ['smiles code', 'active ingredient', 'salt', 'relation', 'assay value', 'target disease',
                 'protein target', 'treatment concentration', 'treatment unit', 'measurement', 'type',
                 'model type', 'description', 'dataset', 'CMPD_ID']

df_final = df_final[final_columns]
df_final.to_csv('set17_result.csv', index=False)