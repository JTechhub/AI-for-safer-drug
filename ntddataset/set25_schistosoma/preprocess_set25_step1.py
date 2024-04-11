import pandas as pd
from rdkit import Chem
"""
Sheet1(REFERENCE.txt) = Reference 논문 데이터
Sheet2(COMPOUND_RECORD.txt) = 스키스토좀병 치료에 효과적일 것으로 예상되는 연구에 사용된 화합물 이름과 번호, CIDX = 화합물 인덱스
Sheet3(ASSAY.txt) = 실험에 사용된 assay 실험 정보
Sheet4(ACTIVITY.txt) = 각 화합물의 결과값(핵심 데이터)
화합물 구조는 sdf파일에서 확인할 수 있다.
"""
datafile = 'Schistoma_Gardner.xlsx'
df = pd.read_excel(datafile,sheet_name='ACTIVITY.txt')

df['smiles code'] = ''
df['active ingredient'] = ''
df['salt'] = ''
df['relation'] = df['RELATION']
df['assay value'] = df['VALUE']
df['target disease'] = 'Schistosomiasis'  #주혈흡층증(말라리아 다음으로 강한 질병)
df['protein target'] = ''
df['treatment concentration'] = ''
df['treatment unit'] = df['UNITS']
df['exposure time (h)'] = ''
df['measurement'] = df['TYPE']
df['type'] = ''
df['model type'] = ''
df['description'] = df['ColumnNames']
df['dataset'] = 'set25'
df['CMPD_ID'] = df['CIDX']

#엑셀에서 CIDX로 구조를 파악할 수 없기 때문에 sdf 파일에서 분자 정보 딕셔너리화
sdf_file = 'COMPOUND_CTAB.sdf'
#key=CIDX와 동일한 데이터, value=분자의 CTAB
def parse_sdf(sdf_file):
    with open(sdf_file, 'r') as f:
        sdf_content = f.read().split('$$$$\n')[:-1]  # 마지막에 오는 빈 항목을 제거하기 위해 슬라이스 사용
    cidx_to_data = {}
    for entry in sdf_content:
        lines = entry.strip().split('\n')
        cidx = lines[0].strip()
        cidx_to_data[cidx] = entry.strip()  # 전체 SDF 엔트리를 저장
    return cidx_to_data

cidx_to_data_full = parse_sdf(sdf_file)

# 파싱된 정보로부터 SMILES 코드 생성
extract_smiles = {}
for cidx, sdf_block in cidx_to_data_full.items():
    mol = Chem.MolFromMolBlock(sdf_block)
    if mol is not None:
        smiles = Chem.MolToSmiles(mol)
        extract_smiles[cidx] = smiles

df['smiles code'] = df['CIDX'].map(extract_smiles)

#염 분리 코드
list_to_remove = []
for index, row in df.iterrows():
    cidx = row['CIDX']
    smi = extract_smiles[cidx]
    mol = Chem.MolFromSmiles(smi)
    if '.' in smi:
        try:
            fragments_list = list(Chem.GetMolFrags(mol, asMols=True))
            fragments_sizes = [frag.GetNumAtoms() for frag in fragments_list]

            largest_i = fragments_sizes.index(max(fragments_sizes))
            largest_frag = fragments_list[largest_i]

            active = Chem.MolToSmiles(largest_frag)

            fragments_list.remove(largest_frag)
            salt_smi_list = [Chem.MolToSmiles(frag) for frag in fragments_list]
            salt = '.'.join(salt_smi_list)

            df.at[index, 'active ingredient'] = active
            df.at[index, 'salt'] = salt
        except:
            print(index, smi)
            list_to_remove.append(index)
    else:
        df.at[index, 'active ingredient'] = smi
        df.at[index, 'salt'] = '.'


# df_merged = df.merge(df_smiles, on='CIDX', how='left')
final_columns = ['smiles code', 'active ingredient', 'salt', 'relation', 'assay value', 'target disease',
                'protein target', 'treatment concentration', 'treatment unit', 'measurement', 'type',
                'model type', 'description', 'dataset', 'CMPD_ID']

df_final = df[final_columns]
df_final.to_csv('set25_result.csv', index=False)