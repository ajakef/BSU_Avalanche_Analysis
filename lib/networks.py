import pandas as pd
LCC_2020 = pd.DataFrame(columns = ['Station', 'CCUBE', 'Datacube'])
LCC_2020.Station = ['LCC0', 'LCC1', 'LCC2']
LCC_2020.CCUBE = ['CC129', 'CC147', 'CC125']
LCC_2020.Datacube = ['ADA', 'BE4', 'AD9']

aftershocks_2020 = pd.DataFrame(columns = ['Station', 'CCUBE', 'Datacube'])
aftershocks_2020.Station = ['IRON']
aftershocks_2020.CCUBE = ['CC147']
aftershocks_2020.Datacube = ['BE4']

Scott_Demo = pd.DataFrame(columns = ['Station', 'CCUBE', 'Datacube'])
Scott_Demo.Station = ['CC125']
Scott_Demo.CCUBE = ['CC125']
Scott_Demo.Datacube = ['AD8']

network = Scott_Demo
