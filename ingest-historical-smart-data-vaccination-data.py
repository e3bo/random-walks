#!/usr/bin/env python
# coding: utf-8

# In[8]:


import boto3
from botocore.config import Config


# In[9]:


oh_config = Config(region_name = 'us-east-2')


# In[10]:


client = boto3.client('dataexchange', config=oh_config)
dsets = client.list_data_sets(Origin='ENTITLED')


# In[24]:


for ds in dsets['DataSets']:
    if ds['Name'] == 'Free: US Covid Vaccinations Data':
        sourceid = ds['Id']


# In[29]:


rev_paginator = client.get_paginator('list_data_set_revisions')
rev_iterator = rev_paginator.paginate(DataSetId=sourceid)


# In[30]:


rids = []
for page in rev_iterator:
    for revision in page['Revisions']:
        rcomments.append(revision['Comment'])
        rupdate.append(revision['UpdatedAt'])
        rids.append(revision['Id'])


# In[37]:


asset_ids = []
asset_names = []
asset_mtimes = []
for rid in rids:
    asset = client.list_revision_assets(DataSetId=sourceid, RevisionId=rid)['Assets'][0]
    asset_ids.append(asset['Id'])
    asset_names.append(asset['Name'])
    asset_mtimes.append(asset['UpdatedAt'])
    


# In[45]:


partnames = [mt.strftime("%Y%m%d") for mt in asset_mtimes]


# In[56]:


keys = []
for (p, c) in zip(partnames, asset_names):
    keys.append('smart-vacc' + '/' + p + '/' + c)


# In[63]:


cj_responses = []
for (a, k, r) in zip(asset_ids, keys, rids):
    resp = client.create_job(
    Details={
        'ExportAssetsToS3': {
            'AssetDestinations': [
                {
                    'AssetId': a,
                    'Bucket': 'covid19-data-aquisitions',
                    'Key': k
                }
            ],
            'DataSetId': sourceid,
            'RevisionId': r
        }},
        Type='EXPORT_ASSETS_TO_S3'
    )
    cj_responses.append(resp)


# In[72]:


[r['ResponseMetadata']['HTTPStatusCode'] for r in start_job_responses]

