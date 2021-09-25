#!/usr/bin/env python
# coding: utf-8

# In[1]:


import boto3
import sagemaker


# In[2]:


print(boto3.__version__)
print(sagemaker.__version__)


# In[5]:


client = boto3.client('dataexchange')


# In[6]:


paginator = client.get_paginator('list_revision_assets')


# In[23]:


response_iterator = paginator.paginate(
    DataSetId='0af53e1b6b368f8ba644927aa540d9b4',
    RevisionId='c24b88baa8545e036e9ea027ae639937')


# In[24]:


target_assets_names = ["delphi-covidcast-covid-19/dataset/csv/jhu-csse/confirmed_incidence_num/day/state.csv",
                "delphi-covidcast-covid-19/dataset/csv/jhu-csse/deaths_incidence_num/day/state.csv",
                 "delphi-covidcast-covid-19/dataset/jsonl/hhs/confirmed_admissions_covid_1d/day/state.jsonl"]


# In[25]:


target_assets = []


# In[26]:


for page in response_iterator:
    for asset in page['Assets']:
        if asset['Name'] in target_assets_names:
            target_assets.append(asset)


# In[27]:


target_assets


# In[28]:


target_assets[0]['Id']


# In[43]:


cj_response = client.create_job(
    Details={
        'ExportAssetsToS3': {
            'AssetDestinations': [
                {
                    'AssetId': target_assets[1]['Id'],
                    'Bucket': 'covid19-data-aquisitions',
                    'Key': target_assets[1]['Name']
                },
                {
                    'AssetId': target_assets[2]['Id'],
                    'Bucket': 'covid19-data-aquisitions',
                    'Key': target_assets[2]['Name']
                },
            ],
            'DataSetId': target_assets[0]['DataSetId'],
            'RevisionId': target_assets[0]['RevisionId']
        }},
    Type='EXPORT_ASSETS_TO_S3'
)


# In[44]:


cj_response


# In[45]:


start_job_response = client.start_job(JobId = cj_response['Id'])


# In[46]:


start_job_response


# In[45]:


for asset in response['Assets']:
    print(asset['Name'])


# In[39]:


s3 = boto3.resource('s3')


# In[40]:


bucket = s3.Bucket('covid19-data-aquisitions')


# In[47]:


for bobject in bucket.objects.all():
    print(bobject)


# # Tutorial

# In[17]:


sqs = boto3.resource('sqs')
queue = sqs.create_queue(QueueName='test', Attributes={'DelaySeconds': '5'})
print(queue.url)
print(queue.attributes.get('DelaySeconds'))


# In[18]:


print(queue.url)


# In[22]:


for queue in sqs.queues.all():
    print(queue.attributes['QueueArn'].split(':')[-1])


# In[23]:


response = queue.send_message(MessageBody='world')


# In[24]:


print(response.get('MessageId'))


# In[25]:


print(response.get('MD5OfMessageBody'))


# In[27]:


queue.send_message(MessageBody='boto3', MessageAttributes={
    'Author': {
        'StringValue' : 'Daniel',
        'DataType': 'String'
    }
})


# In[31]:


response = queue.send_messages(Entries=[
    {
        'Id': '1',
        'MessageBody': 'world'
    },
    {
        'Id': '2',
        'MessageBody': 'boto3',
        'MessageAttributes': {
            'Author': {
                'StringValue': 'Daniel',
                'DataType': 'String'
            }
        }
    }
])

# Print out any failures
print(response.get('Failed'))


# In[37]:


for message in queue.receive_messages(MessageAttributeNames=['Author']):
    author_text = ''
    if message.message_attributes is not None:
        author_name = message.message_attributes.get('Author').get('StringValue')
        if author_name:
            author_text = ' ({0})'.format(author_name)
    
    print('Hello, {0}!{1}'.format(message.body, author_text))
    message.delete()


# In[ ]:





# In[48]:


1 - .9**50


# # Scraps

# In[40]:


response = client.list_revision_assets(
    DataSetId='0af53e1b6b368f8ba644927aa540d9b4',
    MaxResults=200,
    RevisionId='c24b88baa8545e036e9ea027ae639937'
)

