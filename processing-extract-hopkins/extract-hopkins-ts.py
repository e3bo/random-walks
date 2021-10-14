#!/usr/bin/env python
# coding: utf-8

# In[1]:


import boto3
import sagemaker
from sagemaker.processing import Processor, ProcessingInput, ProcessingOutput

region = boto3.session.Session().region_name
account = boto3.client("sts").get_caller_identity().get("Account")
role = 'arn:aws:iam::' + account + ':role/service-role/AmazonSageMaker-ExecutionRole-20210729T125479'


# In[4]:


domain = (
    "amazonaws.com.cn"
    if (region == "cn-north-1" or region == "cn-northwest-1")
    else "amazonaws.com"
)

processor = Processor(
    image_uri="{}.dkr.ecr.{}.{}/cssetransform:latest".format(account, region, domain),
    role=role,
    instance_count=1,
    instance_type="local"
)


# In[23]:


bucket = 'covid19-data-aquisitions'
prefix = 'hopkins/'  
case_file = "time_series_covid19_confirmed_US.csv"
death_file = "time_series_covid19_deaths_US.csv"
reform_file = "time-series-covid19-US.parquet"

client = boto3.client('s3')
s3 = boto3.resource('s3')
bucket_obj = s3.Bucket('covid19-data-aquisitions')
output_prefix = 'hopkins-reformatted/'
output_objs = list(bucket_obj.objects.filter(Prefix=output_prefix))

result = client.list_objects(Bucket=bucket, Prefix=prefix, Delimiter='/')
for o in result.get('CommonPrefixes'):
    stem = "s3://" + bucket + "/" + o.get('Prefix')
    print("Running reformat job for " + stem)
    dsource = stem + death_file
    csource = stem + case_file
    oprefix = o.get('Prefix').replace('hopkins', 'hopkins-reformatted')
    odest = "s3://" + bucket + "/" + oprefix
    
    if any([w.key == oprefix + reform_file for w in output_objs]):
        print("Output already exists, skipping")
    else:
        print("No ouput, running processor")
        processor.run(
            inputs=[ProcessingInput(source=csource, 
                                    destination="/opt/ml/processing/input/cases"),
                    ProcessingInput(source=dsource, 
                                    destination="/opt/ml/processing/input/deaths")],
            outputs=[ProcessingOutput(source="/opt/ml/processing/output",
                                      destination=odest)]
        )

