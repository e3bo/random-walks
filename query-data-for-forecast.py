#!/usr/bin/env python
# coding: utf-8

# In[1]:


import boto3
import datetime

client = boto3.client('athena')


# In[2]:


def gen_query(asof):
    query = """
WITH h AS (SELECT date,
           partition_1 AS fips,
           previous_day_admission_adult_covid_confirmed, 
           previous_day_admission_pediatric_covid_confirmed, 
           previous_day_admission_adult_covid_confirmed_coverage,
           partition_0 AS asof 
           FROM healthdata 
           WHERE partition_1='06' AND partition_0='{asof}'),
d AS (SELECT CAST(DATE_FORMAT(target_end_date, '%Y%m%d') AS BIGINT) AS date, 
      location AS fips,
      value AS deaths, 
      partition_0 AS asof 
      FROM hopkins_reformatted 
      WHERE target_type='day ahead inc death' AND location='06' AND partition_0='{asof}'),
c AS (SELECT CAST(DATE_FORMAT(target_end_date, '%Y%m%d') AS BIGINT) AS date, 
      location AS fips,
      value AS cases, 
      partition_0 AS asof 
      FROM hopkins_reformatted 
      WHERE target_type='day ahead inc case' AND location='06' AND partition_0='{asof}')
SELECT fips,
       asof, 
       date, 
       d.deaths, 
       c.cases,
       h.previous_day_admission_adult_covid_confirmed,
       h.previous_day_admission_pediatric_covid_confirmed,
       h.previous_day_admission_adult_covid_confirmed_coverage
FROM h
FULL OUTER JOIN d
USING (fips, asof, date)
FULL OUTER JOIN c
USING (fips, asof, date)
WHERE date >= 20200301
ORDER BY date;""".format(asof = asof)
    return query


# In[3]:


numweeks = 44
start = datetime.datetime.strptime("2020-06-29", "%Y-%m-%d")
date_list = [start + datetime.timedelta(days=x *7) for x in range(numweeks)]


# In[5]:


responses = []
for date in date_list:
    asof = date.strftime("%Y-%m-%d")
    query = gen_query(asof)
    output_loc = 's3://covid19-data-aquisitions/forecast-inputs/' + asof
    resp = client.start_query_execution(
        QueryString=query,
        QueryExecutionContext={
            'Database': 'covid19'
        },
        ResultConfiguration={
            'OutputLocation': output_loc
        }
    )
    responses.append(resp)

