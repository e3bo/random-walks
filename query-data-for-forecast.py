#!/usr/bin/env python
# coding: utf-8

# In[1]:


import boto3
import datetime
import time

client = boto3.client('athena')


# In[2]:


def gen_query(asof):
    query = """
WITH g AS (
	SELECT DATE_FORMAT(date, '%Y%m%d') AS date2,
		residential_percent_change_from_baseline AS res,
		REPLACE(iso_3166_2_code, 'US-', '') as state,
		partition_0 AS asof
	FROM google_mobility_reports_wayback_parquet
	WHERE partition_0 = '{asof}'
),
g2 AS (
	SELECT g.date2 AS time_value,
		res,
		state,
		asof,
		state_code AS fips
	FROM g
		LEFT JOIN state_fips_lookup USING (state)
),
d AS (
	SELECT DATE_FORMAT(target_end_date, '%Y%m%d') AS time_value,
		location AS fips,
		value AS deaths,
		partition_0 AS asof
	FROM hopkins_reformatted
	WHERE target_type = 'day ahead inc death'
		AND location = '06'
		AND partition_0 = '{asof}'
),
c AS (
	SELECT DATE_FORMAT(target_end_date, '%Y%m%d') AS time_value,
		location AS fips,
		value AS cases,
		partition_0 AS asof
	FROM hopkins_reformatted
	WHERE target_type = 'day ahead inc case'
		AND location = '06'
		AND partition_0 = '{asof}'
),
h AS (
	SELECT CAST(date AS VARCHAR) AS time_value,
		partition_1 AS fips,
		previous_day_admission_adult_covid_confirmed,
		previous_day_admission_pediatric_covid_confirmed,
		previous_day_admission_adult_covid_confirmed_coverage,
		partition_0 AS asof
	FROM healthdata
	WHERE partition_1 = '06'
		AND partition_0 = '{asof}'
)
SELECT time_value,
	asof,
	fips,
	g2.res,
	c.cases,
	d.deaths,
	h.previous_day_admission_adult_covid_confirmed,
	h.previous_day_admission_pediatric_covid_confirmed,
	h.previous_day_admission_adult_covid_confirmed_coverage
FROM d
	FULL JOIN g2 USING (time_value, asof, fips)
	FULL JOIN c USING (time_value, asof, fips)
	FULL JOIN h USING (time_value, asof, fips)
WHERE fips = '06'
	AND time_value > '20200301'
ORDER BY time_value;""".format(asof = asof)
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
    time.sleep(1)

