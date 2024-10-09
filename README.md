# MICB475_24W1_Team_2

**Oct 1 Agenda** 
1. Review collected datasets 
2. Planning action items for the week

Potential Papers:
1 or 4

Aim 1: data wrangling, combining metadata, creating a new column

Aim 2: cross compare across the 6 groups, with indicator, differential abundance, etc. 

*Potentially look at functional differences


**To do**:
Wrangle data this week:
- Tabulate differences between the data, what variable regions, sample size, etc.
- Decide if we want to combine before or after importing
- Create a new metadata table
- Treatment group needs to be common between the two
- Treatment column and treatment detailed
- Origin column (halfverson, FMT)

**P1 Metadata**
- Host_disease includes CD, UC, or n/a (donor)
- Column "Origin" = the dataset the data is from 


**Oct 8 Agenda** 
1. Review differences between the data: discuss potential concerns about sample size (Georgia)
3. Updates on data wrangling
4. Planning action items for the week: demultiplexing data, denoising, clustering, etc.
*Figure out how to distribute tasks


Meeting minutes
- Data Wrangling (ignore dataset overview)
- Subset Halverson: 20 people to 20 IBD to match paper 4: 19 patients and 19 donors

**- Merge afterwards
**
To do:
- Read through the proposal outline
- Wrangle

- Columns:
- Keep diseases as IBD, 

- Paper 4:
- healthy-control, IBD-no_FMT, IBD_FMT-short, IBD_FMT-long


- Halverson: only 
- healthy-control, IBD_no_surgery, IBD_surgery

- Healthy 20 people
- IBD surgery of 3-month time point, 


Subset in the duplier package

Manifest: has sample-id and absolute-file path; in R filter the manifest based off the metadata


