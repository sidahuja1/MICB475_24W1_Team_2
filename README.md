# MICB475_24W1_Team_2

**Oct 1 Agenda**
1. Review collected datasets 
2. Planning action items for the week

Potential Papers:
1 or 4

Aim 1: data wrangling, combining metadata, creating a new column

Aim 2: cross compare across the 6 groups, with indicator, differential abundance, core microbiome, diversity metrics, etc. (5 analysis)

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


**Oct 15 Agenda** 
1. Updates on data-wrangling
2. Discuss Halfvarson (short vs long time points, how we want to subset)
3. Proposal discussion - how to split tasks
- Replace section 8 with the new section
- Use an "ing" word in title -> "comparing" how X and X
- Intro and background: IBD background, explain signficance of research
- Objective: broader question and hypothesis; hypothesis must include a biological rationale
- Experimental aims: 1. wrangling (so categories are can be compared); 
- Approach: analyses
- Data wrangling instead of 

**Halfvarson**: Subsetting numbers
- Healthy-control: n = 9
- IBD_no_surg: n = 19
- IBD_surg_short (timepoint 1): n = 19
- IBD_surg_long (timepoint 2): n = 19

Proposal:

Title: 
- Intro and Background: Victoria and Georgia
- Research Objectives: Alex
- Experimental aims: Lester
- Approach: Sid
- Flowchart: Sid
- Weekly timeframe (Gantt Chart): Alex and Lester
- Data Wrangling: all
- Participation Report: all

**Finish individual parts by Friday**

Oct 22nd 
1. Give updates on data processing
2. Go over what needs to be done in the next week 

- Silva not greengenes

Oct 29:
1. Updates on Wrangling
2. Merge?
3. Proposal revisions


Oct 29th meeting minutes

**Completed previous week: **
- Processing finished  
- Merging finished 

**Notes from meeting:** 
Merge 
- Merge with taxonomy stuff → generate a merged taxonomy file Taxonomy
    - merge rep seqs to make taxonomy files 
- Make diversity metrics on qiime2
- Someone can start making phyloseq object
- Once we have file seq, rest of analysis flows
- Run everything and then we can cut down data after
- Don't do diversity metrics yet on R as we want the crude ones from quiime 2 first
- Generate and send evelyn alpha rarefaction curve before we rarefy 

**Tasks for the week:**
- Alex/sid will do the merging and taxonomy, 
Note: dont bother with taxa bar plots just need taxonomy file and tree file
- Then rareify → send Evelyn the alpha rarefaction curve to make sure we chose correctly 
- We will split up the rest of the work once the taxonomy file is created 

**Next meeting** 
- Have some diversity metrics ready
- Decide if we want to do a functional analysis
- Every one fix their part over the weekend → specifically address the comments 



Nov 5:
1. Troubleshooting merge:
- Our sampling depth options are extremely low, but the rarefaction curve looks fine.
- Discuss how to fix merge and divide up this work 
2. Confirm proposal changes, and ensure all comments have been addressed: everyone read over the changes
- After confirming changes, everyone will complete fixes to references in their section tonight 
3. Distribute analyses between each member to be done by the next meeting
- Reconfirm which analyses we want to do and plan out completion of these


**Nov 5 Meeting Minutes:**
- Issue with the denoising of fmt paper, rerun denoising for FMT

Nov 5th Meeting minutes 

**Merging issue **
- Reads for FMT paper (4) is too small, so we will need to check the table for just the FMT paper
- Looked at FMT demux file 
- Something went wrong with the denoising of the FMT→ 41 samples and only 989 reads 
- Redo from denoising step and then do merging of Mintz and halverson
- Trim Mintz at 220 forward and 220 reverse when denoising again → if it takes longer than  a day cut it → do a detached screen and check every few hours 
- Code seems to be fine, but we should rerun the mintz denoising as this seems to be where the issue came from 
- Try and get the diversity metrics as soon as possible 
**References **
- Use a references program for manuscript so we don't have to keep making changes 

**Things to do **
- Relabel to make files more clear 
- Rerun from demux and denoise again mintz (p4) 
- Finish references for resubmission of proposal 
- Organize Github more→ go through and make sure we have all files and that they are well labeled 

**NOTE: no meeting next week **

Nov 19:
1. Discuss merged results/rarefaction curve
2. Review diversity metrics
3. Distribute analyses to perform before the next meeting


Nov 19 meeting minutes 
Meeting minutes:
Issues with denoising are now worked out
Results 
- Looked at faith_pd → fmt shows less diversity then IBD no surgery
- Ibd reduced diversity
- Weighted unifrac: IBD surg clustered together → FMT clustered with healthy → surgery was more drastic so this makes sense
- FMT gives more of a profile that is closer to the healthy individuals
- Unweighted unifrac: healthy IBD no fmt was significantly dif, with no surgery also significantly different
- Abundance does not appear to matter much
- Same patterns as weight unifrac → good!
- Alpha diversity not super exciting but beta diversity is
- evenness → healthy is sig dif from surgery
- Fmt has profile that is closer to healthy
- Surgery appears to go backwards → resections show a loss in diversity
- Fmt restores diversity (directly building flora to be similar to healthy), but surgery reduces diversity giving you chance to rebuild flora on own
- Consistent results
- Proof of concept we outlined in our research, but this comparison has likely not been done before→ we had never seen how drastic it was when compared to same control 

Next step:
- Individuals analysis→ 5→ diversity, indicator species, core microbiome, deseq→ one person will also build the phyloseq
- Each one of us will do one
- Lester: phyloseq object
- Sid: recreate diversity metric
- Alex: indicator species
- Victoria: core microbiome
- Deseq: Georgia 
R will give us a nicer figure for the manuscript 
Don't worry too much about look

Next week: 
Everyone have their analysis done and we will look at all of them and see which ones to include in manuscript 
