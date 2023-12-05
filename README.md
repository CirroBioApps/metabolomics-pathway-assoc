# metabolomics-pathway-assoc
Analysis of metabolomics data for pathway analysis

[Streamlit](https://streamlit.io/)-based app for the analysis of metabolomics data.

## Data Inputs

1. Table of sample metadata
2. Table of metabolite abundances

## Data Assumptions

1. The first **column** of sample metadata contains the unique sample identifier
2. Arbitrary sample metadata may be included as additional columns
3. The first **row** of the metabolite abundances contains the unique sample identifier
4. The first **column** of the metabolite abundances contains the metabolite name

## Interaction

1. User must indicate which column from the sample metadata table is used to group
together samples from the same person
2. After association tests are performed, a specific metabolite may be selected for
display

## Analysis

Consider the example case where metadata has been provided for:

| sample | patient | timepoint |
| ------ | ------- | --------- |
| p1_pre | p1      | pre       |
| p1_mid | p1      | mid       |
| p1_post | p1     | post     |
| p2_pre | p2      | pre       |
| p2_mid | p2      | mid       |
| p2_post | p2     | post     |
| p3_pre | p3      | pre       |
| p3_mid | p3      | mid       |
| p3_post | p3     | post     |

The user has selected the `patient` column with the 'Group By' menu.

For each of the other metadata columns, for each of the unique pairs
of values in that column, a paired t-test will be run for each metabolite
which compares the relative abundance of each metabolite within each
patient.

The comparisons performed for this dataset would then be:

1. `pre` vs. `mid`
2. `pre` vs. `post`
3. `mid` vs. `post`

The results of the analysis will be presented for each metabolite,
for each of those comparisons.
