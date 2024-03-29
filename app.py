#!/usr/bin/env streamlit

from collections import defaultdict
import json
import requests
from typing import Union
import numpy as np
import streamlit as st
import pandas as pd
from time import sleep
from scipy import stats
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
import plotly.express as px


class Main:

    abund: Union[None, pd.DataFrame]

    def __init__(self):

        # Area for logging messages
        self.log_empty = st.empty()
        self.log_container = self.log_empty.expander("Logs", expanded=True)

        # Area for displays
        self.disp_empty = st.empty()
        self.disp_container = self.disp_empty.container()
        self.disp_stats = self.disp_container.container()
        self.disp_plots = self.disp_container.container()
        self.disp_kegg = self.disp_container.container()
        self.res = None

        # Let the user upload a metadata CSV
        self.metadata_uploader = (
            st
            .sidebar
            .file_uploader(
                "Metadata (CSV)",
                help="Upload a CSV file with the experimental design"
            )
        )

        # Let the user upload an abundances CSV
        self.abund_uploader = (
            st
            .sidebar
            .file_uploader(
                "Abundances (CSV)",
                help="Upload a CSV file with the observed abundances"
            )
        )

        # Let the user upload a KEGG CSV
        self.kegg_uploader = (
            st
            .sidebar
            .file_uploader(
                "KEGG Annotations (CSV)",
                help="Upload a CSV file with the KEGG annotations for each metabolite"
            )
        )

        self.abund = None
        self.metadata = None
        self.kegg = None

        self.run()

    def log(self, msg):
        ct = datetime.datetime.now()
        self.log_container.write(f"[{ct}] {msg}")

    def run(self):

        self.parse_abund()
        self.parse_metadata()

        self.stats()
        self.parse_kegg()
        self.kegg_pathways()

    def stats(self):
        if self.abund is None or self.metadata is None:
            return

        self.disp_stats.markdown("## Statistical Analysis")

        # Pick the column used to compare samples
        compare_by = self.disp_stats.selectbox(
            "Compare Samples By:",
            options=list(self.metadata.columns.values),
            index=list(self.metadata.columns.values).index(
                self.metadata.apply(
                    lambda c: c.unique().shape[0]
                ).sort_values().index.values[0]
            )
        )

        # Get the control and the comparitor values
        groups = self.metadata[compare_by].drop_duplicates().sort_values().tolist()
        if len(groups) == 1:
            self.log(f"Cannot compare - all samples have {compare_by} == {groups[0]}")
            return

        control_group = self.disp_stats.selectbox(
            "Control Group:",
            options=groups
        )
        comparison_group = self.disp_stats.selectbox(
            "Comparison Group:",
            options=[gn for gn in groups if gn != control_group]
        )

        # Let the user decide whether to run a paired or unpaired analysis
        if self.metadata.shape[1] < 2 or self.disp_stats.selectbox(
            "Statistical Test:",
            options=[
                "Paired t-test",
                "Unpaired t-test"
            ]
        ) == "Unpaired t-test":

            # Run an unpaired t-test for the selected set of groups
            self.res = self.unpaired_ttest(compare_by, control_group, comparison_group)

            pair_by = None

        else:
            # If a paired approach was selected
            # Ask the user which column to pair samples by
            pair_by = self.disp_stats.selectbox(
                "Pair By:",
                help="Select the column used to match up pairs of samples for comparison",
                options=[
                    cname for cname in self.metadata.columns.values
                    if cname != compare_by
                ]
            )

            # There can only be a single instance of each value per group
            for (compare_val, pair_val), d in self.metadata.groupby([compare_by, pair_by]):
                if compare_val not in [control_group, comparison_group]:
                    continue
                if d.shape[0] > 1:
                    msg = f"multiple samples found with {compare_by} == {compare_val} and {pair_by} == {pair_val}"
                    self.log(f"Cannot perform paired analysis using {pair_by} -- {msg}")
                    return

            # Run a paired t-test for the selected set of groups, and pairing variable
            self.res = self.paired_ttest(pair_by, compare_by, control_group, comparison_group)

        self.res = (
            self.res
            .assign(
                neg_log_pvalue=-self.res["pvalue"].apply(np.log10)
            )
            .sort_values(by="pvalue")
            .reset_index(drop=True)
            .reindex(
                columns=[
                    "metabolite",
                    "pvalue",
                    "neg_log_pvalue",
                    "tvalue",
                    "log_fold_change",
                    "mean_abund",
                ]
            )
        )

        self.disp_stats.markdown("### Results Table")
        self.disp_stats.dataframe(self.res, hide_index=True)

        self.disp_plots.markdown("### Plot Results")

        kwargs = dict(
            hover_name="metabolite",
            hover_data=["pvalue", "mean_abund"],
            color="neg_log_pvalue",
            color_continuous_scale="bluered",
            labels=dict(
                log_fold_change="Fold Change (log2)",
                neg_log_pvalue="p-value (-log10)",
                mean_abund="Mean Abundance",
                pvalue="p-value"
            )
        )

        # Volcano
        fig = px.scatter(
            data_frame=self.res,
            x="log_fold_change",
            y="neg_log_pvalue",
            title="Volcano Plot",
            **kwargs
        )
        fig.add_vline(x=0, line_width=1, line_color="grey")
        fig.add_hline(y=0, line_width=1, line_color="grey")
        self.disp_plots.plotly_chart(fig)

        # MA Plot
        ma_plot = self.disp_plots.empty()
        ma_log_y = self.disp_plots.checkbox("Log Y Scale", key="logy_maplot")
        fig = px.scatter(
            data_frame=self.res,
            x="log_fold_change",
            y="mean_abund",
            title="M-A Plot",
            log_y=ma_log_y,
            **kwargs
        )
        fig.add_vline(x=0, line_width=1, line_color="grey")
        if not ma_log_y:
            fig.add_hline(y=0, line_width=1, line_color="grey")
        ma_plot.plotly_chart(fig)

        self.disp_plots.markdown("### Plot Single Metabolite")
        # Let the user select one to plot
        metab = self.disp_plots.selectbox(
            "Metabolite to plot",
            self.res['metabolite'].values
        )

        plot_df = (
            self.metadata.assign(**{
                metab: self.abund[metab]
            })
            .reindex(
                columns=(
                    [compare_by, metab]
                    if pair_by is None else
                    [compare_by, metab, pair_by]
                )
            )
            .dropna()
        )
        plot_df = plot_df.loc[
            plot_df[compare_by].isin([control_group, comparison_group])
        ]

        kwargs = dict(
            x=compare_by,
            order=[control_group, comparison_group],
            y=metab
        )

        if pair_by is None:

            plot_type = self.disp_plots.selectbox(
                "Plot Type:",
                options=[
                    "box",
                    "boxen",
                    "violin",
                    "bar",
                    "point",
                    "strip",
                    "swarm"
                ]
            )

            fig, _ = plt.subplots()

            dict(
                box=sns.boxplot,
                boxen=sns.boxenplot,
                violin=sns.violinplot,
                bar=sns.barplot,
                point=sns.pointplot,
                strip=sns.stripplot,
                swarm=sns.swarmplot,
            )[plot_type](
                data=plot_df,
                **kwargs
            )

        else:
            # Only keep the observations with proper pairs
            vc = plot_df[pair_by].value_counts()
            plot_df = plot_df.loc[
                plot_df[pair_by].isin(
                    vc.index.values[vc == 2]
                )
            ]

            fig, _ = plt.subplots()
            sns.pointplot(
                data=plot_df,
                hue=pair_by,
                **kwargs
            )

        if self.disp_plots.checkbox("Log Y Scale"):
            plt.yscale("log")

        plt.show()
        self.disp_plots.pyplot(fig)

        self.disp_plots.dataframe(plot_df)

    def unpaired_ttest(self, compare_by, control_group, comparison_group):

        res_list = []

        for metab in self.abund.columns.values:

            control_vals = self.get_vals(compare_by, control_group, metab)
            comparison_vals = self.get_vals(compare_by, comparison_group, metab)

            res = stats.ttest_ind(control_vals, comparison_vals)

            log_fold_change = np.log2(comparison_vals.mean() / control_vals.mean())

            res_list.append(
                dict(
                    pvalue=res.pvalue,
                    tvalue=res.statistic,
                    mean_abund=np.mean([comparison_vals.mean(), control_vals.mean()]),
                    log_fold_change=log_fold_change,
                    metabolite=metab
                )
            )

        return pd.DataFrame(res_list)

    def paired_ttest(self, pair_by, compare_by, control_group, comparison_group):

        res_list = []

        for metab in self.abund.columns.values:

            # Get the v1 and v2 values, merged across replicates
            comp = pd.DataFrame({
                control_group: self.get_vals(compare_by, control_group, metab, index_by=pair_by),
                comparison_group: self.get_vals(compare_by, comparison_group, metab, index_by=pair_by)
            }).dropna()

            res = stats.ttest_rel(comp[control_group].values, comp[comparison_group].values)

            log_fold_change = np.log2((comp[comparison_group] / comp[control_group]).mean())

            res_list.append(
                dict(
                    pvalue=res.pvalue,
                    tvalue=res.statistic,
                    mean_abund=comp.mean().mean(),
                    log_fold_change=log_fold_change,
                    metabolite=metab
                )
            )

        return pd.DataFrame(res_list)

    def get_vals(self, cname, val, metab, index_by=None):

        vals = (
            self.metadata.query(
                f"{cname} == '{val}'" if isinstance(val, str) else f"{cname} == {val}"
            )
            .assign(
                **{val: self.abund[metab]}
            )
            .drop(columns=[cname])
        )
        if index_by is None:
            return vals[val].dropna()
        else:
            return vals.set_index(index_by)[val].dropna()

    def log_table_size(self, header: str, df: pd.DataFrame):

        nrows = df.shape[0]
        ncols = df.shape[1]
        self.log(f"{header} - {nrows:,} rows x {ncols:,} columns")

    def parse_abund(self):
        if self.abund_uploader:
            try:
                self.abund = pd.read_csv(
                    self.abund_uploader,
                    index_col=0
                )
            except Exception as e:
                st.error(
                    f"Error reading abundances CSV:\n\n{str(e)}"
                )
                return

            self.log_table_size("Abundance Table - Uploaded", self.abund)

            # Drop any columns without top-line header
            self.abund = self.abund.reindex(
                columns=[
                    cname for cname in self.abund.columns.values
                    if not cname.startswith("Unnamed: ")
                ]
            )

            # Convert all values to float, or NaN
            self.abund = self.abund.map(self.makefloat)

            # Drop any rows which are entirely lacking values
            self.abund = self.abund.loc[
                [
                    ix
                    for ix, r in self.abund.iterrows()
                    if not r.isnull().all()
                ]
            ]

            # Find any missing values
            empties = [
                f"Row: {row}\tCol: {col}"
                for row, r in self.abund.iterrows()
                for col, val in r.items()
                if val is None
            ]
            if len(empties) > 0:
                empties = '\n - '.join(empties)
                st.error(
                    f"Found missing values in rows:\n\n - {empties}"
                )
                self.abund = None
                return

            self.log_table_size("Abundance Table - Filtered", self.abund)

            self.abund = self.abund.T.sort_index(axis=1)

        else:
            self.disp_empty.write(
                "Please provide table of abundances"
            )

    @staticmethod
    def makefloat(v):
        try:
            return float(v)
        except Exception as e:
            return None

    def parse_metadata(self):
        # If there are no abundances, stop
        if self.abund is None:
            return

        if self.metadata_uploader:

            try:
                self.metadata = pd.read_csv(
                    self.metadata_uploader,
                    index_col=0
                )
            except Exception as e:
                st.error(
                    f"Error reading metadata CSV:\n\n{str(e)}"
                )
                return

            # Get the list of samples which are in the metadata and abundances
            abund_ix = set(self.abund.index.values)
            metadata_ix = set(self.metadata.index.values)
            overlap_ix = abund_ix & metadata_ix

            self.log(f"Number of samples defined in abundance table: {len(abund_ix):,}")
            self.log(f"Number of samples defined in metadata table: {len(metadata_ix):,}")
            self.log(f"Number of samples defined in both: {len(overlap_ix):,}")

            if len(overlap_ix) == 0:
                st.error("Labels do not match between metadata and abundances")
                self.abund = None
                return

            self.abund = self.abund.reindex(index=list(overlap_ix)).sort_index()
            self.metadata = self.metadata.reindex(index=list(overlap_ix)).sort_index()

        else:
            self.disp_empty.write(
                "Please provide table of metadata"
            )

    def parse_kegg(self):
        # If there is no metadata, stop
        if self.metadata is None:
            return

        self.disp_kegg.write("## KEGG Pathways")

        if self.kegg_uploader:

            try:
                if self.kegg_uploader.name.endswith(".xlsx"):
                    self.kegg = pd.read_excel(
                        self.kegg_uploader,
                        index_col=0,
                        header=None
                    )
                else:
                    self.kegg = pd.read_csv(
                        self.kegg_uploader,
                        index_col=0,
                        header=None
                    )
            except Exception as e:
                st.error(
                    f"Error reading KEGG file:\n\n{str(e)}"
                )
                return

            # Only keep the first column
            self.kegg = self.kegg[
                self.kegg.columns.values[0]
            ].dropna()

            self.log(f"Read in {self.kegg.shape[0]:,} KEGG IDs")

            # Only keep those compounds which are in the abundances
            self.kegg = self.kegg.loc[
                [
                    ix for ix in self.kegg.index.values
                    if ix in self.abund.columns.values
                ]
            ]

            self.log(f"KEGG IDs found in abundance table: {self.kegg.shape[0]:,}")
            self.abund = self.abund.reindex(columns=self.kegg.index.values)
            self.log(f"Only using the {self.abund.shape[1]:,} metabolites listed in the KEGG table")

        else:
            self.disp_kegg.write("Please provide KEGG CSV")

    def kegg_pathways(self):
        """Display the KEGG pathways for any significant genes."""

        if self.kegg is None or self.res is None:
            return

        max_p = self.disp_kegg.number_input(
            "Threshold: Maximum p-value",
            help="Only display pathways for metabolites with p-values less than or equal to this value",
            value=0.05
        )

        filt_res = self.res.query(f"pvalue <= {max_p}")

        self.disp_kegg.write(f"There are {filt_res.shape[0]:,} metabolites with p-value <= {max_p}")

        filt_res = filt_res.assign(
            KEGG_ID=filt_res["metabolite"].apply(self.kegg.get),
            KEGG_NAME=lambda d: d["KEGG_ID"].apply(get_kegg_name)
        ).reindex(
            columns=["KEGG_ID", "KEGG_NAME"] + list(self.res.columns.values)
        )
        self.disp_kegg.dataframe(filt_res, use_container_width=True)

        # Make a table of pathways
        pathways = pd.DataFrame({
            f'{r["KEGG_ID"]} {r["KEGG_NAME"]}': {
                pathway: 1
                for pathway in get_kegg_pathways(
                    r["KEGG_ID"]
                )
                if pathway != "None Found"
            }
            for _, r in filt_res.iterrows()
        }).fillna(0)

        self.disp_kegg.dataframe(pathways)

        min_metabolites = self.disp_kegg.number_input(
            "Show Pathways with >= X metabolites",
            min_value=1,
            value=1,
            step=1
        )

        pathways = pathways.loc[
            pathways.sum(axis=1) >= min_metabolites
        ]
        pathways = pathways.loc[:, pathways.any()]

        pathways = pathways.reindex(
            index=pathways.index.values[
                hierarchy.leaves_list(
                    hierarchy.linkage(pathways)
                )
            ],
            columns=pathways.columns.values[
                hierarchy.leaves_list(
                    hierarchy.linkage(pathways.T)
                )
            ]
        )

        fig, _ = plt.subplots()
        sns.heatmap(
            data=pathways,
            cmap="Blues"
        )
        self.disp_kegg.pyplot(fig)
        plt.close()

        # Select a pathway to display
        spacer = " (# of compounds = "
        pathway_options = [
            f"{pathway}{spacer}{n:,})"
            for pathway, n in (
                pathways
                .sum(axis=1)
                .apply(int)
                .sort_values(ascending=False)
                .items()
            )
        ]

        selected_pathway = self.disp_kegg.selectbox(
            "Display Pathway:",
            pathway_options
        )
        selected_pathway = selected_pathway.split(spacer)[0]

        # Get the list of positive and negatively associated compounds
        assoc_compounds = [
            compound.split(" ")[0] for compound, in_pathway in pathways.items()
            if in_pathway[selected_pathway]
        ]
        assoc_compounds = list(set(assoc_compounds))

        # Get the direction of association
        assoc_compounds = (
            filt_res
            .reindex(columns=["KEGG_ID", "log_fold_change"])
            .dropna()
            .set_index("KEGG_ID")
            .reindex(index=assoc_compounds)
            .dropna()
            ["log_fold_change"] > 0
        )

        map_url = self.format_kegg_map_url(
            selected_pathway.split(" ")[0],
            assoc_compounds
        )

        self.disp_kegg.markdown(
            f"""
            [Display Pathway Map: {selected_pathway}]({map_url})

            - Blue - increased in comparison vs. control group
            - Orange - decreased in comparison vs. control group
            """
        )

    @staticmethod
    def format_kegg_map_url(map_id, assoc_compounds):
        url = f"https://www.kegg.jp/kegg-bin/show_pathway?{map_id}"

        for comp, is_pos in assoc_compounds.items():
            color = "blue" if is_pos else "orange"
            url = f"{url}/{comp}%20{color}"

        return url


@st.cache_resource
def get_kegg_info(ko: str):

    ko = str(ko).strip()

    r = None
    for _ in range(10):
        r = requests.get("http://rest.kegg.jp/get/{}".format(ko))
        if r.status_code == 404:
            sleep(0.1)
        else:
            break

    output = defaultdict(list)
    header = None

    for line in r.text.split("\n"):
        if len(line) == 0:
            continue
        elif line[0] == " " and header is not None:
            output[header].append(line.strip())
        else:
            header = line.strip().split(" ")[0]
            if len(line.strip()) > len(header):
                val = line[len(header):].strip()
                if len(val) > 0:
                    output[header].append(val)

    return output


def get_kegg_name(kegg_id):
    if kegg_id is None:
        return "None Found"
    if get_kegg_info(kegg_id) is not None:
        if len(get_kegg_info(kegg_id)["NAME"]) == 0:
            return "None Found"
        else:
            return get_kegg_info(kegg_id)["NAME"][0]
    else:
        return "None Found"


def get_kegg_pathways(kegg_id):
    if kegg_id is None:
        return ["None Found"]
    if get_kegg_info(kegg_id) is not None:
        return get_kegg_info(kegg_id)["PATHWAY"]
    else:
        return []


if __name__ == "__main__":
    Main()
