#!/usr/bin/env streamlit

import streamlit as st
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
# from statsmodels.stats.multitest import multipletests


class Main:

    def __init__(self):

        # Area for messages
        # Only one message will be written at a time
        self.msg_empty = st.empty()

        # Area for displays
        self.disp_empty = st.empty()

        # Show the user the uploaded data
        self.abund_expander = st.expander("Abundances").empty()
        self.metadata_expander = st.expander("Metadata").empty()
        # self.kegg_expander = st.expander("KEGG").empty()

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

        # # Let the user upload a KEGG CSV
        # self.kegg_uploader = (
        #     st
        #     .sidebar
        #     .file_uploader(
        #         "KEGG Annotations (CSV)",
        #         help="Upload a CSV file with the KEGG annotations for each metabolite"
        #     )
        # )

        self.abund = None
        self.metadata = None
        # self.kegg = None

        self.run()

    def run(self):

        self.parse_abund()
        self.parse_metadata()
        # self.parse_kegg()

        self.stats()

    def stats(self):
        if self.abund is None or self.metadata is None:
            return

        # Pick the variable to group by
        guess_group_by = self.metadata.apply(lambda c: c.unique().shape[0]).sort_values().index.values[-1]

        # Let the user pick the variable to group by
        group_by = st.sidebar.selectbox(
            "Group By",
            options=list(self.metadata.columns.values),
            index=list(self.metadata.columns.values).index(guess_group_by)
        )

        # For the other columns, run a paired t-test
        res = self.run_stats(group_by)

        # Set up an area to display a few things
        cont = self.disp_empty.container()
        cont.dataframe(res)

        # Let the user select one to plot
        metab = cont.selectbox(
            "Metabolite to plot",
            res.index.values
        )

        # Let the user select the metadata of interest
        cname = cont.selectbox(
            "Show across metadata",
            [cn for cn in self.metadata.columns.values if cn != group_by]
        )

        # Get the values for that plot
        plot_df = (
            pd.DataFrame({
                vn: self.get_vals(group_by, cname, vn, metab)
                for vn in self.metadata[cname].unique()
            })
            .dropna()
            .T
        )

        fig, ax = plt.subplots()
        plt.xlabel(cname)
        plt.ylabel("Abundance")
        plt.title(metab)
        plot_df.plot(kind="line", ax=ax)
        cont.pyplot(fig)

    def run_stats(self, group_by):

        return (
            pd.DataFrame({
                f"{v1} vs {v2}": self.paired_ttest(group_by, cname, v1, v2)
                for cname in self.metadata.columns.values
                if cname != group_by
                for v1 in self.metadata[cname].unique()
                for v2 in self.metadata[cname].unique()
                if v1 < v2
            })
            .assign(MIN=lambda d: d.min(axis=1))
            .sort_values(by="MIN")
        )

    def paired_ttest(self, group_by, cname, v1, v2):

        res = []

        for metab in self.abund.columns.values:

            # Get the v1 and v2 values, merged across replicates
            comp = pd.DataFrame({
                v1: self.get_vals(group_by, cname, v1, metab),
                v2: self.get_vals(group_by, cname, v2, metab)
            }).dropna()
            res.append(
                dict(
                    pvalue=stats.ttest_rel(
                        comp[v1].values,
                        comp[v2].values
                    ).pvalue,
                    metabolite=metab
                )
            )

        return pd.DataFrame(res).set_index("metabolite")['pvalue']

    def get_vals(self, group_by, cname, val, metab):

        return (
            self.metadata.query(
                f"{cname} == '{val}'" if isinstance(val, str) else f"{cname} == {val}"
            )
            .assign(
                **{val: self.abund[metab]}
            )
            .drop(columns=[cname])
            .groupby(
                group_by
            )
            [val]
            .mean()
        )

    def parse_abund(self):
        if self.abund_uploader:
            try:
                self.abund = pd.read_csv(
                    self.abund_uploader,
                    index_col=0
                )
            except Exception as e:
                self.msg_empty.write(
                    f"Error reading abundances CSV:\n\n{str(e)}"
                )
                return

            # Drop any columns without top-line header
            self.abund = self.abund.reindex(
                columns=[
                    cname for cname in self.abund.columns.values
                    if not cname.startswith("Unnamed: ")
                ]
            )

            # Convert all values to float, or NaN
            self.abund = self.abund.applymap(self.makefloat)

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
                self.disp_empty.write(
                    f"Found missing values in rows:\n\n - {empties}"
                )
                self.abund = None
                return

            self.abund = self.abund.T.sort_index(axis=1)

            self.abund_expander.write(self.abund)

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
                self.msg_empty.write(
                    f"Error reading metadata CSV:\n\n{str(e)}"
                )
                return

            # Make sure that the metadata matches the abundances
            abund_ix = set(self.abund.index.values)
            metadata_ix = set(self.metadata.index.values)

            if not abund_ix == metadata_ix:

                msg = "Samples in metadata do not match abundances"
                if len(abund_ix - metadata_ix) > 0:
                    for s in list(abund_ix - metadata_ix):
                        msg = f"{msg}\n\nIn Abundances only: {s}"
                if len(metadata_ix - abund_ix) > 0:
                    for s in list(metadata_ix - abund_ix):
                        msg = f"{msg}\n\nIn Metadata only: {s}"
                self.disp_empty.write(msg)
                self.metadata = None
                return

            self.metadata_expander.write(self.metadata)

        else:
            self.disp_empty.write(
                "Please provide table of metadata"
            )

    def parse_kegg(self):
        # If there is no metadata, stop
        if self.metadata is None:
            return

        if self.kegg_uploader:

            try:
                self.kegg = pd.read_csv(
                    self.kegg_uploader,
                    index_col=0
                )
            except Exception as e:
                self.msg_empty.write(
                    f"Error reading KEGG CSV:\n\n{str(e)}"
                )
                return

            # Only keep the first column
            self.kegg = self.kegg[
                self.kegg.columns.values[0]
            ]

            # Only keep those compounds which are in the abundances
            self.kegg = self.kegg.loc[
                [
                    ix for ix in self.kegg.index.values
                    if ix in self.abund.columns.values
                ]
            ]

            self.kegg_expander.write(self.kegg.sort_index())

        else:
            self.disp_empty.write("Please provide KEGG CSV")


if __name__ == "__main__":
    Main()
