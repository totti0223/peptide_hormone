import streamlit as st
import pandas as pd
import plotly.express as px
# fastaファイルをもとにsignalp6.0で結果をgff3で得てoutput.gff3と統合して得られたcsvファイルがスタート

class Condition:
    def __init__(self, condition_func, *args, **kwargs):
        self.condition_func = condition_func
        self.args = args
        self.kwargs = kwargs

    def __call__(self, value):
        return self.condition_func(value, *self.args, **self.kwargs)


def is_shorter_than(val, length):
    return len(val) <= length

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5843760/にフィルタ候補

df = pd.read_csv("annotated_230220.csv", index_col=0)
_df = df.copy()
conditions = {}
condition_no = 0

st.title("Signal Peptide Manual Filter Tool")
col1, col2 = st.columns(2)

with col1:
    st.header("Criteria")
with col2:
    st.header("Summary")
    st.subheader("Start")
    st.write("{} peptide".format(len(_df)))


with col1:
    st.subheader("SignalP6.0")
    key = "{}_condition".format(condition_no)
    conditions[key] = {}
    conditions[key]["name"] = "SignalP6"
    chk_spp = st.checkbox("enable", value=True, key="chk_spp")
    conditions[key]["flag"] = chk_spp
    if chk_spp:
        fil_spscore = st.slider("signal_peptide positive score over than", 0.0, 1.0, 0.95)
        conditions[key]["value"] = fil_spscore
        conditions[key]["type"] = "over"
        fig = px.histogram(_df["SP6_score"], x="SP6_score", nbins=20)
        fig.update_layout(height=180)
        fig.add_vline(x=fil_spscore, line_color="red")
        st.plotly_chart(fig,use_container_width=True,use_container_height=False)
        with col2:
            st.subheader("{} filter: {} {}".format(conditions[key]["name"],
                                               conditions[key]["type"],
                                               conditions[key]["value"]))
            _df = _df[_df["SP6_score"] > conditions[key]["value"]]
            st.write("{} peptide".format(len(_df)))
            st.dataframe(_df,height=200)

condition_no += 1
key = "{}_condition".format(condition_no)

with col1:
    st.subheader('Peptide Length')
    conditions[key] = {}
    conditions[key]["name"] = "PeptideLength"
    chk_len = st.checkbox("enable", value=True, key="chk_len")
    conditions[key]["flag"] = chk_len

    if chk_len:
        # if chk_spp:
        # chk_is_propeptide = st.checkbox("protein whole length", value=True)
        fil_peplen = st.slider('Shorter than (a.a.)', 10, 500, 300)
        conditions[key]["value"] = fil_peplen
        conditions[key]["type"] = "under"
        _df["PeptideLength"] = _df["seq"].apply(lambda x: len(x))

        fig = px.histogram(_df["PeptideLength"], x="PeptideLength", nbins=20)
        fig.update_layout(height=180)
        fig.add_vline(x=fil_peplen, line_color="red")
        st.plotly_chart(fig,use_container_width=True,use_container_height=False)


        with col2:
            st.subheader("{} filter: {} {}".format(conditions[key]["name"],
                                               conditions[key]["type"],
                                               conditions[key]["value"]))
            if conditions[key]["type"] == "over":
                _df = _df[_df[conditions[key]["name"]] > conditions[key]["value"]]
            elif conditions[key]["type"] == "under":
                _df = _df[_df[conditions[key]["name"]] < conditions[key]["value"]]
            st.write("{} peptide".format(len(_df)))
            st.dataframe(_df,height=200)

# with col2:
#     st.write("{} peptide".format(len(_df)))
#     for key in conditions.keys():
#         if conditions[key]["flag"]:
#
#             st.subheader("{} filter: {} {}".format(conditions[key]["name"],
#                                                conditions[key]["type"],
#                                                conditions[key]["value"]))
#
#             st.write("{} peptide".format(len(_df)))
#
#             st.dataframe(_df,height=100)

st.header("Final Result")
st.dataframe(_df)
# # full_pos = read_fasta(path = "Arabidopdis_peptide_full.fasta")
# #
#
# pred_pos = []
# pred_neg = []
# for record in full_pos:
#     # id = record.id
#     annot = record.description
#     seq = record.seq
#     # print(all(c(seq) for c in conditions))
#     if all(c(seq) for c in conditions):
#         print(annot)
#         pred_pos.append(annot)
#     else:
#         pred_neg.append(annot)
#
#     # print('id:', id_part)
#     # print('desc:', desc_part)
#     # print('seq:', seq)
# st.write("Pred_positive: {}".format(len(pred_pos)))
# st.write("Pred_negative: {}".format(len(pred_neg)))
# st.write("done")

st.subheader("Debug")
st.write(conditions)