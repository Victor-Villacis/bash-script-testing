
import pandas as pd
import numpy as np
from glob import glob


filepath = glob("*.xlsx")
if len(filepath)==0:
    filepath = glob("*.xls")

print(filepath)
rawdata_df_list = []
for fn in filepath:
    if fn.find("CNV")>-1:
        temp_df = pd.read_excel(fn,header=73,dtype=str,sheet_name="Results")
        temp_df = temp_df[["Sample Name","Target","CN Predicted"]]
        temp_df.columns = ["Sample Name","SNP Assay Name","Call"]
        temp_df["SNP Assay Name"] = "CYP2D6 CNV"
        temp_df["Call"] = temp_df["Call"].fillna("0")
    else:
        temp_df = pd.read_excel(fn,header=46,dtype=str,sheet_name="Results")
        temp_df = temp_df[["Sample Name","SNP Assay Name","Call"]]
    print(fn)
    rawdata_df_list.append(temp_df)

combined_rawdata_df = pd.concat(rawdata_df_list)
combined_rawdata_df["Call"] = combined_rawdata_df["Call"].str.replace("Heterozygous ","")
combined_rawdata_df["Call"] = combined_rawdata_df["Call"].str.replace("Homozygous ","")
combined_rawdata_df = combined_rawdata_df.drop_duplicates()
sample_list1 = combined_rawdata_df["Sample Name"].unique().tolist()
combined_rawdata_df = combined_rawdata_df[combined_rawdata_df["Call"]!="Undetermined"]
combined_rawdata_df = combined_rawdata_df.dropna()
sample_list2 = combined_rawdata_df["Sample Name"].unique().tolist()
failed_samples = [sample for sample in sample_list1 if sample not in sample_list2]
if len(failed_samples)>0:
    print(f"Samples with no data: {failed_samples}")

combined_rawdata_df["check"] = combined_rawdata_df["Sample Name"] + "_" + combined_rawdata_df["SNP Assay Name"]
counts = combined_rawdata_df[["check","Call"]].groupby("check").count()
counts = pd.DataFrame(counts.to_records())
counts["sample"] = counts["check"].str.split("_",expand=True)[0]
failed_samples = counts[counts["Call"]>1]["sample"].unique().tolist()
del combined_rawdata_df["check"]
if len(failed_samples) > 1:
    print(f"Samples with inconsistent Data: {failed_samples}")
    combined_rawdata_df = combined_rawdata_df[~combined_rawdata_df["Sample Name"].isin(failed_samples)]

#combined_rawdata_df.to_csv("interm_file.csv",header=True,index=False,sep=",")
formatteddata_df = combined_rawdata_df.pivot(index="Sample Name",columns="SNP Assay Name",values="Call")
formatteddata_df = formatteddata_df.fillna("Undetermined")
formatteddata_df.to_csv("formatted_qpcr.csv",header=True,index=True,sep=",")