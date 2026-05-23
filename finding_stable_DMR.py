import pandas as pd

file_dict = {
    'trun_CG_leaf': 'DMR_trun_0.17_WGBS_leaf_CG.txt',
    'trun_CG_callus': 'DMR_trun_0.17_WGBS_callus_CG.txt'}
    # 'trun_CHG_leaf': 'DMR_trun_0.14_WGBS_leaf_CHG.txt',
    # 'trun_CHG_callus': 'DMR_trun_0.14_WGBS_callus_CHG.txt',
    # 'trun_CHH_leaf': 'DMR_trun_0.10_WGBS_leaf_CHH.txt',
    # 'trun_CHH_callus': 'DMR_trun_0.10_WGBS_callus_CHH.txt'}

merge_rules = {
    'trun_CG': ['trun_CG_callus', 'trun_CG_leaf']}
    # 'trun_CHG': ['trun_CHG_callus', 'trun_CHG_leaf'],
    # 'trun_CHH': ['trun_CHH_callus', 'trun_CHH_leaf']}

merged_dataframes = {}

for new_file_name, merge_sources in merge_rules.items():
    merged_df = None
    for source in merge_sources:
        source_df = pd.read_csv(file_dict[source], sep='\t')
        source_df = source_df.drop(["p_value", "Un_mean", "Tr_mean"], 1)

        if merged_df is None:
            merged_df = source_df  
        else:
            merged_df = pd.merge(merged_df, source_df, on=["Chr", "start", "end"], how='inner')
    
    if merged_df is not None:
        merged_df["stable"] = ["Y" if row["type_x"] == row["type_y"] else "N" for _, row in merged_df.iterrows()]
        merged_df = merged_df.rename(columns={"delta_x": "callus_delta_methylation_level", "delta_y": "leaf_delta_methylation_level", "type_x": "callus_type", "type_y": "leaf_type"})

        final_file_name = "stable_newcutoff" + "_0.17" + new_file_name
        merged_dataframes[final_file_name] = merged_df
    
        # Print information about the merged dataframe
        print("Total rows in", final_file_name + ":", len(merged_df))
        print("Rows with stable=Y in", final_file_name + ":", len(merged_df[merged_df["stable"] == "Y"]))

        result_path = final_file_name + ".csv"
        merged_df.to_csv(result_path, sep="\t", index=False)
