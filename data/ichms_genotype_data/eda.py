
import json
import pandas as pd

def eda_dtc_62_json(file_path: str):
    with open(file_path, "r") as f:
        data = json.load(f)
    rcp_no = data["rcpNo"]
    chart_no = data["chartNo"]
    results = data["results"]
    rows = []
    for result in results:
        # parent result
        parent_result = result["parentResult"]
        item_code = parent_result["itemCode"]
        item_name = parent_result["itemName"]
        item_risk_score = parent_result["detailResult"]
        # child results
        child_results = result["childResults"]
        for child_result in child_results:
            rs_id = child_result["itemCode"]
            gene_name = child_result["itemName"]
            genotype = child_result["result"]
            marker_risk_score = child_result["detailResult"]
            rows.append({
                "rcp_no": rcp_no,
                "chart_no": chart_no,
                "item_code": item_code,
                "item_name": item_name,
                "item_risk_score": item_risk_score,
                "rs_id": rs_id,
                "gene_name": gene_name,
                "genotype": genotype,
                "marker_risk_score": marker_risk_score,
            })
    return pd.DataFrame(rows)

def main():
    dtc_62_df = eda_dtc_62_json("data/ichms_genotype_data/dtc_62.json")
    dtc_62_df.to_csv(
        "data/ichms_genotype_data/dtc_62.csv",
        index=False,
        encoding="utf-8-sig",  # Excel 호환을 위해 BOM 포함
    )



if __name__ == "__main__":
    main()