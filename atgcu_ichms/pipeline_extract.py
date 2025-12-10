import logging
from pathlib import Path
from typing import List, Dict, Tuple

import pandas as pd

from .io_utils import load_list_file, load_marker_tsv, sanitize_allele, write_oa, write_json
from .qc import QCReporter
from .vcf_parser import VCFParser


def run_extract(conf: dict, args, logger: logging.Logger) -> int:
    product = args.product
    apt_vcf = conf["apt_vcf"]
    imputed_vcf = conf["imputed_vcf"]
    rsmarker_list = conf["rsmarker_list"]
    sample_list_path = args.sample_list
    plate_barcode = args.plate_barcode
    output_dir = Path(conf["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    output_oa_file = output_dir / f"{product}.oa.txt"
    output_json_file = output_oa_file.with_suffix(".json")
    qc_report_path = conf.get("qc_report") or f"{output_oa_file}.qc.report.txt"

    qc = QCReporter(logger)
    qc.add_section(f"{product} Markers Parser QC Report")
    qc.add_line(f"  실행 시간: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")

    sample_list_raw = load_list_file(sample_list_path, logger)
    sample_analysis = qc.analyze_list_file(sample_list_raw, "Sample List")
    sample_list = sample_analysis["unique_items"]

    marker_df = load_marker_tsv(rsmarker_list, logger)
    marker_list = marker_df["RS_ID"].tolist()
    qc.analyze_list_file(marker_list, "RS Marker List")

    qc.add_section("RS Marker TSV 정보")
    qc.add_key_value("총 마커 수", len(marker_df))
    qc.add_line("")
    qc.add_line("  [마커별 Allele 정보 (처음 10개)]")
    for _, row in marker_df.head(10).iterrows():
        qc.add_line(f"    - {row['RS_ID']}: NORMAL={row['NORMAL_ALLELE']}, RISK={row['RISK_ALLELE']}")
    if len(marker_df) > 10:
        qc.add_line(f"    ... 외 {len(marker_df) - 10}개")

    qc.add_section("APT VCF 파싱")
    qc.add_key_value("파일", apt_vcf)
    apt_parser = VCFParser(apt_vcf, logger)
    apt_parser.parse()
    qc.add_key_value("VCF 내 샘플 수", len(apt_parser.get_samples()))
    qc.add_key_value("VCF 내 variant 수", len(apt_parser.data))

    qc.add_section("APT VCF에서 마커 추출")
    apt_results, marker_positions = apt_parser.extract_markers(marker_df, sample_list)
    qc.add_key_value("요청 마커 수", len(marker_list))
    qc.add_key_value("APT VCF에서 찾은 마커 수", len(apt_results))
    missing_in_apt = [m for m in marker_list if m not in apt_results]
    if missing_in_apt:
        qc.add_line("")
        qc.add_line(f"  [APT VCF에 없는 마커: {len(missing_in_apt)}개]")
        for m in missing_in_apt:
            qc.add_line(f"    - {m}")
        logger.info(f"APT VCF에 없는 마커 목록({len(missing_in_apt)}개): {', '.join(missing_in_apt)}")

    qc.add_section("Imputed VCF 파싱")
    qc.add_key_value("파일", imputed_vcf)
    imputed_parser = VCFParser(imputed_vcf, logger)
    imputed_parser.parse_header_only()
    imputed_samples = imputed_parser.get_samples()
    qc.add_key_value("VCF 내 샘플 수", len(imputed_samples))

    needed_positions = set()
    for _, samples_data in apt_results.items():
        for _, geno_info in samples_data.items():
            if geno_info["is_no_call"]:
                chrom = str(geno_info["chrom"]).replace("chr", "")
                pos = str(geno_info["pos"])
                needed_positions.add((chrom, pos))

    missing_positions = {}
    if {"CHROM", "POS"}.issubset(set(marker_df.columns)):
        for _, row in marker_df[marker_df["RS_ID"].isin(missing_in_apt)].iterrows():
            chrom = str(row["CHROM"]).replace("chr", "")
            pos = str(row["POS"])
            missing_positions[row["RS_ID"]] = (chrom, pos)

    missing_position_set = set(missing_positions.values())
    total_positions = needed_positions | missing_position_set

    qc.add_key_value("Imputed 조회 위치 수", len(total_positions))
    if total_positions:
        qc.add_line("")
        qc.add_line("  [Imputed 조회 위치 목록 (chrom:pos)]")
        for chrom, pos in sorted(
            total_positions, key=lambda x: (str(x[0]), int(x[1]) if str(x[1]).isdigit() else str(x[1]))
        ):
            qc.add_line(f"    - {chrom}:{pos}")
        logger.info(
            "Imputed 조회 위치 목록: "
            + ", ".join([f"{c}:{p}" for c, p in sorted(total_positions, key=lambda x: (str(x[0]), str(x[1])))])
        )

    imputed_parser.index_positions(total_positions, target_samples=set(sample_list), missing_rs_ids=set())
    qc.add_key_value("조회된 variant 수", imputed_parser.indexed_variant_count)

    qc.add_section("No Call 처리")
    no_call_count = 0
    filled_from_imputed = 0
    filled_from_normal = 0

    for rs_id, samples_data in apt_results.items():
        for sample_id, geno_info in samples_data.items():
            if geno_info["is_no_call"]:
                no_call_count += 1
                chrom = geno_info["chrom"]
                pos = geno_info["pos"]
                imputed_variant = imputed_parser.get_variant_by_position(chrom, pos)
                if imputed_variant is not None and sample_id in imputed_samples:
                    gt = imputed_variant[sample_id]
                    ref = imputed_variant["REF"]
                    alt = imputed_variant["ALT"]
                    a1, a2 = imputed_parser.decode_genotype(gt, ref, alt)
                    if a1 is not None and a2 is not None:
                        geno_info["allele1"] = a1
                        geno_info["allele2"] = a2
                        geno_info["is_no_call"] = False
                        geno_info["source"] = "IMPUTED"
                        filled_from_imputed += 1
                        continue
                marker_info = marker_df[marker_df["RS_ID"] == rs_id]
                if len(marker_info) > 0:
                    normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                    geno_info["allele1"] = normal_allele
                    geno_info["allele2"] = normal_allele
                    geno_info["is_no_call"] = False
                    geno_info["source"] = "NORMAL_FILLED"
                    filled_from_normal += 1

    imputed_filled_missing = 0
    for rs_id in missing_in_apt:
        pos_info = missing_positions.get(rs_id)
        imputed_variant = None
        if pos_info:
            chrom, pos = pos_info
            imputed_variant = imputed_parser.get_variant_by_position(chrom, pos)
        if imputed_variant is not None:
            chrom = imputed_variant["CHROM"]
            pos = imputed_variant["POS"]
            ref = imputed_variant["REF"]
            alt = imputed_variant["ALT"]
            apt_results[rs_id] = {}
            for sample_id in sample_list:
                if sample_id in imputed_samples:
                    gt = imputed_variant.get(sample_id, "./.")
                    a1, a2 = imputed_parser.decode_genotype(gt, ref, alt)
                    if a1 is not None and a2 is not None:
                        apt_results[rs_id][sample_id] = {
                            "allele1": a1,
                            "allele2": a2,
                            "chrom": chrom,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                            "gt_raw": gt,
                            "is_no_call": False,
                            "source": "IMPUTED",
                        }
                        imputed_filled_missing += 1
                        continue
                marker_info = marker_df[marker_df["RS_ID"] == rs_id]
                if len(marker_info) > 0:
                    normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                    apt_results[rs_id][sample_id] = {
                        "allele1": normal_allele,
                        "allele2": normal_allele,
                        "chrom": ".",
                        "pos": ".",
                        "ref": normal_allele,
                        "alt": ".",
                        "gt_raw": "FILLED",
                        "is_no_call": False,
                        "source": "NORMAL_FILLED",
                    }
        else:
            marker_info = marker_df[marker_df["RS_ID"] == rs_id]
            if len(marker_info) > 0:
                normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                apt_results[rs_id] = {}
                for sample_id in sample_list:
                    apt_results[rs_id][sample_id] = {
                        "allele1": normal_allele,
                        "allele2": normal_allele,
                        "chrom": ".",
                        "pos": ".",
                        "ref": normal_allele,
                        "alt": ".",
                        "gt_raw": "FILLED",
                        "is_no_call": False,
                        "source": "NORMAL_FILLED",
                    }

    qc.add_key_value("총 No Call 수", no_call_count)
    qc.add_key_value("Imputed VCF에서 채운 수", filled_from_imputed)
    qc.add_key_value("NORMAL_ALLELE로 채운 수", filled_from_normal)
    qc.add_key_value("Imputed로 채운 누락 RS 수", imputed_filled_missing)

    for rs_id in missing_in_apt:
        if rs_id not in apt_results:
            marker_info = marker_df[marker_df["RS_ID"] == rs_id]
            if len(marker_info) > 0:
                normal_allele = marker_info.iloc[0]["NORMAL_ALLELE"]
                apt_results[rs_id] = {}
                for sample_id in sample_list:
                    apt_results[rs_id][sample_id] = {
                        "allele1": normal_allele,
                        "allele2": normal_allele,
                        "chrom": ".",
                        "pos": ".",
                        "ref": normal_allele,
                        "alt": ".",
                        "gt_raw": "FILLED",
                        "is_no_call": False,
                        "source": "NORMAL_FILLED",
                    }

    qc.add_section("OA 형식 출력")
    oa_rows = []
    for rs_id in marker_list:
        if rs_id in apt_results:
            for sample_id in sample_list:
                if sample_id in apt_results[rs_id]:
                    geno_info = apt_results[rs_id][sample_id]
                    source = geno_info.get("source", "APT")
                    oa_row = {
                        "Sample ID": sample_id,
                        "Plate Barcode": plate_barcode,
                        "Gene Symbol": source,
                        "NCBI SNP Reference": rs_id,
                        "Assay Name or ID": rs_id,
                        "Allele 1 Call": geno_info["allele1"],
                        "Allele 2 Call": geno_info["allele2"],
                    }
                    oa_rows.append(oa_row)

    result_df = pd.DataFrame(oa_rows)
    result_df["Sample ID"] = pd.Categorical(result_df["Sample ID"], categories=sample_list, ordered=True)
    result_df["NCBI SNP Reference"] = pd.Categorical(
        result_df["NCBI SNP Reference"], categories=marker_list, ordered=True
    )
    result_df = result_df.sort_values(["Sample ID", "NCBI SNP Reference"]).reset_index(drop=True)

    json_rows = []
    for rs_id in marker_list:
        if rs_id not in apt_results:
            continue
        for sample_id in sample_list:
            if sample_id not in apt_results[rs_id]:
                continue
            geno_info = apt_results[rs_id][sample_id]
            a1 = sanitize_allele(geno_info["allele1"])
            a2 = sanitize_allele(geno_info["allele2"])
            json_rows.append(
                {
                    "sampleId": sample_id,
                    "plateBarcode": plate_barcode,
                    "itemCd": rs_id,
                    "result1": a1,
                    "result2": a2,
                }
            )

    qc.add_key_value("총 출력 row 수", len(result_df))
    qc.add_key_value("샘플 수", len(sample_list))
    qc.add_key_value("마커 수", len(marker_list))

    output_path = output_oa_file
    write_oa(
        output_path,
        meta={
            "Source": f"{product} Markers Parser",
            "Mode": product,
            "APT VCF": apt_vcf,
            "Imputed VCF": imputed_vcf,
            "Marker List": rsmarker_list,
            "Sample List": sample_list_path,
            "Plate Barcode": plate_barcode,
            "Generated Date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        },
        data=result_df,
        logger=logger,
    )
    write_json(output_json_file, json_rows, logger)

    qc.add_section("처리 완료")
    qc.add_key_value("출력 파일", str(output_oa_file))
    qc.add_key_value("JSON 파일", str(output_json_file))
    qc.add_key_value("상태", "SUCCESS")
    qc.print_report()
    qc.write_report(qc_report_path)
    return 0

